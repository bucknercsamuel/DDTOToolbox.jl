# ..:: Top-level Solve Function ::..

function solve(params)
    # ..:: Execute solver sequence ::..
    @time begin
        @time begin
            # ..:: Solve for independently-optimal solutions to each target ::..
            scp_solutions = solve_tree_decoupled(params)
            scp_costs = CVector(zeros(params.n_targs))
            for k = 1:params.n_targs
                scp_costs[k] = scp_solutions.targs[k].cost
            end
            println("\n Solve time for generating optimal solutions to each target:")
        end

        @time begin
            # ..:: Solve for DDTO branching solutions to ALL targets ::..
            (_, ddtoscp_solutions) = solve_tree_ddto(params, scp_costs)
            println("\n Solve time for generating DDTO branch solutions to all targets:")
        end
        println("\n Solve time for the full DDTO solution stack:")
    end

    # ..:: Simulate each target solution from I.C. to T.C.
    @time begin
        dynamics = (t,x,sol) -> dynamics_nonlinear(t,x,optimal_controller(t,sol.t,sol.u,params.disc),params)
        scp_simulations = simulate(scp_solutions, dynamics, params.disc)
        ddtoscp_simulations = simulate(ddtoscp_solutions, dynamics, params.disc)
        println("\n Solve time for RK4 simulation:")
    end

    # ..:: Post-processing (problem-specific) ::..
    @time begin
        scp_solutions_proc       = process_solutions(scp_solutions, params)
        scp_simulations_proc     = process_solutions(scp_simulations, params)
        ddtoscp_solutions_proc   = process_solutions(ddtoscp_solutions, params)
        ddtoscp_simulations_proc = process_solutions(ddtoscp_simulations, params)
        println("\n Solve time for post-processing:")
    end

    return (
        scp_solutions_proc, 
        scp_simulations_proc, 
        ddtoscp_solutions_proc, 
        ddtoscp_simulations_proc)
end

# ..:: DDTO-SCP Solver Functions ::..

function solve_tree_ddto(params, ref_costs::CVector; single_iter=false, ref_trajs=nothing)::Tuple{Bool,DDTOSolution}

    # Set node deferrability allocation (may have already been set, overrides)
    set_deferrability_node_allocation!(params)

    # Obtain initial guess for reference trajectories
    if isnothing(ref_trajs)
        ref_trajs = generate_initial_guess_ddtoscp(params)
    end  

    # SCP Iteration
    feas_status = undef
    t_defer = zeros(params.n_targs)
    solution = ref_trajs
    scp_converged = false
    iteration_cap_reached = true
    VERB_OPT && println("\n=== DDTO-SCP Iteration ===")
    for k = 1:params.scp_iters

        # Solve SCP subproblem
        (solution, feas_status, scp_converged, t_defer) = solve_subproblem_ddto(params, ref_costs, solution, k)

        if single_iter
            iteration_cap_reached = false
            scp_converged = true # flag SCP as converged even if it hasn't
            break # skip all convergence criterion, only going to run a single (potentially-infeasible) iterate!
        end

        if feas_status != MOI.OPTIMAL && feas_status != MOI.ALMOST_OPTIMAL
            iteration_cap_reached = false
            scp_converged = false
            @printf("   ! SCP subproblem is infeasible (MOI status: %s), exiting subproblem iteration.\n", feas_status)
            break
        end
        if scp_converged
            iteration_cap_reached = false
            scp_converged = true
            @printf("   > Convergence condition has been reached, exiting subproblem iteration.\n")
            break
        end
    end

    if iteration_cap_reached
        scp_converged = false
        println("   ! SCP subproblem iteration cap reached, exiting subproblem iteration.")
    end

    # Converged solution data
    println("\nDDTO solution properties:")
    for j = 1:params.n_targs
        ϵ_subopt = (solution.targs[j].cost - ref_costs[j])/ref_costs[j] * 100
        @printf("   Target %i -- %2.1f [s] deferred, % 2.1f [%%] suboptimal.\n", j, t_defer[j], ϵ_subopt)
    end 


    return (scp_converged, solution)
end

function solve_subproblem_ddto(params, ref_costs::CVector, ref_trajs::DDTOSolution, scp_iter::Int)::Tuple{DDTOSolution, MOI.TerminationStatusCode, Bool, Vector}
    # Solve the baseline feasibility problem for DDTO.
    #
    # :in params: The params object
    # :in ref_costs: Optimal costs from `solve_optimal_pdg_all_targets`
    # :out ddto_solution: Contains the DDTO solution for this target/branch point
    # :out feas_status: Feasibility problem solution status code (see MOI.TerminationStatusCode documentation)

    # ..:: Setup ::..
    # Optimizer configuration
    if SOLVER == "ECOS"
        mdl = Model(optimizer_with_attributes(ECOS.Optimizer, "verbose" => 0, "max_iters" => 1000))
    elseif SOLVER == "MOSEK"
        mdl = Model(Mosek.Optimizer)
        JuMP.set_optimizer_attribute(mdl, "LOG", 0) # disable debugging
        JuMP.set_optimizer_attribute(mdl, "MAX_NUM_WARNINGS", 0) # disable warni
    else
        error("SOLVER is invalid, please select either ECOS or MOSEK")
    end

    # Sizing parameters
    n = params.n_targs
    N = params.N
    nx = params.nx
    nu = params.nu
    if params.disc == 0
        error("Zero-order hold not currently supported.")
    elseif params.disc != 0 && params.disc != 1
        error("Please select a valid discretization hold order.")
    end
    τ_max = max(params.τ_targs...)
    τ_lu(j) = params.τ_targs[findfirst(i->i==j, params.λ_targs)] # obtain the deferrability index in the trunk of the j-th target

    # Dynamics functions
    dyn_lin = (t,x,u,p) -> dynamics_linearized(t,x,u,params)
    dyn_nl  = (t,x,u,p) -> dynamics_nonlinear(t,x,u,params)

    # ..:: Optimization variables ::..
    # Unscaled variables
    x_trunk_us  = @variable(mdl, [1:nx,1:τ_max])
    u_trunk_us  = @variable(mdl, [1:nu,1:τ_max])
    x_branch_us = Vector{Matrix{JuMP.VariableRef}}(undef,n)  
    u_branch_us = Vector{Matrix{JuMP.VariableRef}}(undef,n)
    for j = 1:n
        τ = τ_lu(j)
        x_branch_us[j] = @variable(mdl, [1:nx,1:N-τ])
        u_branch_us[j] = @variable(mdl, [1:nu,1:N-τ])
    end

    # Apply affine scaling
    x_trunk = params.Sx*x_trunk_us .+ repeat(params.sx, 1, τ_max)
    u_trunk = params.Su*u_trunk_us .+ repeat(params.su, 1, τ_max)
    x_branch = Vector{Matrix{JuMP.AffExpr}}(undef,n)
    u_branch = Vector{Matrix{JuMP.AffExpr}}(undef,n)
    for j = 1:n
        τ = τ_lu(j)
        x_branch[j] = params.Sx*x_branch_us[j] .+ repeat(params.sx, 1, N-τ)
        u_branch[j] = params.Su*u_branch_us[j] .+ repeat(params.su, 1, N-τ)
    end

    # SCP-specific
    ν_ctrl_trunk = @variable(mdl, [1:nx,1:τ_max-1])
    ν_ctrl_branch = Vector{Matrix{JuMP.VariableRef}}(undef,n) 
    ν_ctrl_stitch = @variable(mdl, [1:nx,1:n])
    η_trunk = @variable(mdl, [1:τ_max])
    η_branch = Vector{Vector{JuMP.VariableRef}}(undef,n)
    for j = 1:n
        τ = τ_lu(j)
        ν_ctrl_branch[j] = @variable(mdl, [1:nx,1:N-τ-1])
        η_branch[j] = @variable(mdl, [1:N-τ])
    end
    @variable(mdl, μ_ctrl) # virtual control slack
    @variable(mdl, μ_buff) # virtual buffer slack
    @variable(mdl, η_s) # trust region slack

    # Convenience functions
    X_trunk(k) = x_trunk[:,k]
    U_trunk(k) = u_trunk[:,k]
    X_branch(k,j) = x_branch[j][:,k]
    U_branch(k,j) = u_branch[j][:,k]

    # ..:: Trunk Constraints ::..
    # Build the trunk
    # Take last deferred target trajectory as the reference
    ref_traj_trunk = copy(ref_trajs.targs[params.λ_targs[end]])
    ref_traj_trunk.t = ref_traj_trunk.t[1:τ_max]
    ref_traj_trunk.x = ref_traj_trunk.x[:,1:τ_max]
    ref_traj_trunk.u = ref_traj_trunk.u[:,1:τ_max]

    # Core constraints
    J_obj_trunk,ν_buff_trunk = core_problem(mdl,x_trunk,u_trunk,params,ref_traj_trunk)

    # Dynamics
    Ak,Bmk,Bpk,_,wk,_,_ = c2d_nonlinear(ref_traj_trunk.t,ref_traj_trunk.x,ref_traj_trunk.u,dyn_nl,dyn_lin,params.disc)
    SxInv = inv(params.Sx)
    SuInv = inv(params.Su)
    if params.disc == 0
        @constraint(mdl, [k=1:τ_max-1], SxInv*X_trunk(k+1) .== SxInv*(Ak[:,:,k]*X_trunk(k) + Bmk[:,:,k]*U_trunk(k) + wk[:,k]) + ν_ctrl_trunk[:,k])
    elseif params.disc == 1
        @constraint(mdl, [k=1:τ_max-1], SxInv*X_trunk(k+1) .== SxInv*(Ak[:,:,k]*X_trunk(k) + Bmk[:,:,k]*U_trunk(k) + Bpk[:,:,k]*U_trunk(k+1) + wk[:,k]) + ν_ctrl_trunk[:,k])
    end

    # Trunk time definition
    s_trunk = u_trunk[end,:]
    t_trunk = time_dilation_control_to_wall_clock_time(s_trunk, ref_traj_trunk.t, params.disc)
    @constraint(mdl, [k=1:τ_max-1], s_trunk[k+1] == s_trunk[k])

    # Trust region
    δXt(k) = SxInv*(X_trunk(k) .- ref_traj_trunk.x[:,k])
    δUt(k) = SuInv*(U_trunk(k) .- ref_traj_trunk.u[:,k])
    @constraint(mdl, [k=1:τ_max], δXt(k)'*δXt(k) + δUt(k)'*δUt(k) <= η_trunk[k])

    # ..:: Branch/Target Constraints ::..
    J_obj_branch = Vector{JuMP.AffExpr}(undef,n)
    ν_buff_branch = Vector{Vector{JuMP.AffExpr}}(undef,n)
    for j = 1:n
        # Take jth reference and build it with last N elements
        τ = τ_lu(j)
        ref_traj_branch = copy(ref_trajs.targs[j])
        ref_traj_branch.t = ref_traj_branch.t[τ+1:end]
        ref_traj_branch.x = ref_traj_branch.x[:,τ+1:end]
        ref_traj_branch.u = ref_traj_branch.u[:,τ+1:end]

        # Core constraints
        J_obj_branch_,ν_buff_branch_ = core_problem(mdl, x_branch[j], u_branch[j], params, ref_traj_branch)
        J_obj_branch[j] = J_obj_branch_
        ν_buff_branch[j] = ν_buff_branch_

        # Dynamics (within branch)
        Ak,Bmk,Bpk,_,wk,_,_ = c2d_nonlinear(ref_traj_branch.t,ref_traj_branch.x,ref_traj_branch.u,dyn_nl,dyn_lin,params.disc)
        if params.disc == 0
            @constraint(mdl, [k=1:N-τ-1], SxInv*X_branch(k+1,j) .== SxInv*(Ak[:,:,k]*X_branch(k,j) + Bmk[:,:,k]*U_branch(k,j) + wk[:,k]) + ν_ctrl_branch[j][:,k])
        elseif params.disc == 1
            @constraint(mdl, [k=1:N-τ-1], SxInv*X_branch(k+1,j) .== SxInv*(Ak[:,:,k]*X_branch(k,j) + Bmk[:,:,k]*U_branch(k,j) + Bpk[:,:,k]*U_branch(k+1,j) + wk[:,k]) + ν_ctrl_branch[j][:,k])
        end

        # Dynamics (stitching to trunk)
        ref_traj_stitch = copy(ref_trajs.targs[j])
        ref_traj_stitch.t = ref_traj_stitch.t[τ:τ+1]
        ref_traj_stitch.x = ref_traj_stitch.x[:,τ:τ+1]
        ref_traj_stitch.u = ref_traj_stitch.u[:,τ:τ+1]
        Ak,Bmk,Bpk,_,wk,_,_ = c2d_nonlinear(ref_traj_stitch.t,ref_traj_stitch.x,ref_traj_stitch.u,dyn_nl,dyn_lin,params.disc)
        if params.disc == 0
            @constraint(mdl, SxInv*X_branch(1,j) .== SxInv*(Ak[:,:,1]*X_trunk(τ) + Bmk[:,:,1]*U_trunk(τ) + wk[:,1]) + ν_ctrl_stitch[:,j])
        elseif params.disc == 1
            @constraint(mdl, SxInv*X_branch(1,j) .== SxInv*(Ak[:,:,1]*X_trunk(τ) + Bmk[:,:,1]*U_trunk(τ) + Bpk[:,:,1]*U_branch(1,j) + wk[:,1]) + ν_ctrl_stitch[:,j])
        end

        # Suboptimality constraint
        @constraint(mdl, sum(J_obj_trunk) + sum(J_obj_branch[j]) <= (1 + params.ϵ_targs[j]) * ref_costs[j])

        # Time dilation constraints (from IC to TC for each target)
        s_branch = u_branch[j][end,:]
        t_target = time_dilation_control_to_wall_clock_time([s_trunk[1:τ];s_branch], ref_trajs.targs[j].t, params.disc)
        Δt_target = diff(t_target)
        @constraint(mdl, [k=1:N-1], params.Δt_min/params.Δt_max <= Δt_target[k]/params.Δt_max <= 1)
        @constraint(mdl, t_target[end]/params.ToF_max <= 1)
        @constraint(mdl, [k=1:N-τ-1], s_branch[k+1] == s_branch[k])

        # Trust region
        δXb(k) = SxInv*(X_branch(k,j) .- ref_traj_branch.x[:,k])
        δUb(k) = SuInv*(U_branch(k,j) .- ref_traj_branch.u[:,k])
        @constraint(mdl, [k=1:N-τ], δXb(k)'*δXb(k) + δUb(k)'*δUb(k) <= η_branch[j][k])
    end

    # ..:: Boundary Conditions ::..
    # Note: inf = no boundary condition to be applied
    # Initial conditions
    for k = 1:nx
        if ~isinf(params.z0[k])
            @constraint(mdl, x_trunk[k,1] == params.z0[k])
        end
    end

    # Terminal conditions
    for j = 1:n
        for k = 1:nx
            if ~isinf(params.zf_targs[k,j])
                @constraint(mdl, x_branch[j][k,end] == params.zf_targs[k,j])
            end
        end
    end

    # ..:: Slack Constraints ::..
    ν_ctrl = [vec(ν_ctrl_trunk); vec.(ν_ctrl_branch)...; vec(ν_ctrl_stitch)]
    ν_buff = [vec(ν_buff_trunk); vec.(ν_buff_branch)...]
    @constraint(mdl, vcat(μ_ctrl, ν_ctrl) in MOI.NormOneCone(length(ν_ctrl)+1))
    @constraint(mdl, vcat(μ_buff, ν_buff) in MOI.NormOneCone(length(ν_buff)+1))
    @constraint(mdl, vcat(η_s, [vec(η_trunk); vec.(η_branch)...]) in SecondOrderCone())
    @constraint(mdl, μ_buff >= 0)
    @constraint(mdl, μ_ctrl >= 0)
    @constraint(mdl, η_s >= 0)

    # ..:: Construct cost function and solve ::..
    J_opt  = -sum([params.α_targs[j]*t_trunk[τ_lu(j)] for j=1:n])/max(params.α_targs...)
    obj_scale = 1/sqrt(max(params.w_obj, params.w_trust, params.w_buff, params.w_ctrl))
    @objective(mdl, Min, 
        (params.w_obj * J_opt 
      + params.w_trust * η_s 
      + params.w_buff * μ_buff 
      + params.w_ctrl * μ_ctrl)*obj_scale)

    optimize!(mdl)
    feas_status = JuMP.termination_status(mdl)

    # ..:: Extract the solution ::..
    x = Vector{CMatrix}(undef,n)
    u = Vector{CMatrix}(undef,n)
    for j = 1:n
        τ = τ_lu(j)
        x[j] = hcat(value.(x_trunk[:,1:τ]), reshape(value.(x_branch[j]),nx,N-τ))
        u[j] = hcat(value.(u_trunk[:,1:τ]), reshape(value.(u_branch[j]),nu,N-τ))
    end
    costs_sol = [sum(value.(J_obj_trunk)) + sum(value.(J_obj_branch[j])) for j = 1:n]

    # ..:: Determine if PTR subproblem has converged ::..
    if feas_status == MOI.OPTIMAL || feas_status == MOI.ALMOST_OPTIMAL
        solve_status = "Feasible"
    else
        solve_status = "Not Feasible"
    end
    
    # Obtain evaluation penalties
    μ_buff_pen = value.(μ_buff)
    μ_ctrl_pen = value.(μ_ctrl)
    η_pen = value.(η_s)

    if (feas_status == MOI.OPTIMAL || feas_status == MOI.ALMOST_OPTIMAL) && (μ_ctrl_pen <= params.ϵ_ctrl) && (μ_buff_pen <= params.ϵ_buff) && (η_pen <= params.ϵ_trust)
        scp_sub_cvged = true
    else
        scp_sub_cvged = false
    end
    @printf("   SCP Iter: %2.i | Status: %s | Cost = %.2e | μ_ctrl_pen = %.2e | μ_buff_pen = %.2e | η_pen = %.2e\n", scp_iter, solve_status, value.(J_opt), μ_ctrl_pen, μ_buff_pen, η_pen)
    flush(stdout)

    # ..:: Package the DDTO Solution ::..
    ddto_solution = EmptyDDTOSolution(n)
    for j = 1:n
        τ = τ_lu(j)
        ddto_solution.targs[j].t = ref_trajs.targs[j].t # maintain reference dilated time
        ddto_solution.targs[j].x = x[j]
        ddto_solution.targs[j].u = u[j]
        ddto_solution.targs[j].cost = costs_sol[j]
    end
    deferrability_times = [value.(t_trunk)[τ_lu(j)] for j=1:n]

    return (ddto_solution, feas_status, scp_sub_cvged, deferrability_times)

end

function set_deferrability_node_allocation!(params)
    # Set deferrability node allocation based on uniform distribution up to N/sqrt(2)
    τ_alloc = round.(CVector(range(1,params.N,Int(round(sqrt(2)*params.n_targs))+1)))[2:2+params.n_targs]
    if length(unique(τ_alloc)) < length(τ_alloc)
        # some targets have the same deferrability index due to rounding
        # attempt to instead space by 1 node
        τ_alloc = range(2,2+params.n_targs,params.n_targs) |> Vector
        if τ_alloc[end] > params.N
            error("There are more targets than number of knot points; adjust parameters accordingly!")
        end
    end
    params.τ_targs = τ_alloc
end