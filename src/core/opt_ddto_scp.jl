# ..:: Top-level Solve Function ::..

function solve(params; single_iter=false, ref_trajs=nothing, simulate_solutions=true, process_the_solutions=true)
    # ..:: Customized problem modification ::..
    # Apply custom scaling (if not already done)
    custom_scaling!(params)
    
    # Modify for extra constraint violation state (if CTCS enabled)
    if params.a.ctcs_enabled
        Sϵ,sϵ = scaling_matrices([0],[params.a.ϵ_ctcs])
        params.a.Sx = [params.a.Sx zeros(params.a.nx,1); zeros(1,params.a.nx) Sϵ]
        params.a.sx = vcat(params.a.sx, sϵ)
        params.a.z0 = vcat(params.a.z0, Inf)
        params.a.zf_targs = vcat(params.a.zf_targs, Inf*ones(1,params.a.n_targs))
        params.a.nx += 1
    end

    # Modify for extra time dilation input
    Δτ = 1/(params.a.N-1)
    Ss,ss = scaling_matrices([params.a.Δt_min/Δτ], [params.a.Δt_max/Δτ])
    params.a.Su = [params.a.Su zeros(params.a.nu,1); zeros(1,params.a.nu) Ss]
    params.a.su = vcat(params.a.su, ss)
    params.a.u0 = vcat(params.a.u0, Inf)
    params.a.uf_targs = vcat(params.a.uf_targs, Inf*ones(1,params.a.n_targs))
    params.a.nu += 1

    # ..:: Customized warmstart method ::..
    ref_trajs_cvx = copy(ref_trajs)
    ref_trajs_ddtocvx = copy(ref_trajs)
    if params.a.warmstart_method == "ddto"
        # Remove augmented terms
        Sx_ = copy(params.a.Sx)
        sx_ = copy(params.a.sx)
        Su_ = copy(params.a.Su)
        su_ = copy(params.a.su)
        if params.a.ctcs_enabled
            params.a.nx -= 1
            params.a.Sx = params.a.Sx[1:params.a.nx,1:params.a.nx]
            params.a.sx = params.a.sx[1:params.a.nx]
        end
        params.a.nu -= 1
        params.a.Su = params.a.Su[1:params.a.nu,1:params.a.nu]
        params.a.su = params.a.su[1:params.a.nu]

        # Compute ddto-cvx solution
        ref_trajs_cvx_,ref_trajs_ddtocvx_ = solve_cvx(params; simulate_solutions=false, process_the_solutions=false)
        
        # Add augmented terms back
        params.a.nx = params.a.ctcs_enabled ? params.a.nx + 1 : params.a.nx
        params.a.nu += 1
        params.a.Sx = Sx_
        params.a.sx = sx_
        params.a.Su = Su_
        params.a.su = su_

        # Generate a full initial guess (w/ augmented terms) and add convex solution elements
        ref_trajs = generate_initial_guess_scp(copy(params))

        # > cvx
        ref_trajs_cvx = copy(ref_trajs)
        o = params.a.ctcs_enabled ? 1 : 0
        for j = 1:params.a.n_targs
            ref_trajs_cvx.targs[j].x[1:end-o,:] = ref_trajs_cvx_.targs[j].x
            ref_trajs_cvx.targs[j].u[1:end-1,:] = ref_trajs_cvx_.targs[j].u
            ref_trajs_cvx.targs[j].u[end,:] = wall_clock_time_to_time_dilation_control(ref_trajs_cvx_.targs[j].t, ref_trajs_cvx.targs[j].t, params.a.disc)
        end

        # > ddto-cvx
        ref_trajs_ddtocvx = copy(ref_trajs)
        for j = 1:params.a.n_targs
            ref_trajs_ddtocvx.targs[j].x[1:end-o,:] = ref_trajs_ddtocvx_.targs[j].x
            ref_trajs_ddtocvx.targs[j].u[1:end-1,:] = ref_trajs_ddtocvx_.targs[j].u
            ref_trajs_ddtocvx.targs[j].u[end,:] = wall_clock_time_to_time_dilation_control(ref_trajs_ddtocvx_.targs[j].t, ref_trajs_ddtocvx.targs[j].t, params.a.disc)
        end
    elseif params.a.warmstart_method == "single"
        # Remove augmented terms
        Sx_ = copy(params.a.Sx)
        sx_ = copy(params.a.sx)
        Su_ = copy(params.a.Su)
        su_ = copy(params.a.su)
        if params.a.ctcs_enabled
            params.a.nx -= 1
            params.a.Sx = params.a.Sx[1:params.a.nx,1:params.a.nx]
            params.a.sx = params.a.sx[1:params.a.nx]
        end
        params.a.nu -= 1
        params.a.Su = params.a.Su[1:params.a.nu,1:params.a.nu]
        params.a.su = params.a.su[1:params.a.nu]

        # Compute cvx solution
        ref_trajs_cvx_ = solve_cvx(params; simulate_solutions=false, process_the_solutions=false, solve_ddto=false)
        
        # Add augmented terms back
        params.a.nx = params.a.ctcs_enabled ? params.a.nx + 1 : params.a.nx
        params.a.nu += 1
        params.a.Sx = Sx_
        params.a.sx = sx_
        params.a.Su = Su_
        params.a.su = su_

        # Generate a full initial guess (w/ augmented terms) and add convex solution elements
        ref_trajs = generate_initial_guess_scp(copy(params))

        # > cvx
        ref_trajs_cvx = copy(ref_trajs)
        o = params.a.ctcs_enabled ? 1 : 0
        for j = 1:params.a.n_targs
            ref_trajs_cvx.targs[j].x[1:end-o,:] = ref_trajs_cvx_.targs[j].x
            ref_trajs_cvx.targs[j].u[1:end-1,:] = ref_trajs_cvx_.targs[j].u
            ref_trajs_cvx.targs[j].u[end,:] = wall_clock_time_to_time_dilation_control(ref_trajs_cvx_.targs[j].t, ref_trajs_cvx.targs[j].t, params.a.disc)
        end
    end

    # ..:: Solve for independently-optimal solutions to each target ::.
    if params.a.warmstart_method == "single"
        ref_trajs_scp = ref_trajs_cvx
    elseif params.a.warmstart_method == "ddto"
        ref_trajs_scp = ref_trajs_cvx
    else
        ref_trajs_scp = generate_initial_guess_scp(params)
    end
    @time begin
        scp_solutions, scp_converged = solve_tree_decoupled(params; single_iter=single_iter, ref_trajs=ref_trajs_scp)
        scp_costs = CVector(zeros(params.a.n_targs))
        for k = 1:params.a.n_targs
            scp_costs[k] = scp_solutions.targs[k].cost
        end
        println("\n Solve time for generating optimal solutions to each target:")
    end

    # ..:: Solve for DDTO branching solutions to ALL targets ::..
    if params.a.warmstart_method == "single"
        ref_trajs_ddtoscp = scp_solutions
    elseif params.a.warmstart_method == "ddto"
        ref_trajs_ddtoscp = ref_trajs_ddtocvx
    else
        ref_trajs_ddtoscp = generate_initial_guess_ddtoscp(params)
    end
    set_deferrability_node_allocation!(params)
    if params.a.n_targs > 1
        @time begin
            ddtoscp_solutions, ddtoscp_converged = solve_tree_ddto(params, scp_costs; single_iter=single_iter, ref_trajs=ref_trajs_ddtoscp)
            println("\n Solve time for generating DDTO branch solutions to all targets:")
        end
        println("\n Solve time for the full DDTO solution stack:")
    else
        ddtoscp_solutions = copy(scp_solutions)
        ddtoscp_converged = copy(scp_converged)
    end

    # ..:: Simulate each target solution from I.C. to T.C.
    if simulate_solutions
        @time begin
            if params.a.ctcs_enabled
                dynamics = (t,x,sol) -> dynamics_nonlinear_ctcs(t,x,optimal_controller(t,sol.t,sol.u,params.a.disc),params,0)
            else
                dynamics = (t,x,sol) -> dynamics_nonlinear(t,x,optimal_controller(t,sol.t,sol.u,params.a.disc),params)
            end
            scp_simulations = simulate(scp_solutions, dynamics, params.a.disc; max_steps=params.a.N_sim)
            ddtoscp_simulations = simulate(ddtoscp_solutions, dynamics, params.a.disc; max_steps=params.a.N_sim)
            println("\n Solve time for RK4 simulation:")
        end
    end

    # ..:: Post-processing (problem-specific) ::..
    if process_the_solutions
        @time begin
            scp_solutions       = process_solutions(scp_solutions, params)
            ddtoscp_solutions   = process_solutions(ddtoscp_solutions, params)
            if simulate_solutions
                scp_simulations     = process_solutions(scp_simulations, params)
                ddtoscp_simulations = process_solutions(ddtoscp_simulations, params)
            end
            println("\n Solve time for post-processing:")
        end
    end

    # Undo dynamic sizing changes
    if params.a.ctcs_enabled
        params.a.nx -= 1
    end
    params.a.nu -= 1

    converged = scp_converged && ddtoscp_converged ? true : false
    if simulate_solutions
        return (
            scp_solutions, 
            scp_simulations, 
            ddtoscp_solutions, 
            ddtoscp_simulations,
            converged)
    else
        return (
            scp_solutions, 
            ddtoscp_solutions,
            converged)
    end
end

# ..:: DDTO-SCP Solver Functions ::..

function solve_tree_ddto(params, ref_costs::CVector; single_iter=false, ref_trajs=nothing)::Tuple{DDTOSolution,Bool}

    # Obtain initial guess for reference trajectories
    if isnothing(ref_trajs)
        ref_trajs = generate_initial_guess_ddtoscp(params)
    end  

    # SCP Iteration
    feas_status = undef
    t_defer = zeros(params.a.n_targs)
    solution = ref_trajs
    scp_converged = false
    iteration_cap_reached = true
    params_ = copy(params)
    VERB_OPT && println("\n=== DDTO-SCP Iteration ===")
    for k = 1:params.a.scp_iters

        # Solve SCP subproblem
        (solution, feas_status, scp_converged, t_defer) = solve_subproblem_ddto(params_, ref_costs, solution, k)

        # Update problem parameters
        param_update_law!(params_)

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
            VERB_DDTO && @printf("   > Convergence condition has been reached, exiting subproblem iteration.\n")
            break
        end
    end

    if iteration_cap_reached
        scp_converged = false
        println("   ! SCP subproblem iteration cap reached, exiting subproblem iteration.")
    end

    # Converged solution data
    println("\nDDTO solution properties:")
    for j = 1:params.a.n_targs
        ϵ_subopt = (solution.targs[j].cost - ref_costs[j])/ref_costs[j] * 100
        @printf("   Target %i -- %2.1f [s] deferred, % 2.1f [%%] suboptimal.\n", j, t_defer[j], ϵ_subopt)
    end 

    return solution, scp_converged
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
    if params.a.ctcs_enabled
        mdl, solver_type = solver_setup(SOLVER_CTCS_ENABLED)
    else
        mdl, solver_type = solver_setup(SOLVER_CTCS_DISABLED)
    end
    trust_region_type = solver_type
    # trust_region_type = "QP"

    # Base parameters
    n = params.a.n_targs
    N = params.a.N
    nx = params.a.nx
    nu = params.a.nu
    if params.a.disc == 0
        error("Zero-order hold not currently supported.")
    elseif params.a.disc != 0 && params.a.disc != 1
        error("Please select a valid discretization hold order.")
    end
    τ_max = max(params.a.τ_targs...)
    τ_lu(j) = params.a.τ_targs[findfirst(i->i==j, params.a.λ_targs)] # obtain the deferrability index in the trunk of the j-th target

    # Dynamics functions (non CTCS)
    if ~params.a.ctcs_enabled
        dyn_lin = (t,x,u,p) -> dynamics_linearized(t,x,u,params)
        dyn_nl  = (t,x,u,p) -> dynamics_nonlinear(t,x,u,params)
    end

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
    x_trunk = params.a.Sx*x_trunk_us .+ repeat(params.a.sx, 1, τ_max)
    u_trunk = params.a.Su*u_trunk_us .+ repeat(params.a.su, 1, τ_max)
    x_branch = Vector{Matrix{JuMP.AffExpr}}(undef,n)
    u_branch = Vector{Matrix{JuMP.AffExpr}}(undef,n)
    for j = 1:n
        τ = τ_lu(j)
        x_branch[j] = params.a.Sx*x_branch_us[j] .+ repeat(params.a.sx, 1, N-τ)
        u_branch[j] = params.a.Su*u_branch_us[j] .+ repeat(params.a.su, 1, N-τ)
    end

    # SCP-specific
    ν_ctrl_trunk = @variable(mdl, [1:nx,1:τ_max-1])
    ν_ctrl_branch = Vector{Matrix{JuMP.VariableRef}}(undef,n) 
    ν_ctrl_stitch = @variable(mdl, [1:nx,1:n])
    if trust_region_type != "QP"
        η_trunk = @variable(mdl, [1:τ_max])
        η_branch = Vector{Vector{JuMP.VariableRef}}(undef,n)
    end
    for j = 1:n
        τ = τ_lu(j)
        ν_ctrl_branch[j] = @variable(mdl, [1:nx,1:N-τ-1])
        if trust_region_type != "QP"
            η_branch[j] = @variable(mdl, [1:N-τ])
        end
    end
    @variable(mdl, μ_ctrl) # virtual control objective slack
    @variable(mdl, μ_buff) # virtual buffer objective slack
    @variable(mdl, η_s) # trust region objective slack

    # Convenience functions
    X_trunk(k) = x_trunk[:,k]
    U_trunk(k) = u_trunk[:,k]
    X_branch(k,j) = x_branch[j][:,k]
    U_branch(k,j) = u_branch[j][:,k]

    # ..:: Trunk Constraints ::..
    # Dynamics
    if params.a.ctcs_enabled
        dyn_lin = (t,x,u,p) -> dynamics_linearized_ctcs(t,x,u,params,0)
        dyn_nl  = (t,x,u,p) -> dynamics_nonlinear_ctcs(t,x,u,params,0)
    end

    # Build the trunk
    # Take last deferred target trajectory as the reference
    ref_traj_trunk = copy(ref_trajs.targs[params.a.λ_targs[end]])
    ref_traj_trunk.t = ref_traj_trunk.t[1:τ_max]
    ref_traj_trunk.x = ref_traj_trunk.x[:,1:τ_max]
    ref_traj_trunk.u = ref_traj_trunk.u[:,1:τ_max]
    ref_traj_trunk.x, ref_traj_trunk.u = remove_ref_zeros(ref_traj_trunk.x, ref_traj_trunk.u)

    # Path constraints (problem-specific)
    J_running_trunk,_ = prob_cost(mdl,x_trunk,u_trunk,params)
    if !params.a.ctcs_enabled
        ν_buff_trunk = prob_constraints(mdl,x_trunk,u_trunk,params,ref_traj_trunk,0)
    else
        ν_buff_trunk = []
    end

    # Dynamics
    Ak,Bmk,Bpk,_,wk,_ = c2d_nonlinear(ref_traj_trunk.t,ref_traj_trunk.x,ref_traj_trunk.u,dyn_nl,dyn_lin,params.a.disc,num_disc_steps=params.a.N_msi)
    SxInv = inv(params.a.Sx)
    SuInv = inv(params.a.Su)
    if params.a.disc == 0
        @constraint(mdl, [k=1:τ_max-1], SxInv*X_trunk(k+1) .== SxInv*(Ak[:,:,k]*X_trunk(k) + Bmk[:,:,k]*U_trunk(k) + wk[:,k]) + ν_ctrl_trunk[:,k])
    elseif params.a.disc == 1
        @constraint(mdl, [k=1:τ_max-1], SxInv*X_trunk(k+1) .== SxInv*(Ak[:,:,k]*X_trunk(k) + Bmk[:,:,k]*U_trunk(k) + Bpk[:,:,k]*U_trunk(k+1) + wk[:,k]) + ν_ctrl_trunk[:,k])
    end

    # Trunk time definition
    s_trunk = u_trunk[end,:]
    t_trunk = time_dilation_control_to_wall_clock_time(s_trunk, ref_traj_trunk.t, params.a.disc)
    @constraint(mdl, [k=1:τ_max], s_trunk[k] >= 0)
    
    # CTCS violation
    if params.a.ctcs_enabled
        @constraint(mdl, [k=1:τ_max], x_trunk[end,k]/params.a.ϵ_ctcs <= 1)
        @constraint(mdl, [k=1:τ_max], x_trunk[end,k] >= 0)
    end

    # Trust region
    δXt(k) = SxInv*(X_trunk(k) .- ref_traj_trunk.x[:,k])
    δUt(k) = SuInv*(U_trunk(k) .- ref_traj_trunk.u[:,k])
    if trust_region_type == "QP"
        η_s = sum([δXt(k)'*δXt(k) + δUt(k)'*δUt(k) for k=1:τ_max])
    else
        @constraint(mdl, [k=1:τ_max], δXt(k)'*δXt(k) + δUt(k)'*δUt(k) <= η_trunk[k])
    end

    # ..:: Branch/Target Constraints ::..
    J_cost = Vector(undef,n)
    ν_buff_branch = Vector{Vector{JuMP.AffExpr}}(undef,n)
    for j = 1:n
        # Dynamics
        if params.a.ctcs_enabled
            dyn_lin = (t,x,u,p) -> dynamics_linearized_ctcs(t,x,u,params,j)
            dyn_nl  = (t,x,u,p) -> dynamics_nonlinear_ctcs(t,x,u,params,j)
        end

        # Take jth reference and build it with last N elements
        τ = τ_lu(j)
        ref_traj_branch = copy(ref_trajs.targs[j])
        ref_traj_branch.t = ref_traj_branch.t[τ+1:end]
        ref_traj_branch.x = ref_traj_branch.x[:,τ+1:end]
        ref_traj_branch.u = ref_traj_branch.u[:,τ+1:end]
        ref_traj_branch.x, ref_traj_branch.u = remove_ref_zeros(ref_traj_branch.x, ref_traj_branch.u)

        # Path constraints (problem-specific)
        J_running_branch,J_term_branch = prob_cost(mdl,x_branch[j],u_branch[j],params)
        if !params.a.ctcs_enabled
            ν_buff_branch_ = prob_constraints(mdl, x_branch[j], u_branch[j], params, ref_traj_branch, j)
        else
            ν_buff_branch_ = []
        end
        ν_buff_branch[j] = ν_buff_branch_

        # Dynamics (within branch)
        Ak,Bmk,Bpk,_,wk,_ = c2d_nonlinear(ref_traj_branch.t,ref_traj_branch.x,ref_traj_branch.u,dyn_nl,dyn_lin,params.a.disc,num_disc_steps=params.a.N_msi)
        if params.a.disc == 0
            @constraint(mdl, [k=1:N-τ-1], SxInv*X_branch(k+1,j) .== SxInv*(Ak[:,:,k]*X_branch(k,j) + Bmk[:,:,k]*U_branch(k,j) + wk[:,k]) + ν_ctrl_branch[j][:,k])
        elseif params.a.disc == 1
            @constraint(mdl, [k=1:N-τ-1], SxInv*X_branch(k+1,j) .== SxInv*(Ak[:,:,k]*X_branch(k,j) + Bmk[:,:,k]*U_branch(k,j) + Bpk[:,:,k]*U_branch(k+1,j) + wk[:,k]) + ν_ctrl_branch[j][:,k])
        end

        # Dynamics (stitching to trunk)
        ref_traj_stitch = copy(ref_trajs.targs[j])
        ref_traj_stitch.t = ref_traj_stitch.t[τ:τ+1]
        ref_traj_stitch.x = ref_traj_stitch.x[:,τ:τ+1]
        ref_traj_stitch.u = ref_traj_stitch.u[:,τ:τ+1]
        ref_traj_stitch.x, ref_traj_stitch.u = remove_ref_zeros(ref_traj_stitch.x, ref_traj_stitch.u)
        Ak,Bmk,Bpk,_,wk,_ = c2d_nonlinear(ref_traj_stitch.t,ref_traj_stitch.x,ref_traj_stitch.u,dyn_nl,dyn_lin,params.a.disc,num_disc_steps=params.a.N_msi)
        if params.a.disc == 0
            @constraint(mdl, SxInv*X_branch(1,j) .== SxInv*(Ak[:,:,1]*X_trunk(τ) + Bmk[:,:,1]*U_trunk(τ) + wk[:,1]) + ν_ctrl_stitch[:,j])
        elseif params.a.disc == 1
            @constraint(mdl, SxInv*X_branch(1,j) .== SxInv*(Ak[:,:,1]*X_trunk(τ) + Bmk[:,:,1]*U_trunk(τ) + Bpk[:,:,1]*U_branch(1,j) + wk[:,1]) + ν_ctrl_stitch[:,j])
        end

        # Suboptimality constraint
        J_cost[j] = sum(J_running_trunk) + sum(J_running_branch) + J_term_branch
        @constraint(mdl, J_cost[j] / ((1 + params.a.ϵ_targs[j]) * ref_costs[j]) <= 1)

        # Time dilation constraints (from IC to TC for each target)
        s_branch = u_branch[j][end,:]
        t_target = time_dilation_control_to_wall_clock_time([s_trunk[1:τ];s_branch], ref_trajs.targs[j].t, params.a.disc)
        Δt_target = diff(t_target)
        @constraint(mdl, [k=1:N-1], params.a.Δt_min/params.a.Δt_max <= Δt_target[k]/params.a.Δt_max <= 1)
        @constraint(mdl, params.a.ToF_min/params.a.ToF_max <= t_target[end]/params.a.ToF_max <= 1)
        @constraint(mdl, [k=1:N-τ], s_branch[k] >= 0)

        # CTCS violation
        if params.a.ctcs_enabled
            @constraint(mdl, [k=1:N-τ], x_branch[j][end,k]/params.a.ϵ_ctcs <= 1)
            @constraint(mdl, [k=1:N-τ], x_branch[j][end,k] >= 0)
        end

        # Trust region
        δXb(k) = SxInv*(X_branch(k,j) .- ref_traj_branch.x[:,k])
        δUb(k) = SuInv*(U_branch(k,j) .- ref_traj_branch.u[:,k])
        if trust_region_type == "QP"
            η_s += sum([δXb(k)'*δXb(k) + δUb(k)'*δUb(k) for k=1:N-τ])
        else
            @constraint(mdl, [k=1:N-τ], δXb(k)'*δXb(k) + δUb(k)'*δUb(k) <= η_branch[j][k])
        end
    end

    # ..:: Boundary Conditions ::..
    # >> State <<
    for k = 1:nx
        if ~isinf(params.a.z0[k])
            @constraint(mdl, SxInv[k,k]*x_trunk[k,1] == SxInv[k,k]*params.a.z0[k])
        end
    end
    for j = 1:n
        for k = 1:nx
            if ~isinf(params.a.zf_targs[k,j])
                @constraint(mdl, SxInv[k,k]*x_branch[j][k,end] == SxInv[k,k]*params.a.zf_targs[k,j])
            end
        end
    end

    # >> Input <<
    for k = 1:nu
        if ~isinf(params.a.u0[k])
            @constraint(mdl, SuInv[k,k]*u_trunk[k,1] == SuInv[k,k]*params.a.u0[k])
        end
    end
    for j = 1:n
        for k = 1:nu
            if ~isinf(params.a.uf_targs[k,j])
                @constraint(mdl, SuInv[k,k]*u_branch[j][k,end] == SuInv[k,k]*params.a.uf_targs[k,j])
            end
        end
    end

    # Trust region constraints
    if trust_region_type != "QP"
        @constraint(mdl, vcat(η_s, [vec(η_trunk); vec.(η_branch)...]) in SecondOrderCone())
        @constraint(mdl, η_s >= 0)
    end
    
    # Virtualization constraints
    ν_ctrl = [vec(ν_ctrl_trunk); vec.(ν_ctrl_branch)...; vec(ν_ctrl_stitch)]
    ν_buff = [vec(ν_buff_trunk); vec.(ν_buff_branch)...]
    @constraint(mdl, vcat(μ_ctrl, ν_ctrl) in MOI.NormOneCone(length(ν_ctrl)+1))
    if length(ν_buff) > 0
        @constraint(mdl, vcat(μ_buff, vec(ν_buff)) in MOI.NormOneCone(length(vec(ν_buff))+1))
        @constraint(mdl, μ_buff >= 0)
    else
        @constraint(mdl, μ_buff == 0)
    end
    @constraint(mdl, μ_ctrl >= 0)

    # ..:: Construct cost function and solve ::..
    α = params.a.α_targs
    λ = params.a.λ_targs
    τ_grid = params.a.τ_targs
    α[λ[n-1]] = (α[λ[n-1]] + α[λ[n]])/2
    J_opt = -sum([params.a.α_targs[j]*t_trunk[τ_lu(j)] for j=1:n])/max(params.a.α_targs...)
    # J_opt = - α[λ[1]]*t_trunk[τ_grid[1]] - sum([α[λ[j]]*(t_trunk[τ_grid[j]]-t_trunk[τ_grid[j-1]]) for j=2:n-1])
    obj_scale = 1/sqrt(max(params.a.w_obj_ddto, params.a.w_trust, params.a.w_buff, params.a.w_ctrl))
    @objective(mdl, Min, 
        (params.a.w_obj_ddto * J_opt
      + params.a.w_trust * η_s 
      + params.a.w_buff * μ_buff 
      + params.a.w_ctrl * μ_ctrl)*obj_scale)

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

    if (feas_status == MOI.OPTIMAL || feas_status == MOI.ALMOST_OPTIMAL) && (μ_ctrl_pen <= params.a.ϵ_ctrl) && (μ_buff_pen <= params.a.ϵ_buff) && (η_pen <= params.a.ϵ_trust)
        scp_sub_cvged = true
    else
        scp_sub_cvged = false
    end
    VERB_DDTO && @printf("   SCP Iter: %2.i | Status: %s | Cost = % .2e | μ_ctrl_pen = %s | μ_buff_pen = %s | η_pen = %s\n",
        scp_iter, 
        convert_to_colored_string(solve_status,("Feasible",)),
        value.(J_opt),
        convert_to_colored_string(μ_ctrl_pen,params.a.ϵ_ctrl),
        convert_to_colored_string(μ_buff_pen,params.a.ϵ_buff), 
        convert_to_colored_string(η_pen,params.a.ϵ_trust))
    flush(stdout)

    # ..:: Package the DDTO Solution ::..
    ddto_solution = EmptyDDTOSolution(n)
    for j = 1:n
        ddto_solution.targs[j].t = ref_trajs.targs[j].t # maintain reference dilated time
        ddto_solution.targs[j].x = x[j]
        ddto_solution.targs[j].u = u[j]
        ddto_solution.targs[j].cost = value.(J_cost[j])
    end
    deferrability_times = [value.(t_trunk)[τ_lu(j)] for j=1:n]

    return (ddto_solution, feas_status, scp_sub_cvged, deferrability_times)

end

function set_deferrability_node_allocation!(params)
    # Set deferrability node allocation based on uniform distribution up to N/sqrt(2)
    if params.a.n_targs > 1
        params.a.τ_targs = round.(CVector(range(2,Int(round(params.a.N/sqrt(2))),params.a.n_targs+2)))[2:end-1]
        params.a.τ_targs[end] = params.a.τ_targs[end-1]
    else
        params.a.τ_targs = [params.a.N]
    end
end