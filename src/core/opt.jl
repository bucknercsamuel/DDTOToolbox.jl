# ..:: Top-level Solve Function ::..

function solve(params)
    # ..:: Execute solver sequence ::..
    @time begin
        @time begin
            # ..:: Solve for independently-optimal solutions to each target ::..
            scp_solutions = solve_tree_decoupled(params)
            scp_costs = CVector(zeros(params.n_targs))
            for k = 1:params.n_targs
                scp_costs[k] = scp_solutions[k].cost
            end
            println("\n Solve time for generating optimal solutions to each target:")
        end

        @time begin
            # ..:: Solve for DDTO branching solutions to ALL targets ::..
            (feas_ddtoscp, ddtoscp_solutions) = solve_tree_ddto(deepcopy(params), scp_costs)
            println("\n Solve time for generating DDTO branch solutions to all targets:")
        end
        println("\n Solve time for the full DDTO solution stack:")
    end

    # Convert DDTO solutions to branch solutions
    ddtoscp_bsolutions, defer_bsolutions = extract_target_trajectories(params, ddtoscp_solutions; SCP=true)

    # Port decoupled SCP solutions to `BranchSolution` objects for type conformance
    scp_bsolutions_ = Vector{BranchSolution}(undef,params.n_targs)
    for j=1:params.n_targs
        scp_bsolutions_[j] = BranchSolution(scp_solutions[j],-1,-1)
    end
    scp_bsolutions = scp_bsolutions_

    # ..:: Simulate each target solution from I.C. to T.C.
    @time begin
        dynamics = (t,x,sol) -> dynamics_nonlinear(t,x,optimal_controller(t,sol.t,sol.u,params.disc),params)
        scp_bsimulations = simulate_branches(scp_bsolutions, dynamics, params.disc)
        ddtoscp_bsimulations = simulate_branches(ddtoscp_bsolutions, dynamics, params.disc)
        defer_bsimulations = simulate_branches(defer_bsolutions, dynamics, params.disc)
        println("\n Solve time for RK4 simulation:")
    end

    # ..:: Post-processing ::..
    @time begin
        scp_solutions_proc       = process_solutions(scp_bsolutions, params)
        scp_simulations_proc     = process_solutions(scp_bsimulations, params)
        ddtoscp_solutions_proc   = process_solutions(ddtoscp_bsolutions, params)
        ddtoscp_simulations_proc = process_solutions(ddtoscp_bsimulations, params)
        defer_solutions_proc     = process_solutions(defer_bsolutions, params)
        defer_simulations_proc   = process_solutions(defer_bsimulations, params)
        println("\n Solve time for post-processing:")
    end

    # Deferrable segment solution/simulation should just be a scalar ProcessedSolution
    defer_bsolutions = defer_bsolutions[1].sol
    defer_bsimulations = defer_bsimulations[1].sol
    defer_solutions_proc = defer_solutions_proc[1].sol
    defer_simulations_proc = defer_simulations_proc[1].sol

    return (
        scp_solutions_proc, 
        scp_simulations_proc, 
        ddtoscp_solutions_proc, 
        ddtoscp_simulations_proc, 
        defer_solutions_proc, 
        defer_simulations_proc)
end

# ..:: Single-Target (Decoupled) Solver Functions ::..

function solve_tree_decoupled(params)::Vector{Solution}
    # Solve the OPC for a given set of params and all targets independently
    #
    # :in params: The params object
    # :out solutions: Vectorized container for all single-target solutions

    # Define container for each `solve_optimal_target` solution
    solutions = Vector{Solution}(undef, params.n_targs)

    # ..:: Define initial guess reference trajectories using linear interpolations ::..
    ref_trajs = Vector{Solution}(undef, params.n_targs)
    for j = 1:params.n_targs
        ref_trajs[j] = generate_initial_guess_scp(params,j)
    end

    # ..:: SCP Iteration ::..
    VERB_OPT && println("\n=== Decoupled SCP solutions for each target ===")
    for j = 1:params.n_targs
        VERB_OPT && @printf("Target: %i\n", params.T_targs[j])
        feas_status = undef
        solution = undef
        scp_converged = false

        for k = 1:params.scp_iters 
            # Solve SCP subproblem
            (solution, feas_status, scp_converged) = solve_subproblem_decoupled(params, ref_trajs[j], j, k)

            if feas_status != MOI.OPTIMAL && feas_status != MOI.ALMOST_OPTIMAL
                @printf("   > SCP subproblem is infeasible (MOI status: %s), exiting subproblem iteration.\n", feas_status)
                break
            else
                # Use solution results for new reference trajectory
                ref_trajs[j] = solution
            end
            if scp_converged
                @printf("   > Convergence condition has been reached, exiting subproblem iteration.\n")
                break
            end
        end
        solutions[j] = solution
        VERB_OPT && @printf("   > Total cost: %.3f\n\n", solution.cost)
    end

    return solutions
end

function solve_subproblem_decoupled(params, ref_traj::Solution, j_targ::Int, scp_iter::Int)::Tuple{Solution, MOI.TerminationStatusCode, Bool}

    # ..:: Setup ::..
    # Optimizer configuration
    if SOLVER == "ECOS"
        mdl = Model(optimizer_with_attributes(ECOS.Optimizer, "verbose" => 0, "max_iters" => 1000))
    elseif SOLVER == "MOSEK"
        mdl = Model(Mosek.Optimizer)
        JuMP.set_optimizer_attribute(mdl, "LOG",  0) # disable debugging
        JuMP.set_optimizer_attribute(mdl, "MAX_NUM_WARNINGS", 0) # disable warnings
    else
        error("SOLVER is invalid, please select either ECOS or MOSEK")
    end

    # Sizing parameters
    nx = params.nx
    nu = params.nu
    N  = params.N
    if params.disc == 0
        N_ctrl = N-1
    elseif params.disc == 1
        N_ctrl = N
    end

    # Param check(s)
    if params.disc != 0 && params.disc != 1
        error("Please select a valid discretization hold order.")
    end
    
    # ..:: Reference trajectory ::..
    x_ref = ref_traj.x
    u_ref = ref_traj.u

    # ..:: Optimization variables ::..
    # Baseline
    x = @variable(mdl, [1:nx,1:N])
    u = @variable(mdl, [1:nu,1:N_ctrl])

    # SCP-specific
    ν_ctrl = @variable(mdl, [1:nx,1:(N-1)]) # virtual control
    η = @variable(mdl, [1:N]) # trust region
    @variable(mdl, μ_ctrl) # virtual control slack
    @variable(mdl, μ_buff) # virtual buffer slack
    @variable(mdl, η_s) # trust region slack

    # Expressions
    Δt = Array{AffExpr}(undef,N-1) # Wall-clock time step

    # ..:: Make the optimization problem ::..

    # Problem-specific constraints
    J_obj,ν_buff,Sx,Su = core_problem(mdl,x,u,params,ref_traj)
    SxInv = inv(Sx)
    SuInv = inv(Su)
    
    # Dynamics
    X(k) = x[:,k]
    U(k) = u[:,k]
    dyn_lin = (t,x,u) -> dynamics_linearized(t,x,u,params)
    dyn_nl  = (t,x,u) -> dynamics_nonlinear(t,x,u,params)
    Ak,Bmk,Bpk,wk,_ = c2d_nonlinear(ref_traj,dyn_nl,dyn_lin,params.disc)
    if params.disc == 0
        @constraint(mdl, [k=1:N-1], SxInv*X(k+1) .== SxInv*(Ak[:,:,k]*X(k) + Bmk[:,:,k]*U(k) + wk[:,k]) + ν_ctrl[:,k])
    elseif params.disc == 1
        @constraint(mdl, [k=1:N-1], SxInv*X(k+1) .== SxInv*(Ak[:,:,k]*X(k) + Bmk[:,:,k]*U(k) + Bpk[:,:,k]*U(k+1) + wk[:,k]) + ν_ctrl[:,k])
    end

    # Time dilation
    s = u[end,:]
    for k=1:(N-1)
        if params.disc == 0
            Δt[k] = @expression(mdl, params.Δτ[k] * s[k])
        elseif params.disc == 1
            Δt[k] = @expression(mdl, (1/2) * params.Δτ[k] * (s[k] + s[k+1]))
        end
    end
    @constraint(mdl, sum(Δt) <= params.ToF_max)
    @constraint(mdl, [k=1:N_ctrl], params.s_min <= s[k] <= params.s_max)

    # Boundary conditions
    z0 = params.z0
    zf = params.zf_targs[:,j_targ]
    for k = 1:nx # inf = no boundary condition to be applied
        if ~isinf(z0[k])
            @constraint(mdl, x[k,1] == z0[k])
        end
        if ~isinf(zf[k])
            @constraint(mdl, x[k,N] == zf[k])
        end
    end

    # Trust region constraints
    δX(k) = SxInv*(X(k) .- x_ref[:,k])
    δU(k) = SuInv*(U(k) .- u_ref[:,k])
    @constraint(mdl, [k=1:N_ctrl], δX(k)'*δX(k) + δU(k)'*δU(k) <= η[k])
    @constraint(mdl, δX(N)'*δX(N) <= η[N])

    # SCP slack constraints
    @constraint(mdl, vcat(μ_ctrl, vec(ν_ctrl)) in MOI.NormOneCone(length(vec(ν_ctrl))+1))
    @constraint(mdl, vcat(μ_buff, vec(ν_buff)) in MOI.NormOneCone(length(vec(ν_buff))+1))
    @constraint(mdl, vcat(η_s, η) in SecondOrderCone())
    @constraint(mdl, μ_ctrl >= 0)
    @constraint(mdl, μ_buff >= 0)
    @constraint(mdl, η_s >= 0)

    # ..:: Solve the problem and save the solution ::..

    # Cost function
    J_opt  = J_obj
    J_ptr  = η_s
    J_buff = μ_buff
    J_ctrl = μ_ctrl
    @objective(mdl, Min, 
        params.w_obj * J_opt 
      + params.w_trust * J_ptr 
      + params.w_buff * J_buff 
      + params.w_ctrl * J_ctrl)

    optimize!(mdl)
    feas_status = JuMP.termination_status(mdl)
    if feas_status == MOI.OPTIMAL || feas_status == MOI.ALMOST_OPTIMAL
        solve_status = "Feasible"
    else
        solve_status = "Not Feasible"
        return (EmptySolution(), feas_status, false)
    end

    # Package the solution
    x = value.(x)
    u = value.(u)
    cost = value.(J_opt)
    sol = Solution(params.τ,x,u,cost)

    # Obtain evaluation penalties
    μ_buff_pen = value.(μ_buff)
    μ_ctrl_pen = value.(μ_ctrl)
    η_pen = value.(η_s)

    # Determine convergence based on SCP penalties
    if (μ_ctrl_pen <= params.ϵ_ctrl) && (μ_buff_pen <= params.ϵ_buff) && (η_pen <= params.ϵ_trust)
        scp_sub_cvged = true
    else
        scp_sub_cvged = false
    end
    @printf("   SCP Iter: %2.i | Status: %s | μ_ctrl_pen = %.2e | μ_buff_pen = %.2e | η_pen = %.2e\n", scp_iter, solve_status, μ_ctrl_pen, μ_buff_pen, η_pen)
    flush(stdout)

    return (sol, feas_status, scp_sub_cvged)
end

# ..:: DDTO SCP Core Solver Functions ::..

function solve_tree_ddto(params, ref_costs::CVector)::Tuple{Bool,Vector{DDTOSolution}}
    # Top-level DDTO solver for all branch points

    # Define container for each DDTO branch solution
    ddto_branch_sols = Vector{DDTOSolution}(undef, params.n_targs)
    for k = 1:(params.n_targs)
        ddto_branch_sols[k] = EmptyDDTOSolution(params.n_targs-k+1)
    end

    # Define running deferred-decision (DD) trajectory segment cost sum
    cost_dd_sum = 0.

    # Constant node allocation for each successive branch
    τ = Int(floor(params.N / params.n_targs))

    # Obtain initial guess for reference trajectories
    ref_trajs = generate_initial_guess_ddtoscp(τ, params)

    # Define first "previous" ddto solution for first branch using optimal solutions
    previous_ddto_solution = EmptyDDTOSolution(params.n_targs)
    for k=1:(params.n_targs)
        previous_ddto_solution.targ_sols[k] = deepcopy(ref_trajs[k])
        previous_ddto_solution.costs_sol[k] = deepcopy(ref_costs[k])
    end

    # Perform staging in the order of preference
    n_targs_total = deepcopy(params.n_targs)
    params_ = deepcopy(params) # Temp object to be mutated through DDTO loop
    pop_idx = 0
    for k = 1:(n_targs_total-1)

        if VERB_DDTO
            specifiers = repeat("%.3f, ", params_.n_targs)
            specifiers = specifiers[1:end-2] # Remove string and comma at end
            format_string = "   Chosen suboptimality tolerances: {"*specifiers*"}\n"

            @printf("\n========= Solving DDTO for Branch #%i =========\n", k)
            @eval @printf($format_string, $params_.ϵ_targs...)
        end

        # Obtain DDTO solution *if* no deferrability could be made (just uses the previous solution)
        if k > 1
            previous_ddto_solution = deepcopy(ddto_branch_sols[k-1])
            previous_ddto_solution.idx_dd = 0
            previous_ddto_solution.cost_dd = 0
        
            # Update previous_ddto_solution to not have the previously-removed target
            deleteat!(previous_ddto_solution.targ_sols, pop_idx)
            deleteat!(previous_ddto_solution.costs_sol, pop_idx)

            # Truncate previous_ddto_solution for previous solution's deferral
            trunc_start = ddto_branch_sols[k-1].idx_dd+1
            for j = 1:params_.n_targs
                previous_ddto_solution.targ_sols[j].t = previous_ddto_solution.targ_sols[j].t[trunc_start:end] .- previous_ddto_solution.targ_sols[j].t[trunc_start]
                previous_ddto_solution.targ_sols[j].x = previous_ddto_solution.targ_sols[j].x[:,trunc_start:end]
                previous_ddto_solution.targ_sols[j].u = previous_ddto_solution.targ_sols[j].u[:,trunc_start:end]
            end
        end        

        # Obtain DDTO solution for this branch
        (ddto_branch_sols[k], ref_trajs) = solve_scp_iteration_ddto(params_, τ, cost_dd_sum, ref_costs, ref_trajs, previous_ddto_solution)

        # Determine target to be removed (first in the current list of λ_targs)
        λ_targ = params_.λ_targs[1]
        deleteat!(params_.λ_targs, 1)
        pop_idx = findfirst(i->i==λ_targ, params_.T_targs)

        # Have to do some slicing magic for matrices
        matrix_slice = collect(1:params_.n_targs)
        deleteat!(matrix_slice, pop_idx)

        # Update params_ target and IC properties for next branch iteration
        params_.n_targs -= 1
        deleteat!(params_.T_targs, pop_idx)
        deleteat!(params_.ϵ_targs, pop_idx)
        params_.N -= ddto_branch_sols[k].idx_dd
        params_.z0 = ddto_branch_sols[k].targ_sols[1].x[:,ddto_branch_sols[k].idx_dd+1]
        params_.zf_targs = params_.zf_targs[:,matrix_slice]

        # Truncate reference trajectories to deferral point for next branch iteration
        for j = 1:params_.n_targs
            trunc_start = ddto_branch_sols[k].idx_dd+1
            ref_trajs[j].t = ref_trajs[j].t[trunc_start:end]
            ref_trajs[j].x = ref_trajs[j].x[:,trunc_start:end]
            ref_trajs[j].u = ref_trajs[j].u[:,trunc_start:end]
        end

        # Update deferred-decision (DD) cost for next branch iteration
        cost_dd_sum += ddto_branch_sols[k].cost_dd

        # Parameter update print statements
        if VERB_DDTO && (k < n_targs_total-1)
            @printf("   Removed target %i for next branch iteration\n", λ_targ)
        end
    end

    if params.n_targs > 1
        # Add a final element to the branch solutions for the final target
        if params.λ_targs[end-1] > params.λ_targs[end]
            final_idx = 1
        else
            final_idx = 2
        end
        ddto_branch_sols[end].targ_sols = Vector{Solution}(undef, params.n_targs)
        ddto_branch_sols[end].costs_sol = [ddto_branch_sols[end-1].costs_sol[final_idx]]
        ddto_branch_sols[end].idx_dd    = 0
        ddto_branch_sols[end].cost_dd   = 0

        # Remove deferred states/controls from previous solution final target
        for j = 1:params.n_targs
            idx_dd = ddto_branch_sols[end-1].idx_dd
            ddto_branch_sols[end].targ_sols[j]   = EmptySolution()
            ddto_branch_sols[end].targ_sols[j].t = ddto_branch_sols[end-1].targ_sols[final_idx].t[idx_dd+1:end] .- ddto_branch_sols[end-1].targ_sols[final_idx].t[idx_dd+1]
            ddto_branch_sols[end].targ_sols[j].x = ddto_branch_sols[end-1].targ_sols[final_idx].x[:,idx_dd+1:end]
            ddto_branch_sols[end].targ_sols[j].u = ddto_branch_sols[end-1].targ_sols[final_idx].u[:,idx_dd+1:end]
        end
    end

    return (true,ddto_branch_sols)
end

function solve_scp_iteration_ddto(params, τ::Int, cost_dd::CReal, ref_costs::CVector, ref_trajs::Vector{Solution}, previous_solution::DDTOSolution)::Tuple{DDTOSolution, Vector{Solution}}

    # SCP subproblem iteration
    feas_status = undef
    solution = undef
    scp_converged = false
    iteration_cap_reached = true
    for k = 1:params.scp_iters

        # Solve SCP subproblem
        (solution, feas_status, scp_converged) = solve_subproblem_ddto(params, τ, ref_costs, cost_dd, ref_trajs, k)

        if feas_status != MOI.OPTIMAL && feas_status != MOI.ALMOST_OPTIMAL
            iteration_cap_reached = false
            @printf("   ! SCP subproblem is infeasible (MOI status: %s), exiting subproblem iteration.\n", feas_status)
            break
        else
            # Use solution results for new reference trajectory
            for j = 1:params.n_targs
                ref_trajs[j] = deepcopy(solution.targ_sols[j])
            end
        end
        if scp_converged
            iteration_cap_reached = false
            @printf("   > Convergence condition has been reached, exiting subproblem iteration.\n")
            break
        end
    end

    if iteration_cap_reached
        println("   ! SCP subproblem iteration cap reached, exiting subproblem iteration.")
    end

    if scp_converged
        solution.idx_dd = τ
        @printf("   > Time deferred: %.3f seconds\n", solution.targ_sols[1].t[τ])
        return (solution, ref_trajs)
    else
        solution.idx_dd = 0
        @printf("   > Time deferred: 0 seconds\n")
        return (previous_solution, ref_trajs)
    end
end

function solve_subproblem_ddto(params, τ::Int, ref_costs::CVector, cost_dd::CReal, reference_targ_trajs::Vector{Solution}, scp_iter::Int)::Tuple{DDTOSolution, MOI.TerminationStatusCode, Bool}
    # Solve the baseline feasibility problem for DDTO.
    #
    # :in params: The params object
    # :in τ: Branch point index
    # :in ref_costs: Optimal costs from `solve_optimal_pdg_all_targets`
    # :in cost_dd: Running cost for decision deferral
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
        N_ctrl = N-1
    elseif params.disc == 1
        N_ctrl = N
    end

    # Param check(s)
    if params.disc != 0 && params.disc != 1
        error("Please select a valid discretization hold order.")
    end

    # Dynamics functions
    dyn_lin = (t,x,u) -> dynamics_linearized(t,x,u,params)
    dyn_nl  = (t,x,u) -> dynamics_nonlinear(t,x,u,params)

    # ..:: Optimization variables ::..
    # Baseline
    x_trunk  = @variable(mdl, [1:nx,1:τ])
    u_trunk  = @variable(mdl, [1:nu,1:τ])
    x_branch = @variable(mdl, [1:nx,1:N-τ,1:n])
    u_branch = @variable(mdl, [1:nu,1:N_ctrl-τ,1:n])

    # SCP-specific
    ν_ctrl_trunk = @variable(mdl, [1:nx,1:τ-1])
    ν_ctrl_branch = @variable(mdl, [1:nx,1:N-τ-1,1:n])
    ν_ctrl_stitch = @variable(mdl, [1:nx,1:n])
    η_trunk = @variable(mdl, [1:τ])
    η_branch = @variable(mdl, [1:N-τ,1:n])
    @variable(mdl, μ_ctrl) # virtual control slack
    @variable(mdl, μ_buff) # virtual buffer slack
    @variable(mdl, η_s) # trust region slack

    # Convenience functions
    X_trunk(k) = x_trunk[:,k]
    U_trunk(k) = u_trunk[:,k]
    X_branch(k,j) = x_branch[:,k,j]
    U_branch(k,j) = u_branch[:,k,j]

    # ..:: Trunk Constraints ::..
    # Build the trunk
    # Take <any> reference and build it with first τ elements
    ref_traj_trunk = copy(reference_targ_trajs[1])
    ref_traj_trunk.x = ref_traj_trunk.x[:,1:τ]
    ref_traj_trunk.u = ref_traj_trunk.u[:,1:τ]

    # Core constraints
    J_obj_trunk,ν_buff_trunk,Sx,Su = core_problem(mdl,x_trunk,u_trunk,params,ref_traj_trunk)
    SxInv = inv(Sx)
    SuInv = inv(Su)

    # Dynamics
    Ak,Bmk,Bpk,wk,_ = c2d_nonlinear(ref_traj_trunk,dyn_nl,dyn_lin,params.disc)
    if params.disc == 0
        @constraint(mdl, [k=1:τ-1], SxInv*X_trunk(k+1) .== SxInv*(Ak[:,:,k]*X_trunk(k) + Bmk[:,:,k]*U_trunk(k) + wk[:,k]) + ν_ctrl_trunk[:,k])
    elseif params.disc == 1
        @constraint(mdl, [k=1:τ-1], SxInv*X_trunk(k+1) .== SxInv*(Ak[:,:,k]*X_trunk(k) + Bmk[:,:,k]*U_trunk(k) + Bpk[:,:,k]*U_trunk(k+1) + wk[:,k]) + ν_ctrl_trunk[:,k])
    end

    # Time dilation
    s_trunk = u_trunk[end,:]
    Δt_trunk = Array{AffExpr}(undef,τ-1)
    for k=1:τ-1
        if params.disc == 0
            Δt_trunk[k] = @expression(mdl, params.Δτ[k] * s_trunk[k])
        elseif params.disc == 1
            Δt_trunk[k] = @expression(mdl, (1/2) * params.Δτ[k] * (s_trunk[k] + s_trunk[k+1]))
        end
    end
    @constraint(mdl, sum(Δt_trunk) <= params.ToF_max)
    @constraint(mdl, [k=1:τ], params.s_min <= s_trunk[k] <= params.s_max)

    # Trust region
    δXt(k) = SxInv*(X_trunk(k) .- ref_traj_trunk.x[:,k])
    δUt(k) = SuInv*(U_trunk(k) .- ref_traj_trunk.u[:,k])
    @constraint(mdl, [k=1:τ], δXt(k)'*δXt(k) + δUt(k)'*δUt(k) <= η_trunk[k])

    # ..:: Branch Constraints ::..
    J_obj_branch = Array{JuMP.VariableRef}(undef,n)
    ν_buff_branch = Array{Vector{JuMP.VariableRef}}(undef,n)
    for j = 1:n
        # Take jth reference and build it with last n-τ elements
        ref_traj_branch_ = copy(reference_targ_trajs[j])
        ref_traj_branch_.x = ref_traj_branch_.x[:,τ+1:end]
        ref_traj_branch_.u = ref_traj_branch_.u[:,τ+1:end]

        # Core constraints
        J_obj_branch_,ν_buff_branch_,_,_ = core_problem(mdl, x_branch[:,:,j], u_branch[:,:,j], params, ref_traj_branch_)
        J_obj_branch[j] = J_obj_branch_
        ν_buff_branch[j] = ν_buff_branch_

        # Dynamics
        Ak,Bmk,Bpk,wk,_ = c2d_nonlinear(ref_traj_branch_,dyn_nl,dyn_lin,params.disc)
        if params.disc == 0
            @constraint(mdl, [k=1:N-τ-1], SxInv*X_branch(k+1,j) .== SxInv*(Ak[:,:,k]*X_branch(k,j) + Bmk[:,:,k]*U_branch(k,j) + wk[:,k]) + ν_ctrl_branch[:,k,j])
        elseif params.disc == 1
            @constraint(mdl, [k=1:N-τ-1], SxInv*X_branch(k+1,j) .== SxInv*(Ak[:,:,k]*X_branch(k,j) + Bmk[:,:,k]*U_branch(k,j) + Bpk[:,:,k]*U_branch(k+1,j) + wk[:,k]) + ν_ctrl_branch[:,k,j])
        end

        # Time dilation
        s_branch = u_branch[end,:,j]
        Δt_branch = Array{AffExpr}(undef,N-τ-1)
        for k=1:N-τ-1
            if params.disc == 0
                Δt_branch[k] = @expression(mdl, params.Δτ[k] * s_branch[k])
            elseif params.disc == 1
                Δt_branch[k] = @expression(mdl, (1/2) * params.Δτ[k] * (s_branch[k] + s_branch[k+1]))
            end
        end
        @constraint(mdl, sum(Δt_branch) <= params.ToF_max)
        @constraint(mdl, [k=1:N_ctrl-τ], params.s_min <= s_branch[k] <= params.s_max)

        # Trust region
        δXb(k) = SxInv*(X_branch(k,j) .- ref_traj_branch_.x[:,k])
        δUb(k) = SuInv*(U_branch(k,j) .- ref_traj_branch_.u[:,k])
        @constraint(mdl, [k=1:N_ctrl-τ], δXb(k)'*δXb(k) + δUb(k)'*δUb(k) <= η_branch[k,j])
    end

    # ..:: DDTOSCP Constraints ::..
    for j = 1:n
        # Apply suboptimality constraint
        @constraint(mdl, sum(J_obj_trunk) + sum(J_obj_branch[j]) + cost_dd <= (1 + params.ϵ_targs[j]) * ref_costs[j])

        # Apply dynamics stitching
        ref_traj_stitch = copy(reference_targ_trajs[j])
        ref_traj_stitch.x = ref_traj_stitch.x[:,τ:τ+1]
        ref_traj_stitch.u = ref_traj_stitch.u[:,τ:τ+1]
        Ak,Bmk,Bpk,wk,_ = c2d_nonlinear(ref_traj_stitch,dyn_nl,dyn_lin,params.disc)
        if params.disc == 0
            @constraint(mdl, SxInv*X_branch(1,j) .== SxInv*(Ak[:,:,1]*X_trunk(τ) + Bmk[:,:,1]*U_trunk(τ) + wk[:,1]) + ν_ctrl_stitch[:,j])
        elseif params.disc == 1
            @constraint(mdl, SxInv*X_branch(1,j) .== SxInv*(Ak[:,:,1]*X_trunk(τ) + Bmk[:,:,1]*U_trunk(τ) + Bpk[:,:,1]*U_branch(1,j) + wk[:,1]) + ν_ctrl_stitch[:,j])
        end
    end

    # Maintain continuity for FOH discretization if we have already deferred by some amount (usually after first branch iteration)
    # (acts as a initial condition on control constrained to the reference initial condition)
    if params.disc == 1 && cost_dd > 0
        u_ref = reference_targ_trajs[1].u # any traj works
        @constraint(mdl, U_trunk(1) .== u_ref[:,1])
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
                @constraint(mdl, x_branch[k,end,j] == params.zf_targs[k,j])
            end
        end
    end

    # ..:: Slack Constraints ::..
    ν_ctrl = [vec(ν_ctrl_trunk); vec(ν_ctrl_branch); vec(ν_ctrl_stitch)]
    ν_buff = [vec(ν_buff_trunk); vec.(ν_buff_branch)...]
    @constraint(mdl, vcat(μ_ctrl, ν_ctrl) in MOI.NormOneCone(length(ν_ctrl)+1))
    @constraint(mdl, vcat(μ_buff, ν_buff) in MOI.NormOneCone(length(ν_buff)+1))
    @constraint(mdl, vcat(η_s, [vec(η_trunk); vec(η_branch)]) in SecondOrderCone())
    @constraint(mdl, μ_buff >= 0)
    @constraint(mdl, μ_ctrl >= 0)
    @constraint(mdl, η_s >= 0)

    # ..:: Construct cost function and solve ::..
    J_opt  = -sum(Δt_trunk)
    J_ptr  = η_s
    J_buff = μ_buff
    J_ctrl = μ_ctrl
    @objective(mdl, Min, 
        params.w_obj * J_opt 
      + params.w_trust * J_ptr 
      + params.w_buff * J_buff 
      + params.w_ctrl * J_ctrl)

    optimize!(mdl)
    feas_status = JuMP.termination_status(mdl)

    # ..:: Extract the solution ::..
    x = zeros(nx, N, n)
    u = zeros(nu, N, n)
    for j = 1:n
        x[:,:,j] = hcat(value.(x_trunk), reshape(value.(x_branch[:,:,j]),nx,N-τ))
        u[:,:,j] = hcat(value.(u_trunk), reshape(value.(u_branch[:,:,j]),nu,N-τ))
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

    if feas_status == MOI.OPTIMAL && (μ_ctrl_pen <= params.ϵ_ctrl) && (μ_buff_pen <= params.ϵ_buff) && (η_pen <= params.ϵ_trust)
        scp_sub_cvged = true
    else
        scp_sub_cvged = false
    end
    @printf("   SCP Iter: %2.i | Status: %s | Cost = %.2e | μ_ctrl_pen = %.2e | μ_buff_pen = %.2e | η_pen = %.2e\n", scp_iter, solve_status, value.(J_opt), μ_ctrl_pen, μ_buff_pen, η_pen)
    flush(stdout)

    # ..:: Determine optimal cost and deferred-decision (DD) cost ::..
    costs_sol = [sum(value.(J_obj_trunk)) + sum(value.(J_obj_branch[j])) for j = 1:n]
    cost_dd = sum(value.(J_obj_trunk))

    # ..:: Package the DDTO Solution ::..
    ddto_solution = EmptyDDTOSolution(n)
    for j = 1:n
        ddto_solution.targ_sols[j].t = params.τ
        ddto_solution.targ_sols[j].x = x[:,:,j]
        ddto_solution.targ_sols[j].u = u[:,:,j]
        ddto_solution.targ_sols[j].cost = costs_sol[j]
    end
    ddto_solution.costs_sol = costs_sol
    ddto_solution.cost_dd   = cost_dd

    return (ddto_solution, feas_status, scp_sub_cvged)

end

# ..:: Numerical Scaling ::..

function scaling_matrices(xmin, xmax)
    make_diagonal(x) = Diagonal(x)
    make_diagonal(x::Number) = x
    s = (xmin + xmax) / 2
    S = make_diagonal(max.(1.0, abs.((xmax - xmin) / 2)))
    return S,s
end

function unscale(xs, xmin, xmax)
    S,s = scaling_matrices(xmin, xmax)
    dims = size(xs)
    if length(dims) == 2
        x = S*xs .+ s
    else
        xs_reshape = reshape(xs, dims[1], prod(dims[2:end]))
        x_reshape = S*xs_reshape .+ s
        x = reshape(x_reshape, dims...)
    end
    return x
end

# ..:: Line Search Optimization ::..

function bisection_search_min_feasible(fun::Function, τ_min::Int, τ_max::Int, ϵ_tol::Int)::Int
    # Use bisection search to find the minimum feasible solution
    # to a function (not to be confused with DDTO bisection search, 
    # which finds the *maximum* feasible solution) 
    #
    # :in fun: Function to be evaluated (must take opt variable τ as input and return cost)
    # :in τ_min: Bracket search minimum bound
    # :in τ_max: Bracket search maximum bound
    # :in ϵ_tol: Suboptimality convergence tolerance

    iter = 1
    while (τ_max - τ_min) > ϵ_tol
        # Update τ
        τ = Int(ceil(0.5*(τ_max + τ_min)))

        # Compute feasible DDTO
        cost = fun(τ)

        # Update τ_max or τ_min based on solution convergence
        if ~isinf(cost)
            τ_max = τ
            solve_status = "Feasible"
        else
            τ_min = τ
            solve_status = "Not Feasible"
        end
        VERB_OPT && @printf("Iteration: %i, τ_min: %i, τ_max: %i -- %s\n", iter, τ_min, τ_max, solve_status)

        # Update iteration count
        iter += 1
    end

    # Set optimal τ
    τ_opt = τ_max

    return τ_opt
end