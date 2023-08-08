#= DDTO SCP Formulation.

Author: Samuel Buckner (UW-ACL)
=#

function solve_decoupled_scp_tree(params)::Vector{Solution}
    # Solve the OPC for a given set of params and all targets independently
    # using `solve_scp_pseudooptimal_target`
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

    # # ..:: Define initial guess reference trajectories using optimal solutions w/out nonconvexities ::..
    # ref_trajs = solve_optimal_tree(params)

    # ..:: SCP Iteration ::..
    VERB_OPT && println("\n=== Decoupled SCP solutions for each target ===")
    for j = 1:params.n_targs
        VERB_OPT && @printf("Target: %i\n", params.T_targs[j])
        feas_status = undef
        solution = undef
        scp_converged = false

        for k = 1:params.scp_iters 
            if params.free_final_time
                N = params.N_fft
            else
                N = params.N_targs[j]
            end

            # Solve SCP subproblem
            (solution, feas_status, scp_converged) = solve_scp_target(params, ref_trajs[j], N, j, k)

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

function solve_ddtoscp_tree(params, ref_costs::CVector, ref_trajs::Vector{Solution})::Tuple{Bool,Vector{DDTOSolution}}
    # Top-level DDTO solver for all branch points

    # Define container for each DDTO branch solution
    ddto_branch_sols = Vector{DDTOSolution}(undef, params.n_targs)
    for k = 1:(params.n_targs)
        ddto_branch_sols[k] = EmptyDDTOSolution(params.n_targs-k+1)
    end

    # Define first "previous" ddto solution for first branch using optimal solutions
    previous_ddto_solution = EmptyDDTOSolution(params.n_targs)
    for k=1:(params.n_targs)
        previous_ddto_solution.targ_sols[k] = deepcopy(ref_trajs[k])
        previous_ddto_solution.costs_sol[k] = deepcopy(ref_costs[k])
    end

    # Define running deferred-decision (DD) trajectory segment cost sum
    cost_dd_sum = 0.

    # Constant node allocation for each successive branch
    τ = Int(floor(params.N_fft / (params.n_targs-1)))

    # Perform branching in the order of preference
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

        # Obtain Bisection-optimal DDTO solution for this branch
        (ddto_branch_sols[k], ref_trajs) = solve_ddtoscp_subproblem(params_, cost_dd_sum, ref_costs, ref_trajs, previous_ddto_solution)

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
        deleteat!(params_.N_targs, pop_idx)
        deleteat!(params_.ϵ_targs, pop_idx)
        params_.N_targs .-= ddto_branch_sols[k].idx_dd
        params_.N_fft -= ddto_branch_sols[k].idx_dd
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
            ddto_branch_sols[end].targ_sols[j].t = ddto_branch_sols[end-1].targ_sols[final_idx].t[idx_dd+1:end]
            ddto_branch_sols[end].targ_sols[j].x = ddto_branch_sols[end-1].targ_sols[final_idx].x[:,idx_dd+1:end]
            ddto_branch_sols[end].targ_sols[j].u = ddto_branch_sols[end-1].targ_sols[final_idx].u[:,idx_dd+1:end]
        end
    end

    return (true,ddto_branch_sols)
end

function solve_ddtoscp_subproblem(params, τ::Int, cost_dd::CReal, ref_costs::CVector, ref_trajs::Vector{Solution}, previous_solution::DDTOSolution)::Tuple{DDTOSolution, Vector{Solution}}

    # SCP subproblem iteration
    feas_status = undef
    solution = undef
    scp_converged = false
    iteration_cap_reached = true
    for k = 1:params.scp_iters

        # Solve SCP subproblem
        (solution, feas_status, scp_converged) = solve_feasible_ddtoscp(params, τ, ref_costs, cost_dd, ref_trajs, k)

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