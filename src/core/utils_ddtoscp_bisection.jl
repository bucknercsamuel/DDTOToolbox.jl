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
        ref_trajs[j] = generate_initial_guess(params,j)
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

            if feas_status != MOI.OPTIMAL
                @printf("   > SCP subproblem is infeasible, exiting subproblem iteration.\n")
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
        (ddto_branch_sols[k], ref_trajs, status_feas) = solve_bisection_ddtoscp(params_, ref_costs, ref_trajs, cost_dd_sum, previous_ddto_solution)
        if ~status_feas
            # If bisection search failed, create and return an empty DDTO solution
            ddto_branch_sols = Vector{DDTOSolution}(undef, params.n_targs)
            for k = 1:(params.n_targs)
                ddto_branch_sols[k] = EmptyDDTOSolution(params.n_targs-k+1)
                for j=1:(params.n_targs-k+1)
                    N_targ = params.N_targs[j]
                    N_targ_ctrl = N_targ - 1
                    ddto_branch_sols[k].targ_sols[j].t    = zeros(N_targ)
                    ddto_branch_sols[k].targ_sols[j].x    = zeros(params.n,N_targ)
                    ddto_branch_sols[k].targ_sols[j].u    = zeros(params.m,N_targ_ctrl)
                    ddto_branch_sols[k].targ_sols[j].cost = 0
                end
            end
            return (false,ddto_branch_sols)
        end

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

function solve_bisection_ddtoscp(params, ref_costs::CVector, ref_trajs::Vector{Solution}, cost_dd::CReal, previous_ddto_solution::DDTOSolution)::Tuple{DDTOSolution, Vector{Solution}, Bool}
    # Uses bisection search to solve quasiconvex optimization problem 
    # to branch to the next-queued target for rejection.
    #
    # :in params: The params object
    # :in ref_costs: Optimal costs
    # :in cost_dd: Running cost for decision deferral
    # :out ddto_solution: Contains the DDTO solution for this target/branch point

    # Initial search bracket
    τ_min = 0
    if params.free_final_time  
        τ_max = min(params.N_fft - 2, params.τ_max)
    else
        τ_max = min(min(params.N_targs...) - 2, params.τ_max)
    end

    # Bisection search to solve quasiconvex (QCvx) optimization problem
    VERB_DDTO && println("=== Bisection Search for QCvx Optimization ===")
    iter = 0
    while (τ_max - τ_min) > 1
        # Update τ
        τ = Int(ceil(0.5*(τ_max + τ_min)))

        # Compute SCP solution
        (~, ref_traj_update, status_feas) = solve_ddtoscp_subproblem(params, τ, ref_costs, cost_dd, ref_trajs)

        # Update τ_min or τ_max based on solution convergence
        if status_feas
            τ_min = τ
            solve_status = "Feasible"
            ref_trajs = ref_traj_update # Update ref traj
        else
            τ_max = τ
            solve_status = "Not Feasible"
        end

        # Update iteration count
        iter += 1
        VERB_DDTO && @printf("Bisection Iteration #%i -- τ_min: %i, τ_max: %i, status: %s\n", iter, τ_min, τ_max, solve_status)
    end
    if iter == 0
        ref_traj_update = deepcopy(ref_trajs)
    end

    # Set optimal τ
    τ_opt = τ_min
    VERB_DDTO && println("Bisection search terminated -- reached convergence condition (τ_max - τ_min) = 1")

    # Compute converged DDTO solution SCP iteration
    # (just re-use previous solution if cannot be deferred)
    if τ_opt == 0
        ddto_solution = deepcopy(previous_ddto_solution)
        status_feas = true
    else
        (ddto_solution, ref_traj_update, status_feas) = solve_ddtoscp_subproblem(params, τ_opt, ref_costs, cost_dd, ref_trajs)
        ddto_solution.idx_dd = τ_opt
    end

    if status_feas
        ref_trajs = deepcopy(ref_traj_update)
        @printf("Bisection search successful -- τ_opt: %i\n", τ_opt)

        VERB_DDTO && println("Updated costs to each remaining target from initial condition:")
        for j = 1:params.n_targs
            VERB_DDTO && @printf("   Target: %i, Cost: %.3f\n", params.T_targs[j], ddto_solution.costs_sol[j] + cost_dd)
        end
    else
        # error("Bisection search unsuccessful. Problem is unsolved.")
        print("Bisection search unsuccessful. Problem is unsolved.")
        ddto_solution = EmptyDDTOSolution(params.n_targs)
    end

    return (ddto_solution, ref_trajs, status_feas)

end

function solve_ddtoscp_subproblem(params, τ::Int, ref_costs::CVector, cost_dd::CReal, ref_trajs::Vector{Solution})::Tuple{DDTOSolution, Vector{Solution}, Bool}

    # SCP subproblem iteration
    feas_status = undef
    ddto_subsolution = undef
    scp_converged = false
    for k = 1:params.scp_iters

        # Solve SCP subproblem
        (ddto_subsolution, feas_status, scp_converged) = solve_feasible_ddtoscp(params, τ, ref_costs, cost_dd, ref_trajs, k)

        if feas_status != MOI.OPTIMAL
            @printf("   > SCP subproblem is infeasible, exiting subproblem iteration.\n")
            break
        else
            # Use solution results for new reference trajectory
            for j = 1:params.n_targs
                ref_trajs[j] = deepcopy(ddto_subsolution.targ_sols[j])
            end
        end
        if scp_converged
            @printf("   > Convergence condition has been reached, exiting subproblem iteration.\n")
            break
        end
    end

    return (ddto_subsolution, ref_trajs, scp_converged)
end