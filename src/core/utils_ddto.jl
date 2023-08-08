#= DDTO for Landing -- DDTO PDG Functions.

Author: Samuel Buckner (UW-ACL)
=#

function solve_optimal_tree(params)::Vector{Solution}
    # Solve the OPC for a given set of params and all targets independently
    # using `solve_optimal_target`
    #
    # :in params: The params object
    # :out solutions: Vectorized container for all single-target solutions

    # Define container for each `solve_optimal_target` solution
    solutions = Vector{Solution}(undef, params.n_targs)

    # Obtain solutions for each target
    VERB_OPT && println("\n=== Optimal solutions for each target ===")
    for j = 1:params.n_targs
        solutions[j] = solve_optimal_target(params, params.N_targs[j], j)
        VERB_OPT && @printf("Target: %i, Cost: %.3f\n", params.T_targs[j], solutions[j].cost)
    end

    return solutions
end

function solve_ddto_tree(params, costs_optimal_0::CVector)::Vector{DDTOSolution}
    # Top-level DDTO solver for all branch points
    #
    # :in params: The params object
    # :in costs_optimal_0: Optimal costs from initial condition
    # :out ddto_branch_sols: Vectorized container for all DDTO branch solutions

    # Define container for each DDTO branch solution
    ddto_branch_sols = Vector{DDTOSolution}(undef, params.n_targs)
    for k = 1:(params.n_targs)
        ddto_branch_sols[k] = EmptyDDTOSolution(params.n_targs-k+1)
    end

    # Define running deferred-decision (DD) trajectory segment cost sum
    cost_dd_sum = 0.

    # Perform branching in the order of preference
    n_targs_total = copy(params.n_targs)
    params_ = copy(params) # Temp object to be mutated through DDTO loop
    for k = 1:(n_targs_total-1)

        if VERB_DDTO
            # specifiers = repeat("%.3f, ", params_.n_targs)
            # specifiers = specifiers[1:end-2] # Remove string and comma at end
            # format_string = "   Chosen suboptimality tolerances: {"*specifiers*"}\n"

            @printf("\n========= Solving DDTO for Branch #%i =========\n", k)
            # @eval @printf($format_string, $params_.ϵ_targs...)
        end

        # Obtain Bisection-optimal DDTO solution for this branch
        ddto_branch_sols[k] = solve_bisection_ddto(params_, costs_optimal_0, cost_dd_sum)

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
        params_.z0 = ddto_branch_sols[k].targ_sols[1].x[:,ddto_branch_sols[k].idx_dd+1]
        params_.zf_targs = params_.zf_targs[:,matrix_slice]

        # Update deferred-decision (DD) cost for next branch iteration
        cost_dd_sum += ddto_branch_sols[k].cost_dd

        # Parameter update print statements
        if VERB_DDTO && (k < n_targs_total-1)
            @printf("   Removed target %i for next branch iteration\n", λ_targ)
        end
    end

    # Add a final element to the branch solutions for the final target
    if params.λ_targs[end-1] > params.λ_targs[end]
        final_idx = 1
    else
        final_idx = 2
    end
    ddto_branch_sols[end].targ_sols = [ddto_branch_sols[end-1].targ_sols[final_idx]]
    ddto_branch_sols[end].costs_sol = [ddto_branch_sols[end-1].costs_sol[final_idx]]
    ddto_branch_sols[end].idx_dd    = 0
    ddto_branch_sols[end].cost_dd   = 0

    return ddto_branch_sols
end

function solve_bisection_ddto(params, costs_optimal::CVector, cost_dd::CReal)::DDTOSolution
    # Uses bisection search to solve quasiconvex optimization problem 
    # to branch to the next-queued target for rejection.
    #
    # :in params: The params object
    # :in costs_optimal: Optimal costs from `solve_optimal_pdg_all_targets`
    # :in cost_dd: Running cost for decision deferral
    # :out ddto_solution: Contains the DDTO solution for this target/branch point

    # Initial search bracket
    τ_min = 0
    τ_max = min(params.N_targs...) - 2

    # Bisection search to solve quasiconvex (QCvx) optimization problem
    VERB_DDTO && println("=== Bisection Search for QCvx Optimization ===")
    iter = 1
    while (τ_max - τ_min) > 1
        # Update τ
        τ = Int(ceil(0.5*(τ_max + τ_min)))

        # Compute feasible DDTO
        (~, status_feas) = solve_feasible_ddto(params, τ, costs_optimal, cost_dd)

        # Update τ_min or τ_max based on solution convergence
        if status_feas == MOI.OPTIMAL
            τ_min = τ
            solve_status = "Feasible"
        else
            τ_max = τ
            solve_status = "Not Feasible"
        end
        VERB_DDTO && @printf("Iteration: %i, τ_min: %i, τ_max: %i -- %s\n", iter, τ_min, τ_max, solve_status)

        # Update iteration count
        iter += 1
    end

    # Set optimal τ
    τ_opt = τ_min
    VERB_DDTO && println("Bisection search terminated -- reached convergence condition (τ_max - τ_min) = 1")

    # Compute converged DDTO solution
    (ddto_solution, status_feas) = solve_feasible_ddto(params, τ_opt, costs_optimal, cost_dd)
    ddto_solution.idx_dd = τ_opt

    # Determine solution convergence
    if status_feas == MOI.OPTIMAL
        VERB_DDTO &&  @printf("Bisection search successful -- τ_opt: %i\n", τ_opt)
    else
        error("Bisection search unsuccessful. Problem is unsolved.")
    end
    VERB_DDTO && println("New costs to each remaining target:")
    for j = 1:params.n_targs
        VERB_DDTO && @printf("   Target: %i, Cost: %.3f\n", params.T_targs[j], ddto_solution.costs_sol[j] + cost_dd)
    end

    # Remove excess state/control nodes from solution
    for j = 1:length(ddto_solution.targ_sols)
        N_targ = params.N_targs[j]
        N_targ_ctrl = N_targ - 1
        ddto_solution.targ_sols[j].t = ddto_solution.targ_sols[j].t[1:N_targ]
        ddto_solution.targ_sols[j].x = ddto_solution.targ_sols[j].x[:,1:N_targ]
        ddto_solution.targ_sols[j].u = ddto_solution.targ_sols[j].u[:,1:N_targ_ctrl]
    end

    return ddto_solution

end

function extract_target_trajectories(params, sols_ddto::Array{DDTOSolution})::Tuple{Vector{BranchSolution},Vector{BranchSolution}}

    # Obtain full solutions to each target
    DDTO_target_solutions = Vector{BranchSolution}(undef, params.n_targs)
    net_deferral_idx = 1
    net_cost_dd = 0
    leading_cost_traj = 0
    T_targs = copy(params.T_targs)
    n_ = size(sols_ddto[1].targ_sols[1].x,1)
    m_ = size(sols_ddto[1].targ_sols[1].u,1)
    t_trunk = CVector(undef,0)
    x_trunk = CMatrix(undef,n_,0)
    u_trunk = CMatrix(undef,m_,0)
    t_offset = 0

    for j in 1:params.n_targs

        # Obtain branch to the desired target
        deferral_idx = sols_ddto[j].idx_dd
        net_deferral_idx += deferral_idx
        λ_targ = params.λ_targs[j]
        rej_idx = findfirst(i->i==λ_targ, T_targs)
        deleteat!(T_targs,rej_idx)
        sol_branch = sols_ddto[j].targ_sols[rej_idx]
        t_branch = deepcopy(sol_branch.t)
        x_branch = deepcopy(sol_branch.x)
        u_branch = deepcopy(sol_branch.u)
        t_offset_  = copy(t_branch[deferral_idx+1]) # - t_branch[1]
        t_branch .+= t_offset
        t_offset  += copy(t_offset_)

        # Compute costs
        if j < params.n_targs
            leading_cost_traj = net_cost_dd
        end
        total_cost = sols_ddto[j].targ_sols[rej_idx].cost + leading_cost_traj
        net_cost_dd += sols_ddto[j].cost_dd

        # Concatenate to create the solution to the given target
        t_target = vcat(t_trunk, t_branch[1:end])
        x_target = hcat(x_trunk, x_branch[:,1:end])
        u_target = hcat(u_trunk, u_branch[:,1:end])
        sol_target = Solution(t_target, x_target, u_target, total_cost)

        # Build the "trunk" to the deferral point for the next solution
        t_trunk = vcat(t_trunk, t_branch[1:deferral_idx])
        x_trunk = hcat(x_trunk, x_branch[:,1:deferral_idx])
        u_trunk = hcat(u_trunk, u_branch[:,1:deferral_idx])

        # Add solution
        def_idx = findfirst(i->i==λ_targ, params.T_targs)
        DDTO_target_solutions[def_idx] = EmptyBranchSolution()
        DDTO_target_solutions[def_idx].sol = sol_target
        DDTO_target_solutions[def_idx].cost_dd = net_cost_dd
        DDTO_target_solutions[def_idx].idx_dd = net_deferral_idx
    end
    DDTO_trunk = Vector{BranchSolution}(undef, 1)
    DDTO_trunk_sol = Solution(t_trunk, x_trunk, u_trunk, net_cost_dd)
    DDTO_trunk[1] = EmptyBranchSolution()
    DDTO_trunk[1].sol = DDTO_trunk_sol
    DDTO_trunk[1].cost_dd = -1
    DDTO_trunk[1].idx_dd = -1

    return (DDTO_target_solutions, DDTO_trunk)
end