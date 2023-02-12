#= DDTO for Landing -- Optimal PDG Functions.

Author: Samuel Buckner (UW-ACL)
=#

function solve_optimal_tree(params::Params)::Vector{BranchSolution}
    # Solve the OPC for a given set of params and all targets independently
    # using `solve_optimal_target`
    # ** (Not DDTO formulation, but used for comparison) **
    #
    # :in params: The params object
    # :out solutions: Vectorized container for all single-target solutions

    # Define container for each `solve_optimal_target` solution
    solutions = Vector{BranchSolution}(undef, params.n_targs)

    # Obtain solutions for each target
    VERB_OPT && println("\n=== Optimal solutions for each target ===")
    for j = 1:params.n_targs
        solution = solve_optimal_target(params, params.N_targs[j], j)
        solutions[j] = BranchSolution(solution,-1,-1)
        VERB_OPT && @printf("Target: %i, Cost: %.3f\n", params.T_targs[j], solution.cost)
    end

    return solutions
end

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
