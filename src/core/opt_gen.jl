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