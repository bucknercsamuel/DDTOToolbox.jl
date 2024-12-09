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

function remove_ref_zeros!(x_ref, u_ref; ϵ_small=1e-6)
    x_ref[x_ref .== 0] .= ϵ_small
    u_ref[u_ref .== 0] .= ϵ_small
end

# ..:: Line Search Optimization ::..

function bisection_search_min_feasible(fun::Function, τ_min::Int, τ_max::Int; ϵ_tol::Int=1, verbose::Bool=true)::Int
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
        verbose && @printf("Iteration: %i, τ_min: %i, τ_max: %i -- %s\n", iter, τ_min, τ_max, solve_status)

        # Update iteration count
        iter += 1
    end

    # Set optimal τ
    τ_opt = τ_max

    return τ_opt
end

function bisection_search_min_feasible(fun::Function, τ_min::Float64, τ_max::Float64; ϵ_tol::Float64=1e-3, verbose::Bool=true)::Float64
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
        τ = 0.5*(τ_max + τ_min)

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
        verbose && @printf("Iteration: %i, τ_min: %.2f, τ_max: %.2f -- %s\n", iter, τ_min, τ_max, solve_status)

        # Update iteration count
        iter += 1
    end

    # Set optimal τ
    τ_opt = τ_max

    return τ_opt
end

function golden_section(f::Function, a::Float64, b::Float64; tol::Float64=1e-3, get_first_feasible::Bool=false, verbose::Bool=true)::Tuple{Float64, Float64, Float64}
    # Golden search for minimizing a unimodal function f(x) on the
    # interval [a,b] to within a prescribed golerance in
    # x. Implementation is based on [1].
    #
    # [1] M. J. Kochenderfer and T. A. Wheeler, Algorithms for
    # Optimization. Cambridge, Massachusetts: The MIT Press, 2019.
    #
    # :in f: oracle with call signature v=f(x) where v is saught to be
    #        minimized.
    # :in a: search domain lower bound.
    # :in b: search domain upper bound.
    # :in tol: tolerance in terms of maximum distance that the
    #          minimizer x∈[a,b] is away from a or b.
    # :in get_first_feasible: Return the first **feasible** solution rather than solving
    #          all the way to the tolerance.
    # :out sol: a tuple where s[1] is the argmin and s[2] is the argmax.
    
    ϕ = (1+√5)/2
    n = ceil(log((b-a)/tol)/log(ϕ)+1)
    ρ = ϕ-1
    d = ρ*b+(1-ρ)*a
    yd = f(d)
    x_sol_last_feas = Inf
    for ~ = 1:n-1
        c = ρ*a+(1-ρ)*b
        yc = f(c)
        if yc < yd
            b,d,yd = d,c,yc
        else
            a,b = b,c
        end
        bracket = sort([a,b,c,d])
        verbose && @printf("Golden Bracket: [%.3f,%.3f,%.3f,%.3f] -- Loss: %.3f\n", bracket..., yc)
        if get_first_feasible && !isinf(yc)
            verbose && println("Feasible solution found, breaking Golden Section Search.")
            break
        end
        if !isinf(yc)
            x_sol_last_feas = b
        end
        flush(stdout)
    end
    x_sol = b
    sol = (x_sol,f(x_sol),x_sol_last_feas)
    return sol
end

# ..:: Other ::..

heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), one(x)) # needed for symbolic maximum differentiation