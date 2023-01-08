#= DDTO for Landing -- Optimal PDG Functions.

Author: Samuel Buckner (UW-ACL)
=#

function solve_optimal_pdg_all_targets(lander::Lander)::Vector{Solution}
    # Solve the optimal landing (PDG) problem for a given lander and all targets
    # using `solve_optimal_pdg_single_target`
    # ** (Not DDTO formulation, but used for comparison) **
    #
    # :in lander: The lander object
    # :out solutions: Vectorized container for all single-target solutions

    # Define container for each `solve_optimal_pdg_single_target` solution
    solutions = Vector{Solution}(undef, lander.n_targs)
    N_optimal = CVector(undef, lander.n_targs)

    # Obtain solutions for each target
    VERB_OPT && println("\n=== Optimal solutions for each target ===")
    for j = 1:lander.n_targs
        solutions[j] = solve_optimal_pdg_single_target(lander, lander.N_targs[j], j)
        VERB_OPT && @printf("Target: %i, Cost: %.3f\n", lander.T_targs[j], solutions[j].cost)
    end

    return solutions
end

function solve_optimal_pdg_single_target(lander::Lander, N::Int, j_targ::Int)::Solution
    # Solve the optimal landing (PDG) problem for a given lander and single target
    # ** (Not DDTO formulation, but used for comparison) **
    #
    # :in lander: The lander object.
    # :in N: Time horizon (used to obtain the optimal time horizon populated into lander.N_targs)
    # :in j_targ: Target index
    # :out sol: Container for solution variables


    # ..:: Discrete time interval ::..
    Δt = lander.Δt
    tf = Δt * (N-1)
    t  = CVector(range(0, stop=tf, length=N))
    N_ctrl = N-1 # Number of nodes to apply control constraints for (N-1 for ZOH)


    # ..:: Make the optimization problem ::..

    # >> Optimizer setup <<
    if SOLVER == "ECOS"
        mdl = Model(optimizer_with_attributes(ECOS.Optimizer, "verbose" => 0))
    elseif SOLVER == "MOSEK"
        mdl = Model(Mosek.Optimizer)
        JuMP.set_optimizer_attribute(mdl, "LOG", 0) # disable debugging
    else
        error("SOLVER is invalid, please select either ECOS or MOSEK")
    end

    # >> Optimization variables <<
    @variable(mdl, r[1:3,1:N])
    @variable(mdl, v[1:3,1:N])
    @variable(mdl, T[1:3,1:N_ctrl])
    @variable(mdl, Γ[1:N_ctrl])

    # >> Convenience functions <<
    X = (k) -> [r[:,k]; v[:,k]] # State at time index k
    U = (k) -> [T[:,k]; Γ[k]]   # Input at time index k

    # >> Cost function <<
    @objective(mdl, Min, Δt*sum(Γ))


    # ..:: Constraints ::..

    # >> Dynamics <<
    A,B,p = c2d_zoh(lander,Δt)
    @constraint(mdl, [k=1:N-1], X(k+1) .==  A*X(k) + B*U(k) + p)

    # >> Constant altitude constraint <<
    @constraint(mdl, [k=1:N-1], r[3,k+1] == r[3,k])

    # >> Thrust bounds <<
    @constraint(mdl, [k=1:N_ctrl], Γ[k] >= lander.ρ_min)
    @constraint(mdl, [k=1:N_ctrl], Γ[k] <= lander.ρ_max)
    @constraint(mdl, [k=1:N_ctrl], vcat(Γ[k], T[:,k]) in MOI.SecondOrderCone(4))

    # >> Attitude pointing constraint <<
    @constraint(mdl, [k=1:N_ctrl], dot(T[:,k],e_z) >= Γ[k]*cos(lander.γ_p))

    # >> Velocity upper bound <<
    # @constraint(mdl, [k=1:N], vcat(lander.v_max_V,v[3,k])   in MOI.SecondOrderCone(2))
    @constraint(mdl, [k=1:N], vcat(lander.v_max_L,v[1:2,k]) in MOI.SecondOrderCone(3))

    # >> Cage bounds <<
    # @constraint(mdl, [k=1:N], r[1,k] >= lander.x_arena_lims[1])
    # @constraint(mdl, [k=1:N], r[1,k] <= lander.x_arena_lims[2])
    # @constraint(mdl, [k=1:N], r[2,k] >= lander.y_arena_lims[1])
    # @constraint(mdl, [k=1:N], r[2,k] <= lander.y_arena_lims[2])
    # @constraint(mdl, [k=1:N], r[3,k] >= lander.z_arena_lims[1])
    # @constraint(mdl, [k=1:N], r[3,k] <= lander.z_arena_lims[2])

    # >> Boundary conditions <<
    @constraint(mdl, r[:,1] .== lander.r0)
    @constraint(mdl, v[:,1] .== lander.v0)
    # @constraint(mdl, T[3,1]  == -dot(e_z, lander.g)*lander.mass) # Vertical orientation constraint
    @constraint(mdl, r[:,N] .== lander.rf_targs[:,j_targ])
    @constraint(mdl, v[:,N] .== lander.vf_targs[:,j_targ])
    # @constraint(mdl, T[3,N_ctrl]  == Γ[N_ctrl]) # Vertical orientation constraint

    # ..:: Solve the problem and save the solution ::..

    optimize!(mdl)
    if termination_status(mdl) != MOI.OPTIMAL
        return EmptySolution()
    end

    # Raw data
    r = value.(r)
    v = value.(v)
    T = value.(T)
    Γ = value.(Γ)
    r0_relax = CVector(undef,0)
    rf_relax = CVector(undef,0)
    cost = objective_value(mdl)

    # Processed data
    T_nrm = CVector([norm(T[:,i],2) for i=1:N_ctrl])
    γ = CVector([acos(dot(T[:,k],e_z)/norm(T[:,k],2)) for k=1:N_ctrl])

    # Package the solution
    sol = Solution(t,r,v,T,Γ,r0_relax,rf_relax,cost,T_nrm,γ)

    return sol
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
