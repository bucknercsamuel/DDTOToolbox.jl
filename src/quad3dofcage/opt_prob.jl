function solve_optimal_target(params::Params, N::Int, j_targ::Int)::Solution
    # Solve the optimal landing (PDG) problem for a given params and single target
    # ** (Not DDTO formulation, but used for comparison) **
    #
    # :in params: The params object.
    # :in N: Time horizon (used to obtain the optimal time horizon populated into params.N_targs)
    # :in j_targ: Target index
    # :out sol: Container for solution variables


    # ..:: Discrete time interval ::..
    Δt = params.Δt
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
    A,B,p = c2d_zoh(params,Δt)
    @constraint(mdl, [k=1:N-1], X(k+1) .==  A*X(k) + B*U(k) + p)

    # >> Constant altitude constraint <<
    @constraint(mdl, [k=1:N-1], r[3,k+1] == r[3,k])

    # >> Thrust bounds <<
    @constraint(mdl, [k=1:N_ctrl], Γ[k] >= params.ρ_min)
    @constraint(mdl, [k=1:N_ctrl], Γ[k] <= params.ρ_max)
    @constraint(mdl, [k=1:N_ctrl], vcat(Γ[k], T[:,k]) in MOI.SecondOrderCone(4))

    # >> Attitude pointing constraint <<
    @constraint(mdl, [k=1:N_ctrl], dot(T[:,k],e_z) >= Γ[k]*cos(params.γ_p))

    # >> Velocity upper bound <<
    # @constraint(mdl, [k=1:N], vcat(params.v_max_V,v[3,k])   in MOI.SecondOrderCone(2))
    @constraint(mdl, [k=1:N], vcat(params.v_max_L,v[1:2,k]) in MOI.SecondOrderCone(3))

    # >> Cage bounds <<
    # @constraint(mdl, [k=1:N], r[1,k] >= params.x_arena_lims[1])
    # @constraint(mdl, [k=1:N], r[1,k] <= params.x_arena_lims[2])
    # @constraint(mdl, [k=1:N], r[2,k] >= params.y_arena_lims[1])
    # @constraint(mdl, [k=1:N], r[2,k] <= params.y_arena_lims[2])
    # @constraint(mdl, [k=1:N], r[3,k] >= params.z_arena_lims[1])
    # @constraint(mdl, [k=1:N], r[3,k] <= params.z_arena_lims[2])

    # >> Boundary conditions <<
    @constraint(mdl, r[:,1] .== params.r0)
    @constraint(mdl, v[:,1] .== params.v0)
    # @constraint(mdl, T[3,1]  == -dot(e_z, params.g)*params.mass) # Vertical orientation constraint
    @constraint(mdl, r[:,N] .== params.rf_targs[:,j_targ])
    @constraint(mdl, v[:,N] .== params.vf_targs[:,j_targ])
    # @constraint(mdl, T[3,N_ctrl]  == Γ[N_ctrl]) # Vertical orientation constraint

    # ..:: Solve the problem and save the solution ::..

    optimize!(mdl)
    if termination_status(mdl) != MOI.OPTIMAL
        return EmptySolution()
    end

    # Obtain optimized decision variables
    cost = objective_value(mdl)
    r = value.(r)
    v = value.(v)
    T = value.(T)
    Γ = value.(Γ)
    x = vcat(r,v)
    u = vcat(T,reshape(Γ,1,N_ctrl))

    # Package the solution
    sol = Solution(t,x,u,cost)

    return sol
end