function solve_feasible_ddto(params::Quad3DoFCageParams, τ::Int, costs_optimal::CVector, cost_dd::CReal)::Tuple{DDTOSolution, MOI.TerminationStatusCode}
    # Solve the baseline feasibility problem for DDTO.
    #
    # :in params: The params object
    # :in τ: Branch point index
    # :in costs_optimal: Optimal costs from `solve_optimal_pdg_all_targets`
    # :in cost_dd: Running cost for decision deferral
    # :out ddto_solution: Contains the DDTO solution for this target/branch point
    # :out feas_status: Feasibility problem solution status code (see MOI.TerminationStatusCode documentation)

    # ..:: Discrete time interval ::..
    if params.disc != 0 && params.disc != 1
        error("Please select a valid discretization hold order.")
    end

    N  = max(params.N_targs...)
    n  = params.n_targs
    Δt = params.Δt
    tf = Δt * (N-1)
    if params.disc == 0
        N_ctrl = N-1
        A,B,p = c2d_LTI_affine_zoh(params,Δt)
    elseif paramd.disc == 1
        error("Have not implemented FOH yet...")
        N_ctrl = N
    end


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
    @variable(mdl, r[1:3,1:N,1:n])
    @variable(mdl, v[1:3,1:N,1:n])
    @variable(mdl, T[1:3,1:N,1:n])
    @variable(mdl, Γ[1:N,1:n])

    # >> Expression holders <<
    subopt = Array{QuadExpr}(undef, N_ctrl, n)

    # >> Convenience functions <<
    X = (k,j) -> [r[:,k,j]; v[:,k,j]] # State at time index k and target j
    U = (k,j) -> [T[:,k,j]; Γ[k,j]]    # Input at time index k and target j


    # ..:: Constraints ::..

    # >> Iterate through targets <<
    for j = 1:n

        # Target N
        N_targ = params.N_targs[j]
        if params.disc == 0
            N_targ_ctrl = N-1
        elseif paramd.disc == 1
            error("Have not implemented FOH yet...")
            N_targ_ctrl = N
        end

        # Slice indexing to n without current target j
        J = collect(1:n)
        deleteat!(J, j)

        # >> Dynamics <<
        @constraint(mdl, [k=1:N_targ-1], X(k+1,j) .== A*X(k,j) + B*U(k,j) + p)

        # >> Constant altitude constraint <<
        @constraint(mdl, [k=1:N_targ-1], r[3,k+1,j] == r[3,k,j])

        # >> Thrust bounds <<
        @constraint(mdl, [k=1:N_targ_ctrl], Γ[k,j] >= params.ρ_min)
        @constraint(mdl, [k=1:N_targ_ctrl], Γ[k,j] <= params.ρ_max)
        @constraint(mdl, [k=1:N_targ_ctrl], vcat(Γ[k,j], T[:,k,j]) in MOI.SecondOrderCone(4))

        # >> Attitude pointing constraint <<
        @constraint(mdl, [k=1:N_targ_ctrl], dot(T[:,k,j],e_z) >= Γ[k,j]*cos(params.γ_p))

        # >> Velocity upper bound <<
        # @constraint(mdl, [k=1:N_targ], vcat(params.v_max_V,v[3,k,j])   in MOI.SecondOrderCone(2))
        @constraint(mdl, [k=1:N_targ], vcat(params.v_max_L,v[1:2,k,j]) in MOI.SecondOrderCone(3))

        # >> Cage bounds <<
        @constraint(mdl, [k=1:N_targ], r[1,k,j] >= params.x_arena_lims[1])
        @constraint(mdl, [k=1:N_targ], r[1,k,j] <= params.x_arena_lims[2])
        @constraint(mdl, [k=1:N_targ], r[2,k,j] >= params.y_arena_lims[1])
        @constraint(mdl, [k=1:N_targ], r[2,k,j] <= params.y_arena_lims[2])
        @constraint(mdl, [k=1:N_targ], r[3,k,j] >= params.z_arena_lims[1])
        @constraint(mdl, [k=1:N_targ], r[3,k,j] <= params.z_arena_lims[2])

        # >> Identicality << 
        for k = 1:N_targ_ctrl
            if τ > 0
                if k <= τ
                    for l = 1:n-1
                        @constraint(mdl, U(k,j) .== U(k,J[l]))
                    end
                end
            end
        end

        # >> Suboptimality <<
        for k = 1:N_ctrl
            if k <= N_targ_ctrl
                subopt[k,j] = @expression(mdl, Δt*Γ[k,j])
            else
                subopt[k,j] = @expression(mdl, 0.0)
            end
        end

        # >> Zero out state/control nodes from N_targ+1 to N <<
        @constraint(mdl, [k=N_targ+1:N],           X(k,j) .== zeros(params.n,1))
        @constraint(mdl, [k=N_targ_ctrl+1:N_ctrl], U(k,j) .== zeros(params.m,1))

        # >> Boundary conditions << 
        @constraint(mdl, X(1,j)      .== params.z0)
        @constraint(mdl, X(N_targ,j) .== params.zf_targs[:,j])

        # >> Sub-optimality <<
        @constraint(mdl, sum(subopt[:,j]) + cost_dd <= (1 + params.ϵ_targs[j]) * costs_optimal[j])
    end

    # >> Cost function <<
    @objective(mdl, Min, sum(subopt))
    # @objective(mdl, Min, 0)

    optimize!(mdl)
    feas_status = JuMP.termination_status(mdl)

    r = value.(r)
    v = value.(v)
    T = value.(T)
    Γ = value.(Γ)
    x = vcat(r,v)
    u = vcat(T,reshape(Γ,1,N,n))


    # ..:: Determine optimal cost and deferred-decision (DD) cost ::..

    costs_sol = CVector(zeros(n))
    cost_dd  = 0
    for j = 1:n
        N_targ = params.N_targs[j]
        for k = 1:N_targ-1
            costs_sol[j] += Δt*Γ[k,j]
            if k==τ && j==1
                cost_dd = costs_sol[j]
            end
        end
    end


    # ..:: Package the DDTO Solution ::..

    ddto_solution = EmptyDDTOSolution(n)
    for j = 1:n
        ddto_solution.targ_sols[j].t = CVector(range(0, stop=tf, length=params.N_targs[j]))
        ddto_solution.targ_sols[j].x = x[:,:,j]
        ddto_solution.targ_sols[j].u = u[:,:,j]
        ddto_solution.targ_sols[j].cost = costs_sol[j]
    end
    ddto_solution.costs_sol = costs_sol
    ddto_solution.cost_dd   = cost_dd

    return (ddto_solution, feas_status)

end

function solve_optimal_target(params::Quad3DoFCageParams, N::Int, j_targ::Int)::Solution
    # Solve the optimal landing (PDG) problem for a given params and single target
    # ** (Not DDTO formulation, but used for comparison) **
    #
    # :in params: The params object.
    # :in N: Time horizon (used to obtain the optimal time horizon populated into params.N_targs)
    # :in j_targ: Target index
    # :out sol: Container for solution variables


    # ..:: Discrete time interval ::..
    if params.disc != 0 && params.disc != 1
        error("Please select a valid discretization hold order.")
    end
    
    Δt = params.Δt
    tf = Δt * (N-1)
    t  = CVector(range(0, stop=tf, length=N))
    if params.disc == 0
        N_ctrl = N-1
    elseif paramd.disc == 1
        error("Have not implemented FOH yet...")
        N_ctrl = N
    end


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
    A,B,p = c2d_LTI_affine_zoh(params,Δt)
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
    @constraint(mdl, X(1) .== params.z0)
    @constraint(mdl, X(N) .== params.zf_targs[:,j_targ])

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