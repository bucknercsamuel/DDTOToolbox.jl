function solve_feasible_ddtoscp(params::Params, τ::Int, ref_costs::CVector, cost_dd::CReal, reference_targ_trajs::Vector{Solution}, scp_iter::Int)::Tuple{DDTOSolution, MOI.TerminationStatusCode, Bool}
    # Solve the baseline feasibility problem for DDTO.
    #
    # :in params: The params object
    # :in τ: Branch point index
    # :in ref_costs: Optimal costs from `solve_optimal_pdg_all_targets`
    # :in cost_dd: Running cost for decision deferral
    # :out ddto_solution: Contains the DDTO solution for this target/branch point
    # :out feas_status: Feasibility problem solution status code (see MOI.TerminationStatusCode documentation)

    # ..:: Discrete time interval ::..

    N  = max(params.N_targs...)
    n  = params.n_targs
    Δt = params.Δt
    tf = Δt * (N-1)
    N_ctrl = N-1 # Number of nodes to apply control constraints for (N-1 for ZOH)
    A,B,p = c2d_zoh(params,Δt)

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

    # >> Base optimization variables <<
    @variable(mdl, r[1:3,1:N,1:n])
    @variable(mdl, v[1:3,1:N,1:n])
    @variable(mdl, T[1:3,1:N,1:n])
    @variable(mdl, Γ[1:N,1:n])

    # >> SCP variables <<
    # Boundary condition relaxations
    # @variable(mdl, r0_relax[1:3])
    # @variable(mdl, rf_relax[1:3,1:n])

    # Virtual buffers
    @variable(mdl, ν[1:params.n_obstacles,1:N,1:n])
    @variable(mdl, μ[1:params.n_obstacles,1:N,1:n])

    # Trust region variables
    @variable(mdl, η_x[1:N])
    @variable(mdl, η_u[1:N_ctrl])

    # Slack variables for objective function
    @variable(mdl, μ_s)
    # @variable(mdl, η_s)
    @variable(mdl, η_x_s)
    @variable(mdl, η_u_s)
    # @variable(mdl, r0_relax_s)
    # @variable(mdl, rf_relax_s)

    # >> Expression holders <<
    subopt = Array{AffExpr}(undef, N_ctrl, n)

    # >> Convenience functions <<
    X = (k,j) -> [r[:,k,j]; v[:,k,j]] # State at time index k and target j
    U = (k,j) -> [T[:,k,j]; Γ[k,j]]   # Input at time index k and target j

    # ..:: Constraints ::..

    # >> Iterate through targets <<
    for j = 1:n

        # Target N
        N_targ = params.N_targs[j]
        N_targ_ctrl = N_targ - 1

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
        # @constraint(mdl, [k=1:N_targ], r[1,k,j] >= params.x_arena_lims[1])
        # @constraint(mdl, [k=1:N_targ], r[1,k,j] <= params.x_arena_lims[2])
        # @constraint(mdl, [k=1:N_targ], r[2,k,j] >= params.y_arena_lims[1])
        # @constraint(mdl, [k=1:N_targ], r[2,k,j] <= params.y_arena_lims[2])
        # @constraint(mdl, [k=1:N_targ], r[3,k,j] >= params.z_arena_lims[1])
        # @constraint(mdl, [k=1:N_targ], r[3,k,j] <= params.z_arena_lims[2])

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
        # @constraint(mdl, r[:,1,j]      .== params.r0 + r0_relax)
        @constraint(mdl, r[:,1,j]      .== params.r0)
        @constraint(mdl, v[:,1,j]      .== params.v0)
        # @constraint(mdl, r[:,N_targ,j] .== params.rf_targs[:,j] + rf_relax[:,j])
        @constraint(mdl, r[:,N_targ,j] .== params.rf_targs[:,j])
        @constraint(mdl, v[:,N_targ,j] .== params.vf_targs[:,j])

        # >> Suboptimality <<
        @constraint(mdl, sum(subopt[:,j]) + cost_dd <= (1 + params.ϵ_targs[j]) * ref_costs[j])

        # >> SCP constraints <<
        # Extract reference trajectory for target j
        x_ref = reference_targ_trajs[j].x
        u_ref = reference_targ_trajs[j].u
        r_ref = x_ref[1:3,:]
        v_ref = x_ref[4:6,:]
        T_ref = u_ref[1:3,:]
        Γ_ref = u_ref[4,:]

        # Linearization constraints
        for o = 1:params.n_obstacles
            H = params.H_obstacles[o]
            for k = 1:N
                Δr = r_ref[:,k] - params.p_obstacles[:,o]
                δr = r[:,k,j] - r_ref[:,k]
                ξ  = norm(H*Δr,2)
                ζ  = transpose(H)*H*Δr / ξ
                @constraint(mdl, ξ + dot(ζ,δr) >= params.R_obstacles[o] + ν[o,k,j])
                @constraint(mdl, vcat(μ[o,k,j], ν[o,k,j]) in MOI.NormOneCone(2))
            end
        end

        # Trust region constraints
        @constraint(mdl, [k=1:N_targ],      vcat(η_x[k], X(k,j) - x_ref[:,k]) in MOI.SecondOrderCone(params.n+1))
        @constraint(mdl, [k=1:N_targ_ctrl], vcat(η_u[k], U(k,j) - u_ref[:,k]) in MOI.SecondOrderCone(params.m+1))
    end

    # Cost function slack constraints
    @constraint(mdl, vcat(μ_s, vec(μ)) in MOI.SecondOrderCone(params.n_obstacles*N*n+1))
    @constraint(mdl, vcat(η_x_s, η_x) in MOI.SecondOrderCone(N+1))
    @constraint(mdl, vcat(η_u_s, η_u) in MOI.SecondOrderCone(N_ctrl+1))
    # @constraint(mdl, vcat(r0_relax_s, r0_relax) in MOI.SecondOrderCone(3+1))
    # @constraint(mdl, vcat(rf_relax_s, vec(rf_relax)) in MOI.SecondOrderCone(3*n+1))
    # @constraint(mdl, sum(μ.^2) <= μ_s)
    # @constraint(mdl, sum(η_x.^2) + sum(η_u.^2) <= η_s)
    # @constraint(mdl, sum(r0_relax.^2) <= r0_relax_s)
    # @constraint(mdl, sum(rf_relax.^2) <= rf_relax_s)

    # >> Cost function <<
    # @objective(mdl, Min, 
    #         # sum(subopt) + 
    #         params.w_buff * sum(μ.^2) + 
    #         params.w_trust * (sum(η_x.^2) + sum(η_u.^2)) +
    #         params.w_r0 * sum(r0_relax.^2) + 
    #         params.w_rf * sum(rf_relax.^2))
    # @objective(mdl, Min, 0)
    @objective(mdl, Min, 
        sum(subopt) + 
        params.w_buff * μ_s + 
        params.w_trust * (η_x_s + η_u_s))
        # params.w_r0 * r0_relax_s + 
        # params.w_rf * rf_relax_s)

    optimize!(mdl)
    feas_status = JuMP.termination_status(mdl)

    r = value.(r)
    v = value.(v)
    T = value.(T)
    Γ = value.(Γ)
    ν = value.(ν)
    μ = value.(μ)
    η_x = value.(η_x)
    η_u = value.(η_u)
    η = [η_x;η_u]
    x = vcat(r,v)
    u = vcat(T,reshape(Γ,1,N,n))

    # ..:: Determine if PTR subproblem has converged ::..
    if feas_status == MOI.OPTIMAL
        solve_status = "Feasible"
    else
        solve_status = "Not Feasible"
    end
    μ_max_nodes = []
    for k = 1:N
        append!(μ_max_nodes, max(μ[:,k,:]...,0))
    end
    μ_pen = sum(μ_max_nodes)
    η_pen = norm(η,2)

    @printf("   SCP Iter: %2.i | Status: %s | μ_pen = %.2e | η_pen = %.2e\n", scp_iter, solve_status, μ_pen, η_pen)
    if (μ_pen <= params.ϵ_cvg) && (η_pen <= params.ϵ_cvg)
        scp_sub_cvged = true
    else
        scp_sub_cvged = false
    end

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

    return (ddto_solution, feas_status, scp_sub_cvged)

end