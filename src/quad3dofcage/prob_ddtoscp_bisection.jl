function solve_ddtoscp(params::Quad3DoFCageParams)

    # Require free-final-time formulation for current main method
    # (Time dilation necessary for deferral)
    if !(params.free_final_time)
        error("Please select free-final-time formulation for DDTO-SCP.")
    end

    # Augment boundary states if using free-final-time (do not want affine form)
    if params.free_final_time
        params.z0 = vcat(params.z0, 1)
        params.zf_targs = vcat(params.zf_targs, ones(1,params.n_targs))
    end

    # ..:: Execute solver sequence ::..
    @time begin
        @time begin
            # ..:: Solve for independently-optimal solutions to each target ::..
            scp_solutions = solve_decoupled_scp_tree(params)
            scp_costs = CVector(zeros(params.n_targs))
            for k = 1:params.n_targs
                scp_costs[k] = scp_solutions[k].cost
            end
            println("\n Solve time for generating optimal solutions to each target:")
        end

        @time begin
            # ..:: Solve for DDTO branching solutions to ALL targets ::..
            (feas_ddtoscp, ddtoscp_solutions) = solve_ddtoscp_tree(deepcopy(params), scp_costs, deepcopy(scp_solutions))
            println("\n Solve time for generating DDTO branch solutions to all targets:")
        end
        println("\n Solve time for the full DDTO solution stack:")
    end

    # Convert DDTO solutions to branch solutions
    ddtoscp_bsolutions, defer_bsolutions = extract_target_trajectories(params, ddtoscp_solutions)

    # Port decoupled SCP solutions to `BranchSolution` objects for type conformance
    scp_bsolutions_ = Vector{BranchSolution}(undef,params.n_targs)
    for j=1:params.n_targs
        scp_bsolutions_[j] = BranchSolution(scp_solutions[j],-1,-1)
    end
    scp_bsolutions = scp_bsolutions_

    # ..:: Simulate each target solution from I.C. to T.C.
    @time begin
        if params.free_final_time
            dynamics = (t,x,sol) -> dyn_nl(t,x,optimal_controller(t,sol.t,sol.u,params.disc),params)
        else
            dynamics = (t,x,sol) -> params.A_c*x + params.B_c*optimal_controller(t,sol.t,sol.u,params.disc) + params.p_c
        end
        scp_bsimulations = simulate_branches(scp_bsolutions, dynamics, params.disc)
        ddtoscp_bsimulations = simulate_branches(ddtoscp_bsolutions, dynamics, params.disc)
        defer_bsimulations = simulate_branches(defer_bsolutions, dynamics, params.disc)
        println("\n Solve time for RK4 simulation:")
    end

    # ..:: Post-processing ::..
    @time begin
        scp_solutions_proc       = process_solutions(scp_bsolutions, params)
        scp_simulations_proc     = process_solutions(scp_bsimulations, params)
        ddtoscp_solutions_proc   = process_solutions(ddtoscp_bsolutions, params)
        ddtoscp_simulations_proc = process_solutions(ddtoscp_bsimulations, params)
        defer_solutions_proc     = process_solutions(defer_bsolutions, params)
        defer_simulations_proc   = process_solutions(defer_bsimulations, params)
        println("\n Solve time for post-processing:")
    end

    # Deferrable segment solution/simulation should just be a scalar ProcessedSolution
    defer_bsolutions = defer_bsolutions[1].sol
    defer_bsimulations = defer_bsimulations[1].sol
    defer_solutions_proc = defer_solutions_proc[1].sol
    defer_simulations_proc = defer_simulations_proc[1].sol

    return (scp_solutions_proc, scp_simulations_proc, ddtoscp_solutions_proc, ddtoscp_simulations_proc, defer_solutions_proc, defer_simulations_proc)
end

function solve_feasible_ddtoscp(params::Quad3DoFCageParams, τ::Int, ref_costs::CVector, cost_dd::CReal, reference_targ_trajs::Vector{Solution}, scp_iter::Int)::Tuple{DDTOSolution, MOI.TerminationStatusCode, Bool}
    # Solve the baseline feasibility problem for DDTO.
    #
    # :in params: The params object
    # :in τ: Branch point index
    # :in ref_costs: Optimal costs from `solve_optimal_pdg_all_targets`
    # :in cost_dd: Running cost for decision deferral
    # :out ddto_solution: Contains the DDTO solution for this target/branch point
    # :out feas_status: Feasibility problem solution status code (see MOI.TerminationStatusCode documentation)

    # ..:: Discrete time interval ::..
    if params.disc != 0 && params.disc != 1
        error("Please select a valid discretization hold order.")
    end

    n = params.n_targs
    if params.free_final_time
        N = params.N_fft
        n_ = params.n+1
        m_ = params.m+1
    else
        N  = max(params.N_targs...)
        Δt = params.Δt
        tf = Δt * (N-1)
        n_ = params.n
        m_ = params.m
    end
    if params.disc == 0
        N_ctrl = N-1
        if !params.free_final_time
            A,B,p = c2d_LTI_affine_zoh(params.A_c,params.B_c,params.p_c,Δt)
        end
    elseif params.disc == 1
        N_ctrl = N
        if !params.free_final_time
            A,Bm,Bp,p = c2d_LTI_affine_foh(params.A_c, params.B_c, params.p_c, Δt)
        end
    end

    # ..:: Make the optimization problem ::..

    # >> Optimizer setup <<
    if SOLVER == "ECOS"
        mdl = Model(optimizer_with_attributes(ECOS.Optimizer, "verbose" => 0))
    elseif SOLVER == "MOSEK"
        mdl = Model(Mosek.Optimizer)
        JuMP.set_optimizer_attribute(mdl, "LOG", 0) # disable debugging
        JuMP.set_optimizer_attribute(mdl, "MAX_NUM_WARNINGS", 0) # disable warnings
    else
        error("SOLVER is invalid, please select either ECOS or MOSEK")
    end

    # >> Base optimization variables <<
    @variable(mdl, r[1:3,1:N,1:n])
    @variable(mdl, v[1:3,1:N,1:n])
    @variable(mdl, T[1:3,1:N_ctrl,1:n])
    @variable(mdl, Γ[1:N_ctrl,1:n])
    if params.free_final_time
        @variable(mdl, s[1:N_ctrl,1:n])
        Δt = Array{AffExpr}(undef,N-1,n)
    end
    
    # >> SCP variables <<

    # Virtual buffers
    @variable(mdl, ν_obs[1:params.n_obstacles,1:N,1:n])
    @variable(mdl, μ_obs[1:params.n_obstacles,1:N,1:n])

    # Virtual control
    @variable(mdl, ν_ctrl[1:(n_),1:(N-1),1:n])
    @variable(mdl, μ_ctrl[1:(n_),1:(N-1),1:n])

    # Trust region variables
    @variable(mdl, η_x[1:N])
    @variable(mdl, η_u[1:N_ctrl])

    # Slack variables for objective function
    @variable(mdl, μ_obs_s)
    @variable(mdl, μ_ctrl_s)
    @variable(mdl, η_x_s)
    @variable(mdl, η_u_s)

    # >> Expression holders <<
    subopt = Array{AffExpr}(undef, N_ctrl, n)

    # >> Convenience functions <<
    if !params.free_final_time
        X = (k,j) -> [r[:,k,j]; v[:,k,j]] # State at time index k and target j
        U = (k,j) -> T[:,k,j]             # Input at time index k and target j
    else
        X = (k,j) -> [r[:,k,j]; v[:,k,j]; 1] # Augmented state (to bring in affine term)
        U = (k,j) -> [T[:,k,j]; s[k,j]]      # Augmented control (with time dilation term)
    end

    # ..:: Constraints ::..

    # Dynamics functions (same for all targets)
    dyn_lin_ = (t,x,u) -> dyn_lin(t,x,u,params)
    dyn_nl_  = (t,x,u) -> dyn_nl(t,x,u,params)

    # >> Iterate through targets <<
    for j = 1:n

        # Target N
        if params.free_final_time
            N_targ = params.N_fft
        else
            N_targ = params.N_targs[j]
        end
        if params.disc == 0
            N_targ_ctrl = N_targ - 1
        elseif params.disc == 1
            N_targ_ctrl = N_targ
        end

        # Slice indexing to n without current target j
        J = collect(1:n)
        deleteat!(J, j)

        # >> Convex State & Control Constraints <<

        # Dynamics (fixed-final-time)
        if !params.free_final_time
            if params.disc == 0
                @constraint(mdl, [k=1:N_targ-1], X(k+1,j) .== A*X(k,j) + B*U(k,j) + p)
            elseif params.disc == 1
                @constraint(mdl, [k=1:N_targ-1], X(k+1,j) .== A*X(k,j) + Bm*U(k,j) + Bp*U(k+1,j) + p)
            end
        end

        # Constant altitude
        @constraint(mdl, [k=1:N_targ-1], r[3,k+1,j] == r[3,k,j])

        # Thrust bounds
        # @constraint(mdl, [k=1:N_targ_ctrl], Γ[k,j] >= params.ρ_min)
        # @constraint(mdl, [k=1:N_targ_ctrl], Γ[k,j] <= params.ρ_max)
        # @constraint(mdl, [k=1:N_targ_ctrl], vcat(Γ[k,j], T[:,k,j]) in MOI.SecondOrderCone(4))
        @constraint(mdl, [k=1:N_targ_ctrl], vcat(params.ρ_max, T[:,k,j]) in MOI.SecondOrderCone(4))
        @constraint(mdl, [k=1:N_targ_ctrl], vcat(Γ[k,j], T[:,k,j]) in MOI.SecondOrderCone(4))

        # Attitude pointing
        # @constraint(mdl, [k=1:N_targ_ctrl], dot(T[:,k,j],e_z) >= Γ[k,j]*cos(params.γ_p))
        @constraint(mdl, [k=1:N_targ_ctrl], vcat(dot(T[:,k,j],e_z)/cos(params.γ_p), T[:,k,j]) in MOI.SecondOrderCone(4))

        # Velocity upper bound
        # @constraint(mdl, [k=1:N_targ], vcat(params.v_max_V,v[3,k,j])   in MOI.SecondOrderCone(2))
        @constraint(mdl, [k=1:N_targ], vcat(params.v_max_L,v[1:2,k,j]) in MOI.SecondOrderCone(3))

        # Cage bounds
        # @constraint(mdl, [k=1:N_targ], r[1,k,j] >= params.x_arena_lims[1])
        # @constraint(mdl, [k=1:N_targ], r[1,k,j] <= params.x_arena_lims[2])
        # @constraint(mdl, [k=1:N_targ], r[2,k,j] >= params.y_arena_lims[1])
        # @constraint(mdl, [k=1:N_targ], r[2,k,j] <= params.y_arena_lims[2])
        # @constraint(mdl, [k=1:N_targ], r[3,k,j] >= params.z_arena_lims[1])
        # @constraint(mdl, [k=1:N_targ], r[3,k,j] <= params.z_arena_lims[2])

        # Identicality
        for k = 1:N_targ_ctrl
            if τ > 0
                if k <= τ
                    for l = 1:n-1
                        @constraint(mdl, U(k,j) .== U(k,J[l]))
                    end
                end
            end
        end

        # Suboptimality
        for k = 1:N_ctrl
            if k <= N_targ_ctrl
                subopt[k,j] = @expression(mdl, Γ[k,j])
            else
                subopt[k,j] = @expression(mdl, 0.0)
            end
        end

        # Zero out state/control nodes from N_targ+1 to N
        @constraint(mdl, [k=N_targ+1:N],           X(k,j) .== zeros(n_,1))
        @constraint(mdl, [k=N_targ_ctrl+1:N_ctrl], U(k,j) .== zeros(m_,1))

        # Boundary conditions
        z0 = params.z0
        zf = params.zf_targs[:,j]
        @constraint(mdl, X(1,j)      .== z0)
        @constraint(mdl, X(N_targ,j) .== zf)

        # Suboptimality
        @constraint(mdl, sum(subopt[:,j]) + cost_dd <= (1 + params.ϵ_targs[j]) * ref_costs[j])

        # Time dilation
        if params.free_final_time
            for k=1:(N-1)
                if params.disc == 0
                    Δt[k,j] = @expression(mdl, params.Δτ[k] * s[k,j])
                elseif params.disc == 1
                    Δt[k,j] = @expression(mdl, (1/2) * params.Δτ[k] * (s[k,j] + s[k+1,j]))
                end
            end
            @constraint(mdl, sum(Δt[:,j]) <= params.ToF_max)
            @constraint(mdl, [k=1:N-1], params.Δt_min <= Δt[k,j] <= params.Δt_max)
            @constraint(mdl, [k=1:N_ctrl], params.s_min <= s[k,j] <= params.s_max)
        end

        # >> PTR constraints <<
        
        # Extract reference trajectory for target j
        ref_traj = reference_targ_trajs[j]
        t_ref = ref_traj.t
        x_ref = ref_traj.x
        u_ref = ref_traj.u
        r_ref = x_ref[1:3,:]
        v_ref = x_ref[4:6,:]
        T_ref = u_ref[1:3,:]
        # Γ_ref = u_ref[4,:]

        # Dynamics (free-final-time)
        if params.free_final_time
            # Obtain approximate LTV discrete-time dynamics
            if params.disc == 0
                Ak,Bk,_,wk,_ = c2d_nonlinear(ref_traj,dyn_nl_,dyn_lin_,params.disc)
            elseif params.disc == 1
                Ak,Bmk,Bpk,wk,_ = c2d_nonlinear(ref_traj,dyn_nl_,dyn_lin_,params.disc)
            end

            # Apply constraints
            if params.disc == 0
                @constraint(mdl, [k=1:N-1], X(k+1,j) .== Ak[:,:,k]*X(k,j) + Bk[:,:,k]*U(k,j) + wk[:,k] + ν_ctrl[:,k,j])
            elseif params.disc == 1
                @constraint(mdl, [k=1:N-1], X(k+1,j) .== Ak[:,:,k]*X(k,j) + Bmk[:,:,k]*U(k,j) + Bpk[:,:,k]*U(k+1,j) + wk[:,k] + ν_ctrl[:,k,j])
            end
            @constraint(mdl, [k=1:N-1,c=1:n_], vcat(μ_ctrl[c,k,j], ν_ctrl[c,k,j]) in MOI.NormOneCone(2))
        end

        # Ellipsoidal obstacle constraints
        for o = 1:params.n_obstacles
            H = params.H_obstacles[o]
            for k = 1:N
                Δr = r_ref[:,k] - params.p_obstacles[:,o]
                δr = r[:,k,j] - r_ref[:,k]
                ξ  = max(norm(H*Δr,2),1e-4)
                ζ  = transpose(H)*H*Δr / ξ
                @constraint(mdl, ξ + dot(ζ,δr) >= params.R_obstacles[o] + ν_obs[o,k,j])
                @constraint(mdl, vcat(μ_obs[o,k,j], ν_obs[o,k,j]) in MOI.NormOneCone(2))
            end
        end

        # Trust region constraints
        @constraint(mdl, [k=1:N_targ],      vcat(η_x[k], X(k,j) - x_ref[:,k]) in MOI.SecondOrderCone(n_+1))
        @constraint(mdl, [k=1:N_targ_ctrl], vcat(η_u[k], U(k,j) - u_ref[:,k]) in MOI.SecondOrderCone(m_+1))

        # Maintain continuity for FOH discretization if we have already deferred by some amount (usually after first branch iteration)
        # (acts as a initial condition on control constrained to the reference initial condition)
        if params.disc == 1 && cost_dd > 0
            @constraint(mdl, U(1,j) .== u_ref[:,1])
        end
    end

    # Cost function slack constraints
    @constraint(mdl, vcat(μ_obs_s, vec(μ_obs)) in MOI.SecondOrderCone(params.n_obstacles*N*n+1))
    @constraint(mdl, vcat(μ_ctrl_s, vec(μ_ctrl)) in MOI.SecondOrderCone((n_)*(N-1)*n+1))
    @constraint(mdl, vcat(η_x_s, η_x) in MOI.SecondOrderCone(N+1))
    @constraint(mdl, vcat(η_u_s, η_u) in MOI.SecondOrderCone(N_ctrl+1))

    # >> Cost function <<
    J_opt  = sum(subopt)
    J_ptr  = params.w_trust * (η_x_s + η_u_s)
    J_buff = params.w_buff * μ_obs_s
    J_ctrl = params.w_ctrl * μ_ctrl_s
    if !params.free_final_time
        J_ctrl = 0
    end
    if params.n_obstacles == 0
        J_buff = 0
    end
    @objective(mdl, Min, J_opt + J_ptr + J_buff + J_ctrl)
        
    optimize!(mdl)
    feas_status = JuMP.termination_status(mdl)

    r = value.(r)
    v = value.(v)
    T = value.(T)
    Γ = value.(Γ)
    ν_obs = value.(ν_obs)
    μ_obs = value.(μ_obs)
    η_x = value.(η_x)
    η_u = value.(η_u)
    η = [η_x;η_u]
    if params.free_final_time
        s = value.(s)
        x = vcat(r,v,ones(1,N,n))
        u = vcat(T,reshape(s,1,N_ctrl,n))
    else
        x = vcat(r,v)
        u = vcat(T)
    end

    # ..:: Determine if PTR subproblem has converged ::..
    if feas_status == MOI.OPTIMAL
        solve_status = "Feasible"
    else
        solve_status = "Not Feasible"
    end
    μ_obs_max_nodes = []
    for k = 1:N
        append!(μ_obs_max_nodes, max(μ_obs[:,k,:]...,0))
    end
    μ_obs_pen = sum(μ_obs_max_nodes)
    η_pen = norm(η,2)
    μ_ctrl_pen = value.(μ_ctrl_s)

    if feas_status == MOI.OPTIMAL && (μ_ctrl_pen <= params.ϵ_ctrl) && (μ_obs_pen <= params.ϵ_buff) && (η_pen <= params.ϵ_trust)
        scp_sub_cvged = true
    else
        scp_sub_cvged = false
    end
    @printf("   SCP Iter: %2.i | Status: %s | Cost = %.2e | μ_ctrl_pen = %.2e | μ_obs_pen = %.2e | η_pen = %.2e\n", scp_iter, solve_status, value.(J_opt), μ_ctrl_pen, μ_obs_pen, η_pen)
    flush(stdout)

    # ..:: Determine optimal cost and deferred-decision (DD) cost ::..

    costs_sol = CVector(zeros(n))
    cost_dd  = 0
    for j = 1:n
        if params.free_final_time
            N_targ = params.N_fft
        else
            N_targ = params.N_targs[j]
        end
        N_targ_ctrl = N_targ - 1
        for k = 1:N_targ_ctrl
            costs_sol[j] += Γ[k,j]
            if k==τ && j==1
                cost_dd = costs_sol[j]
            end
        end
    end

    # ..:: Package the DDTO Solution ::..

    ddto_solution = EmptyDDTOSolution(n)
    for j = 1:n
        if params.free_final_time
            Δtj = value.(Δt[:,j])
            t = vcat(0,cumsum(Δtj))
        else
            t = CVector(range(0, stop=tf, length=params.N_targs[j]))
        end

        ddto_solution.targ_sols[j].t = t
        ddto_solution.targ_sols[j].x = x[:,:,j]
        ddto_solution.targ_sols[j].u = u[:,:,j]
        ddto_solution.targ_sols[j].cost = costs_sol[j]
    end
    ddto_solution.costs_sol = costs_sol
    ddto_solution.cost_dd   = cost_dd

    return (ddto_solution, feas_status, scp_sub_cvged)

end