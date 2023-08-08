function solve_ddtoscp(params::Quad3DoFCageParams)

    # Require free-final-time formulation for current main method
    # (Time dilation necessary for deferral)
    if !(params.free_final_time)
        error("Please select free-final-time formulation for DDTO-SCP.")
    end

    # Augment boundary states do not want affine form)
    params.z0 = vcat(params.z0, 1)
    params.zf_targs = vcat(params.zf_targs, ones(1,params.n_targs))

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
            initial_guess = generate_initial_guess_ddtoscp(params)
            (feas_ddtoscp, ddtoscp_solutions) = solve_ddtoscp_tree(deepcopy(params), scp_costs, initial_guess)
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
        dynamics = (t,x,sol) -> dyn_nl(t,x,optimal_controller(t,sol.t,sol.u,params.disc),params)
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
    N = params.N_fft
    nx = params.n+1
    nu = params.m+1
    if params.disc == 0
        N_ctrl = N-1
    elseif params.disc == 1
        N_ctrl = N
    end

    # ..:: Make the optimization problem ::..

    # >> Optimizer setup <<
    if SOLVER == "ECOS"
        mdl = Model(optimizer_with_attributes(ECOS.Optimizer, "verbose" => 0))
    elseif SOLVER == "MOSEK"
        mdl = Model(Mosek.Optimizer)
        JuMP.set_optimizer_attribute(mdl, "LOG", 0) # disable debugging
        JuMP.set_optimizer_attribute(mdl, "MAX_NUM_WARNINGS", 0) # disable warni
    else
        error("SOLVER is invalid, please select either ECOS or MOSEK")
    end

    # >> Scaled optimization variables <<
    @variable(mdl, r_s[1:3,1:N,1:n+1])
    @variable(mdl, v_s[1:3,1:N,1:n+1])
    @variable(mdl, T_s[1:3,1:N_ctrl,1:n+1])
    @variable(mdl, s_s[1:N_ctrl,1:n+1])
    
    # >> Unscaling <<
    rmin = [params.x_arena_lims[1]; params.y_arena_lims[1]; params.z_arena_lims[1]]
    rmax = [params.x_arena_lims[2]; params.y_arena_lims[2]; params.z_arena_lims[2]]
    r = unscale(r_s, rmin, rmax)
    v = unscale(v_s, -params.v_max_L, params.v_max_L)
    T = unscale(T_s, -params.ρ_max, params.ρ_max)
    s = unscale(s_s, params.s_min, params.s_max)
    Sx,_ = scaling_matrices([rmin; -params.v_max_L*ones(3);1], [rmax; params.v_max_L*ones(3);1])
    Su,_ = scaling_matrices([-params.ρ_max*ones(3); params.s_min], [params.ρ_max*ones(3); params.s_max])
    SxInv, SuInv = inv(Sx), inv(Su)

    # >> SCP variables <<zzz
    @variable(mdl, ν_obs[1:params.n_obstacles,1:N,1:n+1])
    @variable(mdl, μ_obs)
    @variable(mdl, ν_ctrl[1:nx,1:N-1,1:n+1])
    @variable(mdl, μ_ctrl)
    @variable(mdl, η[1:N,1:n+1])
    @variable(mdl, η_s)

    # >> Expressions <<
    Δt = Array{AffExpr}(undef,N-1,n+1)
    subopt_trunk = Array{QuadExpr}(undef, N_ctrl, 1)
    subopt_branch = Array{QuadExpr}(undef, N_ctrl, n)

    # Indexing
    J_branch = 1:n
    J_trunk = n+1

    # >> Convenience functions <<
    X(k,j) = [r[:,k,j]; v[:,k,j]; 1] # Augmented state (to bring in affine term)
    U(k,j) = [T[:,k,j]; s[k,j]]      # Augmented control (with time dilation term)

    # ..:: Constraints for all segments ::..

    # Dynamics functions
    dyn_lin_ = (t,x,u) -> dyn_lin(t,x,u,params)
    dyn_nl_  = (t,x,u) -> dyn_nl(t,x,u,params)

    for j = 1:n+1

        # >> Reference trajectory <<
        if j == J_trunk
            ref_traj = reference_targ_trajs[1] # Any solution works
            t_ref = ref_traj.t[1:N]
            x_ref = ref_traj.x[:,1:N]
            u_ref = ref_traj.u[:,1:N]
            r_ref = x_ref[1:3,:]
            v_ref = x_ref[4:6,:]
            T_ref = u_ref[1:3,:]
        else
            ref_traj = reference_targ_trajs[j]
            t_ref = ref_traj.t[N+1:end]
            x_ref = ref_traj.x[:,N+1:end]
            u_ref = ref_traj.u[:,N+1:end]
            r_ref = x_ref[1:3,:]
            v_ref = x_ref[4:6,:]
            T_ref = u_ref[1:3,:]
        end

        # >> Dynamics <<
        Ak,Bmk,Bpk,wk,_ = c2d_nonlinear(ref_traj,dyn_nl_,dyn_lin_,params.disc)
        if params.disc == 0
            @constraint(mdl, [k=1:N-1], SxInv*X(k+1,j) .== SxInv*(Ak[:,:,k]*X(k,j) + Bk[:,:,k]*U(k,j) + wk[:,k]) + ν_ctrl[:,k,j])
        elseif params.disc == 1
            @constraint(mdl, [k=1:N-1], SxInv*X(k+1,j) .== SxInv*(Ak[:,:,k]*X(k,j) + Bmk[:,:,k]*U(k,j) + Bpk[:,:,k]*U(k+1,j) + wk[:,k]) + ν_ctrl[:,k,j])
        end

        # >> Global State & Control Constraints <<
        # Constant altitude
        @constraint(mdl, [k=1:N-1], r[3,k+1,j] == r[3,k,j])

        # Thrust bounds
        # Χ(k) = normalize(T_ref[:,k])
        @constraint(mdl, [k=1:N_ctrl], vcat(params.ρ_max, T[:,k,j]) in SecondOrderCone())
        # @constraint(mdl, [k=1:N_ctrl], dot(Χ(k),T[:,k,j]) >= params.ρ_min)

        # Attitude pointing
        @constraint(mdl, [k=1:N_ctrl], vcat(dot(T[:,k,j],e_z)/cos(params.γ_p), T[:,k,j]) in SecondOrderCone())

        # Velocity upper bound
        @constraint(mdl, [k=1:N], vcat(params.v_max_L,v[1:2,k,j]) in SecondOrderCone())

        # Cage bounds
        @constraint(mdl, [k=1:N], r[1,k,j] >= params.x_arena_lims[1])
        @constraint(mdl, [k=1:N], r[1,k,j] <= params.x_arena_lims[2])
        @constraint(mdl, [k=1:N], r[2,k,j] >= params.y_arena_lims[1])
        @constraint(mdl, [k=1:N], r[2,k,j] <= params.y_arena_lims[2])
        @constraint(mdl, [k=1:N], r[3,k,j] >= params.z_arena_lims[1])
        @constraint(mdl, [k=1:N], r[3,k,j] <= params.z_arena_lims[2])

        # Time dilation
        for k=1:(N-1)
            if params.disc == 0
                Δt[k,j] = @expression(mdl, params.Δτ[k] * s[k,j])
            elseif params.disc == 1
                Δt[k,j] = @expression(mdl, (1/2) * params.Δτ[k] * (s[k,j] + s[k+1,j]))
            end
        end
        @constraint(mdl, sum(Δt[:,j]) <= params.ToF_max)
        @constraint(mdl, [k=1:N_ctrl], params.s_min <= s[k,j] <= params.s_max)

        # Ellipsoidal obstacle constraints
        for o = 1:params.n_obstacles
            H = params.H_obstacles[o]
            for k = 1:N
                Δr = r_ref[:,k] - params.p_obstacles[:,o]
                δr = r[:,k,j] - r_ref[:,k]
                ξ  = max(norm(H*Δr,2),1e-4)
                ζ  = transpose(H)*H*Δr / ξ
                @constraint(mdl, ξ + dot(ζ,δr) >= params.R_obstacles[o] + ν_obs[o,k,j])
            end
        end

        # Trust region constraints
        δX(k) = SxInv*(X(k,j) .- x_ref[:,k])
        δU(k) = SuInv*(U(k,j) .- u_ref[:,k])
        @constraint(mdl, [k=1:N_ctrl], δX(k)'*δX(k) + δU(k)'*δU(k) <= η[k,j])
        @constraint(mdl, δX(N)'*δX(N) <= η[N,j])
    end

    # Maintain continuity for FOH discretization if we have already deferred by some amount (usually after first branch iteration)
    # (acts as a initial condition on control constrained to the reference initial condition)
    if params.disc == 1 && cost_dd > 0
        u_ref = reference_targ_trajs[1].u # any traj works
        @constraint(mdl, U(1,J_trunk) .== u_ref[:,1])
    end

    # ..:: DDTO Stitching Constraints ::..
    # >> Trunk <<
    # Add trunk suboptimality
    for k = 1:N_ctrl
        subopt_trunk[k] = @expression(mdl, -sum(T[:,k,J_trunk]))
    end

    # Initial condition
    @constraint(mdl, X(1,J_trunk) .== params.z0)

    # >> Branches <<
    for j in J_branch
        # Add branch suboptimality
        for k = 1:N_ctrl
            subopt_branch[k,j] = @expression(mdl, -sum(T[:,k,j]))
        end

        # Apply suboptimality constraint
        @constraint(mdl, sum(subopt_branch[:,j]) + sum(subopt_trunk) + cost_dd <= (1 + params.ϵ_targs[j]) * ref_costs[j])

        # Apply boundary conditions
        @constraint(mdl, X(1,j) .== X(N,J_trunk))
        @constraint(mdl, X(N,j) .== params.zf_targs[:,j])
    end

    # ..:: PTR Constraints ::..
    @constraint(mdl, vcat(μ_obs, vec(ν_obs)) in MOI.NormOneCone(params.n_obstacles*N*(n+1)+1))
    @constraint(mdl, vcat(μ_ctrl, vec(ν_ctrl)) in MOI.NormOneCone(nx*(N-1)*(n+1)+1))
    @constraint(mdl, vcat(η_s, vec(η)) in SecondOrderCone())
    @constraint(mdl, μ_obs >= 0)
    @constraint(mdl, μ_ctrl >= 0)
    @constraint(mdl, η_s >= 0)

    # >> Cost function <<
    J_opt  = -sum(s[:,J_trunk])
    J_ptr  = η_s
    J_buff = μ_obs
    J_ctrl = μ_ctrl
    if params.n_obstacles == 0
        J_buff = 0
    end
    @objective(mdl, Min, J_opt + params.w_trust * J_ptr + params.w_buff * J_buff + params.w_ctrl * J_ctrl)
        
    optimize!(mdl)
    feas_status = JuMP.termination_status(mdl)

    # display(sum(value.(subopt_trunk)))
    # display(value.(ref_costs))
    # display([sum(value.(subopt_branch[:,j])) for j = 1:n])

    r = value.(r)
    v = value.(v)
    T = value.(T)
    s = value.(s)
    ν_obs = value.(ν_obs)
    μ_obs = value.(μ_obs)

    r_trunk = r[:,:,J_trunk]
    v_trunk = v[:,:,J_trunk]
    T_trunk = T[:,:,J_trunk]
    s_trunk = s[:,J_trunk]
    r_branch = r[:,:,J_branch]
    v_branch = v[:,:,J_branch]
    T_branch = T[:,:,J_branch]
    s_branch = s[:,J_branch]

    x = zeros(nx, 2*N, n)
    u = zeros(nu, N+N_ctrl, n)

    for j in J_branch
        x[:,1:N,j] = vcat(r_trunk, v_trunk, ones(1,N))
        u[:,1:N,j] = vcat(T_trunk, reshape(s_trunk,1,N_ctrl))
        x[:,N+1:end,j] = vcat(r_branch[:,:,j], v_branch[:,:,j], ones(1,N))
        u[:,N+1:end,j] = vcat(T_branch[:,:,j], reshape(s_branch[:,j],1,N_ctrl))
    end

    # ..:: Determine if PTR subproblem has converged ::..
    if feas_status == MOI.OPTIMAL || feas_status == MOI.ALMOST_OPTIMAL
        solve_status = "Feasible"
    else
        solve_status = "Not Feasible"
    end
    
    # Obtain evaluation penalties
    μ_obs_pen = value.(μ_obs)
    μ_ctrl_pen = value.(μ_ctrl)
    η_pen = value.(η_s)

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
        N_ctrl = N - 1
        for k = 1:N_ctrl
            costs_sol[j] += sum(T[:,k,j].^2)
            if k==τ && j==1
                cost_dd = costs_sol[j]
            end
        end
    end

    # ..:: Package the DDTO Solution ::..

    ddto_solution = EmptyDDTOSolution(n)
    Δt_trunk = value.(Δt[:,J_trunk])
    for j = 1:n
        Δtj = value.(Δt[:,J_branch[j]])
        t = vcat(0,cumsum([Δt_trunk;0;Δtj]))
        ddto_solution.targ_sols[j].t = t
        ddto_solution.targ_sols[j].x = x[:,:,j]
        ddto_solution.targ_sols[j].u = u[:,:,j]
        ddto_solution.targ_sols[j].cost = costs_sol[j]
    end
    ddto_solution.costs_sol = costs_sol
    ddto_solution.cost_dd   = cost_dd

    return (ddto_solution, feas_status, scp_sub_cvged)

end

function solve_scp_target(params::Quad3DoFCageParams, ref_traj::Solution, N::Int, j_targ::Int, scp_iter::Int)::Tuple{Solution, MOI.TerminationStatusCode, Bool}

    # ..:: Discrete time interval ::..
    if params.disc != 0 && params.disc != 1
        error("Please select a valid discretization hold order.")
    end
    
    nx = params.n+1
    nu = params.m+1
    if params.disc == 0
        N_ctrl = N-1
    elseif params.disc == 1
        N_ctrl = N
    end


    # ..:: Make the optimization problem ::..

    # >> Optimizer setup <<
    if SOLVER == "ECOS"
        mdl = Model(optimizer_with_attributes(ECOS.Optimizer, "verbose" => 0))
    elseif SOLVER == "MOSEK"
        mdl = Model(Mosek.Optimizer)
        JuMP.set_optimizer_attribute(mdl, "LOG",  0) # disable debugging
        JuMP.set_optimizer_attribute(mdl, "MAX_NUM_WARNINGS", 0) # disable warnings
    else
        error("SOLVER is invalid, please select either ECOS or MOSEK")
    end

    # >> Scaled Optimization variables <<
    @variable(mdl, r_s[1:3,1:N])
    @variable(mdl, v_s[1:3,1:N])
    @variable(mdl, T_s[1:3,1:N_ctrl])
    @variable(mdl, s_s[1:N_ctrl])
    Δt = Array{AffExpr}(undef,N-1)

    # >> Unscaling <<
    rmin = [params.x_arena_lims[1]; params.y_arena_lims[1]; params.z_arena_lims[1]]
    rmax = [params.x_arena_lims[2]; params.y_arena_lims[2]; params.z_arena_lims[2]]
    r = unscale(r_s, rmin, rmax)
    v = unscale(v_s, -params.v_max_L, params.v_max_L)
    T = unscale(T_s, -params.ρ_max, params.ρ_max)
    s = unscale(s_s, params.s_min, params.s_max)
    Sx,_ = scaling_matrices([rmin; -params.v_max_L*ones(3);1], [rmax; params.v_max_L*ones(3);1])
    Su,_ = scaling_matrices([-params.ρ_max*ones(3); params.s_min], [params.ρ_max*ones(3); params.s_max])
    SxInv, SuInv = inv(Sx), inv(Su)

    # >> SCP variables <<
    # Virtual buffers
    @variable(mdl, ν_obs[1:params.n_obstacles,1:N])
    @variable(mdl, μ_obs)
    @variable(mdl, ν_ctrl[1:nx,1:(N-1)])
    @variable(mdl, μ_ctrl)
    @variable(mdl, η[1:N])
    @variable(mdl, η_s)

    # >> Convenience functions <<
    X = (k) -> [r[:,k]; v[:,k]; 1] # Augmented state (to bring in affine term)
    U = (k) -> [T[:,k]; s[k]] # Augmented control (with time dilation term)

    # Extract reference trajectory
    t_ref = ref_traj.t
    x_ref = ref_traj.x
    u_ref = ref_traj.u
    r_ref = x_ref[1:3,:]
    v_ref = x_ref[4:6,:]
    T_ref = u_ref[1:3,:]

    # ..:: Constraints ::..

    # >> Dynamics <<
    dyn_lin_ = (t,x,u) -> dyn_lin(t,x,u,params)
    dyn_nl_  = (t,x,u) -> dyn_nl(t,x,u,params)
    Ak,Bmk,Bpk,wk,_ = c2d_nonlinear(ref_traj,dyn_nl_,dyn_lin_,params.disc)
    if params.disc == 0
        @constraint(mdl, [k=1:N-1], SxInv*X(k+1) .== SxInv*(Ak[:,:,k]*X(k) + Bmk[:,:,k]*U(k) + wk[:,k]) + ν_ctrl[:,k])
    elseif params.disc == 1
        @constraint(mdl, [k=1:N-1], SxInv*X(k+1) .== SxInv*(Ak[:,:,k]*X(k) + Bmk[:,:,k]*U(k) + Bpk[:,:,k]*U(k+1) + wk[:,k]) + ν_ctrl[:,k])
    end

    # >> Convex State & Control Constraints <<
    # Constant altitude constraint
    @constraint(mdl, [k=1:N-1], r[3,k+1] == r[3,k])

    # Thrust bounds
    # Χ(k) = normalize(T_ref[:,k])
    @constraint(mdl, [k=1:N_ctrl], vcat(params.ρ_max, T[:,k]) in SecondOrderCone())
    # @constraint(mdl, [k=1:N_ctrl], dot(Χ(k),T[:,k]) >= params.ρ_min)

    # Attitude pointing constraint
    @constraint(mdl, [k=1:N_ctrl], vcat(dot(T[:,k],e_z)/cos(params.γ_p), T[:,k]) in SecondOrderCone())

    # Velocity upper bound
    @constraint(mdl, [k=1:N], vcat(params.v_max_L,v[1:2,k]) in SecondOrderCone())

    # Cage bounds
    @constraint(mdl, [k=1:N], r[1,k] >= params.x_arena_lims[1])
    @constraint(mdl, [k=1:N], r[1,k] <= params.x_arena_lims[2])
    @constraint(mdl, [k=1:N], r[2,k] >= params.y_arena_lims[1])
    @constraint(mdl, [k=1:N], r[2,k] <= params.y_arena_lims[2])
    @constraint(mdl, [k=1:N], r[3,k] >= params.z_arena_lims[1])
    @constraint(mdl, [k=1:N], r[3,k] <= params.z_arena_lims[2])

    # Time dilation
    for k=1:(N-1)
        if params.disc == 0
            Δt[k] = @expression(mdl, params.Δτ[k] * s[k])
        elseif params.disc == 1
            Δt[k] = @expression(mdl, (1/2) * params.Δτ[k] * (s[k] + s[k+1]))
        end
    end
    @constraint(mdl, sum(Δt) <= params.ToF_max)
    @constraint(mdl, [k=1:N_ctrl], params.s_min <= s[k] <= params.s_max)

    # Obstacle constraints
    for o = 1:params.n_obstacles
        H = params.H_obstacles[o]
        for k = 1:N
            Δr = r_ref[:,k] - params.p_obstacles[:,o]
            δr = r[:,k] - r_ref[:,k]
            ξ  = max(norm(H*Δr,2),1e-4)
            ζ  = transpose(H)*H*Δr / ξ
            @constraint(mdl, ξ + dot(ζ,δr) >= params.R_obstacles[o] + ν_obs[o,k])
        end
    end

    # >> Boundary conditions <<
    z0 = params.z0
    zf = params.zf_targs[:,j_targ]
    @constraint(mdl, X(1) .== z0)
    @constraint(mdl, X(N) .== zf)

    # Trust region constraints
    δX(k) = SxInv*(X(k) .- x_ref[:,k])
    δU(k) = SuInv*(U(k) .- u_ref[:,k])
    @constraint(mdl, [k=1:N_ctrl], δX(k)'*δX(k) + δU(k)'*δU(k) <= η[k])
    @constraint(mdl, δX(N)'*δX(N) <= η[N])

    # Cost function slack constraints
    @constraint(mdl, vcat(μ_obs, vec(ν_obs)) in MOI.NormOneCone(params.n_obstacles*N+1)) 
    @constraint(mdl, vcat(μ_ctrl, vec(ν_ctrl)) in MOI.NormOneCone(nx*(N-1)+1))
    @constraint(mdl, vcat(η_s, η) in SecondOrderCone())
    @constraint(mdl, μ_obs >= 0)
    @constraint(mdl, μ_ctrl >= 0)
    @constraint(mdl, η_s >= 0)

    # ..:: Solve the problem and save the solution ::..

    # >> Cost function <<
    J_opt  = sum(vec(T).^2)
    J_ptr  = η_s
    J_buff = μ_obs
    J_ctrl = μ_ctrl
    if params.n_obstacles == 0
        J_buff = 0
    end
    @objective(mdl, Min, J_opt + params.w_trust * J_ptr + params.w_buff * J_buff + params.w_ctrl * J_ctrl)

    optimize!(mdl)
    feas_status = JuMP.termination_status(mdl)
    if feas_status == MOI.OPTIMAL || feas_status == MOI.ALMOST_OPTIMAL
        solve_status = "Feasible"
    else
        solve_status = "Not Feasible"
        return (EmptySolution(), feas_status, false)
    end

    # Obtain optimized decision variables
    cost = objective_value(mdl)
    r = value.(r)
    v = value.(v)
    T = value.(T)
    μ_obs = value.(μ_obs)
    s = value.(s)
    x = vcat(r,v,ones(1,N))
    u = vcat(T,reshape(s,1,N_ctrl))
    Δt = value.(Δt)
    t = vcat(0,cumsum(Δt))

    # Package the solution
    sol = Solution(t,x,u,cost)

    # Obtain evaluation penalties
    μ_obs_pen = value.(μ_obs)
    μ_ctrl_pen = value.(μ_ctrl)
    η_pen = value.(η_s)

    # Determine convergence based on SCP penalties
    if (μ_ctrl_pen <= params.ϵ_ctrl) && (μ_obs_pen <= params.ϵ_buff) && (η_pen <= params.ϵ_trust)
        scp_sub_cvged = true
    else
        scp_sub_cvged = false
    end
    @printf("   SCP Iter: %2.i | Status: %s | μ_ctrl_pen = %.2e | μ_obs_pen = %.2e | η_pen = %.2e\n", scp_iter, solve_status, μ_ctrl_pen, μ_obs_pen, η_pen)
    flush(stdout)

    return (sol, feas_status, scp_sub_cvged)
end