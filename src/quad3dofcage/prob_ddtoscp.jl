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
            (feas_ddtoscp, ddtoscp_solutions) = solve_ddtoscp_tree(deepcopy(params), scp_costs)
            println("\n Solve time for generating DDTO branch solutions to all targets:")
        end
        println("\n Solve time for the full DDTO solution stack:")
    end

    # Convert DDTO solutions to branch solutions
    ddtoscp_bsolutions, defer_bsolutions = extract_target_trajectories(params, ddtoscp_solutions; SCP=true)

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

    # Dynamics functions
    dyn_lin_ = (t,x,u) -> dyn_lin(t,x,u,params)
    dyn_nl_  = (t,x,u) -> dyn_nl(t,x,u,params)

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

    # Extra variables
    @variable(mdl, ν_ctrl_stitch[1:nx,1:n])

    # >> Build the trunk <<
    # Take <any> reference and build it with first τ elements
    ref_traj = copy(reference_targ_trajs[1])
    ref_traj.x = ref_traj.x[:,1:τ]
    ref_traj.u = ref_traj.u[:,1:τ]

    # Build base SCP problem with trunk
    r_trunk,v_trunk,T_trunk,s_trunk,ν_trunk_ctrl,ν_trunk_obs,η_trunk,Δt_trunk,SxInv = joint_problem(mdl, τ, params, ref_traj)

    # >> Build each branch
    r_branch = Array{AffExpr}(undef, 3, N-τ, n)
    v_branch = Array{AffExpr}(undef, 3, N-τ, n)
    T_branch = Array{AffExpr}(undef, 3, N_ctrl-τ, n)
    s_branch = Array{AffExpr}(undef, N_ctrl-τ, n)
    ν_branch_ctrl = Array{VariableRef}(undef, nx, N-τ-1, n)
    ν_branch_obs = Array{VariableRef}(undef, params.n_obstacles, N-τ, n)
    η_branch = Array{VariableRef}(undef, N-τ, n)
    Δt_branch = Array{AffExpr}(undef, N-τ-1, n)
    for j = 1:n
        # Take jth reference and build it with last n-τ elements
        ref_traj = copy(reference_targ_trajs[j])
        ref_traj.x = ref_traj.x[:,τ+1:end]
        ref_traj.u = ref_traj.u[:,τ+1:end]

        # Build branch SCP problem
        r_branch[:,:,j],v_branch[:,:,j],T_branch[:,:,j],s_branch[:,j],ν_branch_ctrl[:,:,j],ν_branch_obs[:,:,j],η_branch[:,j],Δt_branch[:,j],_ = joint_problem(mdl, N-τ, params, ref_traj)
    end

    # ..:: Segment Stitching Constraints ::..
    # >> Convenience functions <<
    X_trunk(k) = [r_trunk[:,k]; v_trunk[:,k]; 1]           # Augmented state (to bring in affine term)
    U_trunk(k) = [T_trunk[:,k]; s_trunk[k]]                # Augmented control (with time dilation term)
    X_branch(k,j) = [r_branch[:,k,j]; v_branch[:,k,j]; 1]  # Augmented state (to bring in affine term)
    U_branch(k,j) = [T_branch[:,k,j]; s_branch[k,j]]       # Augmented control (with time dilation term)

    # >> Trunk <<
    # Initial condition
    @constraint(mdl, X_trunk(1) .== params.z0)

    # >> Branches <<
    for j = 1:n
        # Apply suboptimality constraint
        @constraint(mdl, sum(Δt_trunk) + sum(Δt_branch[:,j]) + cost_dd <= (1 + params.ϵ_targs[j]) * ref_costs[j])

        # Apply dynamics stitching
        ref_traj_stitch = copy(reference_targ_trajs[j])
        ref_traj_stitch.x = ref_traj_stitch.x[:,τ:τ+1]
        ref_traj_stitch.u = ref_traj_stitch.u[:,τ:τ+1]
        dyn_lin_ = (t,x,u) -> dyn_lin(t,x,u,params)
        dyn_nl_  = (t,x,u) -> dyn_nl(t,x,u,params)
        Ak,Bmk,Bpk,wk,_ = c2d_nonlinear(ref_traj_stitch,dyn_nl_,dyn_lin_,params.disc)
        if params.disc == 0
            @constraint(mdl, SxInv*X_branch(1,j) .== SxInv*(Ak[:,:,1]*X_trunk(τ) + Bmk[:,:,1]*U_trunk(τ) + wk[:,1]) + ν_ctrl_stitch[:,j])
        elseif params.disc == 1
            @constraint(mdl, SxInv*X_branch(1,j) .== SxInv*(Ak[:,:,1]*X_trunk(τ) + Bmk[:,:,1]*U_trunk(τ) + Bpk[:,:,1]*U_branch(1,j) + wk[:,1]) + ν_ctrl_stitch[:,j])
        end

        # # Apply boundary conditions
        # @constraint(mdl, X_branch(1,j) .== X_trunk(τ))
        # @constraint(mdl, X_branch(N-τ,j) .== params.zf_targs[:,j])

        # # Maintain continuity in control if using FOH
        # if params.disc == 1
        #     @constraint(mdl, U_branch(1,j) .== U_trunk(τ))
        # end
    end

    # Maintain continuity for FOH discretization if we have already deferred by some amount (usually after first branch iteration)
    # (acts as a initial condition on control constrained to the reference initial condition)
    if params.disc == 1 && cost_dd > 0
        u_ref = reference_targ_trajs[1].u # any traj works
        @constraint(mdl, U_trunk(1) .== u_ref[:,1])
    end

    # ..:: PTR Constraints ::..
    ν_ctrl = [vec(ν_trunk_ctrl); vec(ν_branch_ctrl); vec(ν_ctrl_stitch)]
    ν_obs = [vec(ν_trunk_obs); vec(ν_branch_obs)]
    @variable(mdl, μ_obs)
    @variable(mdl, μ_ctrl)
    @variable(mdl, η_s)
    @constraint(mdl, vcat(μ_ctrl, ν_ctrl) in MOI.NormOneCone(length(ν_ctrl)+1))
    @constraint(mdl, vcat(μ_obs, ν_obs) in MOI.NormOneCone(length(ν_obs)+1))
    @constraint(mdl, vcat(η_s, [vec(η_trunk); vec(η_branch)]) in SecondOrderCone())
    @constraint(mdl, μ_obs >= 0)
    @constraint(mdl, μ_ctrl >= 0)
    @constraint(mdl, η_s >= 0)

    # >> Cost function <<
    J_opt  = -sum(Δt_trunk)
    J_ptr  = η_s
    J_buff = μ_obs
    J_ctrl = μ_ctrl
    if params.n_obstacles == 0
        J_buff = 0
    end
    @objective(mdl, Min, 
        params.w_obj * J_opt 
      + params.w_trust * J_ptr 
      + params.w_buff * J_buff 
      + params.w_ctrl * J_ctrl)
        
    optimize!(mdl)
    feas_status = JuMP.termination_status(mdl)

    # ..:: Extract the solution ::..
    r_trunk = value.(r_trunk)
    v_trunk = value.(v_trunk)
    T_trunk = value.(T_trunk)
    s_trunk = value.(s_trunk)
    r_branch = value.(r_branch)
    v_branch = value.(v_branch)
    T_branch = value.(T_branch)
    s_branch = value.(s_branch)

    x = zeros(nx, N, n)
    u = zeros(nu, N, n)

    for j = 1:n
        x[:,1:τ,j] = vcat(r_trunk, v_trunk, ones(1,τ))
        u[:,1:τ,j] = vcat(T_trunk, reshape(s_trunk,1,τ))
        x[:,τ+1:end,j] = vcat(r_branch[:,:,j], v_branch[:,:,j], ones(1,N-τ))
        u[:,τ+1:end,j] = vcat(T_branch[:,:,j], reshape(s_branch[:,j],1,N_ctrl-τ))
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
    Δt_trunk = value.(Δt_trunk)
    Δt_branch = value.(Δt_branch)
    costs_sol = [sum(Δt_trunk) + sum(Δt_branch[:,j]) for j = 1:n]
    cost_dd = sum(Δt_trunk)

    # ..:: Package the DDTO Solution ::..
    ddto_solution = EmptyDDTOSolution(n)
    for j = 1:n
        ddto_solution.targ_sols[j].t = params.τ
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

    # >> Convenience functions <<
    X = (k) -> [r[:,k]; v[:,k]; 1] # Augmented state (to bring in affine term)
    U = (k) -> [T[:,k]; s[k]] # Augmented control (with time dilation term)

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

    # >> Joint problem <<
    r,v,T,s,ν_ctrl,ν_obs,η,Δt,_ = joint_problem(mdl, N, params, ref_traj)

    # >> Boundary conditions <<
    z0 = params.z0
    zf = params.zf_targs[:,j_targ]
    @constraint(mdl, X(1) .== z0)
    @constraint(mdl, X(N) .== zf)

    # Cost function slack constraints
    @variable(mdl, μ_obs)
    @variable(mdl, μ_ctrl)
    @variable(mdl, η_s)
    @constraint(mdl, vcat(μ_obs, vec(ν_obs)) in MOI.NormOneCone(params.n_obstacles*N+1)) 
    @constraint(mdl, vcat(μ_ctrl, vec(ν_ctrl)) in MOI.NormOneCone(nx*(N-1)+1))
    @constraint(mdl, vcat(η_s, η) in SecondOrderCone())
    @constraint(mdl, μ_obs >= 0)
    @constraint(mdl, μ_ctrl >= 0)
    @constraint(mdl, η_s >= 0)

    # ..:: Solve the problem and save the solution ::..

    # >> Cost function <<
    J_opt  = sum(Δt)
    J_ptr  = η_s
    J_buff = μ_obs
    J_ctrl = μ_ctrl
    if params.n_obstacles == 0
        J_buff = 0
    end
    @objective(mdl, Min, 
        params.w_obj * J_opt 
      + params.w_trust * J_ptr 
      + params.w_buff * J_buff 
      + params.w_ctrl * J_ctrl)

    optimize!(mdl)
    feas_status = JuMP.termination_status(mdl)
    if feas_status == MOI.OPTIMAL || feas_status == MOI.ALMOST_OPTIMAL
        solve_status = "Feasible"
    else
        solve_status = "Not Feasible"
        return (EmptySolution(), feas_status, false)
    end

    # Obtain optimized decision variables
    cost = value.(J_opt)
    r = value.(r)
    v = value.(v)
    T = value.(T)
    μ_obs = value.(μ_obs)
    s = value.(s)
    x = vcat(r,v,ones(1,N))
    u = vcat(T,reshape(s,1,N_ctrl))

    # Package the solution
    sol = Solution(params.τ,x,u,cost)

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

function joint_problem(mdl::JuMP.Model, N::Int, params::Quad3DoFCageParams, ref_traj::Solution)
    """
    NOTE: This function contains the joint problem formulation shared by `solve_scp_target` and `solve_feasible_ddtoscp`
    """
    nx = params.n+1
    nu = params.m+1
    if params.disc == 0
        N_ctrl = N-1
    elseif params.disc == 1
        N_ctrl = N
    end

    # >> Scaled Optimization variables <<
    r_s = @variable(mdl, [1:3,1:N])
    v_s = @variable(mdl, [1:3,1:N])
    T_s = @variable(mdl, [1:3,1:N_ctrl])
    s_s = @variable(mdl, [1:N_ctrl])
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
    ν_obs = @variable(mdl, [1:params.n_obstacles,1:N])
    ν_ctrl = @variable(mdl, [1:nx,1:(N-1)])
    η = @variable(mdl, [1:N])

    # >> Convenience functions <<
    X = (k) -> [r[:,k]; v[:,k]; 1] # Augmented state (to bring in affine term)
    U = (k) -> [T[:,k]; s[k]] # Augmented control (with time dilation term)

    # Extract reference trajectory elements
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
    Χ(k) = normalize(T_ref[:,k])
    @constraint(mdl, [k=1:N_ctrl], vcat(params.ρ_max, T[:,k]) in SecondOrderCone())
    # @constraint(mdl, [k=1:N_ctrl], dot(Χ(k),T[:,k]) >= params.ρ_min)

    # Attitude pointing constraint
    @constraint(mdl, [k=1:N_ctrl], vcat(dot(T[:,k],e_z)/cos(params.γ_p), T[:,k]) in SecondOrderCone())

    # Velocity upper bound
    @constraint(mdl, [k=1:N], vcat(params.v_max_L,v[1:2,k]) in SecondOrderCone())

    # Cage bounds
    # @constraint(mdl, [k=1:N], r[1,k] >= params.x_arena_lims[1])
    # @constraint(mdl, [k=1:N], r[1,k] <= params.x_arena_lims[2])
    # @constraint(mdl, [k=1:N], r[2,k] >= params.y_arena_lims[1])
    # @constraint(mdl, [k=1:N], r[2,k] <= params.y_arena_lims[2])
    # @constraint(mdl, [k=1:N], r[3,k] >= params.z_arena_lims[1])
    # @constraint(mdl, [k=1:N], r[3,k] <= params.z_arena_lims[2])

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

    # Trust region constraints
    δX(k) = SxInv*(X(k) .- x_ref[:,k])
    δU(k) = SuInv*(U(k) .- u_ref[:,k])
    @constraint(mdl, [k=1:N_ctrl], δX(k)'*δX(k) + δU(k)'*δU(k) <= η[k])
    @constraint(mdl, δX(N)'*δX(N) <= η[N])

    # Return created optimization variables
    return r,v,T,s,ν_ctrl,ν_obs,η,Δt,SxInv
end