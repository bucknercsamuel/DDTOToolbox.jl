function solve_scp_target(params::Params, ref_traj::Solution, N::Int, j_targ::Int, scp_iter::Int)::Tuple{Solution, MOI.TerminationStatusCode, Bool}

    # ..:: Discrete time interval ::..
    if params.disc != 0 && params.disc != 1
        error("Please select a valid discretization hold order.")
    end
    
    if !params.free_final_time
        Δt = params.Δt
        tf = Δt * (N-1)
        t  = CVector(range(0, stop=tf, length=N))
        n_ = params.n
        m_ = params.m
    else
        n_ = params.n+1
        m_ = params.m+1
    end
    if params.disc == 0
        N_ctrl = N-1
        if !params.free_final_time
            A,B,p = c2d_LTI_affine_zoh(params.A_c, params.B_c, params.p_c, Δt)
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
        JuMP.set_optimizer_attribute(mdl, "LOG",  0) # disable debugging
        JuMP.set_optimizer_attribute(mdl, "MAX_NUM_WARNINGS", 0) # disable warnings
    else
        error("SOLVER is invalid, please select either ECOS or MOSEK")
    end

    # >> Optimization variables <<
    @variable(mdl, r[1:3,1:N])
    @variable(mdl, v[1:3,1:N])
    @variable(mdl, T[1:3,1:N_ctrl])
    @variable(mdl, Γ[1:N_ctrl])
    if params.free_final_time
        @variable(mdl, s[1:N_ctrl])
        Δt = Array{AffExpr}(undef,N-1)
    end

    # >> SCP variables <<
    # Virtual buffers
    @variable(mdl, ν_obs[1:params.n_obstacles,1:N])
    @variable(mdl, μ_obs[1:params.n_obstacles,1:N])

    # Virtual control
    @variable(mdl, ν_ctrl[1:n_,1:(N-1)])
    @variable(mdl, μ_ctrl[1:n_,1:(N-1)])

    # Trust region variables
    @variable(mdl, η_x[1:N])
    @variable(mdl, η_u[1:N_ctrl])

    # Slack variables for objective function
    @variable(mdl, μ_obs_s)
    @variable(mdl, μ_ctrl_s)
    @variable(mdl, η_x_s)
    @variable(mdl, η_u_s)

    # >> Convenience functions <<
    if !params.free_final_time
        X = (k) -> [r[:,k]; v[:,k]] # State at time index k
        U = (k) -> T[:,k]   # Input at time index k
    else
        X = (k) -> [r[:,k]; v[:,k]; 1] # Augmented state (to bring in affine term)
        U = (k) -> [T[:,k]; s[k]] # Augmented control (with time dilation term)
    end

    # ..:: Constraints ::..

    # >> Convex State & Control Constraints <<

    # >> Dynamics (convex if fixed-final-time)
    if !params.free_final_time
        if params.disc == 0
            @constraint(mdl, [k=1:N-1], X(k+1) .== A*X(k) + B*U(k) + p)
        elseif params.disc == 1
            @constraint(mdl, [k=1:N-1], X(k+1) .== A*X(k) + Bm*U(k) + Bp*U(k+1) + p)
        end
    end

    # >> Constant altitude constraint <<
    @constraint(mdl, [k=1:N-1], r[3,k+1] == r[3,k])

    # >> Thrust bounds <<
    # @constraint(mdl, [k=1:N_ctrl], Γ[k] >= params.ρ_min)
    # @constraint(mdl, [k=1:N_ctrl], Γ[k] <= params.ρ_max)
    # @constraint(mdl, [k=1:N_ctrl], vcat(Γ[k], T[:,k]) in MOI.SecondOrderCone(4))
    @constraint(mdl, [k=1:N_ctrl], vcat(params.ρ_max, T[:,k]) in MOI.SecondOrderCone(4))
    @constraint(mdl, [k=1:N_ctrl], vcat(Γ[k], T[:,k]) in MOI.SecondOrderCone(4))

    # >> Attitude pointing constraint <<
    # @constraint(mdl, [k=1:N_ctrl], dot(T[:,k],e_z) >= norm(T[:,k])*cos(params.γ_p))
    @constraint(mdl, [k=1:N_ctrl], vcat(dot(T[:,k],e_z)/cos(params.γ_p), T[:,k]) in MOI.SecondOrderCone(4))

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

    # Time dilation
    if params.free_final_time
        for k=1:(N-1)
            if params.disc == 0
                Δt[k] = @expression(mdl, params.Δτ[k] * s[k])
            elseif params.disc == 1
                Δt[k] = @expression(mdl, (1/2) * params.Δτ[k] * (s[k] + s[k+1]))
            end
        end
        @constraint(mdl, sum(Δt) <= params.ToF_max)
        @constraint(mdl, [k=1:(N-1)], params.Δt_min <= Δt[k] <= params.Δt_max)
        @constraint(mdl, [k=1:N_ctrl], params.s_min <= s[k] <= params.s_max)
    end

    # >> SCP constraints <<
    # Extract reference trajectory for target j
    t_ref = ref_traj.t
    x_ref = ref_traj.x
    u_ref = ref_traj.u
    r_ref = x_ref[1:3,:]
    v_ref = x_ref[4:6,:]
    T_ref = u_ref[1:3,:]

    # Dynamics (free-final-time)
    if params.free_final_time
        dyn_lin_ = (t,x,u) -> dyn_lin(t,x,u,params)
        dyn_nl_  = (t,x,u) -> dyn_nl(t,x,u,params)

        # Obtain approximate LTV discrete-time dynamics
        if params.disc == 0
            Ak,Bk,_,wk,_ = c2d_nonlinear(ref_traj,dyn_nl_,dyn_lin_,params.disc)
        elseif params.disc == 1
            Ak,Bmk,Bpk,wk,_ = c2d_nonlinear(ref_traj,dyn_nl_,dyn_lin_,params.disc)
        end

        # Apply constraints
        if params.disc == 0
            @constraint(mdl, [k=1:N-1], X(k+1) .== Ak[:,:,k]*X(k) + Bk[:,:,k]*U(k) + wk[:,k] + ν_ctrl[:,k])
        elseif params.disc == 1
            @constraint(mdl, [k=1:N-1], X(k+1) .== Ak[:,:,k]*X(k) + Bmk[:,:,k]*U(k) + Bpk[:,:,k]*U(k+1) + wk[:,k] + ν_ctrl[:,k])
        end
        @constraint(mdl, [k=1:N-1,j=1:n_], vcat(μ_ctrl[j,k], ν_ctrl[j,k]) in MOI.NormOneCone(2))
    end

    # Linearization constraints
    for o = 1:params.n_obstacles
        H = params.H_obstacles[o]
        for k = 1:N
            Δr = r_ref[:,k] - params.p_obstacles[:,o]
            δr = r[:,k] - r_ref[:,k]
            ξ  = max(norm(H*Δr,2),1e-4)
            ζ  = transpose(H)*H*Δr / ξ
            @constraint(mdl, ξ + dot(ζ,δr) >= params.R_obstacles[o] + ν_obs[o,k])
            @constraint(mdl, vcat(μ_obs[o,k], ν_obs[o,k]) in MOI.NormOneCone(2))
        end
    end

    # Trust region constraints
    @constraint(mdl, [k=1:N],      vcat(η_x[k], X(k) - x_ref[:,k]) in MOI.SecondOrderCone(n_+1))
    @constraint(mdl, [k=1:N_ctrl], vcat(η_u[k], U(k) - u_ref[:,k]) in MOI.SecondOrderCone(m_+1))

    # Cost function slack constraints
    @constraint(mdl, vcat(μ_obs_s, vec(μ_obs)) in MOI.SecondOrderCone(params.n_obstacles*N+1))
    @constraint(mdl, vcat(μ_ctrl_s, vec(μ_ctrl)) in MOI.SecondOrderCone((n_)*(N-1)+1))
    @constraint(mdl, vcat(η_x_s, η_x) in MOI.SecondOrderCone(N+1))
    @constraint(mdl, vcat(η_u_s, η_u) in MOI.SecondOrderCone(N_ctrl+1))

    # >> Boundary conditions <<
    z0 = params.z0
    zf = params.zf_targs[:,j_targ]
    @constraint(mdl, X(1) .== z0)
    @constraint(mdl, X(N) .== zf)

    # ..:: Solve the problem and save the solution ::..

    # >> Cost function <<
    J_opt  = sum(Γ)
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
    if feas_status != MOI.OPTIMAL
        return (EmptySolution(), feas_status, false)
    end

    # Obtain optimized decision variables
    cost = objective_value(mdl)
    r = value.(r)
    v = value.(v)
    T = value.(T)
    Γ = value.(Γ)
    μ_obs = value.(μ_obs)
    η_x = value.(η_x)
    η_u = value.(η_u)
    η = [η_x;η_u]
    if params.free_final_time
        s = value.(s)
        x = vcat(r,v,ones(1,N))
        u = vcat(T,reshape(s,1,N_ctrl))
    else
        x = vcat(r,v)
        u = vcat(T)
    end

    # Obtain physical time "t" if using free-final-time formulation
    if params.free_final_time
        Δt = value.(Δt)
        t = vcat(0,cumsum(Δt))
    end
    
    # Package the solution
    sol = Solution(t,x,u,cost)

     # ..:: Determine if PTR subproblem has converged ::..
     if feas_status == MOI.OPTIMAL
        solve_status = "Feasible"
    else
        solve_status = "Not Feasible"
    end

    # Obtain evaluation penalty for virtual buffers
    μ_obs_max_nodes = []
    for k = 1:N
        append!(μ_obs_max_nodes, max(μ_obs[:,k]...,0))
    end
    μ_obs_pen = sum(μ_obs_max_nodes)

    # Obtain evaluation penalty for virtual control
    μ_ctrl_pen = value.(μ_ctrl_s)

    # Obtain evaluation penalty for trust region
    η_pen = norm(η,2)

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