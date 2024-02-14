Base.@ccallable function skyenet_ddtoscp_interface(
        num_targs::UInt32,
        r0_ptr::Ptr{Cdouble}, r0_size::Cint,
        v0_ptr::Ptr{Cdouble}, v0_size::Cint,
        rf_ptr::Ptr{Cdouble}, rf_size::Cint,
        vf_ptr::Ptr{Cdouble}, vf_size::Cint,
        K::UInt32,
        tf::CReal, # not in use currently, TODO: implement fixed-final-time
        tf_max::CReal,
        dt_min::CReal,
        dt_max::CReal,
        a_min::CReal,
        a_max::CReal,
        v_max::CReal,
        theta_max::CReal,
        ri_relax::CReal,
        rf_relax::CReal,
        subopt_tol::CReal,
        w_obj::CReal,
        w_trust::CReal,
        w_buff::CReal,
        scp_iters::UInt32,
        tau_max::UInt32,
        eps_cvg::CReal,
        n::UInt32,
        interp_ref::Bool,
        MAX_HORIZON::UInt32,
        MAX_TARGETS::UInt32,
        MAX_OBS::UInt32,
        c_x_ptr::Ptr{Cdouble}, c_x_size::Cint,
        c_y_ptr::Ptr{Cdouble}, c_y_size::Cint,
        R_ptr::Ptr{Cdouble}, R_size::Cint,
        M0_ptr::Ptr{Cdouble}, M0_size::Cint,
        M1_ptr::Ptr{Cdouble}, M1_size::Cint,
        t_out_ptr::Ptr{Cdouble}, t_out_size::Cint,
        r_out_ptr::Ptr{Cdouble}, r_out_size::Cint,
        v_out_ptr::Ptr{Cdouble}, v_out_size::Cint,
        a_out_ptr::Ptr{Cdouble}, a_out_size::Cint,
        r0_relax_out_ptr::Ptr{Cdouble}, r0_relax_out_size::Cint,
        rf_relax_out_ptr::Ptr{Cdouble}, rf_relax_out_size::Cint
    )::Cvoid

    # Set up logging
    LOGGING = false
    if LOGGING
        log_path = "/Users/samuelbuckner/Documents/ACL/Code/ONR_DDTO_Demo_2022/ddtoscp_debugging_log.txt" # Currently hardcoded path for debugging purposes only!
        io = open(log_path, "a")
        write(io, "\n\n ======= RUNNING PARSERS.JL =======\n\n")
        logger = SimpleLogger(io)
    end

    ## Unwrap input array pointer/sizes into arrays that Julia can use
    r0 = unsafe_wrap(Array, r0_ptr, r0_size, own=false)
    v0 = unsafe_wrap(Array, v0_ptr, v0_size, own=false)
    rf = unsafe_wrap(Array, rf_ptr, rf_size, own=false)
    vf = unsafe_wrap(Array, vf_ptr, vf_size, own=false)
    c_x = unsafe_wrap(Array, c_x_ptr, c_x_size, own=false)
    c_y = unsafe_wrap(Array, c_y_ptr, c_y_size, own=false)
    R = unsafe_wrap(Array, R_ptr, R_size, own=false)
    M0 = unsafe_wrap(Array, M0_ptr, M0_size, own=false)
    M1 = unsafe_wrap(Array, M1_ptr, M1_size, own=false)

    ## Unwrap output array pointers/sizes into arrays that Julia can use (will be overwritten at end)
    t_out = unsafe_wrap(Array, t_out_ptr, t_out_size, own=false)
    r_out = unsafe_wrap(Array, r_out_ptr, r_out_size, own=false)
    v_out = unsafe_wrap(Array, v_out_ptr, v_out_size, own=false)
    a_out = unsafe_wrap(Array, a_out_ptr, a_out_size, own=false)
    r0_relax_out = unsafe_wrap(Array, r0_relax_out_ptr, r0_relax_out_size, own=false)
    rf_relax_out = unsafe_wrap(Array, rf_relax_out_ptr, rf_relax_out_size, own=false)

    ## Define the base params and scenario params
    params = Quad3DoFCageParams()

    # >> Vehicle parameters <<
    params.ρ_min = params.mass * a_min
    params.ρ_max = params.mass * a_max

    # >> Constraint parameters <<
    params.γ_p = theta_max
    params.v_max_L = v_max

    # >> Obstacle parameters <<
    params.n_obstacles = n
    params.R_obstacles = R[1:n]
    params.p_obstacles = vcat(
        reshape(c_x[1:n],1,params.n_obstacles),
        reshape(c_y[1:n],1,params.n_obstacles),
        zeros(1,n)
    )
    params.H_obstacles = repeat([zeros(3,3)],n)
    for j = 1:n
        H = zeros(3,3)
        for i = 1:2
            ind = j + MAX_OBS*(i-1)
            H[i,1] = M0[ind]
            H[i,2] = M1[ind]
        end
        params.H_obstacles[j] = H
    end

    # >> Boundary conditions <<
    params.a.z0 = zeros(params.a.nx)
    params.a.z0[1:3] = r0
    params.a.z0[4:6] = v0
    params.a.z0[7] = 0
    params.a.zf_targs = zeros(params.a.nx,num_targs)
    for j = 1:num_targs
        for i = 1:3
            ind = j + MAX_TARGETS*(i-1)
            params.a.zf_targs[i,j] = rf[ind]
            params.a.zf_targs[i+3,j] = vf[ind]
        end
        params.a.zf_targs[7,j] = Inf
    end

    # >> Target conditions <<
    params.a.n_targs = num_targs
    params.a.λ_targs = collect(1:num_targs)
    params.a.T_targs = collect(1:num_targs)
    params.a.α_targs = ones(num_targs)
    params.a.ϵ_targs = fill(subopt_tol, num_targs)

    # >> Time dilation & discretization <<
    params.a.N = K
    params.a.Δt_min = dt_min
    params.a.Δt_max = dt_max
    params.a.ToF_max = tf_max

    # >> SCP Params <<
    params.w_obj = w_obj
    params.a.w_trust = w_trust
    params.a.w_ctrl = w_buff # keep same to minimize parameters
    params.a.w_buff = w_buff
    params.a.scp_iters = scp_iters
    params.a.ϵ_ctrl = eps_cvg
    params.a.ϵ_buff = eps_cvg
    params.a.ϵ_trust = eps_cvg

    # >> Reference trajectory extraction <<
    # Obtain from current values of {t_out, r_out, v_out, a_out}
    t_bar = zeros(K,params.a.n_targs)
    r_bar = zeros(3,K,params.a.n_targs)
    v_bar = zeros(3,K,params.a.n_targs)
    a_bar = zeros(3,K,params.a.n_targs)
    for c = 1:3
        for k = 1:K
            for t = 1:num_targs
                ind = t + MAX_TARGETS*(k-1) + MAX_HORIZON*MAX_TARGETS*(c-1)
                r_bar[c,k,t] = r_out[ind]
                v_bar[c,k,t] = v_out[ind]
                a_bar[c,k,t] = a_out[ind]
                if c == 1
                    t_bar[k,t] = t_out[ind]
                end
            end
        end
    end

    ref_trajs = generate_initial_guess_ddtoscp(params)
    if !interp_ref
        for j = 1:params.a.n_targs
            s_bar = wall_clock_time_to_time_dilation_control(t_bar[:,j], ref_trajs.targs[j].t, params.a.disc) #TODO: fix this
            Δt_bar = diff(t_bar[:,j])
            ∫T_bar = CVector(cumsum([params.a.z0[7],[norm(Δt_bar[k]*a_bar[:,k,j]*params.mass) for k=1:length(s_bar)-1]...]))
            ref_trajs.targs[j] = EmptySolution()
            ref_trajs.targs[j].t = t_bar[:,j]
            ref_trajs.targs[j].x = vcat(r_bar[:,:,j], v_bar[:,:,j], reshape(∫T_bar, 1, length(∫T_bar)))
            ref_trajs.targs[j].u = vcat(a_bar[:,:,j] * params.mass, reshape(s_bar, 1, length(s_bar)))
        end
    end

    ## Call DDTO
    DDTO_target_solutions = solve_skyenet(params, ref_trajs)

    # Write outputs to memory
    for c = 1:3
        for k = 1:K
            for t = 1:num_targs
                ind = t + MAX_TARGETS*(k-1) + MAX_HORIZON*MAX_TARGETS*(c-1)
                r_out[ind] = DDTO_target_solutions.targs[t].r[c,k]
                v_out[ind] = DDTO_target_solutions.targs[t].v[c,k]
                a_out[ind] = DDTO_target_solutions.targs[t].T[c,k] / params.mass
                if c == 1
                    t_out[ind] = DDTO_target_solutions.targs[t].t[k]
                end
            end
        end
    end
    # for c = 1:3
    #     r0_relax_out[c] = DDTO_target_solutions[1].sol.r0_relax[c]
    #     for t = 1:num_targs
    #         ind = t + MAX_TARGETS*(c-1)
    #         # ind = c + 3*(t-1)
    #         rf_relax_out[ind] = DDTO_target_solutions[t].sol.rf_relax[c]
    #     end
    # end
end

function solve_skyenet(params::Quad3DoFCageParams, ref_trajs::DDTOSolution)::Quad3DoFCageDDTOSolution
    ddtoscp_solutions_unprocessed = EmptyDDTOSolution(params.a.n_targs)
    # try
        @time begin
            @time begin
                # ..:: Solve for independently-optimal solutions to each target ::..
                scp_solutions = solve_tree_decoupled(params; single_iter=false, ref_trajs=ref_trajs)
                scp_costs = CVector(zeros(params.a.n_targs))
                for k = 1:params.a.n_targs
                    scp_costs[k] = scp_solutions.targs[k].cost
                end
                println("\n Solve time for generating optimal solutions to each target:")
            end
            
            @time begin
                # ..:: Solve for DDTO branching solutions to ALL targets ::..
                (feas_ddtoscp, ddtoscp_solutions_unprocessed) = solve_tree_ddto(params, scp_costs; single_iter=false, ref_trajs=ref_trajs)
                println("\n Solve time for generating DDTO branch solutions to all targets:")
            end
            println("\n Solve time for the full DDTO solution stack:")
        end
    # catch e
    #     println("!! Error thrown during DDTO solve:")
    #     try
    #         println(e)
    #     catch
    #         println("No error message available...")
    #     end
    #     ddtoscp_solutions_unprocessed = generate_initial_guess_ddtoscp(params)
    # end
    
    # Process and output solutions
    ddtoscp_solutions = process_solutions(ddtoscp_solutions_unprocessed, params)
    return ddtoscp_solutions
end