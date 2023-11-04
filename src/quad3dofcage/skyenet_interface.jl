Base.@ccallable function skyenet_ddtoscp_interface(
        num_targs::UInt32,
        r0_ptr::Ptr{Cdouble}, r0_size::Cint,
        v0_ptr::Ptr{Cdouble}, v0_size::Cint,
        rf_ptr::Ptr{Cdouble}, rf_size::Cint,
        vf_ptr::Ptr{Cdouble}, vf_size::Cint,
        K::UInt32,
        tf::CReal, # not in use currently, TODO: implement fixed-final-time
        a_min::CReal,
        a_max::CReal,
        v_max::CReal,
        theta_max::CReal,
        ri_relax::CReal,
        rf_relax::CReal,
        subopt_tol::CReal,
        w_buff::CReal,
        w_trust::CReal,
        scp_iters::UInt32,
        tau_max::UInt32,
        eps_cvg::CReal,
        n::UInt32,
        MAX_HORIZON::UInt32,
        MAX_TARGETS::UInt32,
        MAX_OBS::UInt32,
        c_x_ptr::Ptr{Cdouble}, c_x_size::Cint,
        c_y_ptr::Ptr{Cdouble}, c_y_size::Cint,
        R_ptr::Ptr{Cdouble}, R_size::Cint,
        M0_ptr::Ptr{Cdouble}, M0_size::Cint,
        M1_ptr::Ptr{Cdouble}, M1_size::Cint,
        t_out_ptr::Ptr{Cdouble}, t_out_size::Cint,
        s_out_ptr::Ptr{Cdouble}, s_out_size::Cint,
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
    s_out = unsafe_wrap(Array, s_out_ptr, s_out_size, own=false)
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
    params.z0 = zeros(params.nx)
    params.z0[1:3] = r0
    params.z0[4:6] = v0
    params.z0[7] = 0
    params.zf_targs = zeros(params.nx,num_targs)
    for j = 1:num_targs
        for i = 1:3
            ind = j + MAX_TARGETS*(i-1)
            params.zf_targs[i,j] = rf[ind]
            params.zf_targs[i+3,j] = vf[ind]
        end
        params.zf_targs[7,j] = Inf
    end

    # >> Target conditions <<
    params.n_targs = num_targs
    params.λ_targs = collect(1:num_targs)
    params.T_targs = collect(1:num_targs)
    params.ϵ_targs = fill(subopt_tol, num_targs)

    # >> Discretization <<
    params.N = K
    params.τ = CVector(range(0, stop=1, length=params.N))
    params.Δτ = diff(params.τ)

    # >> SCP Params <<
    params.w_obj = 1
    params.w_ctrl = w_buff # keep same to minimize parameters
    params.w_buff = w_buff
    params.w_trust = w_trust
    params.scp_iters = scp_iters
    params.ϵ_ctrl = eps_cvg
    params.ϵ_buff = eps_cvg
    params.ϵ_trust = eps_cvg

    # >> Other <<
    params.τ_max = tau_max

    ## Call DDTO
    DDTO_target_solutions = solve_skyenet(params)

    # Write outputs to memory
    for k = 1:K
        t_out[k] = DDTO_target_solutions[1].sol.t[k]
    end
    for c = 1:3
        for k = 1:K
            for t = 1:num_targs
                ind = t + MAX_TARGETS*(k-1) + MAX_HORIZON*MAX_TARGETS*(c-1)
                # ind = c + 3*(k-1) + 3*MAX_HORIZON*(t-1)
                r_out[ind] = DDTO_target_solutions[t].sol.r[c,k]
                v_out[ind] = DDTO_target_solutions[t].sol.v[c,k]
                if k < K
                    a_out[ind] = DDTO_target_solutions[t].sol.T[c,k] / params.mass
                else
                    a_out[ind] = 0
                end
            end
        end
    end
    for c = 1:3
        r0_relax_out[c] = DDTO_target_solutions[1].sol.r0_relax[c]
        for t = 1:num_targs
            ind = t + MAX_TARGETS*(c-1)
            # ind = c + 3*(t-1)
            rf_relax_out[ind] = DDTO_target_solutions[t].sol.rf_relax[c]
        end
    end
end

function solve_skyenet(params::Quad3DoFCageParams)::Vector{BranchSolution}
    ddtoscp_solutions = Vector{DDTOSolution}(undef, params.n_targs)
    try
        @time begin
            @time begin
                # ..:: Solve for independently-optimal solutions to each target ::..
                scp_solutions = solve_tree_decoupled(params)
                scp_costs = CVector(zeros(params.n_targs))
                for k = 1:params.n_targs
                    scp_costs[k] = scp_solutions[k].cost
                end
                println("\n Solve time for generating optimal solutions to each target:")
            end
    
            @time begin
                # ..:: Solve for DDTO branching solutions to ALL targets ::..
                (feas_ddtoscp, ddtoscp_solutions) = solve_tree_ddto(deepcopy(params), scp_costs)
                println("\n Solve time for generating DDTO branch solutions to all targets:")
            end
            println("\n Solve time for the full DDTO solution stack:")
        end
    catch
        for k = 1:(params.n_targs)
            ddtoscp_solutions[k] = EmptyDDTOSolution(params.n_targs-k+1)
            N = params.N
            N_ctrl = N - 1
            for j=1:(params.n_targs-k+1)
                ddtoscp_solutions[k].targ_sols[j].t        = zeros(N)
                ddtoscp_solutions[k].targ_sols[j].r        = zeros(3,N)
                ddtoscp_solutions[k].targ_sols[j].v        = zeros(3,N)
                ddtoscp_solutions[k].targ_sols[j].T        = zeros(3,N_ctrl)
                ddtoscp_solutions[k].targ_sols[j].Γ        = zeros(N_ctrl)
                ddtoscp_solutions[k].targ_sols[j].r0_relax = zeros(3)
                ddtoscp_solutions[k].targ_sols[j].rf_relax = zeros(3)
                ddtoscp_solutions[k].targ_sols[j].cost     = 0
                ddtoscp_solutions[k].targ_sols[j].T_nrm    = zeros(N_ctrl)
                ddtoscp_solutions[k].targ_sols[j].γ        = zeros(N_ctrl)
            end
        end
    end
    
    # Convert DDTO solutions to branch solutions
    ddtoscp_branch_solutions,~ = extract_target_trajectories(params, ddtoscp_solutions; SCP=true)

    return ddtoscp_branch_solutions
end