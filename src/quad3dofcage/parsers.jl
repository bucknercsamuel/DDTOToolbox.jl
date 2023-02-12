Base.@ccallable function skyenet_ddtoscp_interface(
        num_targs::UInt32,
        r0_ptr::Ptr{Cdouble}, r0_size::Cint,
        v0_ptr::Ptr{Cdouble}, v0_size::Cint,
        rf_ptr::Ptr{Cdouble}, rf_size::Cint,
        vf_ptr::Ptr{Cdouble}, vf_size::Cint,
        K::UInt32,
        tf::CReal,
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
    params = Params()

    # >> Vehicle parameters <<
    params.ρ_min = params.mass * a_min
    params.ρ_max = params.mass * a_max

    # >> Constraint parameters <<
    params.γ_p = theta_max
    params.v_max_L = v_max

    # >> Dynamics <<
    params.Δt = tf / (K-1)

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
    params.r0 = r0
    params.v0 = v0
    params.rf_targs = zeros(3,num_targs)
    params.vf_targs = zeros(3,num_targs)
    for i = 1:3
        for j = 1:num_targs
            ind = j + MAX_TARGETS*(i-1)
            # ind = i + 3*(j-1)
            params.rf_targs[i,j] = rf[ind]
            params.vf_targs[i,j] = vf[ind]
        end
    end

    # >> Target conditions <<
    params.n_targs = num_targs
    params.N_targs = fill(K, num_targs)
    params.λ_targs = collect(1:num_targs)
    params.T_targs = collect(1:num_targs)
    params.ϵ_targs = fill(subopt_tol, num_targs)

    # >> SCP Params <<
    params.w_buff = w_buff
    params.w_trust = w_trust
    params.w_r0 = ri_relax
    params.w_rf = rf_relax
    params.sub_iters = scp_iters
    params.ϵ_cvg = eps_cvg

    # >> Other <<
    params.τ_max = tau_max
    method = "SCP"

    ## Call DDTO
    ~, DDTO_target_solutions = execute_ddto_solution(params, method)

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

function standalone_interface(scenario::String="default", method::String="SCP")
    ## Define the params params and scenario parameters
    params = Params()

    if scenario=="default"
        ~ # Do nothing
    elseif scenario=="toy1"
        scenario_toy1!(params)
    elseif scenario=="onr_demo"
        scenario_onr_demo!(params)
    end

    ## Call DDTO
    sols_optimal, DDTO_target_solutions = execute_ddto_solution(params, method)

    ## Plot solutions
    # Setup
    set_fonts()
    PyPlot.close("all")
    pygui(false)

    # Plot functionality
    plot_parametric_optimal_trajectories(params, sols_optimal)
    plot_parametric_ddto_trajectories(params, DDTO_target_solutions)
    plot_states(params, DDTO_target_solutions)
end

function execute_ddto_solution(params::Params, method::String)::Tuple{Vector{Solution},Vector{BranchSolution}}
    if method == "Baseline"
        @time begin
            @time begin
                # ..:: Solve for independently-optimal solutions to each target ::..
                (sols_optimal) = solve_optimal_tree(params)
                costs_optimal = CVector(zeros(params.n_targs))
                for k = 1:params.n_targs
                    costs_optimal[k] = sols_optimal[k].cost
                end
                println("\n Solve time for generating optimal solutions to each target:")
            end
    
            @time begin
                # ..:: Solve for DDTO branching solutions to ALL targets ::..
                sols_ddto = solve_ddto_tree(params, costs_optimal)
                println("\n Solve time for generating DDTO branch solutions to all targets:")
            end
            println("\n Solve time for the full DDTO solution stack:")
        end
        DDTO_target_solutions = extract_target_trajectories(params, sols_ddto)
    
    elseif method == "SCP"
        try
            @time begin
                @time begin
                    # ..:: Solve for independently-optimal solutions to each target ::..
                    (sols_optimal) = solve_optimal_tree(params)
                    costs_optimal = CVector(zeros(params.n_targs))
                    for k = 1:params.n_targs
                    costs_optimal[k] = sols_optimal[k].cost
                    end
                    println("\n Solve time for generating optimal solutions to each target:")
                end
        
                @time begin
                    # ..:: Solve for DDTO branching solutions to ALL targets ::..
                    sols_ddto = solve_ddtoscp_tree(params, costs_optimal, deepcopy(sols_optimal))
                    println("\n Solve time for generating DDTO branch solutions to all targets:")
                end
                println("\n Solve time for the full DDTO solution stack:")
            end
        catch
            sols_ddto = Vector{DDTOSolution}(undef, params.n_targs)
            for k = 1:(params.n_targs)
                sols_ddto[k] = EmptyDDTOSolution(params.n_targs-k+1)
                for j=1:(params.n_targs-k+1)
                    N_targ = params.N_targs[j]
                    N_targ_ctrl = N_targ - 1
                    sols_ddto[k].targ_sols[j].t        = zeros(N_targ)
                    sols_ddto[k].targ_sols[j].r        = zeros(3,N_targ)
                    sols_ddto[k].targ_sols[j].v        = zeros(3,N_targ)
                    sols_ddto[k].targ_sols[j].T        = zeros(3,N_targ_ctrl)
                    sols_ddto[k].targ_sols[j].Γ        = zeros(N_targ_ctrl)
                    sols_ddto[k].targ_sols[j].r0_relax = zeros(3)
                    sols_ddto[k].targ_sols[j].rf_relax = zeros(3)
                    sols_ddto[k].targ_sols[j].cost     = 0
                    sols_ddto[k].targ_sols[j].T_nrm    = zeros(N_targ_ctrl)
                    sols_ddto[k].targ_sols[j].γ        = zeros(N_targ_ctrl)
                end
            end
        end
        DDTO_target_solutions = extract_target_trajectories(params, sols_ddto)
        
    else
        println("Selection invalid.")
    end
    
    return sols_optimal, DDTO_target_solutions
end