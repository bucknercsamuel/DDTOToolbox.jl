#=
C-callable FFI bridge exposing DDTO-SCP for the Skyenet / external C++ stack:
unwraps pointer buffers into a cage-scenario problem, solves, and writes
trajectories back into preallocated output arrays.
=#

"""
    skyenet_ddtoscp_interface(...) -> Cvoid

`@ccallable` entry point for external callers. Packs C pointer buffers into a
`Quad3DoFCageParams` problem, runs `solve`, and writes nodal /
simulated position-velocity-acceleration trajectories into the provided output
arrays. Sets `ddto_converged` on success.

# Arguments
**Problem sizing**
- `num_targs`: number of landing targets.
- `K`: number of discretization nodes per target trajectory.
- `n`: number of cylindrical obstacles.
- `scp_iters`: maximum SCP/PTR iterations.
- `sim_steps`: RK4 simulation steps per inter-node interval.

**Initial / terminal conditions (inputs)**
- `r0_ptr`, `r0_size`: initial position `[m]` (length 3).
- `v0_ptr`, `v0_size`: initial velocity `[m/s]` (length 3).
- `rf_ptr`, `rf_size`: terminal positions packed as `rf[t + MAX_TARGETS*(i-1)]` for
  component `i ∈ {1,2,3}` and target `t`.
- `vf_ptr`, `vf_size`: terminal velocities with the same packing as `rf`.

**Obstacle geometry (inputs)**
- `c_x_ptr`, `c_x_size`: obstacle center x-coordinates (length `n`).
- `c_y_ptr`, `c_y_size`: obstacle center y-coordinates (length `n`).
- `R_ptr`, `R_size`: obstacle radii (length `n`).
- `M0_ptr`, `M0_size`: ellipse shape matrix entries `M0[j + MAX_OBS*(i-1)]`.
- `M1_ptr`, `M1_size`: ellipse shape matrix entries `M1[j + MAX_OBS*(i-1)]`.

**Vehicle / constraint limits (inputs)**
- `tf`: reserved fixed-final-time parameter (currently unused).
- `tf_max`: maximum time of flight `[s]`.
- `dt_min`, `dt_max`: minimum and maximum segment durations `[s]`.
- `a_min`, `a_max`: minimum and maximum specific thrust/acceleration `[m/s²]`.
- `v_max`: maximum lateral speed `[m/s]`.
- `theta_max`: maximum pointing angle `[rad]`.
- `ri_relax`, `rf_relax`: reserved boundary-relaxation parameters (outputs allocated).
- `subopt_tol`: per-target DDTO suboptimality tolerance fraction.
- `w_obj`: single-target objective weight.
- `w_trust`, `w_buff`: trust-region and virtual-buffer objective weights.
- `eps_cvg`: SCP convergence tolerances for control, buffer, and trust penalties.
- `eps_ctcs`: CTCS violation tolerance.

**Solver flags (inputs)**
- `ctcs_enabled`: enable CTCS constraint augmentation.
- `autogen_init_guess`: if `true`, ignore warmstart buffers and auto-generate guesses.
- `ddto_init_guess`: select DDTO warmstart behavior passed to the solver stack.

**Buffer layout constants (inputs)**
- `MAX_HORIZON`: allocation count for trajectory discretization nodes in output arrays.
- `MAX_TARGETS`: allocation count for per-target packing in output arrays.
- `MAX_OBS`: allocation count for obstacle ellipse parameter packing.
- `MAX_SIM_NODES`: allocation count for simulated trajectory step node packing.

**Warmstart / output trajectory buffers**
- `t_out_ptr`, `t_out_size`: nodal wall-clock times (input warmstart, output solution).
- `r_out_ptr`, `r_out_size`: nodal positions (input warmstart, output solution).
- `v_out_ptr`, `v_out_size`: nodal velocities (input warmstart, output solution).
- `a_out_ptr`, `a_out_size`: nodal accelerations/thrust (input warmstart, output solution).
- `t_sim_out_ptr`, `t_sim_out_size`: simulated time vector (output).
- `r_sim_out_ptr`, `r_sim_out_size`: simulated positions (output).
- `v_sim_out_ptr`, `v_sim_out_size`: simulated velocities (output).
- `a_sim_out_ptr`, `a_sim_out_size`: simulated accelerations (output).
- `r0_relax_out_ptr`, `r0_relax_out_size`: reserved initial-relaxation output buffer.
- `rf_relax_out_ptr`, `rf_relax_out_size`: reserved terminal-relaxation output buffer.
- `ddto_converged_ptr`: scalar boolean written to `true` when `solve` converges.

# Returns
- none; results are written in place to the output pointer buffers.
"""
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
        sim_steps::UInt32,
        eps_cvg::CReal,
        eps_ctcs::CReal,
        n::UInt32,
        ctcs_enabled::Bool,
        autogen_init_guess::Bool,
        ddto_init_guess::Bool,
        MAX_HORIZON::UInt32,
        MAX_TARGETS::UInt32,
        MAX_OBS::UInt32,
        MAX_SIM_NODES::UInt32,
        c_x_ptr::Ptr{Cdouble}, c_x_size::Cint,
        c_y_ptr::Ptr{Cdouble}, c_y_size::Cint,
        R_ptr::Ptr{Cdouble}, R_size::Cint,
        M0_ptr::Ptr{Cdouble}, M0_size::Cint,
        M1_ptr::Ptr{Cdouble}, M1_size::Cint,
        t_out_ptr::Ptr{Cdouble}, t_out_size::Cint,
        r_out_ptr::Ptr{Cdouble}, r_out_size::Cint,
        v_out_ptr::Ptr{Cdouble}, v_out_size::Cint,
        a_out_ptr::Ptr{Cdouble}, a_out_size::Cint,
        t_sim_out_ptr::Ptr{Cdouble}, t_sim_out_size::Cint,
        r_sim_out_ptr::Ptr{Cdouble}, r_sim_out_size::Cint,
        v_sim_out_ptr::Ptr{Cdouble}, v_sim_out_size::Cint,
        a_sim_out_ptr::Ptr{Cdouble}, a_sim_out_size::Cint,
        r0_relax_out_ptr::Ptr{Cdouble}, r0_relax_out_size::Cint,
        rf_relax_out_ptr::Ptr{Cdouble}, rf_relax_out_size::Cint,
        ddto_converged_ptr::Ptr{Bool}
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
    t_sim_out = unsafe_wrap(Array, t_sim_out_ptr, t_sim_out_size, own=false)
    r_sim_out = unsafe_wrap(Array, r_sim_out_ptr, r_sim_out_size, own=false)
    v_sim_out = unsafe_wrap(Array, v_sim_out_ptr, v_sim_out_size, own=false)
    a_sim_out = unsafe_wrap(Array, a_sim_out_ptr, a_sim_out_size, own=false)
    r0_relax_out = unsafe_wrap(Array, r0_relax_out_ptr, r0_relax_out_size, own=false)
    rf_relax_out = unsafe_wrap(Array, rf_relax_out_ptr, rf_relax_out_size, own=false)
    ddto_converged = unsafe_wrap(Array, ddto_converged_ptr, 1, own=false)

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
    params.a.uf_targs = repeat(params.a.u0,1,num_targs) # repeat initial input cond
    for j = 1:num_targs
        for i = 1:3
            ind = j + MAX_TARGETS*(i-1)
            params.a.zf_targs[i,j] = rf[ind]
            params.a.zf_targs[i+3,j] = vf[ind]
        end
        params.a.zf_targs[7,j] = Inf
    end
    params.h_constant = params.a.z0[3]

    # >> Target conditions <<
    params.a.n_targs = num_targs
    params.a.λ_targs = collect(1:num_targs)
    params.a.J_targs = collect(1:num_targs)
    params.a.α_targs = ones(num_targs)
    params.a.ϵ_targs = fill(subopt_tol, num_targs)

    # >> Time dilation & discretization <<
    params.a.N = K
    params.a.Δt_min = dt_min
    params.a.Δt_max = dt_max
    params.a.ToF_max = tf_max

    # >> SCP Params <<
    params.a.ctcs_enabled = ctcs_enabled
    params.a.ddto_warmstart = ddto_init_guess
    params.a.w_obj_sing = w_obj
    params.a.w_obj_ddto = w_obj/num_targs
    params.a.w_trust = w_trust
    params.a.w_ctrl = w_buff # keep same to minimize parameters
    params.a.w_buff = w_buff
    params.a.scp_iters = scp_iters
    params.a.sim_steps = sim_steps
    params.a.ϵ_ctrl = eps_cvg
    params.a.ϵ_buff = eps_cvg
    params.a.ϵ_trust = eps_cvg
    params.a.ϵ_ctcs = eps_ctcs

    # >> Build custom scaling matrices <<
    custom_scaling!(params)
    
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

    if autogen_init_guess
        ref_trajs = nothing
    else
        ref_trajs = EmptyDDTOSolution(num_targs)
        for j = 1:params.a.n_targs
            τ_bar = range(0, stop=1, length=K) |> CVector
            s_bar = wall_clock_time_to_time_dilation_control(t_bar[:,j], τ_bar, params.a.disc)
            Δt_bar = diff(t_bar[:,j])
            ∫T_bar = CVector(cumsum([params.a.z0[7],[norm(Δt_bar[k]*a_bar[:,k,j]*params.mass) for k=1:length(s_bar)-1]...]))
            ref_trajs.targs[j] = EmptySolution()
            ref_trajs.targs[j].t = τ_bar
            ref_trajs.targs[j].x = vcat(r_bar[:,:,j], v_bar[:,:,j], reshape(∫T_bar, 1, length(∫T_bar)), zeros(1,length(∫T_bar)))
            ref_trajs.targs[j].u = vcat(a_bar[:,:,j] * params.mass, reshape(s_bar, 1, length(s_bar)))
        end
    end

    ## Call DDTO
    error_thrown = false
    DDTO_target_solutions = EmptyDDTOSolution(num_targs)
    DDTO_target_simulations = EmptyDDTOSolution(num_targs)
    try
        _,_,DDTO_target_solutions,DDTO_target_simulations,converged = solve(params; ref_trajs=ref_trajs)
        ddto_converged[1] = converged
    catch e
        println("!! Error thrown during DDTO solve:")
        try
            println(e)
        catch
            println("No error message available...")
        end
        error_thrown = true
    end

    # Write outputs to memory
    if !error_thrown
        for c = 1:3
            # Solution outputs
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
            # Simulation outputs
            for k = 1:sim_steps*(K-1)
                for t = 1:num_targs
                    ind = t + MAX_TARGETS*(k-1) + MAX_SIM_NODES*MAX_TARGETS*(c-1)
                    r_sim_out[ind] = DDTO_target_simulations.targs[t].r[c,k]
                    v_sim_out[ind] = DDTO_target_simulations.targs[t].v[c,k]
                    a_sim_out[ind] = DDTO_target_simulations.targs[t].T[c,k] / params.mass
                    if c == 1
                        t_sim_out[ind] = DDTO_target_simulations.targs[t].t[k]
                    end
                end
            end
        end
    end
end