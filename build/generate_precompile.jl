# for testing, cd into DDTOSCP folder and run as follows:
# julia --startup-file=no --project=. build/generate_precompile.jl
# Note: to see help for Julia flags use julia -h, or julia --help-hidden.

using DDTOSCP

# enable debug printing
# ENV["JULIA_DEBUG"]=DDTOSCP

begin
    # Set test variables based on default Params
    params = DDTOSCP.Quad3DoFCageSampleScenario()

    # Define maximum sizing parameters (since C++ wants static array sizing)
    MAX_HORIZON = UInt32(50) # Set as arbitrarily large number for test
    MAX_TARGETS = UInt32(50) # Set as arbitrarily large number for test
    MAX_OBS     = UInt32(50) # Set as arbitrarily large number for test
    MAX_SIM_STEPS = UInt32(50) # Set as arbitrarily large number for test
    MAX_SIM_NODES = UInt32(MAX_SIM_STEPS * (MAX_HORIZON - 1))

    # Define empty arrays with static sizing
    r0 = zeros(3)
    v0 = zeros(3)
    rf = zeros(MAX_TARGETS,3)
    vf = zeros(MAX_TARGETS,3)
    c_x = zeros(MAX_OBS)
    c_y = zeros(MAX_OBS)
    R = zeros(MAX_OBS)
    M0 = zeros(MAX_OBS,2)
    M1 = zeros(MAX_OBS,2)
    t_out = zeros(MAX_TARGETS,MAX_HORIZON)
    r_out = zeros(MAX_TARGETS,MAX_HORIZON,3)
    v_out = zeros(MAX_TARGETS,MAX_HORIZON,3)
    a_out = zeros(MAX_TARGETS,MAX_HORIZON,3)
    r_sim_out = zeros(MAX_TARGETS,MAX_SIM_NODES,3)
    r0_relax_out = zeros(3)
    rf_relax_out = zeros(MAX_TARGETS,3)

    # >> Params parameters <<
    a_min = Float64(params.ρ_min / params.mass)
    a_max = Float64(params.ρ_max / params.mass)

    # >> Constraint parameters <<
    theta_max = Float64(params.γ_p)
    v_max = Float64(params.v_max_L)

    # >> Obstacle parameters <<
    n = UInt32(params.n_obstacles)
    for j = 1:n
        R[j] = params.R_obstacles[j]
        c_x[j] = params.p_obstacles[1,j]
        c_y[j] = params.p_obstacles[2,j]
        M0[j,1] = params.H_obstacles[j][1,1]
        M0[j,2] = params.H_obstacles[j][2,1]
        M1[j,1] = params.H_obstacles[j][1,2]
        M1[j,2] = params.H_obstacles[j][2,2]
    end

    # >> Target conditions <<
    num_targs = UInt32(params.a.n_targs)
    subopt_tol = Float64(params.a.ϵ_targs[1])

    # >> Boundary conditions <<
    for i = 1:3
        r0[i] = params.a.z0[i]
        v0[i] = params.a.z0[i+3]
        for j = 1:num_targs
            rf[j,i] = params.a.zf_targs[i,j]
            vf[j,i] = params.a.zf_targs[i+3,j]
        end
    end

    # >> Time dilation & discretization <<
    K = UInt32(params.a.N)
    tf = Float64(0) # free-final-time only currently
    tf_max = Float64(params.a.ToF_max)
    dt_min = Float64(params.a.Δt_min)
    dt_max = Float64(params.a.Δt_max)

    # >> SCP Params <<
    w_obj = Float64(params.a.w_obj_sing)
    w_buff = Float64(params.a.w_buff)
    w_trust = Float64(params.a.w_trust)
    ri_relax = Float64(0)
    rf_relax = Float64(0)
    scp_iters = UInt32(params.a.scp_iters)
    sim_steps = UInt32(params.a.sim_steps)
    eps_cvg = Float64(params.a.ϵ_trust)
    eps_ctcs = Float64(params.a.ϵ_ctcs)
    ctcs_enabled = true
    interp_ref = false

    # >> Populate _out variables with a sample reference trajectory (initial guess generated) <<
    ref_trajs = DDTOSCP.generate_initial_guess_ddtoscp(params)
    for c = 1:3
        for k = 1:K
            for t = 1:num_targs
                ind = t + MAX_TARGETS*(k-1) + MAX_HORIZON*MAX_TARGETS*(c-1)
                r_out[ind] = ref_trajs.targs[t].x[c,k]
                v_out[ind] = ref_trajs.targs[t].x[c+3,k]
                a_out[ind] = ref_trajs.targs[t].u[c,k] / params.mass
                if c == 1
                    t_out[ind] = ref_trajs.targs[t].t[k]
                end
            end
        end
    end
    
    ## Convert all arrays to pointer/size combinations
    # Pointers
    r0_ptr = Ptr{Cdouble}(pointer(r0))
    v0_ptr = Ptr{Cdouble}(pointer(v0))
    rf_ptr = Ptr{Cdouble}(pointer(rf))
    vf_ptr = Ptr{Cdouble}(pointer(vf))
    c_x_ptr = Ptr{Cdouble}(pointer(c_x))
    c_y_ptr = Ptr{Cdouble}(pointer(c_y))
    R_ptr = Ptr{Cdouble}(pointer(R))
    M0_ptr = Ptr{Cdouble}(pointer(M0))
    M1_ptr = Ptr{Cdouble}(pointer(M1))
    t_out_ptr = Ptr{Cdouble}(pointer(t_out))
    r_out_ptr = Ptr{Cdouble}(pointer(r_out))
    v_out_ptr = Ptr{Cdouble}(pointer(v_out))
    a_out_ptr = Ptr{Cdouble}(pointer(a_out))
    r_sim_out_ptr = Ptr{Cdouble}(pointer(r_sim_out))
    r0_relax_out_ptr = Ptr{Cdouble}(pointer(r0_relax_out))
    rf_relax_out_ptr = Ptr{Cdouble}(pointer(rf_relax_out))

    # Sizes
    r0_size = Base.cconvert(Cint, length(r0))
    v0_size = Base.cconvert(Cint, length(v0))
    rf_size = Base.cconvert(Cint, length(rf))
    vf_size = Base.cconvert(Cint, length(vf))
    c_x_size = Base.cconvert(Cint, length(c_x))
    c_y_size = Base.cconvert(Cint, length(c_y))
    R_size = Base.cconvert(Cint, length(R))
    M0_size = Base.cconvert(Cint, length(M0))
    M1_size = Base.cconvert(Cint, length(M1))
    t_out_size = Base.cconvert(Cint, length(t_out))
    r_out_size = Base.cconvert(Cint, length(r_out))
    v_out_size = Base.cconvert(Cint, length(v_out))
    a_out_size = Base.cconvert(Cint, length(a_out))
    r_sim_out_size = Base.cconvert(Cint, length(r_sim_out))
    r0_relax_out_size = Base.cconvert(Cint, length(r0_relax_out))
    rf_relax_out_size = Base.cconvert(Cint, length(rf_relax_out))

    GC.@preserve r0 v0 rf vf c_x c_y R M0 M1 t_out r_out v_out a_out r0_relax_out rf_relax_out begin
        # Call skyenet_interface
        DDTOSCP.skyenet_ddtoscp_interface(
            num_targs,
            r0_ptr, r0_size,
            v0_ptr, v0_size,
            rf_ptr, rf_size,
            vf_ptr, vf_size,
            K,
            tf,
            tf_max,
            dt_min,
            dt_max,
            a_min,
            a_max,
            v_max,
            theta_max,
            ri_relax,
            rf_relax,
            subopt_tol,
            w_obj,
            w_trust,
            w_buff,
            scp_iters,
            sim_steps,
            eps_cvg,
            eps_ctcs,
            n,
            ctcs_enabled,
            interp_ref,
            MAX_HORIZON,
            MAX_TARGETS,
            MAX_OBS,
            MAX_SIM_NODES,
            c_x_ptr, c_x_size,
            c_y_ptr, c_y_size,
            R_ptr, R_size,
            M0_ptr, M0_size,
            M1_ptr, M1_size,
            t_out_ptr, t_out_size,
            r_out_ptr, r_out_size,
            v_out_ptr, v_out_size,
            a_out_ptr, a_out_size,
            r_sim_out_ptr, r_sim_out_size,
            r0_relax_out_ptr, r0_relax_out_size,
            rf_relax_out_ptr, rf_relax_out_size
        )
    end
end