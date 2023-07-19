function default_params()
    # >> Environmental parameters <<
    g = -9.81*e_z
    ρ = 1.225

    # >> Params parameters <<
    n_rotor = 4
    mass = 0.35
    ρ_min = 1.0
    ρ_max = 7.0
    x_arena_lims = CVector([-4.5,+4.5])
    y_arena_lims = CVector([-2.5,+2.5])
    z_arena_lims = CVector([-2,+0])

    # >> Constraint parameters <<
    γ_p = 45 * DEG_2_RAD
    v_max_V = 0
    v_max_L = 5

    # >> Dynamics <<
    # (pre-augmentation in free-final-time case)
    A_c = CMatrix([
        zeros(3,3) I(3);
        zeros(3,3) zeros(3,3)
    ])
    B_c = CMatrix([
        zeros(3,3);
        I(3)/mass
    ])
    # B_c = CMatrix([
    #     zeros(3,3) zeros(3);
    #     I(3)/mass  zeros(3)
    # ])
    p_c = CVector(vcat(zeros(3),g))
    n,m = size(B_c)

    # Default scenario parameters
    # >> Obstacle parameters <<
    R_obstacles = [
        0.3,
        0.5,
        0.6,
        0.2,
        0.2
    ] # Radii of all circular obstacles
    n_obstacles = length(R_obstacles) # Number of obstacles
    p_obstacles = hcat( # Positions of circular obstacless
       -1.25*e_x + 0.5*e_y - 1*e_z,
        0*e_x    + 0*e_y   - 1*e_z,
        1*e_x    + 1*e_y   - 1*e_z,
        1*e_x    - 0.7*e_y - 1*e_z,
        1.8*e_x  + 0.3*e_y - 1*e_z,
    )
    H_obstacles = repeat([I(3)],n_obstacles)

    # >> Initial condition state <<
    r0 = -3*e_x + 1*e_y - 1*e_z
    v0 =  0*e_x + 0*e_y + 0*e_z
    z0 = [r0;v0]

    # >> Target conditions <<
    n_targs = 3
    rf_targs = hcat(
        0.0*e_x - 1.5*e_y - 1*e_z, # Target 1
        2.0*e_x - 1.0*e_y - 1*e_z, # Target 2
        3.0*e_x + 1*e_y   - 1*e_z, # Target 3
    ) 
    vf_targs = hcat(
        0*e_x + 0*e_y + 0*e_z, # Target 1
        0*e_x + 0*e_y + 0*e_z, # Target 2
        0*e_x + 0*e_y + 0*e_z, # Target 3
    )
    zf_targs = vcat(rf_targs,vf_targs)
    λ_targs = [1, 2, 3]
    T_targs = 1:n_targs
    ϵ_targs = CVector([0.2, 0.2, 0.2])

    # >> SCP Params <<
    w_ctrl = 1e7
    w_buff = 1e-2
    w_trust = 1e3
    ϵ_ctrl = 1e-2
    ϵ_buff = 1e-2
    ϵ_trust = 1e-2
    scp_iters = 10

    # >> Time Discretization <<
    free_final_time = true
    disc = 1

    # Fixed-final-time
    Δt = 0.5
    N_targs = [21, 21, 21]

    # Free-final-time
    N_fft = 21
    τ = CVector(range(0, stop=1, length=N_fft))
    Δτ = diff(τ)
    Δt_min = 0.01
    Δt_max = 2
    s_min = 0.01
    s_max = 3
    ToF_max = 10

    # >> Other <<
    τ_max = 1e10

    # >> Make params object <<
    params = Quad3DoFCageParams(
        g,
        ρ,
        n_rotor,
        mass,
        ρ_min,
        ρ_max,
        γ_p,
        v_max_V,
        v_max_L,
        n_obstacles,
        R_obstacles,
        p_obstacles,
        H_obstacles,
        x_arena_lims,
        y_arena_lims,
        z_arena_lims,
        z0,
        n_targs,
        zf_targs,
        λ_targs,
        T_targs,
        ϵ_targs,
        n,
        m,
        A_c,
        B_c,
        p_c,
        w_ctrl,
        w_buff,
        w_trust,
        ϵ_ctrl,
        ϵ_buff,
        ϵ_trust,
        scp_iters,
        free_final_time,
        disc,
        Δt,
        N_targs,
        N_fft,
        τ,
        Δτ,
        Δt_min,
        Δt_max,
        s_min,
        s_max,
        ToF_max,
        τ_max,
    )

    return params
end

function scenario_toy1()

    """
    SCENARIO OBJECTIVE:
    To test performance with many obstacles and tight spaces
    """

    # Load default params first
    params = default_params()

    # High-level settings
    eps = 0.2  # Accepted level of suboptimality
    obs_rad = 0.6 # [m] Radius of all cylindrical obstacles
    height = 1 # [m] Height of the maneuver

    # >> Obstacle parameters <<
    params.n_obstacles = 14 # Number of obstacles
    params.R_obstacles = fill(obs_rad, params.n_obstacles) # Radii of all circular obstacles
    params.p_obstacles = hcat( # Positions of circular obstacless
       -3*e_x + 1.5*e_y - height*e_z,
       -1*e_x + 1.5*e_y - height*e_z,
       +1*e_x + 1.5*e_y - height*e_z,
       +3*e_x + 1.5*e_y - height*e_z,
       -2*e_x + 0.5*e_y - height*e_z,
       +0*e_x + 0.5*e_y - height*e_z,
       +2*e_x + 0.5*e_y - height*e_z,
       -3*e_x - 0.5*e_y - height*e_z,
       -1*e_x - 0.5*e_y - height*e_z,
       +1*e_x - 0.5*e_y - height*e_z,
       +3*e_x - 0.5*e_y - height*e_z,
       -2*e_x - 1.5*e_y - height*e_z,
       +0*e_x - 1.5*e_y - height*e_z,
       +2*e_x - 1.5*e_y - height*e_z,
    )
    params.H_obstacles = repeat([I(3)],params.n_obstacles)

    # >> Dynamics <<
    params.Δt = 0.2

    # >> Initial condition state <<
    r0 = -3*e_x + 0.5*e_y - height*e_z
    v0 =  0*e_x + 0*e_y + 0*e_z
    params.z0 = [r0;v0]

    # >> Target conditions <<
    params.n_targs = 4
    rf_targs = hcat(
        -1*e_x - 1.5*e_y - height*e_z,
        +3*e_x - 1.5*e_y - height*e_z,
        +3*e_x + 0.5*e_y - height*e_z,
        +0*e_x + 1.5*e_y - height*e_z,
    )
    vf_targs = zeros(3,params.n_targs)
    params.zf_targs = vcat(rf_targs,vf_targs)
    params.N_targs = fill(21, params.n_targs)
    params.λ_targs = [3, 2, 4, 1]
    params.T_targs = 1:params.n_targs
    params.ϵ_targs = fill(eps, params.n_targs)
    params.τ_max   = 1e10 # arbitrarily-large value to disable

    # >> SCP Params <<
    params.w_ctrl = 1e5
    params.w_buff = 1e4
    params.w_trust = 1e3

    return params
end

function scenario_toy2()

    """
    SCENARIO OBJECTIVE:
    To test varying-dilation free-final-time formulation with extremely-separated target states
    """

    # Load default params first
    params = default_params()

    # High-level settings
    eps = 0.2  # Accepted level of suboptimality
    height = 1 # [m] Height of the maneuver

    # >> Obstacle parameters <<
    params.n_obstacles = 0 # No obstacles

    # >> Dynamics <<
    params.Δt = 0.2 # not used for free-final-time!

    # >> Initial condition state <<
    r0 =  0*e_x + 0*e_y + height*e_z
    v0 =  0*e_x + 0*e_y + 0*e_z
    params.z0 = [r0;v0]

    # >> Target conditions <<
    params.n_targs = 3
    rf_targs = hcat(
        +10*e_x + 0*e_y + height*e_z,
        +5*e_x  + 3*e_y + height*e_z,
        +2*e_x  - 3*e_y + height*e_z,
    )
    vf_targs = zeros(3,params.n_targs)
    params.zf_targs = vcat(rf_targs,vf_targs)
    params.N_targs = fill(21, params.n_targs) # not used for free-final-time!
    params.λ_targs = [3, 2, 1]
    params.T_targs = 1:params.n_targs
    params.ϵ_targs = fill(eps, params.n_targs)
    params.τ_max   = 1e10 # arbitrarily-large value to disable

    # >> SCP Params <<
    params.w_ctrl = 1e4
    params.w_buff = 0
    params.w_trust = 1e-0

    return params
end