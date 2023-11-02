function default_params()
    # >> Environmental parameters <<
    g = -9.81*e_z
    ρ = 1.225

    # >> Vehicle parameters <<
    n_rotor = 4
    mass = 0.35
    ρ_min = 1.0
    ρ_max = 7.0
    x_arena_lims = CVector([-4.5,+4.5])
    y_arena_lims = CVector([-2.5,+2.5])
    z_arena_lims = CVector([-2,+0])

    # >> Constraint parameters <<
    γ_p = 45 * DEG_2_RAD
    v_max_V = 0.
    v_max_L = 5.
    nx = 6 # (position, velocity)
    nu = 4 # (thrust, time dilation)

    # Obstacle and boundary parameters 
    # (defaults to empty, scenario-specific)
    n_obstacles = -1
    R_obstacles = CVector(undef,0)
    p_obstacles = CMatrix(undef,0,0)
    H_obstacles = Vector(undef,0)
    n_targs = -1
    z0 = CVector(undef,0)
    zf_targs = CMatrix(undef,0,0)
    λ_targs = Array{Int}(undef,0)
    T_targs = Array{Int}(undef,0)
    ϵ_targs = CVector(undef,0)

    # >> SCP Params <<
    w_obj = 1
    w_ctrl = 1e7
    w_buff = 1e-2
    w_trust = 1e3
    ϵ_ctrl = 1e-2
    ϵ_buff = 1e-2
    ϵ_trust = 1e-2
    scp_iters = 10

    # >> Time dilation & discretization <<
    N = 11
    τ = CVector(range(0, stop=1, length=N))
    Δτ = diff(τ)
    Δt_min = 0.01
    Δt_max = 2.
    s_min = 0.01
    s_max = 3.
    ToF_max = 10.
    disc = 1

    # >> Other <<
    τ_max = 1000

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
        nx,
        nu,
        n_targs,
        zf_targs,
        λ_targs,
        T_targs,
        ϵ_targs,
        w_obj,
        w_ctrl,
        w_buff,
        w_trust,
        ϵ_ctrl,
        ϵ_buff,
        ϵ_trust,
        scp_iters,
        N,
        τ,
        Δτ,
        Δt_min,
        Δt_max,
        s_min,
        s_max,
        ToF_max,
        disc,
        τ_max
    )

    return params
end

function scenario_obstacles_hard()

    """
    SCENARIO OBJECTIVE:
    To test performance with many obstacles and tight spaces
    """

    # Load default params first
    params = default_params()

    # High-level settings
    eps = 0.1  # Accepted level of suboptimality
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
    params.λ_targs = [3, 2, 4, 1]
    params.T_targs = 1:params.n_targs
    params.ϵ_targs = fill(eps, params.n_targs)

    # >> SCP Params <<
    params.w_obj = 1e3
    params.w_ctrl = 1e4
    params.w_buff = 1e4
    params.w_trust = 1e3
    params.ϵ_ctrl = 1e-2
    params.ϵ_buff = 1e-2
    params.ϵ_trust = 1e-2
    params.scp_iters = 10

    # >> Time dilation & discretization <<
    params.N = 21
    params.τ = CVector(range(0, stop=1, length=params.N))
    params.Δτ = diff(params.τ)
    Δt_min = 0.01
    Δt_max = .5
    params.s_min = Δt_min / min(params.Δτ...)
    params.s_max = Δt_max / min(params.Δτ...)
    params.ToF_max = 10

    return params
end

function scenario_obstacles_easy()

    """
    SCENARIO OBJECTIVE:
    To test performance with many obstacles and tight spaces
    """

    # Load default params first
    params = default_params()

    # High-level settings
    eps = .1  # Accepted level of suboptimality
    obs_rad = 0.6 # [m] Radius of all cylindrical obstacles
    height = 1 # [m] Height of the maneuver

    # >> Obstacle parameters <<
    params.n_obstacles = 3 # Number of obstacles
    params.R_obstacles = fill(obs_rad, params.n_obstacles) # Radii of all circular obstacles
    params.p_obstacles = hcat( # Positions of circular obstacless
       +2*e_x + 0.5*e_y - height*e_z,
       -2*e_x + 0.5*e_y - height*e_z,
       +0*e_x - 0.5*e_y - height*e_z,
    )
    params.H_obstacles = repeat([I(3)],params.n_obstacles)

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
    params.λ_targs = [3, 2, 4, 1]
    params.T_targs = 1:params.n_targs
    params.ϵ_targs = fill(eps, params.n_targs)

    # >> SCP Params <<
    params.w_obj = 1e0
    params.w_ctrl = 1e5
    params.w_buff = 1e4
    params.w_trust = 1e3
    params.ϵ_ctrl = 1e-2
    params.ϵ_buff = 1e-2
    params.ϵ_trust = 1e-2
    params.scp_iters = 10

    # >> Time dilation & discretization <<
    params.N = 21
    params.τ = CVector(range(0, stop=1, length=params.N))
    params.Δτ = diff(params.τ)
    params.Δt_min = 0.001
    params.Δt_max = 0.1
    params.s_min = params.Δt_min / min(params.Δτ...)
    params.s_max = params.Δt_max / min(params.Δτ...)
    params.ToF_max = 10

    return params
end

function scenario_no_obstacles()

    """
    SCENARIO OBJECTIVE:
    To test varying-dilation free-final-time formulation with extremely-separated target states
    """

    # Load default params first
    params = default_params()

    # High-level settings
    eps = 0.2 # Accepted level of suboptimality
    height = 1 # [m] Height of the maneuver

    # >> Obstacle parameters <<
    params.n_obstacles = 0 # No obstacles

    # >> Initial condition state <<
    r0 =  0*e_x + 0*e_y - height*e_z
    v0 =  0*e_x + 0*e_y + 0*e_z
    params.z0 = [r0;v0]

    # >> Target conditions <<
    params.n_targs = 3
    rf_targs = hcat(
        +10*e_x + 0*e_y - height*e_z,
        +5*e_x  + 3*e_y - height*e_z,
        +2*e_x  - 3*e_y - height*e_z,
    )
    vf_targs = zeros(3,params.n_targs)
    params.zf_targs = vcat(rf_targs,vf_targs)
    params.λ_targs = [3, 2, 1]
    params.T_targs = 1:params.n_targs
    params.ϵ_targs = fill(eps, params.n_targs)

    # >> SCP Params <<
    params.w_obj = 1e1
    params.w_ctrl = 1e6
    params.w_buff = 1e4
    params.w_trust = 1e4
    params.ϵ_ctrl = 1e-2
    params.ϵ_buff = 1e-2
    params.ϵ_trust = 1e-2
    params.scp_iters = 10

    # >> Time dilation & discretization <<
    params.N = 21
    params.τ = CVector(range(0, stop=1, length=params.N))
    params.Δτ = diff(params.τ)
    params.Δt_min = 1e-6
    params.Δt_max = 0.5
    params.s_min = params.Δt_min / min(params.Δτ...)
    params.s_max = params.Δt_max / min(params.Δτ...)
    params.ToF_max = 10

    return params
end