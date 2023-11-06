function scenario_obstacles_hard()

    """
    SCENARIO OBJECTIVE:
    To test performance with many obstacles and tight spaces
    """

    # Load default params first
    params = Quad3DoFCageParams()

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
    params.z0 = [r0;v0;0]

    # >> Target conditions <<
    params.n_targs = 4
    rf_targs = hcat(
        -1*e_x - 1.5*e_y - height*e_z,
        +3*e_x - 1.5*e_y - height*e_z,
        +3*e_x + 0.5*e_y - height*e_z,
        +0*e_x + 1.5*e_y - height*e_z,
    )
    vf_targs = zeros(3,params.n_targs)
    params.zf_targs = vcat(rf_targs,vf_targs,Inf*ones(1,params.n_targs)) # Inf: not constraining this state
    params.λ_targs = [3, 2, 4, 1]
    params.T_targs = 1:params.n_targs
    params.ϵ_targs = fill(eps, params.n_targs)

    # >> SCP Params <<
    params.w_obj = 1e3
    params.w_ctrl = 1e6
    params.w_buff = 1e4
    params.w_trust = 1e3
    params.ϵ_ctrl = 1e-2
    params.ϵ_buff = 1e-2
    params.ϵ_trust = 1e-2
    params.scp_iters = 15

    # >> Time dilation & discretization <<
    params.N = 10
    params.τ = CVector(range(0, stop=1, length=params.N))
    params.Δτ = params.τ[2]-params.τ[1]
    Δt_min = 0.001
    Δt_max = 1.
    params.s_min = Δt_min / params.Δτ
    params.s_max = Δt_max / params.Δτ
    params.ToF_max = 20.

    return params
end

function scenario_obstacles_easy()

    """
    SCENARIO OBJECTIVE:
    To test performance with many obstacles and tight spaces
    """

    # Load default params first
    params = Quad3DoFCageParams()

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
    params.z0 = [r0;v0;0]

    # >> Target conditions <<
    params.n_targs = 4
    rf_targs = hcat(
        -1*e_x - 1.5*e_y - height*e_z,
        +3*e_x - 1.5*e_y - height*e_z,
        +3*e_x + 0.5*e_y - height*e_z,
        +0*e_x + 1.5*e_y - height*e_z,
    )
    vf_targs = zeros(3,params.n_targs)
    params.zf_targs = vcat(rf_targs,vf_targs,Inf*ones(1,params.n_targs)) # Inf: not constraining this state
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
    params.N = 20
    params.τ = CVector(range(0, stop=1, length=params.N))
    params.Δτ = params.τ[2]-params.τ[1]
    params.Δt_min = 0.001
    params.Δt_max = 2
    params.s_min = params.Δt_min / params.Δτ
    params.s_max = params.Δt_max / params.Δτ
    params.ToF_max = 10

    return params
end

function scenario_no_obstacles()

    """
    SCENARIO OBJECTIVE:
    To test varying-dilation free-final-time formulation with extremely-separated target states
    """

    # Load default params first
    params = Quad3DoFCageParams()

    # High-level settings
    eps = 0.5 # Accepted level of suboptimality
    height = 1 # [m] Height of the maneuver

    # >> Obstacle parameters <<
    params.n_obstacles = 0 # No obstacles

    # >> Initial condition state <<
    r0 =  0*e_x + 0*e_y - height*e_z
    v0 =  0*e_x + 0*e_y + 0*e_z
    params.z0 = [r0;v0;0]

    # >> Target conditions <<
    params.n_targs = 3
    rf_targs = hcat(
        +10*e_x + 0*e_y - height*e_z,
        +5*e_x  + 3*e_y - height*e_z,
        +2*e_x  - 3*e_y - height*e_z,
    )
    vf_targs = zeros(3,params.n_targs)
    params.zf_targs = vcat(rf_targs,vf_targs,Inf*ones(1,params.n_targs)) # Inf: not constraining this state
    params.λ_targs = [3, 2, 1]
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
    params.scp_iters = 15

    # >> Time dilation & discretization <<
    params.N = 10
    params.τ = CVector(range(0, stop=1, length=params.N))
    params.Δτ = params.τ[2]-params.τ[1]
    params.Δt_min = 0.005
    params.Δt_max = 2.
    params.s_min = params.Δt_min / params.Δτ
    params.s_max = params.Δt_max / params.Δτ
    params.ToF_max = 20.

    return params
end