function scenario_toy1!(params::Params)

    """
    SCENARIO OBJECTIVE:
    To test performance with many obstacles and tight spaces
    """

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
    params.w_buff = 10
    params.w_trust = 0.01
end

function scenario_toy2!(params::Params)

    """
    SCENARIO OBJECTIVE:
    To test varying-dilation free-final-time formulation with extremely-separated target states
    """

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
    params.λ_targs = [3, 2, 4, 1]
    params.T_targs = 1:params.n_targs
    params.ϵ_targs = fill(eps, params.n_targs)
    params.τ_max   = 1e10 # arbitrarily-large value to disable

    # >> SCP Params <<
    params.w_ctrl = 1e7
    params.w_buff = 0
    params.w_trust = 1e3
end