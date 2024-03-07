function scenario_obstacles_hard()

    """
    SCENARIO OBJECTIVE:
    To test performance with many obstacles and tight spaces
    """

    # Load default params first
    params = Quad3DoFCageParams()

    # High-level settings
    eps = .01  # Accepted level of suboptimality
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
    params.a.z0 = [r0;v0;0]
    params.h_constant = params.a.z0[3]

    # >> Target conditions <<
    params.a.n_targs = 4
    rf_targs = hcat(
        -1*e_x - 1.5*e_y - height*e_z,
        +3*e_x - 1.5*e_y - height*e_z,
        +3*e_x + 0.5*e_y - height*e_z,
        +0*e_x + 1.5*e_y - height*e_z,
    )
    vf_targs = zeros(3,params.a.n_targs)
    params.a.zf_targs = vcat(rf_targs,vf_targs,Inf*ones(1,params.a.n_targs)) # Inf: not constraining this state
    params.a.uf_targs = repeat(params.a.u0,1,params.a.n_targs) # repeat initial input cond
    params.a.λ_targs = [1,4,2,3]
    params.a.T_targs = 1:params.a.n_targs
    params.a.α_targs = ones(params.a.n_targs)
    params.a.ϵ_targs = fill(eps, params.a.n_targs)

    # >> SCP Params <<
    params.a.ctcs_enabled = true
    params.a.ddto_warmstart = false
    params.a.w_obj_sing = .01
    params.a.w_obj_ddto = params.a.w_obj_sing/params.a.n_targs
    params.a.w_ctrl = 50
    params.a.w_buff = params.a.w_ctrl
    params.a.w_trust = 2
    params.a.ϵ_ctrl = 1e-4
    params.a.ϵ_buff = 1e-4
    params.a.ϵ_trust = 1e-4
    params.a.scp_iters = 100

    # >> Time dilation & discretization <<
    params.a.N = 10
    params.a.Δt_min = 0.2
    params.a.Δt_max = 0.7
    params.a.ToF_max = 10.

    # >> Build custom scaling matrices <<
    custom_scaling!(params)

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
    eps = 0.1  # Accepted level of suboptimality
    obs_rad = 0.4 # [m] Radius of all cylindrical obstacles
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
    params.a.z0 = [r0;v0;0]
    params.h_constant = params.a.z0[3]

    # >> Target conditions <<
    params.a.n_targs = 4
    rf_targs = hcat(
        -1*e_x - 1.5*e_y - height*e_z,
        +3*e_x - 1.5*e_y - height*e_z,
        +3*e_x + 0.5*e_y - height*e_z,
        +0*e_x + 1.5*e_y - height*e_z,
    )
    vf_targs = zeros(3,params.a.n_targs)
    params.a.zf_targs = vcat(rf_targs,vf_targs,Inf*ones(1,params.a.n_targs)) # Inf: not constraining this state
    params.a.uf_targs = repeat(params.a.u0,1,params.a.n_targs) # repeat initial input cond
    params.a.λ_targs = [1,2,3,4]
    params.a.T_targs = 1:params.a.n_targs
    params.a.α_targs = ones(params.a.n_targs)
    params.a.ϵ_targs = fill(eps, params.a.n_targs)

    # >> SCP Params <<
    params.a.ctcs_enabled = true
    params.a.w_obj_sing = 1e0
    params.a.w_obj_ddto = 5e-1
    params.a.w_ctrl = 1e4
    params.a.w_buff = params.a.w_ctrl
    params.a.w_trust = 1e3
    params.a.ϵ_ctrl = 1e-4
    params.a.ϵ_buff = 1e-4
    params.a.ϵ_trust = 1e-4
    params.a.scp_iters = 50

    # >> Time dilation & discretization <<
    params.a.N = 12
    params.a.Δt_min = 0.001
    params.a.Δt_max = 2
    params.a.ToF_max = 20

    # >> Build custom scaling matrices <<
    custom_scaling!(params)

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
    params.a.z0 = [r0;v0;0]
    params.h_constant = params.a.z0[3]
    params.cage_bounds_enabled = false

    # >> Target conditions <<
    params.a.n_targs = 3
    rf_targs = hcat(
        +10*e_x + 0*e_y - height*e_z,
        +5*e_x  + 3*e_y - height*e_z,
        +2*e_x  - 3*e_y - height*e_z,
    )
    vf_targs = zeros(3,params.a.n_targs)
    params.a.zf_targs = vcat(rf_targs,vf_targs,Inf*ones(1,params.a.n_targs)) # Inf: not constraining this state
    params.a.uf_targs = repeat(params.a.u0,1,params.a.n_targs) # repeat initial input cond
    params.a.λ_targs = [3,2,1]
    params.a.T_targs = 1:params.a.n_targs
    params.a.α_targs = [0,0,1]
    params.a.ϵ_targs = fill(eps, params.a.n_targs)

    # >> SCP Params <<
    params.a.w_obj_sing = 1e1
    params.a.w_obj_ddto = 1e1
    params.a.w_ctrl = 1e3
    params.a.w_buff = params.a.w_ctrl
    params.a.w_trust = 1e2
    params.a.ϵ_ctrl = 1e-4
    params.a.ϵ_buff = 1e-4
    params.a.ϵ_trust = 1e-4
    params.a.scp_iters = 10

    # >> Time dilation & discretization <<
    params.a.N = 11
    params.a.Δt_min = 0.005
    params.a.Δt_max = 2.
    params.a.ToF_max = 10.

    # >> Build custom scaling matrices <<
    custom_scaling!(params)

    return params
end