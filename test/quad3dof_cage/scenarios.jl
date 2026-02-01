using Distributions

function scenario_obstacles_hard(lex::Bool=false)

    """
    SCENARIO OBJECTIVE:
    To test performance with many obstacles and tight spaces
    """

    # Load default params first
    params = Quad3DoFCageParams()

    # High-level settings
    eps = .2  # Accepted level of suboptimality
    obs_rad = 0.6 # [m] Radius of all cylindrical obstacles
    height = 1 # [m] Height of the maneuver plane above ground level

    # >> Obstacle parameters <<
    params.n_obstacles = 14 # Number of obstacles
    params.R_obstacles = fill(obs_rad, params.n_obstacles) # Radii of all circular obstacles
    params.p_obstacles = hcat( # Positions of circular obstacless
       -3*e_y + 1.5*e_x - height*e_z,
       -1*e_y + 1.5*e_x - height*e_z,
       +1*e_y + 1.5*e_x - height*e_z,
       +3*e_y + 1.5*e_x - height*e_z,
       -2*e_y + 0.5*e_x - height*e_z,
       +0*e_y + 0.5*e_x - height*e_z,
       +2*e_y + 0.5*e_x - height*e_z,
       -3*e_y - 0.5*e_x - height*e_z,
       -1*e_y - 0.5*e_x - height*e_z,
       +1*e_y - 0.5*e_x - height*e_z,
       +3*e_y - 0.5*e_x - height*e_z,
       -2*e_y - 1.5*e_x - height*e_z,
       +0*e_y - 1.5*e_x - height*e_z,
       +2*e_y - 1.5*e_x - height*e_z,
    )
    params.H_obstacles = repeat([I(3)],params.n_obstacles)

    # >> Initial condition state <<
    r0 = -3*e_y + 0.5*e_x - height*e_z
    v0 =  0*e_y + 0*e_x + 0*e_z
    params.a.z0 = [r0;v0;0]
    params.h_constant = params.a.z0[3]

    # >> Target conditions <<
    params.a.n_targs = 4
    rf_targs = hcat(
        -1*e_y - 1.5*e_x - height*e_z,
        +3*e_y - 1.5*e_x - height*e_z,
        +3*e_y + 0.5*e_x - height*e_z,
        +0*e_y + 1.5*e_x - height*e_z,
    )
    vf_targs = zeros(3,params.a.n_targs)
    params.a.zf_targs = vcat(rf_targs,vf_targs,Inf*ones(1,params.a.n_targs)) # Inf: not constraining this state
    params.a.uf_targs = repeat(params.a.u0,1,params.a.n_targs) # repeat initial input cond
    params.a.λ_targs = [1,4,2,3]
    params.a.J_targs = 1:params.a.n_targs
    params.a.α_targs = ones(params.a.n_targs)
    params.a.ϵ_targs = fill(eps, params.a.n_targs)

    # >> SCP Params <<
    params.a.ctcs_enabled = true
    params.a.warmstart_method = "single"
    params.a.w_obj_sing = 1e-1
    params.a.w_obj_ddto = 1e0
    params.a.w_ctrl = 5e1
    params.a.w_buff = params.a.w_ctrl
    params.a.w_trust = 2e0
    params.a.ϵ_ctrl = 5e-4
    params.a.ϵ_buff = 5e-4
    params.a.ϵ_trust = 5e-4
    params.a.scp_iters = 100

    # >> Time dilation & discretization <<
    params.a.N = 15
    params.a.Δt_min = 0.01
    params.a.Δt_max = 1.
    params.a.ToF_max = 20.

    # >> Update some settings for DDTO-LEX specifically <<
    if lex
        # >> SCP Params <<
        params.a.ctcs_enabled = true
        params.a.warmstart_method = "single"
        params.a.w_obj_sing = 1e-1
        params.a.w_obj_ddto = Inf # not used for DDTO-LEX
        params.a.w_ctrl = 5e1
        params.a.w_buff = params.a.w_ctrl
        params.a.w_trust = 2e0
        params.a.ϵ_ctrl = 5e-4
        params.a.ϵ_buff = 5e-4
        params.a.ϵ_trust = 5e-4
        params.a.scp_iters = 200

        # >> Time dilation & discretization <<
        params.a.N = 10
        params.a.Δt_min = 0.01
        params.a.Δt_max = 1.
        params.a.ToF_max = 20.
    end

    # >> Build custom scaling matrices <<
    custom_scaling!(params)

    return params
end

# Write me a new scenario which takes in the definition of scenario_obstacles_hard and randomly generates n targets that are within the arena bounds and at least x distance away from the nearest obstacle
function scenario_obstacles_hard_random_targets(;n_targets::Int=4, min_distance_from_obstacle::Float64=0.01, lex::Bool=false)
    # Instantiation of scenario for n targets
    params = scenario_obstacles_hard(lex)
    params.a.n_targs = n_targets
    params.a.uf_targs = Inf*ones(3,n_targets)
    params.a.λ_targs = collect(1:n_targets)
    params.a.J_targs = collect(1:n_targets)
    params.a.α_targs = ones(n_targets)
    params.a.ϵ_targs = params.a.ϵ_targs[1]*ones(n_targets)

    # Initialize zf_targs
    rf_targs = zeros(3,params.a.n_targs)
    vf_targs = zeros(3,params.a.n_targs)
    params.a.zf_targs = vcat(rf_targs,vf_targs,Inf*ones(1,params.a.n_targs))

    # Generate random distributions 
    x_dist = Uniform(params.x_arena_lims[1], params.x_arena_lims[2])
    y_dist = Uniform(params.y_arena_lims[1], params.y_arena_lims[2])
    z_setpoint = params.a.z0[3] # Keep altitude constant

    # Randomly acquire targets that meet the criteria
    targets_found = 0
    rf_targs = zeros(3,n_targets)
    while targets_found < n_targets
        # Random sample in distributions
        x = rand(x_dist)
        y = rand(y_dist)
        z = z_setpoint
        target = CVector([x;y;z])
        
        # Iterate through all obstacles and ensure the target is at least x distance away from the nearest obstacle
        obstacle_too_close = false
        for i = 1:params.n_obstacles
            if norm(target - params.H_obstacles[i]*params.p_obstacles[:,i]) < (params.R_obstacles[i] + min_distance_from_obstacle)
                obstacle_too_close = true
                break
            end
        end
        if obstacle_too_close
            continue
        end

        targets_found += 1
        rf_targs[:,targets_found] = target
    end

    # Sort the targets by distance from the initial condition
    # first target to defer should be closest, and so on
    dist_from_initial = [norm(rf_targs[:,k] - params.a.z0[1:3]) for k = 1:n_targets]
    idx_sort = sortperm(dist_from_initial)
    rf_targs = rf_targs[:,idx_sort]
    params.a.zf_targs[1:3,:] = rf_targs

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
    eps = 0.1 # Accepted level of suboptimality
    height = 1 # [m] Height of the maneuver

    # >> Obstacle parameters <<
    params.n_obstacles = 0 # No obstacles

    # >> Initial condition state <<
    r0 =  0*e_y + 0*e_x - height*e_z
    v0 =  0*e_y + 0*e_x + 0*e_z
    params.a.z0 = [r0;v0;0]
    params.h_constant = params.a.z0[3]
    params.cage_bounds_enabled = false

    # >> Target conditions <<
    params.a.n_targs = 3
    rf_targs = hcat(
        +10*e_y + 0*e_x - height*e_z,
        +5*e_y  + 3*e_x - height*e_z,
        +2*e_y  - 3*e_x - height*e_z,
    )
    vf_targs = zeros(3,params.a.n_targs)
    params.a.zf_targs = vcat(rf_targs,vf_targs,Inf*ones(1,params.a.n_targs)) # Inf: not constraining this state
    params.a.uf_targs = repeat(params.a.u0,1,params.a.n_targs) # repeat initial input cond
    params.a.λ_targs = [3,2,1]
    params.a.J_targs = 1:params.a.n_targs
    params.a.α_targs = [1,1,1]
    params.a.ϵ_targs = fill(eps, params.a.n_targs)

    # >> SCP Params <<
    params.a.w_obj_sing = 1e-2
    params.a.w_obj_ddto = params.a.w_obj_sing
    params.a.w_ctrl = 5e1
    params.a.w_buff = params.a.w_ctrl
    params.a.w_trust = 1e0
    params.a.ϵ_ctrl = 1e-3
    params.a.ϵ_buff = 1e-3
    params.a.ϵ_trust = 1e-3
    params.a.scp_iters = 100

    # >> Time dilation & discretization <<
    params.a.N = 12
    params.a.Δt_min = 0.01
    params.a.Δt_max = 2.
    params.a.ToF_max = 20.

    # >> Build custom scaling matrices <<
    custom_scaling!(params)

    return params
end