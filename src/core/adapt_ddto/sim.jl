#=
Adaptive-DDTO perception-simulation utilities: synthetic landing-target pools,
noisy radius updates, target refresh / allocation, and obstacle generation.
=#

"""
    LandingTarget

Simulated landing site from a perception stack.

# Fields
- `id::Int`: unique target identifier
- `R::CReal`: estimated safe landing radius ``[m]``
- `rf::CVector`: inertial landing position ``[m]``
"""
mutable struct LandingTarget
    id::Int     # Target ID
    R::CReal    # Target radius
    rf::CVector # Target position
end

"""
    sim_build_target_pool(num_targets, R_landing_region; min_radius=1., starting_radius_fac=2.) -> Vector{LandingTarget}

Build a pool of landing targets placed equidistantly on a circle of radius
`R_landing_region`, each initialized with radius `min_radius * starting_radius_fac`.

# Arguments
- `num_targets`: number of synthetic landing sites to create.
- `R_landing_region`: horizontal placement radius `[m]` about the origin.
- `min_radius`: minimum landing-disk radius used as the base scale `[m]` (default `1.`).
- `starting_radius_fac`: multiplier applied to `min_radius` for initial radii (default `2.`).

# Returns
- `target_pool`: vector of [`LandingTarget`](@ref) instances with IDs, radii, and positions.
"""
function sim_build_target_pool(num_targets::Int, R_landing_region::CReal; min_radius::CReal=1., starting_radius_fac::CReal=2.)::Vector{LandingTarget}
    # Initialize target radii to all be the same value 
    R_targs = fill(min_radius*starting_radius_fac, num_targets)

    # Place targets equidistant along the landing region radius
    rf_targs = CMatrix(undef, 3, num_targets)
    θ_sep = 2 * pi / num_targets
    θ_targ = 0.
    for j = 1:num_targets
        r_targ = R_landing_region
        rf_targs[:,j] = r_targ*cos(θ_targ)*e_x + r_targ*sin(θ_targ)*e_y + 1*e_z
        θ_targ += θ_sep
    end

    # Create the target pool
    target_pool = [LandingTarget(j, R_targs[j], rf_targs[:,j]) for j=1:num_targets]
    return target_pool
end

"""
    sim_update_targets!(params, target_pool; noise_std=0.2, crossweight=0.05)

Simulate a perception update of landing radii with Gaussian noise and a
cross-fade filter; sync allocated targets into `params.R_targs`.

# Arguments
- `params`: HALO parameters; allocated targets matched by `params.a.ID_targs`.
- `target_pool`: mutable pool of [`LandingTarget`](@ref) radii to perturb.
- `noise_std`: standard deviation of Gaussian radius noise (default `0.2`).
- `crossweight`: cross-fade weight retaining the previous radius (default `0.05`).

# Returns
- none

# Notes
Mutates `target_pool` and `params`.
"""
function sim_update_targets!(params, target_pool::Vector{LandingTarget}; noise_std::CReal=.2, crossweight::CReal=.05)
    # Update bounding radii of currently locked targets
    # (Not updating any other parameters currently)
    for k = 1:length(target_pool)
        R_ = add_gauss([target_pool[k].R], noise_std, 0.0, clip=false)[1]
        if target_pool[k].R > params.R_targs_min
            if R_ > target_pool[k].R
                R_ = crossweight * target_pool[k].R + (1-crossweight) * R_
            end
        else
            if R_ < target_pool[k].R
                R_ = crossweight * target_pool[k].R + (1-crossweight) * R_
            end
        end
        target_pool[k].R = max(R_,0.)
        if target_pool[k].id in params.a.ID_targs
            idx = findfirst(params.a.ID_targs .== target_pool[k].id)
            params.R_targs[idx] = target_pool[k].R
        end
    end
end

"""
    sim_refresh_targets!(params, target_pool)

Fill `params` up to `n_targs_max` targets from `target_pool`, keeping currently
allocated sites, computing desirability scores, and setting `λ_targs`.

# Arguments
- `params`: HALO parameters resized to `n_targs_max` active targets.
- `target_pool`: candidate [`LandingTarget`](@ref) pool to draw new sites from.

# Returns
- none

# Notes
Mutates `params`.
"""
function sim_refresh_targets!(params, target_pool::Vector{LandingTarget})
    # Require that we have more targets to choose from than is needed
    if length(target_pool) < params.n_targs_max
        error("Not enough targets to choose from.")
    end

    # Copy current remaining targets as the old targets
    R_targs_old = copy(params.R_targs)
    ID_targs_old = copy(params.a.ID_targs)
    zf_targs_old = copy(params.a.zf_targs)
    if params.n_targs_max == 1 # remove the remaining target anyways
        zf_targs_old = []
    end
    if !isempty(zf_targs_old)
        rf_targs_old = zf_targs_old[1:3,:]
    else
        rf_targs_old = zeros()
    end

    # Acquire `n_missing` best available targets from the available pool
    n_missing = params.n_targs_max - params.a.n_targs # (N_max - N_current)
    rf_targs_new = zeros(3, n_missing)
    R_targs_new = zeros(n_missing)
    ID_targs_new = zeros(n_missing)
    pool_pref_indices = sortperm([target_pool[j].R for j=1:length(target_pool)], rev=true)
    cur_pool_idx = 1
    for k = 1:n_missing
        while true
            # Randomly select a target from the target pool
            new_targ = target_pool[pool_pref_indices[cur_pool_idx]]

            # Check if the target is already in the old targets by target ID (facilitates target reuse)
            if !isempty(ID_targs_old)
                if any(new_targ.id .== ID_targs_old)
                    cur_pool_idx += 1
                    continue
                end
            end

            # If the target is unique, add it to the new targets
            ID_targs_new[k] = new_targ.id
            R_targs_new[k] = new_targ.R
            rf_targs_new[:,k] = new_targ.rf
            cur_pool_idx += 1
            break
        end
    end
    
    # Concatenate to create list of all targets
    R_targs = vcat(R_targs_old, R_targs_new)
    ID_targs = vcat(ID_targs_old, ID_targs_new)
    rf_targs = !isempty(zf_targs_old) ? hcat(rf_targs_old, rf_targs_new) : rf_targs_new
    vf_targs = zeros(3, params.n_targs_max)
    ∫Tf_targs = Inf * ones(1,params.n_targs_max)
    zf_targs = vcat(rf_targs, vf_targs, ∫Tf_targs)

    # Compute parameters of interest
    # Point cloud density
    pcd_targs = fill(10, params.n_targs_max) # Only simulated in AirSim

    # 99th-percentile uncertainty
    µ_99_targs = fill(0, params.n_targs_max) # Only simulated in AirSim

    # Proximity to vehicle
    prox_veh_targs = [norm(rf_targs[1:2,i] - params.a.z0[1:2], 2) for i=1:params.n_targs_max]

    # Proximity to cluster
    prox_clust_targs = zeros(params.n_targs_max)
    for i = 1:params.n_targs_max
        sum = 0
        for j = 1:params.n_targs_max
            sum += norm(rf_targs[1:2,i] - rf_targs[1:2,j], 2)
        end
        prox_clust_targs[i] = sum
    end
    prox_clust_targs ./= sum(prox_clust_targs)

    # Compute confidence score for each target
    des_score = zeros(params.n_targs_max)
    for j = 1:params.n_targs_max
        des_score[j] = 
            pcd_targs[j] * params.w_des[1] +
            prox_veh_targs[j] * params.w_des[2] +
            prox_clust_targs[j] * params.w_des[3] +
            µ_99_targs[j] * params.w_des[4] +
            R_targs[j] * params.w_des[5]
    end
    
    # Sort target preference order by confidence score
    λ_targs = sortperm(des_score)

    # Add all values to the params
    params.a.n_targs    = params.n_targs_max
    params.a.zf_targs   = zf_targs
    params.a.uf_targs   = Inf * ones(params.a.nu,params.a.n_targs)
    params.a.λ_targs    = λ_targs
    params.a.ID_targs   = ID_targs
    params.a.J_targs    = 1:params.n_targs_max
    params.a.τ_targs    = zeros(params.a.n_targs)
    params.a.α_targs    = ones(params.a.n_targs)
    params.a.ϵ_targs    = fill(params.ϵ_subopt, params.a.n_targs)
    params.a.w_obj_ddto = params.a.w_obj_sing / params.a.n_targs
    params.R_targs      = R_targs
    params.p_targs["pcd"] = pcd_targs
    params.p_targs["prox_veh"] = prox_veh_targs
    params.p_targs["prox_clust"] = prox_clust_targs
    params.p_targs["µ_99"] = µ_99_targs
end

"""
    generate_obstacles!(params, n_obstacles, obs_rad_position, obs_rad)

Place `n_obstacles` equally spaced cylindrical obstacles of radius `obs_rad`
on a circle of radius `obs_rad_position` about the origin.

# Arguments
- `params`: HALO parameters whose obstacle fields are populated.
- `n_obstacles`: number of cylindrical obstacles to place.
- `obs_rad_position`: horizontal placement radius `[m]` for obstacle centers.
- `obs_rad`: obstacle cylinder radius `[m]`.

# Returns
- none

# Notes
Mutates `params` obstacle fields.
"""
function generate_obstacles!(params, n_obstacles, obs_rad_position, obs_rad)
    params.n_obstacles = n_obstacles
    params.R_obstacles = fill(obs_rad, params.n_obstacles)
    θ_sep = 2 * pi / params.n_obstacles
    θ_obs = 0.
    params.p_obstacles = CMatrix(undef, 3, params.n_obstacles)
    for j = 1:params.n_obstacles
        r_obs = obs_rad_position
        params.p_obstacles[:,j] = r_obs*cos(θ_obs)*e_x + r_obs*sin(θ_obs)*e_y + 1*e_z
        θ_obs += θ_sep
    end
    obs_shape = 1.0I(3)
    obs_shape[3,3] = 0
    params.H_obstacles = repeat([obs_shape], params.n_obstacles)
end