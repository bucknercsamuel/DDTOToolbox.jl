#= Adaptive-DDTO -- Utility functions to simulate perception interface.
Author: Samuel Buckner (UW-ACL)
=#

mutable struct LandingTarget
    id::Int     # Target ID
    R::CReal    # Target radius
    rf::CVector # Target position
end

function sim_build_target_pool(num_targets::Int, R_landing_region::CReal; min_radius::CReal, max_radius::CReal)::Vector{LandingTarget}
    """
    Generate a set of random targets in a landing region.

    Args:
        num_targets (Int): Number of targets to generate.
        R_landing_region (CReal): Radius of the landing region.
        min_radius (CReal): Minimum radius of the targets.
        max_radius (CReal): Maximum radius of the targets.

    Returns:
        target_pool (Vector{LandingTarget}): A vector of LandingTarget objects.
    """

    # Sample target radii randomly between min and max radius
    rand_vals = rand(CReal, num_targets)
    R_targs = CVector(min_radius .+ (max_radius-min_radius) * rand_vals)

    # Get random target positions uniformly dispersed in a radius of `R_landing_region`
    rf_targs = CMatrix(undef, 3, num_targets)
    for j = 1:num_targets
        r_targ = R_landing_region * rand(CReal)
        θ_targ = 2 * pi * rand(CReal)
        rf_targs[:,j] = r_targ*cos(θ_targ)*e_x + r_targ*sin(θ_targ)*e_y + 1*e_z
    end

    # Create the target pool
    target_pool = [LandingTarget(j, R_targs[j], rf_targs[:,j]) for j=1:num_targets]
    return target_pool
end

function sim_refresh_targets!(params, target_pool::Vector{LandingTarget})
    """
    Acquire params.n_targs_max targets from the available pool of targets while maintaining old ones.
    * NOTE: This function will modify the params object.

    Args:
        params (any): The parameter object.
        R_landing_region (CReal): The radius of interest.
    """

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

function sim_update_targets!(params, target_pool::Vector{LandingTarget}; noise_std::CReal=0.2)
    """
    Simulate the update of locked target parameters from the perception stack.
    * NOTE: This function will modify the target_pool and params objects.

    Args:
        params (any): The parameter object.
        target_pool::Vector{LandingTarget}: target pool
    """

    # Update bounding radii of currently locked targets
    # (Not updating any other parameters currently)
    for k = 1:length(target_pool)
        R_ = add_gauss([target_pool[k].R], noise_std, 0.0, clip=false)[1]
        target_pool[k].R = max(R_,0.)
        if target_pool[k].id in params.a.ID_targs
            idx = findfirst(params.a.ID_targs .== target_pool[k].id)
            params.R_targs[idx] = target_pool[k].R
        end
    end
end