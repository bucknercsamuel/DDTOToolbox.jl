#= Adaptive-DDTO -- Utility functions to simulate perception interface.
Author: Samuel Buckner (UW-ACL)
=#

function sim_acquire_new_targets!(params, R_landing_region::CReal)
    """
    Simulate the acquisition of new targets from the perception stack.
    * NOTE: This function will modify the params object.

    Args:
        params (any): The parameter object.
        R_landing_region (CReal): The radius of interest.
    """

    # Copy current remaining targets as the old targets
    R_targs_old = copy(params.R_targs)
    zf_targs_old = copy(params.a.zf_targs)
    if !isempty(zf_targs_old)
        rf_targs_old = zf_targs_old[1:3,:]
    else
        rf_targs_old = zeros()
    end

    # Generate `n_missing` new targets
    n_missing = params.n_targs_max - params.a.n_targs # (N_max - N_current)
    (R_targs_new, rf_targs_new) = sim_generate_random_targets(params, n_missing, R_landing_region)

    # Concatenate to create list of all targets
    R_targs = vcat(R_targs_old, R_targs_new)
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
    params.a.T_targs    = 1:params.a.n_targs
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

function sim_update_locked_targets!(params; noise_std::CReal=0.2)
    """
    Simulate the update of locked target parameters from the perception stack.
    * NOTE: This function will modify the params object.

    Args:
        params (any): The params object.
    """

    # Update bounding radii of currently locked targets
    # (Not updating any other parameters currently)
    params.R_targs = add_gauss(params.R_targs, noise_std, 0.0, clip=false)
    for k = 1:length(params.R_targs)
        params.R_targs[k] = max(params.R_targs[k], 0)
    end
end

function sim_generate_random_targets(params, N::Int, R_landing_region::CReal)::Tuple{CVector,CMatrix}
    """
    Generate a set of random targets in a landing region.

    Args:
        params (any): The params object.
        N (Int): Number of targets to generate.
        R_landing_region (CReal): Radius of the landing region.

    Returns:
        R_targs (CVector): Radii of each generated target.
        rf_targs (CMatrix): Position of each generated target.
    """

    # Get random target bounding radii at some amount larger than the minimum radius threshold
    rand_vals = rand(CReal, N)
    R_targs = CVector(params.R_targs_min .+ 5 * rand_vals)

    # Get random target positions uniformly dispersed in a radius of `R_landing_region`
    rf_targs = CMatrix(undef, 3, N)
    for j = 1:N
        r_targ = R_landing_region * rand(CReal)
        θ_targ = 2 * pi * rand(CReal)
        rf_targs[:,j] = r_targ*cos(θ_targ)*e_x + r_targ*sin(θ_targ)*e_y + 1*e_z
    end

    return (R_targs, rf_targs)
end