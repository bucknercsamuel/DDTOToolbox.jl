#= Adaptive-DDTO -- Utility functions specific to ADDTO.
Author: Samuel Buckner (UW-ACL)
=#

function extract_trunk_segment(params, sol::Quad3DoFDDTOSolution)::Quad3DoFSolution
    """
    Extract the trunk (deferrable) segment of a full DDTO solution.

    Args:
        params (params): The params object.
        sol (DDTOSolution): The DDTO solution bundle.

    Returns:
        sol_trunk (Solution): Container for trunk (deferrable segment) solution.
    """
    τ_cutoff = params.a.τ_targs[findfirst(i->i==params.a.λ_targs[end-1], params.a.λ_targs)] # obtain the deferrability index of the j-th target (solution)

    τ_trunk = sol.targs[params.a.λ_targs[end-1]].τ[1:τ_cutoff]
    t_trunk = sol.targs[params.a.λ_targs[end-1]].t[1:τ_cutoff]
    x_trunk = sol.targs[params.a.λ_targs[end-1]].x[:,1:τ_cutoff]
    u_trunk = sol.targs[params.a.λ_targs[end-1]].u[:,1:τ_cutoff]
    r_trunk = sol.targs[params.a.λ_targs[end-1]].r[:,1:τ_cutoff]
    v_trunk = sol.targs[params.a.λ_targs[end-1]].v[:,1:τ_cutoff]
    T_trunk = sol.targs[params.a.λ_targs[end-1]].T[:,1:τ_cutoff]
    s_trunk = sol.targs[params.a.λ_targs[end-1]].s[1:τ_cutoff]
    T_nrm_trunk = sol.targs[params.a.λ_targs[end-1]].T_nrm[1:τ_cutoff]
    ∫T_trunk = sol.targs[params.a.λ_targs[end-1]].∫T[1:τ_cutoff]
    γ_trunk = sol.targs[params.a.λ_targs[end-1]].γ[1:τ_cutoff]
    cost_trunk = -1 # TODO: should we find a way to set the cost of the trunk?

    sol_trunk = Quad3DoFSolution(
        τ_trunk,
        t_trunk, 
        x_trunk, 
        u_trunk,
        r_trunk,
        v_trunk,
        T_trunk,
        s_trunk,
        T_nrm_trunk,
        ∫T_trunk,
        γ_trunk,
        cost_trunk)
    return sol_trunk
end

function extract_guid_lock_segment(sol_ddto::Quad3DoFDDTOSolution, defer_targ::Int, λ_targs_org::Vector{Int})::Quad3DoFSolution
    """
    Extract the guidance-locked segment of a full DDTO solution.

    Args:
        params (params): The params object.
        ddto_sol (DDTOSolution): Vectorized container for all DDTO solutions.
        λ_defer_idx (Int): λ (preference order) vector index for the target at which the deferral occurs.
        λ_targs_org (Vector{Int}): Original target preference order when DDTO was originally computed.

    Returns:
        sol_guid (Solution): Container for guidance-locked solution.
    """
    λ_defer_idx = findfirst(i->i==defer_targ, λ_targs_org)
    return sol_ddto.targs[λ_targs_org[λ_defer_idx]]
end

function remove_ddto_target!(params, T_targ::Int)
    """
    Remove a target from the vehicle parameters.
    * NOTE: This function will modify the params object.

    Args:
        params (any): The parameter object.
        T_targ (Int): Target tag for the target to be removed.
    """

    # Determine indices for removal
    pop_idx_T = findfirst(i->i==T_targ, params.a.T_targs)
    pop_idx_λ = findfirst(i->i==T_targ, params.a.λ_targs)
    slice_T = collect(1:params.a.n_targs)
    slice_λ = collect(1:params.a.n_targs)
    deleteat!(slice_T, pop_idx_T)
    deleteat!(slice_λ, pop_idx_λ)

    # Parameter updates for removing the target
    params.a.n_targs -= 1
    params.a.λ_targs = params.a.λ_targs[slice_λ]
    params.a.T_targs = params.a.T_targs[slice_T]
    params.a.ϵ_targs = params.a.ϵ_targs[slice_T]
    params.a.zf_targs = params.a.zf_targs[:,slice_T]
    params.R_targs = params.R_targs[slice_T]
    for (key,~) in params.p_targs
        params.p_targs[key] = params.p_targs[key][slice_T]
    end
end

function switch_decision(params, branch_targ::Int)::Bool
    """
    Determine the switch decision at a branch point.

    Args:
        params (any): The parameter object.
        branch_targ (Int): Target tag for the target being considered for deferral at the branch point.

    Returns:
        switch_decision (Bool): True if deferral to branch_targ should take place, false otherwise.
    """

    # Find the other targs that aren't the branch targ
    other_targs = copy(params.a.T_targs)
    deleteat!(other_targs, findfirst(i->i==branch_targ, other_targs))

    # Get indices for each target
    other_targ_idx = Array{Int}(undef, length(other_targs))
    for j = 1:length(other_targ_idx)
        other_targ_idx[j] = findfirst(i->i==other_targs[j], params.a.T_targs)
    end
    branch_targ_idx = findfirst(i->i==branch_targ, params.a.T_targs)

    # Check if desirability score of branch targ is greater than that of all other targs
    des_score = zeros(params.a.n_targs)
    for j = 1 : params.a.n_targs
        des_score[j] = 
            params.p_targs["pcd"][j] * params.w_des[1] + 
            params.p_targs["prox_veh"][j] * params.w_des[2] + 
            params.p_targs["prox_clust"][j] * params.w_des[3] + 
            params.p_targs["µ_99"][j] * params.w_des[4] + 
            params.R_targs[j] * params.w_des[5]
    end

    if des_score[branch_targ_idx] > maximum(des_score[other_targ_idx])
        switch = true
    else
        switch = false
    end

    return switch
end

function log_results!(params, results::Dict, guid::Dict, flags::Dict, sim_cur_state::Vector{Float64}, sim_cur_time::Float64)
    """
    Logs results
    """
    # Obtain current control for logging
    if flags["ddto_converged"]
        sim_cur_control = optimal_controller(guid["cur_time"], guid["cur_traj"].t, guid["cur_traj"].u, params.a.disc)
    else
        sim_cur_control = zeros(params.a.nu+1)
    end
    
    # Log continuous sim results
    results["sim_state"]   = hcat(results["sim_state"], sim_cur_state)
    results["sim_control"] = hcat(results["sim_control"], sim_cur_control)
    append!(results["sim_time"], sim_cur_time)

    # Log current target radii (if a target index is unallocated, insert -Inf)
    sim_cur_radii = fill(-Inf, params.n_targs_max)
    sim_cur_radii[params.a.T_targs] = params.R_targs
    results["targs_radii"] = hcat(results["targs_radii"], sim_cur_radii)
    
    # Log current target positions (if a target index is unallocated, insert -Inf)
    sim_cur_targ_pos = -Inf * ones(3, params.n_targs_max)
    sim_cur_targ_pos[:,params.a.T_targs] = params.a.zf_targs[1:3,:]
    append!(results["targs_positions"], [sim_cur_targ_pos])

    # Log target status (1 = valid, 0 = lost)
    targs_status = zeros(params.n_targs_max)
    for k=1:params.n_targs_max
        if k in params.a.T_targs
            targs_status[k] = 1
        end
    end
    results["targs_status"] = hcat(results["targs_status"], targs_status)
    
    # Log conditional sim results (DDTO)
    if flags["log_ddto_results"]
        append!(results["guid_update_ddto_bundles"], [guid["cur_ddto"]])
        append!(results["guid_update_trajs"], [guid["cur_traj"]])
        append!(results["guid_update_time"], sim_cur_time)
        flags["log_ddto_results"] = false
    end
end