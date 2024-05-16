#= Adaptive-DDTO -- Utility functions specific to ADDTO.
Author: Samuel Buckner (UW-ACL)
=#

function extract_trunk_segment(params, sol::DDTOSolution)::Solution
    """
    Extract the trunk (deferrable) segment of a full DDTO solution.

    Args:
        params (params): The params object.
        sol (DDTOSolution): The DDTO solution bundle.

    Returns:
        sol_trunk (Solution): Container for trunk (deferrable segment) solution.
    """
    τ_cutoff = params.a.τ_targs[findfirst(i->i==params.a.λ_targs[end-1], params.a.λ_targs)] # obtain the deferrability index of the j-th target (solution)

    t_trunk = sol.targs[[params.a.λ_targs[end-1]]].t[1:τ_cutoff]
    x_trunk = sol.targs[[params.a.λ_targs[end-1]]].x[1:τ_cutoff]
    u_trunk = sol.targs[[params.a.λ_targs[end-1]]].u[1:τ_cutoff]
    cost_trunk = nothing # TODO: should we find a way to set the cost of the trunk?

    sol_trunk = Solution(t_trunk, x_trunk, u_trunk, cost_trunk)
    return sol_trunk
end

function extract_guid_lock_traj(params, sol_ddto::DDTOSolution, λ_defer_idx, λ_targs_org::Vector{Int})::Solution
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
    pop_idx_T = findfirst(i->i==T_targ, params.T_targs)
    pop_idx_λ = findfirst(i->i==T_targ, params.λ_targs)
    slice_T = collect(1:params.n_targs)
    slice_λ = collect(1:params.n_targs)
    deleteat!(slice_T, pop_idx_T)
    deleteat!(slice_λ, pop_idx_λ)

    # Parameter updates for removing the target
    params.n_targs -= 1
    params.λ_targs = params.λ_targs[slice_λ]
    params.T_targs = params.T_targs[slice_T]
    params.N_targs = params.N_targs[slice_T]
    params.R_targs = params.R_targs[slice_T]
    params.ϵ_targs = params.ϵ_targs[slice_T]
    params.rf_targs = params.rf_targs[:,slice_T]
    params.vf_targs = params.vf_targs[:,slice_T]
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
    other_targs = copy(params.T_targs)
    deleteat!(other_targs, findfirst(i->i==branch_targ, other_targs))

    # Get indices for each target
    other_targ_idx = Array{Int}(undef, length(other_targs))
    for j = 1:length(other_targ_idx)
        other_targ_idx[j] = findfirst(i->i==other_targs[j], params.T_targs)
    end
    branch_targ_idx = findfirst(i->i==branch_targ, params.T_targs)

    # Check if desirability score of branch targ is greater than that of all other targs
    des_score = zeros(params.n_targs)
    for j = 1 : params.n_targs
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