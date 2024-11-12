#= Adaptive-DDTO -- Utility functions specific to ADDTO.
Author: Samuel Buckner (UW-ACL)
=#

function extract_trunk_segment(params, sol::Quad3DoFDDTOSolution; sim::Bool=false)::Quad3DoFSolution
    """
    Extract the trunk (deferrable) segment of a full DDTO solution.

    Args:
        params (params): The params object.
        sol (DDTOSolution): The DDTO solution bundle.

    Returns:
        sol_trunk (Solution): Container for trunk (deferrable segment) solution.

    TODO: Figure out what to do with this, deprecated function at this point.
    """

    if params.a.n_targs > 1
        τ_cutoff = params.a.τ_targs[findfirst(i->i==params.a.λ_targs[end-1], params.a.λ_targs)]
        idx = params.a.λ_targs[end-1]
    else
        τ_cutoff = params.a.τ_targs[1]
        idx = 1
    end
    if sim
        τ_cutoff = (τ_cutoff - 1) * params.a.N_sim + 1
    end

    τ_trunk = sol.targs[idx].τ[1:τ_cutoff]
    t_trunk = sol.targs[idx].t[1:τ_cutoff]
    x_trunk = sol.targs[idx].x[:,1:τ_cutoff]
    u_trunk = sol.targs[idx].u[:,1:τ_cutoff]
    r_trunk = sol.targs[idx].r[:,1:τ_cutoff]
    v_trunk = sol.targs[idx].v[:,1:τ_cutoff]
    T_trunk = sol.targs[idx].T[:,1:τ_cutoff]
    s_trunk = sol.targs[idx].s[1:τ_cutoff]
    T_nrm_trunk = sol.targs[idx].T_nrm[1:τ_cutoff]
    ∫T_trunk = sol.targs[idx].∫T[1:τ_cutoff]
    γ_trunk = sol.targs[idx].γ[1:τ_cutoff]
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

function extract_segment(sol_ddto::Quad3DoFDDTOSolution, defer_targ::Int, λ_targs_org::Vector{Int})::Quad3DoFSolution
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

    if params.a.n_targs > 0
        # Determine indices for removal
        pop_idx_J = findfirst(i->i==T_targ, params.a.J_targs)
        pop_idx_λ = findfirst(i->i==T_targ, params.a.λ_targs)
        slice_J = collect(1:params.a.n_targs)
        slice_λ = collect(1:params.a.n_targs)
        deleteat!(slice_J, pop_idx_J)
        deleteat!(slice_λ, pop_idx_λ)

        # Parameter updates for removing the target
        params.a.n_targs -= 1
        params.a.λ_targs  = params.a.λ_targs[slice_λ]
        params.a.ID_targs = params.a.ID_targs[slice_J]
        params.a.J_targs  = params.a.J_targs[slice_J]
        params.a.ϵ_targs  = params.a.ϵ_targs[slice_J]
        params.a.τ_targs  = params.a.τ_targs[slice_λ]
        params.a.α_targs  = params.a.α_targs[slice_J]
        params.a.zf_targs = params.a.zf_targs[:,slice_J]
        params.a.uf_targs = params.a.uf_targs[:,slice_J]
        params.a.w_obj_ddto = params.a.w_obj_sing / params.a.n_targs
        params.R_targs = params.R_targs[slice_J]
        for (key,~) in params.p_targs
            params.p_targs[key] = params.p_targs[key][slice_J]
        end
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
    other_targs = copy(params.a.J_targs)
    deleteat!(other_targs, findfirst(i->i==branch_targ, other_targs))

    # Get indices for each target
    other_targ_idx = Array{Int}(undef, length(other_targs))
    for j = 1:length(other_targ_idx)
        other_targ_idx[j] = findfirst(i->i==other_targs[j], params.a.J_targs)
    end
    branch_targ_idx = findfirst(i->i==branch_targ, params.a.J_targs)

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

function setup_addto_dicts(params)
    # Guidance
    guid = Dict()
    guid["cur_opt"]      = EmptyQuad3DoFDDTOSolution(params.a.n_targs) # Most recently-computed optimal solution set
    guid["cur_ddto"]     = EmptyQuad3DoFDDTOSolution(params.a.n_targs) # Most recently-computed DDTO solution set
    guid["cur_ddto_sim"] = EmptyQuad3DoFDDTOSolution(params.a.n_targs) # Most recently-computed DDTO simulation set
    guid["cur_traj"]     = EmptyQuad3DoFSolution() # Current guidance solution to track
    guid["cur_traj_sim"] = EmptyQuad3DoFSolution() # Current guidance solution to track
    guid["cur_time"]     = 0.0 # Current time in guidance solution
    guid["defer_targ"]   = -1 # Next deferred target in consideration (tag number)
    guid["defer_time"]   = 1.e6 # Time until branch point to next deferred target
    guid["lock_time"]    = 1.e6 # Time at which guidance lock was activated
    guid["λ_targs_org"]  = params.a.λ_targs # Stores initial preference ordering
    guid["comp_params"]  = Quad3DoFHaloParams()

    # Flags
    flags = Dict()
    flags["update_ddto"]           = true
    flags["ddto_converged"]        = false
    flags["log_ddto_results"]      = false # If set to true, log DDTO results
    flags["guid_lock_activated"]   = false # If set to true, Adaptive-DDTO will be disabled and guidance will fix to the best target at the current time
    flags["descent_complete"]      = false # If set to true, signals the end of the simulation/descent phase
    flags["guid_lock_staged"]      = false # If set to true, stage a guidance lock
    flags["guid_recently_updated"] = false

    # Results (to be logged)
    results = Dict()
    results["guid_update_ddto_params"]       = Array{Quad3DoFHaloParams}(undef,0)
    results["guid_update_ddto_bundles"]      = Array{Quad3DoFDDTOSolution}(undef,0)
    results["guid_update_ddto_bundles_sims"] = Array{Quad3DoFDDTOSolution}(undef,0)
    results["guid_update_trajs"]             = Array{Quad3DoFSolution}(undef, 0)
    results["guid_update_trajs_sims"]        = Array{Quad3DoFSolution}(undef, 0)
    results["guid_update_time"]              = CVector(undef, 0)
    results["sim_time"]                      = CVector(undef, 0)
    results["sim_state"]                     = CMatrix(undef, params.a.nx, 0)
    results["sim_control"]                   = CMatrix(undef, params.a.nu, 0)
    results["targs_ID"]                      = Matrix{Int}(undef, params.n_targs_max, 0)
    results["targs_radii"]                   = CMatrix(undef, params.n_targs_max, 0)
    results["targs_status"]                  = Matrix{Bool}(undef, params.n_targs_max, 0)
    results["targs_positions"]               = Array{CMatrix}(undef, 0)
    results["targpool_ID"]                   = Vector{Int}(undef, 0) # filled only once
    results["targpool_positions"]            = CMatrix(undef, 0, 0) # filled only once
    results["targpool_radii"]                = CMatrix(undef, 0, 0) 
    results["targpool_allocated"]            = Matrix{Bool}(undef, 0, 0)

    return guid,flags,results
end

function log_results!(params, results::Dict, guid::Dict, flags::Dict, sim_cur_state::Vector{Float64}, sim_cur_control::Vector{Float64}, sim_cur_time::Float64; target_pool::Vector=[])
    """
    Logs results
    """
    hcat_c = (A,x) -> length(A) == 0 ? x : hcat(A,x)

    # Log continuous sim results
    results["sim_state"]   = hcat_c(results["sim_state"], sim_cur_state)
    results["sim_control"] = hcat_c(results["sim_control"], sim_cur_control)
    append!(results["sim_time"], sim_cur_time)

    # Log current target ID (if a target index is unallocated, insert 0)
    sim_cur_targ_id = zeros(Int,params.n_targs_max)
    sim_cur_targ_id[params.a.J_targs] = params.a.ID_targs
    results["targs_ID"] = hcat_c(results["targs_ID"], sim_cur_targ_id)

    # Log current target radii (if a target index is unallocated, insert -Inf)
    sim_cur_radii = fill(-Inf, params.n_targs_max)
    sim_cur_radii[params.a.J_targs] = params.R_targs
    results["targs_radii"] = hcat_c(results["targs_radii"], sim_cur_radii)
    
    # Log current target positions (if a target index is unallocated, insert -Inf)
    sim_cur_targ_pos = -Inf * ones(3, params.n_targs_max)
    sim_cur_targ_pos[:,params.a.J_targs] = params.a.zf_targs[1:3,:]
    append!(results["targs_positions"], [sim_cur_targ_pos])

    # Log target status (1 = valid, 0 = lost)
    targs_status = zeros(Bool,params.n_targs_max)
    for k=1:params.n_targs_max
        if k in params.a.J_targs
            targs_status[k] = 1
        end
    end
    results["targs_status"] = hcat_c(results["targs_status"], targs_status)
    
    # Log conditional sim results (DDTO)
    if flags["log_ddto_results"]
        append!(results["guid_update_ddto_params"], [guid["comp_params"]])
        append!(results["guid_update_ddto_bundles"], [guid["cur_ddto"]])
        append!(results["guid_update_ddto_bundles_sims"], [guid["cur_ddto_sim"]])
        append!(results["guid_update_trajs"], [guid["cur_traj"]])
        append!(results["guid_update_trajs_sims"], [guid["cur_traj_sim"]])
        append!(results["guid_update_time"], sim_cur_time)
        flags["log_ddto_results"] = false
    end

    # Log target pool status if pool is provided
    if target_pool != []
        if length(results["targpool_ID"]) == 0 # only write first time
            results["targpool_ID"] = [t.id for t in target_pool]
        end
        if length(results["targpool_positions"]) == 0 # only write first time
            results["targpool_positions"] = mapreduce(permutedims, vcat, [t.rf for t in target_pool])'
        end
        results["targpool_radii"] = hcat_c(results["targpool_radii"], [t.R for t in target_pool])
        results["targpool_allocated"] = hcat_c(results["targpool_allocated"], [t.id in params.a.ID_targs for t in target_pool])
    end

    return results, flags
end