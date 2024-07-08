#= Adaptive-DDTO -- Core algorithm functions.
Author: Samuel Buckner (UW-ACL)
=#

function compute_ddto_guidance!(params, guid::Dict, flags::Dict, sim_cur_state::Vector{Float64}, sim_cur_control::Vector{Float64}, sim_cur_time::Float64)
    """
    Compute new DDTO guidance tree (if staged to do so)
    """

    # Set guidance initial state as current sim state
    for k = 1:6
        params.a.z0[k] = sim_cur_state[k]
    end

    # Saturate velocities to satisfy constraints
    vel_lat_idx = 4:5
    vel_vert_idx = 6
    if norm(params.a.z0[vel_lat_idx]) > params.v_max_L
        params.a.z0[vel_lat_idx] = params.a.z0[vel_lat_idx] / norm(params.a.z0[vel_lat_idx]) * params.v_max_L
        display("Warning: velocity saturated!")
    end
    params.a.z0[vel_vert_idx] = max(min(params.a.z0[vel_vert_idx], params.v_max_V), -params.v_max_V)

    # Set guidance initial control to current sim control
    for k = 1:3
        params.a.u0[k] = sim_cur_control[k]
    end

    # Saturate thrusts to satisfy constraints
    thrust_idx = 1:3
    if norm(params.a.u0[thrust_idx]) > params.ρ_max && ~isinf(norm(params.a.u0[thrust_idx]))
        params.a.u0[thrust_idx] = params.a.u0[thrust_idx] / norm(params.a.u0[thrust_idx]) * params.ρ_max
    elseif norm(params.a.u0[thrust_idx]) < params.ρ_min
        params.a.u0[thrust_idx] = params.a.u0[thrust_idx] / norm(params.a.u0[thrust_idx]) * params.ρ_min
    end

    # Guidance solving
    flags["ddto_converged"] = false
    # try
        _,_,guid["cur_ddto"],guid["cur_ddto_sim"],flags["ddto_converged"] = solve(params) # Compute DDTO solution
        guid["comp_params"] = copy(params)
        flags["ddto_converged"] = true # TODO: remove once DDTO converges properly
    # catch e
    #     @printf("  -> DDTO ERROR [%.2f s]: %s\n", sim_cur_time, e)
    # end
    if !flags["ddto_converged"]
        @printf("  -> UPDATE [%.2f s]: Guidance lock staged [DDTO computation unsuccessful -- contingency activated!]\n", sim_cur_time)
        flags["guid_lock_staged"] = true
    end

    if !flags["guid_lock_staged"]
        # guid["cur_traj"] = extract_trunk_segment(params, guid["cur_ddto"]) # Track the trunk of DDTO by default
        # guid["cur_traj_sim"] = extract_trunk_segment(params, guid["cur_ddto_sim"], sim=true) # Track the trunk of DDTO by default
        guid["cur_traj"] = extract_segment(guid["cur_ddto"], params.a.λ_targs[end], params.a.λ_targs)
        guid["cur_traj_sim"] = extract_segment(guid["cur_ddto_sim"], params.a.λ_targs[end], params.a.λ_targs)
        @printf("  -> UPDATE [%.2f s]: DDTO solution successfully recomputed [tracking trunk segment]\n", sim_cur_time)

        # Parameter updates
        guid["cur_time"] = 0.0 # Reset guidance time to zero
        guid["defer_targ"] = params.a.λ_targs[1]
        guid["defer_time"] = guid["cur_ddto"].targs[guid["defer_targ"]].t[params.a.τ_targs[1]]
        guid["defer_state"] = guid["cur_ddto"].targs[guid["defer_targ"]].x[:,params.a.τ_targs[1]]
        guid["λ_targs_org"] = params.a.λ_targs

        # If trunk segment has zero length (no deferring could take place),
        # lock guidance to the best target at the current point in time (last index of last DDTO branch solution)
        # as a contingency measure
        if length(guid["cur_traj"].t) == 0
            @printf("  -> UPDATE [%.2f s]: Guidance lock staged [DDTO deferral was not possible -- contingency activated!]\n", sim_cur_time)
            flags["guid_lock_staged"] = true
        end
    end

    # Flag updates
    flags["update_ddto"] = false
    flags["log_ddto_results"] = true

    return guid, flags
end

function check_unsafe_targets!(params, guid::Dict, flags::Dict, sim_cur_time::Float64)
    """
    Check for unsafe targets (radii check)
    """
    cur_targs = copy(params.a.T_targs)
    for targ in cur_targs
        targ_idx = findfirst(i->i==targ, params.a.T_targs)

        # Remove target if unsafe
        if params.R_targs[targ_idx] <= params.R_targs_min
            @printf("  -> UPDATE [%.2f s]: Removing target %i [bounding radius below the minimum threshold]\n", sim_cur_time, targ)
            remove_ddto_target!(params, targ)

            # If this target was queued for deferral, move to next target for deferral
            if targ == guid["defer_targ"] && params.a.n_targs > 0
                guid["defer_targ"] = params.a.λ_targs[1] # Add the next target in the queue to consideration for deferral
                guid["defer_time"] = guid["cur_ddto"].targs[guid["defer_targ"]].t[params.a.τ_targs[1]]
                guid["defer_state"] = guid["cur_ddto"].targs[guid["defer_targ"]].x[:,params.a.τ_targs[1]]
            end
        end

        # Reached minimum target threshold
        if params.a.n_targs < params.n_targs_min
            ~flags["guid_lock_staged"] && @printf("  -> UPDATE [%.2f s]: DDTO recomputation staged [target set count below the minimum threshold]\n", sim_cur_time)
            flags["update_ddto"] = true
            break
        end
    end

    return guid, flags
end

function check_branch_switch!(params, guid::Dict, flags::Dict, sim_cur_state::Vector{Float64}, sim_cur_time::Float64; criteria::String="time")
    """
    Check for branch switching decision
    """
    if criteria == "time"
        criterion = guid["cur_time"] >= guid["defer_time"]
    elseif criteria == "altitude"
        criterion = sim_cur_state[3] <= guid["defer_state"][3]
    end

    if criterion
        while criterion
            # Determine if we should switch or not
            switch_branch = switch_decision(params, guid["defer_targ"])

            # Engage switch by staging DDTO update
            if switch_branch
                ~flags["guid_lock_staged"] && @printf("  -> UPDATE [%.2f s]: DDTO recomputation staged [chose to defer to target %i]\n", sim_cur_time, guid["defer_targ"])

                # Remove all targets except for switch target (`guid["defer_targ"]`)
                other_targs = copy(params.a.T_targs)
                deleteat!(other_targs, findfirst(i->i==guid["defer_targ"], other_targs))
                for targ ∈ other_targs
                    remove_ddto_target!(params, targ)
                end

                flags["update_ddto"] = true
                break

            # Remove the current target for deferral and go to the next one
            else
                @printf("  -> UPDATE [%.2f s]: Removing target %i [chose to stay on trunk segment]\n", sim_cur_time, guid["defer_targ"])
                remove_ddto_target!(params, guid["defer_targ"]) # Remove the target that was in consideration for deferral
                guid["defer_targ"] = params.a.λ_targs[1] # Add the next target in the queue to consideration for deferral
                guid["defer_time"] = guid["cur_ddto"].targs[guid["defer_targ"]].t[params.a.τ_targs[1]]
                guid["defer_state"] = guid["cur_ddto"].targs[guid["defer_targ"]].x[:,params.a.τ_targs[1]]
            end

            # Reached minimum target threshold
            if params.a.n_targs < params.n_targs_min
                ~flags["guid_lock_staged"] && @printf("  -> UPDATE [%.2f s]: DDTO recomputation staged [target set count below the minimum threshold]\n", sim_cur_time)
                flags["update_ddto"] = true
                break
            end

            # Reassess criterion
            if criteria == "time"
                criterion = guid["cur_time"] >= guid["defer_time"]
            elseif criteria == "altitude"
                criterion = sim_cur_state[3] <= guid["defer_state"][3]
            end
        end
    end

    return guid, flags
end

function check_cutoff_altitude!(sim_cur_state::Vector{Float64}, sim_cur_time::Float64, cutoff_altitude::Float64, flags::Dict)
    """
    Check if we have reached the cutoff altitude
    """
    cur_altitude = sim_cur_state[3]
    if cur_altitude <= cutoff_altitude
        @printf("  -> UPDATE [%.2f s]: Guidance locked [Cutoff altitude reached!]\n", sim_cur_time)
        flags["guid_lock_staged"] = true
    end

    return flags
end

function activate_guidance_lock!(params, guid::Dict, flags::Dict, sim_cur_time::Float64)
    """
    Lock guidance to best current target if necessary
    """
    if params.a.n_targs == 1 # Wait to lock until we have only one target remaining
        # Determine the current "best" target in terms of radius and obtain the corresponding trajectory
        targ_best_idx = argmax(params.R_targs)
        targ_best = params.a.T_targs[targ_best_idx]
        guid["cur_traj"] = extract_segment(guid["cur_ddto"], targ_best, guid["λ_targs_org"])
        guid["cur_traj_sim"] = extract_segment(guid["cur_ddto_sim"], targ_best, guid["λ_targs_org"])
        
        # Parameter updates
        guid["defer_targ"] = targ_best
        guid["defer_time"] = 1e6
        if ~flags["guid_lock_activated"]
            flags["log_ddto_results"] = true
        end
        flags["guid_lock_activated"] = true
        flags["guid_lock_staged"] = false
        guid["lock_time"] = sim_cur_time
        @printf("  -> UPDATE [%.2f s]: Guidance locked to target %i\n", sim_cur_time, guid["defer_targ"])

        # Remove all targets except for locked target (`guid["defer_targ"]`)
        other_targs = copy(params.a.T_targs)
        deleteat!(other_targs, findfirst(i->i==guid["defer_targ"], other_targs))
        for targ in other_targs
            remove_ddto_target!(params, targ)
        end
    end

    return guid, flags
end