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
    params.a.z0[7] = 0.0 # Set initial cum thrust to zero

    # Saturate velocities to satisfy constraints
    vel_lat_idx = 4:5
    vel_vert_idx = 6
    eps = 1e-3
    if norm(params.a.z0[vel_lat_idx]) > params.v_max_L - eps
        params.a.z0[vel_lat_idx] = params.a.z0[vel_lat_idx] / norm(params.a.z0[vel_lat_idx]) * (params.v_max_L - eps)
        display("Warning: velocity saturated!")
    end
    params.a.z0[vel_vert_idx] = max(min(params.a.z0[vel_vert_idx], params.v_max_V-eps), params.v_min_V+eps)

    # Set guidance initial control to current sim control
    for k = 1:3
        params.a.u0[k] = sim_cur_control[k]
    end

    # Saturate thrusts to satisfy constraints
    thrust_idx = 1:3
    eps = 0.1 * (params.ρ_max - params.ρ_min) # buffer to both satisfy constraint and give some room for guidance
    if norm(params.a.u0[thrust_idx]) > params.ρ_max - eps && ~isinf(norm(params.a.u0[thrust_idx]))
        params.a.u0[thrust_idx] = params.a.u0[thrust_idx] / norm(params.a.u0[thrust_idx]) * (params.ρ_max - eps)
    elseif norm(params.a.u0[thrust_idx]) < params.ρ_min + eps
        params.a.u0[thrust_idx] = params.a.u0[thrust_idx] / norm(params.a.u0[thrust_idx]) * (params.ρ_min + eps)
    end

    # Guidance solving
    flags["ddto_converged"] = false
    # try
        _,_,guid["cur_ddto"],guid["cur_ddto_sim"],flags["ddto_converged"] = solve(params) # Compute DDTO solution
        # guid["cur_ddto"],guid["cur_ddto_sim"],_,_,flags["ddto_converged"] = solve(params) # Compute DDTO solution
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
        guid["cur_traj"] = guid["cur_ddto"].targs[params.a.λ_targs[end]]
        guid["cur_traj_sim"] = guid["cur_ddto_sim"].targs[params.a.λ_targs[end]]
        @printf("  -> UPDATE [%.2f s]: DDTO solution successfully recomputed [tracking trunk segment]\n", sim_cur_time)

        # Parameter updates
        guid["cur_time"] = 0.0 # Reset guidance time to zero
        guid["defer_targ"] = params.a.λ_targs[1]
        guid["defer_time"] = guid["cur_ddto"].targs[guid["defer_targ"]].t[params.a.τ_targs[1]]
        guid["defer_state"] = guid["cur_ddto"].targs[guid["defer_targ"]].x[:,params.a.τ_targs[1]]
        guid["λ_targs_org"] = params.a.λ_targs
    end

    # Flag updates
    flags["update_ddto"] = false
    flags["log_ddto_results"] = true

    return guid, flags
end

function check_unsafe_targets!(params, guid::Dict, flags::Dict, sim_cur_time::Float64)
    """
    Check for unsafe targets (radii check) -- if too many unsafe targets, remove and schedule update
    """
    # Find all targets that are unsafe
    cur_targs = copy(params.a.J_targs)
    unsafe_targs = []
    num_targs_unsafe = 0
    for targ in cur_targs
        targ_idx = findfirst(i->i==targ, params.a.J_targs)
        if params.R_targs[targ_idx] <= params.R_targs_min
            append!(unsafe_targs, targ)
            num_targs_unsafe += 1
        end
    end

    # Check if we have too many unsafe targets; if so, remove them and schedule DDTO update
    if (params.a.n_targs - num_targs_unsafe) < params.n_targs_min
        for targ in unsafe_targs
            remove_ddto_target!(params, targ)
        end
        ~flags["guid_lock_staged"] && @printf("  -> UPDATE [%.2f s]: DDTO recomputation staged [target set count below the minimum threshold]\n", sim_cur_time)
        flags["update_ddto"] = true
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
                other_targs = copy(params.a.J_targs)
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

function check_cutoff_altitude!(sim_cur_time::Float64, altitude::Float64, cutoff_altitude::Float64, flags::Dict)
    """
    Check if we have reached the cutoff altitude
    """
    if altitude <= cutoff_altitude
        @printf("  -> UPDATE [%.2f s]: Guidance locked [Cutoff altitude reached!]\n", sim_cur_time)
        flags["guid_lock_staged"] = true
    end

    return flags
end

function activate_guidance_lock!(params, guid::Dict, flags::Dict, sim_cur_time::Float64)
    """
    Lock guidance to best current target if necessary
    """
    # Wait to lock until we have only one target remaining to fully lock the guidance
    if params.a.n_targs == 1
        # Determine the current "best" target in terms of radius and obtain the corresponding trajectory
        # (unnecessary if we use the n_targs==1 condition above, but kept for consistency)
        targ_best_idx = argmax(params.R_targs)
        targ_best = params.a.J_targs[targ_best_idx]
        guid["cur_traj"] = guid["cur_ddto"].targs[targ_best]
        guid["cur_traj_sim"] = guid["cur_ddto_sim"].targs[targ_best]
        
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
    end

    return guid, flags
end

function step_halo_sim(params, sim_cur_time::CReal, sim_cur_state::CVector, guid::Dict; Δt_sim::Float64=0.01)
    feedforward_controller = (t) -> optimal_controller(t, guid["cur_traj_sim"].t, guid["cur_traj_sim"].u, params.a.disc)
    dynamics = (t,x) -> dynamics_nonlinear_nondilated(t, x, feedforward_controller(t), params)
    sim_cur_state = rk4_step(sim_cur_state, dynamics, guid["cur_time"], Δt_sim)
    guid["cur_time"] += Δt_sim
    sim_cur_time     += Δt_sim
    return sim_cur_time, sim_cur_state
end