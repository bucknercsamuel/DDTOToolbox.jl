function compute_cum_thrust(data)
    cum_thrust = 0.
    for k = 1:length(data["sim_time"])-1
        delta_time = data["sim_time"][k+1] - data["sim_time"][k]
        cur_thrust = data["sim_control"][:,k]
        cum_thrust += norm(cur_thrust) * delta_time
    end
    return cum_thrust
end
function compute_induced_energy(data; rho=1.225, S_A=0.18*0.11)
    induced_energy = 0.
    for k = 1:length(data["sim_time"])-1
        delta_time = data["sim_time"][k+1] - data["sim_time"][k]
        cur_thrust = data["sim_control"][:,k]
        induced_energy += norm(cur_thrust)^(3/2) * delta_time / sqrt(2*rho*S_A)
    end
    return induced_energy
end
function compute_mechanical_energy(data)
    mechanical_energy = 0.
    for k = 1:length(data["sim_time"])-1
        delta_time = data["sim_time"][k+1] - data["sim_time"][k]
        cur_thrust = data["sim_control"][:,k]
        cur_velocity = data["sim_state"][4:6,k]
        mechanical_energy += dot(cur_thrust, cur_velocity) * delta_time
    end
    return mechanical_energy
end
function compute_ate(data; track_ahead_alt=1., downsample=10)
    # Compute ATE
    ate = 0.
    iter = 0
    for k = 1:downsample:length(data["sim_time"])
        cur_time = data["sim_time"][k]
        cur_guid_idx = findlast(x -> x <= cur_time, data["guid_update_times"])
        guid_defer_idx = data["guid_prefer_vecs"][cur_guid_idx][end]
        cur_guid = data["guid_ddto_trajs_sims"][cur_guid_idx,guid_defer_idx]
        cur_state = data["sim_state"][1:3,k]
        cur_state_shifted = reshape(cur_state + [0,0,-track_ahead_alt],3,1)
        dists = [norm(cur_state_shifted - cur_guid[1:3,kk]) for kk = 1:size(cur_guid,2)]
        track_idx = argmin(dists)
        ate += dists[track_idx]
        iter += 1
    end
    return ate/iter
end
function compute_radius_at_cutoff(data; window_size=10)
    # Find the largest radius at the cutoff time
    # we use window_size because the software may not capture radii timing updates precisely
    cutoff_time = data["guid_update_times"][end]
    cutoff_idx = findlast(x -> x <= cutoff_time, data["sim_time"])
    if cutoff_idx == length(data["sim_time"])
        display(data)
    end
    largest_radius = maximum(data["sim_targs_radii"][:,cutoff_idx-window_size:cutoff_idx])

    return largest_radius
end
function get_altitude_at_cutoff(data)
    # Cutoff time/index, matching compute_radius_at_cutoff
    cutoff_time = data["guid_update_times"][end]
    cutoff_idx = findlast(x -> x <= cutoff_time, data["sim_time"])

    # Check altitude (drone z minus terrain z at the cutoff x,y) against the expected cutoff altitude
    x_cutoff = data["sim_state"][1, cutoff_idx]
    y_cutoff = data["sim_state"][2, cutoff_idx]
    z_cutoff = data["sim_state"][3, cutoff_idx]
    alt_at_cutoff = nothing
    try
        z_terrain = map_data["zlookup"][Int(round(y_cutoff)), Int(round(x_cutoff))] # swapped due to NED to ENU conversion
        alt_at_cutoff = z_cutoff - z_terrain
    catch e
        alt_at_cutoff = -1 # invalid if we can't use map lookup, return negative altitude to indicate invalid
    end
    return alt_at_cutoff, cutoff_idx
end
function compute_safety_of_run(data; R_min = 1., E_max = 200.)
    # Make sure induced energy falls below P_max
    # induced_energy = compute_induced_energy(data)
    induced_energy = compute_cum_thrust(data)
    if induced_energy > E_max
        return false
    end

    # Make sure the largest radius at the cutoff time is greater than R_min
    largest_radius = compute_radius_at_cutoff(data)
    if largest_radius < R_min
        return false
    end

    # Make sure the error code is 1
    if data["error_code"] != 1
        return false
    end

    return true
end
function validate_run(data, map_data; min_allowable_altitude=20.0, max_allowable_altitude=Inf)
    # Check the altitude at the cutoff time
    alt_at_cutoff,idx_at_cutoff = get_altitude_at_cutoff(data)
    if alt_at_cutoff < min_allowable_altitude || alt_at_cutoff > max_allowable_altitude
        @warn "Cutoff altitude check failed" alt_at_cutoff min_allowable_altitude max_allowable_altitude
        return false
    end

    # Ensure that at the cutoff idx onwards, only one radius is non-zero
    cutoff_radii = data["sim_targs_radii"][:,idx_at_cutoff:end]
    test1 = ~all([sum(cutoff_radii[:,k] .>= 0) <= 1 for k = 1:size(cutoff_radii,2)])
    test2 = ~(sum(cutoff_radii[:,1] .>= 0) == 1)
    if test1 || test2
        @warn "Cutoff radius check failed" cutoff_radii
        return false
    end

    # Get the index corresponding to the only non-zero radius at cutoff idx
    cutoff_idx = findfirst(x -> x >= 0, cutoff_radii[:,1])

    # End state of the cutoff trajectory (last guidance update, preferred deferred branch)
    guid_defer_idx = data["guid_prefer_vecs"][end][end]
    cutoff_traj = data["guid_ddto_trajs_sims"][end, cutoff_idx]
    cutoff_radius = data["sim_targs_radii"][cutoff_idx, idx_at_cutoff]
    cutoff_end_pos = cutoff_traj[1:3, end]

    # # Check the final simulated position against the cutoff trajectory's end position
    # final_pos = data["sim_state"][1:3, end]
    # pos_err = norm(final_pos - cutoff_end_pos)
    # if pos_err > cutoff_radius
    #     @warn "Final target position check failed" final_pos cutoff_end_pos pos_err cutoff_radius
    #     return false
    # end

    return true
end