# Scoring split:
#   1. usable_log        — pickle/IC integrity (algorithm-agnostic)
#   2. landed            — error_code == 1
#   3. operational_success — landed + committed + site radius + energy budget
#
# Cutoff/commitment is an *outcome*, not a reason to drop the trial.

# ---------------------------------------------------------------------------
# Integrals (use gap-filled state/control when available)
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Terrain / allocation helpers
# ---------------------------------------------------------------------------
function _map_zlookup(map_data)
    map_data === nothing && return nothing
    return get(map_data, "zlookup", nothing)
end

function n_allocated_targets(radii_col)
    count(x -> isfinite(x) && x >= 0, radii_col)
end

function altitude_agl(data, idx, map_data)
    zlookup = _map_zlookup(map_data)
    zlookup === nothing && return nothing
    x = data["sim_state"][1, idx]
    y = data["sim_state"][2, idx]
    z = data["sim_state"][3, idx]
    try
        z_terrain = zlookup[Int(round(y)), Int(round(x))]
        return z - z_terrain
    catch
        return nothing
    end
end

"""
Algorithm-agnostic log integrity. Does *not* inspect cutoff shape or lock.
Start pose is a shared IC, so a failed start lookup is a trial-level defect.
"""
function usable_log(data, map_data=nothing)
    for key in ("sim_time", "sim_state", "sim_control", "sim_targs_radii", "error_code")
        haskey(data, key) || return false
    end
    t = data["sim_time"]
    n = length(t)
    n < 2 && return false
    size(data["sim_state"], 2) == n || return false
    size(data["sim_control"], 2) == n || return false
    size(data["sim_targs_radii"], 2) == n || return false
    any(isfinite, t) || return false
    agl0 = altitude_agl(data, 1, map_data)
    agl0 === nothing && return false
    return true
end

function landed(data)
    try
        return Int(data["error_code"]) == 1
    catch
        return false
    end
end

"""
Commit index:
- Graph-DDTO: first sample of the terminal suffix where exactly one target
  stays allocated through the end (the N→1 lock).
- Greedy (alloc==1 for the whole flight): first sample with AGL ≤ h_cut,
  i.e. the actual cutoff-altitude policy rather than t=0.
Returns `nothing` if the vehicle never locked / never reached h_cut.
"""
function find_commit_index(data, map_data=nothing; h_cut=50.0)
    radii = data["sim_targs_radii"]
    n = size(radii, 2)
    n == 0 && return nothing
    alloc = [n_allocated_targets(view(radii, :, k)) for k in 1:n]
    alloc[end] != 1 && return nothing
    k = n
    while k > 1 && alloc[k - 1] == 1
        k -= 1
    end
    if k > 1
        return k
    end
    for i in 1:n
        agl = altitude_agl(data, i, map_data)
        if agl !== nothing && agl <= h_cut
            return i
        end
    end
    return nothing
end

function radius_at_index(data, idx; window_size=10)
    radii = data["sim_targs_radii"]
    n = size(radii, 2)
    i0 = max(1, idx - window_size)
    i1 = min(n, idx)
    return maximum(radii[:, i0:i1])
end

function compute_radius_at_cutoff(data, map_data=nothing; window_size=10, h_cut=50.0)
    idx = find_commit_index(data, map_data; h_cut=h_cut)
    idx === nothing && return NaN
    return radius_at_index(data, idx; window_size=window_size)
end

function get_altitude_at_cutoff(data, map_data=nothing; h_cut=50.0)
    idx = find_commit_index(data, map_data; h_cut=h_cut)
    idx === nothing && return (NaN, 0)
    agl = altitude_agl(data, idx, map_data)
    agl === nothing && return (NaN, idx)
    return (agl, idx)
end

"""
Operational success. Commitment/altitude/radius are outcomes of *this*
algorithm, never a reason to delete the paired trial.

Defaults:
- R_min = 1 m matches HALO `R_targs_min`
- E_max = 2500 J induced energy (see analysis notes; analog of the old 200 N·s impulse cap)
- min_commit_agl = 20 m: must lock with enough altitude left to track the branch
"""
function compute_safety_of_run(data, map_data=nothing;
        R_min=1.0, E_max=2500.0, min_commit_agl=20.0, h_cut=50.0)
    if haskey(data, "landed") && haskey(data, "committed") && haskey(data, "radius_at_cutoff")
        return operational_success_from_fields(data; R_min=R_min, E_max=E_max, min_commit_agl=min_commit_agl)
    end
    landed(data) || return false
    idx = find_commit_index(data, map_data; h_cut=h_cut)
    idx === nothing && return false
    radius_at_index(data, idx) < R_min && return false
    agl = altitude_agl(data, idx, map_data)
    if agl !== nothing && agl < min_commit_agl
        return false
    end
    energy = haskey(data, "induced_energy") ? data["induced_energy"] : compute_induced_energy(data)
    energy > E_max && return false
    return true
end

function operational_success_from_fields(data; R_min=1.0, E_max=2500.0, min_commit_agl=20.0)
    data["landed"] || return false
    data["committed"] || return false
    !(isfinite(data["radius_at_cutoff"]) && data["radius_at_cutoff"] >= R_min) && return false
    agl = data["altitude_at_cutoff"]
    isfinite(agl) && agl < min_commit_agl && return false
    energy = get(data, "induced_energy", Inf)
    energy > E_max && return false
    return true
end

# Deprecated name: previously mixed cutoff-shape into data quality.
# Now only log integrity; keep the name so older scripts keep compiling.
function validate_run(data, map_data=nothing; kwargs...)
    return usable_log(data, map_data)
end

"""
Score a run *before* gap-filling (events) and optionally refresh integrals after.
Does not mutate `error_code`.
"""
function score_run!(data, map_data=nothing; h_cut=50.0, R_min=1.0, E_max=2500.0, min_commit_agl=20.0)
    data["usable"] = usable_log(data, map_data)
    data["landed"] = landed(data)
    idx = find_commit_index(data, map_data; h_cut=h_cut)
    data["committed"] = idx !== nothing
    data["commit_idx"] = idx === nothing ? 0 : idx
    data["radius_at_cutoff"] = idx === nothing ? NaN : radius_at_index(data, idx)
    agl, _ = get_altitude_at_cutoff(data, map_data; h_cut=h_cut)
    data["altitude_at_cutoff"] = agl
    guid_t = haskey(data, "guid_update_times_orig") ? data["guid_update_times_orig"] :
             get(data, "guid_update_times", Float64[])
    data["num_recomputations"] = max(length(guid_t) - 1, 0)
    data["cum_thrust"] = compute_cum_thrust(data)
    data["induced_energy"] = compute_induced_energy(data)
    data["mechanical_energy"] = compute_mechanical_energy(data)
    data["ATE"] = 0.0
    data["safe_run"] = compute_safety_of_run(data, map_data;
        R_min=R_min, E_max=E_max, min_commit_agl=min_commit_agl, h_cut=h_cut)
    return data
end
