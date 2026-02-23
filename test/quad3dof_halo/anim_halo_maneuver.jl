"""
Animated version of plot_paper_demo_traj_history: the drone body moves in real time
along the trajectory instead of showing a timelapse of multiple drone poses.
The most recent DDTO trajectory bundle (from guid_ddto_trajs_sims) is shown and
updated whenever a new solution is computed (according to guid_update_times).
"""
# using GLMakie
using LinearAlgebra
include("plots.jl")

"""Index into guid_ddto_trajs_sims for the bundle active at time t."""
function current_guid_bundle_index(guid_update_times, t)
    idx = findlast(x -> x <= t, guid_update_times)
    return idx === nothing ? 1 : idx
end

"""Extract (x, y, z) vectors for lines! from a single trajectory (3×n or similar)."""
function traj_to_xyz(traj)
    m = collect(traj)
    if ndims(m) == 1 && length(m) == 1
        m = collect(m[1])
    end
    if ndims(m) == 2 && size(m, 1) >= 3
        return (m[1, :], m[2, :], m[3, :])
    end
    return (Float64[], Float64[], Float64[])
end

function animate_paper_demo_traj_history(
        run_data,
        map_data;
        azel = (pi/4, pi/4),
        fps = 30,
        playback_speed = 1.0,
        loop = true,
        scale = 15,
        show_time_label = true
    )
    sim_state = run_data["sim_state"]
    sim_control = run_data["sim_control"]
    sim_time = run_data["sim_time"]
    n_steps = size(sim_state, 2)

    guid_ddto_trajs_sims = run_data["guid_ddto_trajs_sims"]
    guid_update_times = run_data["guid_update_times"]
    n_bundle_trajs = ndims(guid_ddto_trajs_sims) == 2 ? size(guid_ddto_trajs_sims, 2) : length(guid_ddto_trajs_sims[1])
    sim_duration = sim_time[end] - sim_time[1]

    # Build the same static scene as plot_paper_demo_traj_history (no timelapse drones)
    f = Figure(size = (800, 800))
    ax_cell = show_time_label ? (f[2, 1]) : (f[1, 1])
    ax = Axis3(
        ax_cell,
        xlabel = "",
        ylabel = "",
        zlabel = "",
        aspect = :equal,
        azimuth = azel[1],
        elevation = azel[2],
        xgridvisible = false,
        ygridvisible = false,
        zgridvisible = false,
        xticklabelsvisible = false,
        yticklabelsvisible = false,
        zticklabelsvisible = false,
        xgridcolor = :transparent,
        ygridcolor = :transparent,
        zgridcolor = :transparent,
        xticksvisible = false,
        yticksvisible = false,
        zticksvisible = false,
    )
    hidespines!(ax)

    rmin = [min(sim_state[k, :]...) for k ∈ 1:3]
    rmax = [max(sim_state[k, :]...) for k ∈ 1:3]
    xLims, yLims, zLims = get_equal_3d_lims(rmin, rmax)
    # Lock axis limits so the view does not jump when the DDTO bundle updates (avoids "teleporting" drone)
    xlims!(ax, xLims...)
    ylims!(ax, yLims...)
    zlims!(ax, zLims...)

    # Terrain (same as plot_paper_demo_traj_history)
    xlims_terrain = yLims
    ylims_terrain = xLims
    xs = LinRange(xlims_terrain[1], xlims_terrain[2], Int(round(xlims_terrain[2] - xlims_terrain[1])))
    ys = LinRange(ylims_terrain[1], ylims_terrain[2], Int(round(ylims_terrain[2] - ylims_terrain[1])))
    zs = [map_data["zlookup"][Int(round(x)), Int(round(y))] for x in xs, y in ys]
    xs, ys, zs = ys, xs, transpose(zs)

    surface!(ax,
        xs, ys, zs;
        colormap = range(parse(Colorant, "orangered4"), stop = parse(Colorant, "darkorange2"), length = 100),
        alpha = 0.8)

    # Drone trajectory trace: only where the drone has been (updated each frame)
    positions = sim_state[1:3, :]
    traj_x_obs = Observable(positions[1, 1:1])
    traj_y_obs = Observable(positions[2, 1:1])
    traj_z_obs = Observable(positions[3, 1:1])
    lines!(ax, traj_x_obs, traj_y_obs, traj_z_obs;
        style3D_ct..., :alpha => 1.0, :color => :black)

    # DDTO bundle: most recent bundle, updated when guid_update_times is crossed
    # Colors: gist_rainbow (matplotlib-like) linearly interpolated 0→1 over n trajectories
    bundle_colors = try
        cgrad(:gist_rainbow, n_bundle_trajs, categorical = true)
    catch
        cgrad(:rainbow, n_bundle_trajs, categorical = true)  # fallback if :gist_rainbow unavailable
    end
    style_ddto = merge(copy(style3D_ct), Dict(:linewidth => 3, :alpha => 1.0))
    x_obs_list = [Observable(Float64[]) for _ in 1:n_bundle_trajs]
    y_obs_list = [Observable(Float64[]) for _ in 1:n_bundle_trajs]
    z_obs_list = [Observable(Float64[]) for _ in 1:n_bundle_trajs]
    for j in 1:n_bundle_trajs
        lines!(ax, x_obs_list[j], y_obs_list[j], z_obs_list[j];
            style_ddto..., color = bundle_colors[j])
    end

    function update_bundle_lines!(guid_idx)
        if ndims(guid_ddto_trajs_sims) == 2
            for j in 1:n_bundle_trajs
                traj = guid_ddto_trajs_sims[guid_idx, j]
                x, y, z = traj_to_xyz(traj)
                x_obs_list[j][] = x
                y_obs_list[j][] = y
                z_obs_list[j][] = z
            end
        else
            bundle = guid_ddto_trajs_sims[guid_idx]
            n = length(bundle)
            for j in 1:n_bundle_trajs
                if j <= n
                    x, y, z = traj_to_xyz(bundle[j])
                    x_obs_list[j][] = x
                    y_obs_list[j][] = y
                    z_obs_list[j][] = z
                else
                    x_obs_list[j][] = Float64[]
                    y_obs_list[j][] = Float64[]
                    z_obs_list[j][] = Float64[]
                end
            end
        end
    end

    current_guid_idx = Ref(current_guid_bundle_index(guid_update_times, sim_time[1]))
    update_bundle_lines!(current_guid_idx[])

    # Thrust direction helper: avoid zero norm
    function safe_thrust_direction(control_col)
        t = control_col
        n = norm(t)
        if n < 1e-10
            return [1.0, 0.0, 0.0]  # fallback
        end
        return t / n
    end

    # Observables for the single moving drone
    position_obs = Observable(sim_state[1:3, 1])
    thrust_direction_obs = Observable(safe_thrust_direction(sim_control[1:3, 1]))

    plot_drone_observable(ax, position_obs, thrust_direction_obs; scale = scale)

    # Optional time label at top center of figure (e.g. "Sim Time: 5.04s (5x)")
    time_label_obs = Observable("Sim Time: 0.00s (1x)")
    if show_time_label
        f[1, 1] = Label(f, time_label_obs; fontsize = 18, halign = :center)
        rowsize!(f.layout, 1, Auto(-1))  # minimal height for label row
    end

    # Update all observables for a given frame index i (used by both timer and record).
    function update_frame(i)
        i = min(max(1, i), n_steps)
        position_obs[] = sim_state[1:3, i]
        thrust_direction_obs[] = safe_thrust_direction(sim_control[1:3, i])
        # traj_x_obs[] = positions[1, 1:i]
        # traj_y_obs[] = positions[2, 1:i]
        # traj_z_obs[] = positions[3, 1:i]
        # if show_time_label
        #     speed_str = playback_speed == 1 ? "1x" : "$(playback_speed)x"
        #     time_label_obs[] = "Sim Time: $(round(sim_time[i], digits=2))s ($speed_str)"
        # end
        # new_guid_idx = current_guid_bundle_index(guid_update_times, sim_time[i])
        # if new_guid_idx != current_guid_idx[]
        #     current_guid_idx[] = new_guid_idx
        #     update_bundle_lines!(new_guid_idx)
        # end
    end

    # Animation driven by wall clock: playback_speed = 5 means 5 sim seconds per 1 real second
    real_start = Ref{Float64}(NaN)
    timer_interval = 1.0 / fps

    function update_drone(tim)
        if isnan(real_start[])
            real_start[] = time()
        end
        elapsed_real = time() - real_start[]
        target_sim = sim_time[1] + (loop ? mod(playback_speed * elapsed_real, sim_duration) : playback_speed * elapsed_real)
        if !loop && target_sim >= sim_time[end]
            update_frame(n_steps)
            close(tim)
            return
        end
        # Use last index with time <= target_sim, then first index with that same time.
        # This avoids "teleporting" when sim_time has duplicate entries at guidance updates
        # (e.g. pre- and post-update state logged at the same t); we show the earlier state.
        i = findlast(j -> sim_time[j] <= target_sim, 1:n_steps)
        i = i === nothing ? 1 : findfirst(j -> sim_time[j] == sim_time[i], 1:n_steps)
        update_frame(i)
    end

    tim = Timer(tim -> update_drone(tim), 0.0, interval = timer_interval)

    screen = GLMakie.Screen()
    display(screen, f)
    return (; screen, timer = tim, position_obs, thrust_direction_obs)
end

# Example usage (uncomment and run from repo root or test/quad3dof_halo):
# using DDTOToolbox
# using Pandas
# include("anim_halo_maneuver.jl")
# local_path = abspath(@__DIR__)
# run_data = read_pickle(joinpath(local_path, "data\\paper_demo_data.pkl"))
# map_data = Dict("zlookup" => read_pickle(joinpath(local_path, "map_lookups\\maps\\dunes_test_hard\\lookup_table.pkl")))
# with_theme(theme3d) do
#     # Live animation only:
#     result = animate_paper_demo_traj_history(run_data, map_data; fps=30, playback_speed=1.0, loop=true)
#     # To stop: close(result.timer); close(result.screen)
# end
