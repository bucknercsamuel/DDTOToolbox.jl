using CairoMakie
using Colors
using LinearAlgebra
include("../../utils/plot_utils.jl")
include("plot_defaults.jl")
include("plot_drone.jl")

# ----------------------------------------------------------------------------
# Shared scene setup (terrain + axis styling)
# ----------------------------------------------------------------------------

function _paper_demo_axis3(f, ax_idx; azel=(pi/4, pi/4))
    ax = Axis3(
        f[ax_idx...],
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
    return ax
end

"""Plot terrain surface; `xLims`/`yLims` are ENU limits (swapped to NED for lookup)."""
function _paper_demo_terrain!(ax, map_data, xLims, yLims; terrain_alpha=0.8)
    xlims_terrain = yLims
    ylims_terrain = xLims
    xs = LinRange(xlims_terrain[1], xlims_terrain[2], Int(round(xlims_terrain[2] - xlims_terrain[1])))
    ys = LinRange(ylims_terrain[1], ylims_terrain[2], Int(round(ylims_terrain[2] - ylims_terrain[1])))
    zs = [map_data["zlookup"][Int(round(x)), Int(round(y))] for x in xs, y in ys]
    xs, ys, zs = ys, xs, transpose(zs)
    surface!(ax, xs, ys, zs;
        colormap = range(parse(Colorant, "orangered4"), stop=parse(Colorant, "darkorange2"), length = 100),
        alpha = terrain_alpha)
end

function _safe_thrust_direction(thrust)
    thrust = vec(collect(Float64, thrust))
    n = norm(thrust)
    return n < 1e-10 ? [1.0, 0.0, 0.0] : thrust / n
end

"""
Pose at the start of a DDTO bundle: position from the first point of the plan
(all branches share this trunk origin); thrust from logged sim control at the
corresponding guidance-update time.
"""
function _bundle_start_pose(run_data, bundle_idx, trajs)
    pos = nothing
    for traj in trajs
        m = collect(traj)
        if size(m, 1) >= 3 && size(m, 2) >= 1
            pos = vec(Float64.(m[1:3, 1]))
            break
        end
    end

    sim_time = vec(collect(Float64, run_data["sim_time"]))
    sim_state = collect(run_data["sim_state"])
    sim_control = collect(run_data["sim_control"])
    guid_t = Float64(run_data["guid_update_times"][bundle_idx])
    i = findlast(j -> sim_time[j] <= guid_t, eachindex(sim_time))
    i = something(i, 1)

    if pos === nothing
        pos = vec(Float64.(sim_state[1:3, i]))
    end
    thrust = vec(Float64.(sim_control[1:3, i]))
    return pos, thrust
end

"""Drone `scale` matched to scene size (see `plot_3d_trajs`, which uses `scale=15`)."""
function _auto_drone_scale(xLims, yLims, zLims; fraction = 0.04)
    span = maximum([xLims[2] - xLims[1], yLims[2] - yLims[1], zLims[2] - zLims[1]])
    return max(span * fraction, 1.0)
end

"""Collect target trajectories for one guidance bundle (matrix or nested layout)."""
function _guid_bundle_trajs(guid_ddto_trajs, bundle_idx)
    if ndims(guid_ddto_trajs) == 2
        n = size(guid_ddto_trajs, 2)
        return [collect(guid_ddto_trajs[bundle_idx, j]) for j in 1:n]
    else
        bundle = guid_ddto_trajs[bundle_idx]
        return [collect(bundle[j]) for j in 1:length(bundle)]
    end
end

function _lims_from_positions(positions::AbstractVector{<:AbstractVector})
    rmin = [minimum(p[k] for p in positions) for k in 1:3]
    rmax = [maximum(p[k] for p in positions) for k in 1:3]
    return get_equal_3d_lims(rmin, rmax)
end

function _finalize_paper_demo_plot(f, interactive, save_path, default_name)
    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    else
        out = save_path === nothing ? joinpath(fig_path, default_name * fig_ext) : save_path
        CairoMakie.save(out, f, px_per_unit = 4)
        return f
    end
end

# ----------------------------------------------------------------------------
# Full history figure (all bundles + sim trace)
# ----------------------------------------------------------------------------

function plot_paper_demo_traj_history(
        run_data,
        map_data;
        interactive=true,
        azel=(pi/4,pi/4),
        save_path=nothing,
        show_drone_timelapse=true
    )
    f = Figure(size = (800, 800))
    ax = _paper_demo_axis3(f, (1, 1); azel = azel)

    # Results parsing
    sim_time = run_data["sim_time"]
    sim_state = run_data["sim_state"]
    # sim_control = run_data["sim_control"]
    guid_ddto_trajs = run_data["guid_ddto_trajs_sims"]
    guid_defer_vecs = run_data["guid_defer_vecs"]
    guid_update_times = run_data["guid_update_times"]
    n_ddto_sols = length(guid_update_times)

    # # Clean part of sim state where the vehicle resets to above its position (should not be plotted)
    # idx_reset = findfirst(dh->dh>0, diff(sim_state[3,:]))
    # sim_state = sim_state[:,1:idx_reset]

    # Results parsing
    rmin = [min(sim_state[k,:]...) for k∈1:3]
    rmax = [max(sim_state[k,:]...) for k∈1:3]
    xLims,yLims,zLims = get_equal_3d_lims(rmin, rmax)
    # xlims!(ax, xLims...)
    # ylims!(ax, yLims...)
    # zlims!(ax, zLims...)

    _paper_demo_terrain!(ax, map_data, xLims, yLims)

    # Plot DDTO trajectory bundles, color-coded by temporal progress along the
    # trajectory (guid_update_times[k] mapped onto [sim_time[1], sim_time[end]]
    # via the viridis colormap).
    sim_t0, sim_tf = sim_time[1], sim_time[end]
    Δsim = sim_tf - sim_t0
    bundle_cmap = cgrad(:viridis)
    is_matrix_layout = ndims(guid_ddto_trajs) == 2
    for k = 1:n_ddto_sols
        frac = Δsim > 0 ? clamp((guid_update_times[k] - sim_t0) / Δsim, 0.0, 1.0) : 0.0
        color_k = bundle_cmap[frac]

        # Number of targets per bundle is data-driven (no reliance on params)
        n_targs = is_matrix_layout ? size(guid_ddto_trajs, 2) : length(guid_ddto_trajs[k])

        for j = 1:n_targs
            traj = collect(is_matrix_layout ? guid_ddto_trajs[k, j] : guid_ddto_trajs[k][j])
            (size(traj, 1) >= 3 && size(traj, 2) >= 1) || continue
            lines!(ax,
                traj[1,:], traj[2,:], traj[3,:];
                style3D_ct_ddto..., :alpha=>0.9, :color=>color_k)
        end
    end

    # Plot the recorded sim state
    positions = sim_state[1:3,:]
    lines!(ax,
        positions[1,:], positions[2,:], positions[3,:];
        style3D_ct..., :alpha=>.5, :color=>:black)

    # Plot the vehicle body timelapse
    if show_drone_timelapse
        num_frames = 10
        idx_frames = round.(Int, LinRange(1, size(sim_state,2), num_frames))
        for idx in idx_frames
            position = sim_state[1:3,idx]
            thrust = run_data["sim_control"][1:3,idx]
            thrust_direction = thrust / norm(thrust)
            plot_drone(ax, position, thrust_direction; scale=.2)
        end
    end

    return _finalize_paper_demo_plot(f, interactive, save_path, "paper_demo_traj")
end

# ----------------------------------------------------------------------------
# Single-bundle figure (initial drone + one DDTO bundle, custom colors)
# Requires `using DDTOToolbox` in the calling script (for Quad3DoFHaloParams).
# ----------------------------------------------------------------------------

"""
    plot_paper_demo_traj_bundle(run_data, map_data; bundle_idx=1, kwargs...)

Like `plot_paper_demo_traj_history`, but plots only one DDTO trajectory bundle
(`bundle_idx` = 1 is the first guidance update) and omits the simulated flight
history. The vehicle is drawn once at the start of the selected bundle.

Deferrable (trunk) and target (branch) segments are colored separately when
`guid_defer_vecs` / `guid_prefer_vecs` are present in `run_data` (via `plot_bundle`).
Otherwise each target trajectory is drawn in full using `branch_colors`.

# Color kwargs
- `defer_color`: trunk / deferrable segment (default `:gray`).
- `branch_colors`: per-target branch colors (default `[:red, :blue, :green, :orange]`).
- `show_ddto_split`: split trunk vs branches; default `true` when defer/prefer vecs exist.
- `drone_scale`: Makie scale for `plot_drone`; default `nothing` auto-sizes from scene extent (~4% of span).
"""
function plot_paper_demo_traj_bundle(
        run_data,
        map_data;
        bundle_idx = 1,
        interactive = true,
        azel = (pi/4, pi/4),
        save_path = nothing,
        defer_color = :gray,
        branch_colors = [:red, :blue, :green, :orange],
        show_ddto_split = nothing,
        show_defer_nodes = false,
        show_sol_nodes = false,
        drone_scale = nothing,
        terrain_alpha = 0.8,
        traj_alpha = 0.9,
    )
    guid_ddto_sims = run_data["guid_ddto_trajs_sims"]
    guid_ddto_sols = get(run_data, "guid_ddto_trajs_sols", nothing)
    guid_defer_vecs = get(run_data, "guid_defer_vecs", nothing)
    guid_prefer_vecs = get(run_data, "guid_prefer_vecs", nothing)
    n_bundles = length(run_data["guid_update_times"])

    bundle_idx in 1:n_bundles ||
        error("bundle_idx=$bundle_idx out of range 1:$n_bundles")

    trajs = _guid_bundle_trajs(guid_ddto_sims, bundle_idx)
    n_targs = length(trajs)
    n_targs >= 1 || error("bundle $bundle_idx has no target trajectories")

    # Axis limits from bundle origin + selected bundle trajectories
    bundle_start, thrust_start = _bundle_start_pose(run_data, bundle_idx, trajs)
    traj_pts = Vector{Vector{Float64}}()
    push!(traj_pts, bundle_start)
    for traj in trajs
        m = collect(traj)
        (size(m, 1) >= 3 && size(m, 2) >= 1) || continue
        for col in 1:size(m, 2)
            push!(traj_pts, m[1:3, col])
        end
    end
    xLims, yLims, zLims = _lims_from_positions(traj_pts)

    f = Figure(size = (800, 800))
    ax = _paper_demo_axis3(f, (1, 1); azel = azel)
    xlims!(ax, xLims...)
    ylims!(ax, yLims...)
    zlims!(ax, zLims...)
    _paper_demo_terrain!(ax, map_data, xLims, yLims; terrain_alpha = terrain_alpha)

    color_branch = j -> branch_colors[mod1(j, length(branch_colors))]
    proj_idxs = [1, 2, 3]
    data_sims = [[trajs[j][c, :] for j in 1:n_targs] for c in proj_idxs]
    data_sols = if guid_ddto_sols !== nothing
        sol_trajs = _guid_bundle_trajs(guid_ddto_sols, bundle_idx)
        [[sol_trajs[j][c, :] for j in 1:n_targs] for c in proj_idxs]
    else
        data_sims
    end

    can_split = guid_defer_vecs !== nothing && guid_prefer_vecs !== nothing &&
        bundle_idx <= length(guid_defer_vecs) && bundle_idx <= length(guid_prefer_vecs)
    do_split = show_ddto_split === nothing ? can_split : show_ddto_split

    if do_split && can_split
        virt_params = Quad3DoFHaloParams()
        virt_params.a.τ_targs = guid_defer_vecs[bundle_idx]
        virt_params.a.λ_targs = guid_prefer_vecs[bundle_idx]
        virt_params.a.n_targs = n_targs
        plot_bundle(ax,
            data_sols,
            data_sims,
            virt_params,
            style3D_ct_ddto,
            style3D_dt;
            color_trunk = defer_color,
            color_branch = color_branch,
            show_sol_nodes = show_sol_nodes,
            show_defer_nodes = show_defer_nodes,
            show_ddto_split = true,
            alpha = traj_alpha,
        )
    else
        for j in 1:n_targs
            m = collect(trajs[j])
            (size(m, 1) >= 3 && size(m, 2) >= 1) || continue
            lines!(ax, m[1, :], m[2, :], m[3, :];
                style3D_ct_ddto..., :alpha => traj_alpha, :color => color_branch(j))
        end
    end

    # scale_drone = drone_scale === nothing ? _auto_drone_scale(xLims, yLims, zLims) : drone_scale
    plot_drone(ax, bundle_start, _safe_thrust_direction(thrust_start);
        scale = 15, overdraw = true)

    default_name = "paper_demo_bundle$(bundle_idx)"
    return _finalize_paper_demo_plot(f, interactive, save_path, default_name)
end
