using CairoMakie
using Colors
using LinearAlgebra
using Statistics
include("../../utils/plot_utils.jl")
include("plot_defaults.jl")

"""Index into guid_ddto_trajs_sims for the bundle active at time t."""
function _current_guid_bundle_index(guid_update_times, t)
    idx = findlast(x -> x <= t, guid_update_times)
    return idx === nothing ? 1 : idx
end

"""Closest-point distance from state position (3-vector) to trajectory (3×n matrix of points)."""
function _guidance_error(state_pos, traj)
    m = collect(traj)
    size(m, 1) >= 3 || return 0.0
    n_pts = size(m, 2)
    n_pts == 0 && return 0.0
    p = state_pos isa AbstractVector ? state_pos : vec(state_pos)
    dists = [norm(p - m[1:3, k]) for k in 1:n_pts]
    return minimum(dists)
end

"""
Compute the guidance tracking error history for a run, using the same logic as
`animate_paper_demo_traj_history`: for each sim step, find the active bundle,
the tracked trajectory within that bundle (via `guid_prefer_vecs`, with a
final-update lock-on to the closest-endpoint trajectory), and return the
minimum distance from the drone position to any point on that trajectory.
"""
function compute_guidance_error_history(run_data)
    sim_state = run_data["sim_state"]
    sim_time = run_data["sim_time"]
    guid_ddto_trajs_sims = run_data["guid_ddto_trajs_sims"]
    guid_update_times = run_data["guid_update_times"]
    guid_prefer_vecs = get(run_data, "guid_prefer_vecs", nothing)
    n_steps = size(sim_state, 2)
    n_bundle_trajs = ndims(guid_ddto_trajs_sims) == 2 ?
        size(guid_ddto_trajs_sims, 2) : length(guid_ddto_trajs_sims[1])
    n_guid = length(guid_update_times)

    # Final guidance update locks onto one trajectory: pick the one whose end-state
    # is closest to the final drone position.
    final_lockon_traj_idx = 1
    if n_guid >= 1
        final_state = sim_state[1:3, end]
        dists = Float64[]
        for k in 1:n_bundle_trajs
            traj_k = ndims(guid_ddto_trajs_sims) == 2 ?
                collect(guid_ddto_trajs_sims[n_guid, k]) :
                collect(guid_ddto_trajs_sims[n_guid][k])
            if size(traj_k, 1) < 3 || size(traj_k, 2) < 1
                push!(dists, Inf)
                continue
            end
            push!(dists, norm(final_state .- traj_k[1:3, end]))
        end
        final_lockon_traj_idx = isempty(dists) ? 1 : argmin(dists)
    end

    function get_tracked_traj(guid_idx)
        traj_idx = if guid_idx == n_guid
            final_lockon_traj_idx
        else
            (guid_prefer_vecs !== nothing && guid_idx <= length(guid_prefer_vecs)) ?
                Int(guid_prefer_vecs[guid_idx][end]) : 1
        end
        traj_idx = max(1, min(traj_idx, n_bundle_trajs))
        if ndims(guid_ddto_trajs_sims) == 2
            return collect(guid_ddto_trajs_sims[guid_idx, traj_idx])
        else
            return collect(guid_ddto_trajs_sims[guid_idx][traj_idx])
        end
    end

    error_history = Float64[]
    for j in 1:n_steps
        guid_idx = _current_guid_bundle_index(guid_update_times, sim_time[j])
        traj = get_tracked_traj(guid_idx)
        push!(error_history, _guidance_error(sim_state[1:3, j], traj))
    end
    return error_history
end

"""
Plot guidance tracking error vs. sim time for one or more flight cases.

`cases` is a Vector of NamedTuples `(run_data, label, color)`. Each `run_data`
is the contents of a flight-results pickle (must include `sim_state`,
`sim_time`, `guid_ddto_trajs_sims`, `guid_update_times`, optionally
`guid_prefer_vecs`).

`align_start=true` shifts each run's time so it starts at 0 (using
`minimum(t)`, robust to non-monotonic sim_time vectors). Set to false to plot
against raw sim_time.
"""
function plot_flight_tracking_error(
        cases;
        interactive = true,
        save_path = nothing,
        title = "Guidance Tracking Error",
        align_start = true,
    )
    f = Figure(size = (900, 500))
    ax = Axis(
        f[1, 1],
        xlabel = align_start ? "Time since maneuver start [s]" : "Sim time [s]",
        ylabel = "Tracking error [m]",
        title = title,
        xautolimitmargin = (0, 0),
        topspinevisible = true,
        rightspinevisible = true,
        xgridvisible = true,
        ygridvisible = true,
    )

    for case in cases
        run_data = case.run_data
        label = case.label
        color = case.color
        err = compute_guidance_error_history(run_data)
        t_raw = collect(run_data["sim_time"])
        t = vec(t_raw)
        is_monotonic = all(diff(t) .>= 0)
        @info "case stats" label length(t) first_t=t[1] last_t=t[end] min_t=minimum(t) max_t=maximum(t) monotonic=is_monotonic n_err=length(err)
        if align_start
            t = t .- minimum(t)
        end
        lines!(ax, t, err; color = color, linewidth = 2.5, label = label)
    end

    axislegend(ax; position = :rt, framevisible = true)

    if save_path !== nothing
        CairoMakie.save(save_path, f)
        println("Saved: ", save_path)
    end
    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    end
    return f
end


"""
Plot Z position (sim_state[3, :]) vs raw sim_time for one or more flight cases.
Intended for diagnostics — no time shifting/alignment is applied.

`show_indices = true` (default) adds a small text label next to each scatter
point with that sample's original 1-based index, so outlier indices can be
read off the plot interactively (zoom in with GLMakie to read).

`show_resampled = true` overlays each case's `resampled` field (another data
dict with `sim_time` and `sim_state`) as a dashed line in the same color.
Cases without a `resampled` field are silently skipped.

`cases` is the same `(run_data, label, color)` vector used by
`plot_flight_tracking_error`, optionally with an additional `resampled` field.
"""
function plot_flight_z_position(
        cases;
        interactive = true,
        save_path = nothing,
        title = "Z position vs sim_time (raw)",
        show_indices = true,
        index_fontsize = 8,
        show_resampled = false,
    )
    f = Figure(size = (1200, 700))
    ax = Axis(
        f[1, 1],
        xlabel = "sim_time (raw) [s]",
        ylabel = "Z position [m]",
        title = title,
        xautolimitmargin = (0, 0),
        topspinevisible = true,
        rightspinevisible = true,
        xgridvisible = true,
        ygridvisible = true,
    )

    for case in cases
        run_data = case.run_data
        label = case.label
        color = case.color
        sim_state = run_data["sim_state"]
        sim_time = run_data["sim_time"]
        t = vec(collect(sim_time))
        z = vec(collect(sim_state[3, :]))
        is_monotonic = all(diff(t) .>= 0)
        @info "z-pos case stats" label length_t=length(t) length_z=length(z) first_t=t[1] last_t=t[end] min_t=minimum(t) max_t=maximum(t) monotonic=is_monotonic min_z=minimum(z) max_z=maximum(z)
        lines!(ax, t, z; color = color, linewidth = 2.5, label = label)
        scatter!(ax, t, z; color = color, markersize = 5, strokewidth = 0.5, strokecolor = :black)
        if show_indices
            positions = [Point2f(t[i], z[i]) for i in 1:length(t)]
            text!(ax,
                positions;
                text = string.(1:length(t)),
                fontsize = index_fontsize,
                color = :black,
                offset = (4, 4),
                align = (:left, :bottom),
            )
        end

        if show_resampled && hasproperty(case, :resampled) && case.resampled !== nothing
            r = case.resampled
            t_r = vec(collect(r["sim_time"]))
            z_r = vec(collect(r["sim_state"][3, :]))
            lines!(ax, t_r, z_r; color = color, linewidth = 1.5, linestyle = :dash, label = "$label (resampled)")
        end
    end

    axislegend(ax; position = :rt, framevisible = true)

    if save_path !== nothing
        CairoMakie.save(save_path, f)
        println("Saved: ", save_path)
    end
    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    end
    return f
end


"""
Numerical time derivative of `y` with respect to `t`, returning a vector of the
same length. Endpoints use one-sided differences; interior points use
second-order central differences. `t` need not be uniformly spaced.
"""
function _ddt(y::AbstractVector, t::AbstractVector)
    n = length(y)
    @assert length(t) == n "_ddt: y and t must be the same length"
    dy = similar(y, Float64)
    if n == 1
        dy[1] = 0.0
        return dy
    end
    dy[1] = (y[2] - y[1]) / (t[2] - t[1])
    dy[end] = (y[end] - y[end-1]) / (t[end] - t[end-1])
    for i in 2:n-1
        dy[i] = (y[i+1] - y[i-1]) / (t[i+1] - t[i-1])
    end
    return dy
end


"""
Hampel-style local outlier detection. For each sample, compute the median and
MAD of a `window`-sized neighborhood (`half = window ÷ 2` on each side). A
point is flagged as an outlier when its deviation from the local median exceeds
`k × 1.4826 × MAD` (1.4826 makes MAD a consistent estimator of σ for Gaussian
noise). Returns the keep-mask. Robust to bursts of outliers up to ~half the
window size.
"""
function _outlier_mask(y::AbstractVector; window::Int = 11, k::Real = 4.0)
    n = length(y)
    keep = trues(n)
    half = window ÷ 2
    for i in 1:n
        lo = max(1, i - half)
        hi = min(n, i + half)
        w = @view y[lo:hi]
        med = median(w)
        mad_val = median(abs.(w .- med))
        threshold = k * 1.4826 * mad_val
        if threshold > 0 && abs(y[i] - med) > threshold
            keep[i] = false
        end
    end
    return keep
end

"""
Final paper figure: two vertically stacked subplots sharing the time axis.
Top — Vertical velocity, computed numerically as `d/dt(sim_state[3, :])`
      (central differences with one-sided endpoints), with local Hampel-style
      outlier rejection applied to remove spike samples that numerical
      differentiation amplifies. The sign follows whatever Z convention
      `sim_state` uses (typically positive-up, so descent is negative).
Bot — Guidance tracking error (via `compute_guidance_error_history`) vs sim_time.

For each case, the **resampled** data is used for both subplots. If a case has
a `resampled` field (as in `cases_z`), that field is preferred; otherwise
`run_data` is used directly (as in `cases_resampled`, where `run_data` already
holds the resampled dict).

Outlier rejection on the velocity uses a local median/MAD test with
`outlier_window` samples and threshold `outlier_k × 1.4826 × MAD`. Set
`reject_velocity_outliers = false` to disable.

If a case carries a `velocity_limit` field (e.g. `-0.2`), a transparent
constraint band is drawn in the case's color from a floor below the data up
to the limit, with a dotted horizontal line at the limit value. Cases without
this field get no band. `constraint_alpha` controls the band opacity.
"""
function plot_TRO_flightdata_figure(
        cases;
        interactive = true,
        save_path = nothing,
        title = nothing,
        reject_velocity_outliers = true,
        outlier_window::Int = 100,
        outlier_k::Real = 2.0,
        constraint_alpha::Real = 0.15,
    )
    f = Figure(size = (800, 600))
    ax_vel = Axis(
        f[1, 1],
        ylabel = "Descent velocity [m/s]",
        title = isnothing(title) ? "" : title,
        xautolimitmargin = (0, 0),
        topspinevisible = true,
        rightspinevisible = true,
        xgridvisible = true,
        ygridvisible = true,
        xticklabelsvisible = false,
    )
    ax_err = Axis(
        f[2, 1],
        xlabel = "Maneuver time [s]",
        ylabel = "Tracking error [m]",
        xautolimitmargin = (0, 0),
        topspinevisible = true,
        rightspinevisible = true,
        xgridvisible = true,
        ygridvisible = true,
    )
    linkxaxes!(ax_vel, ax_err)
    rowgap!(f.layout, 8)

    # First pass: assemble per-case plot data so we know the y-extent before
    # drawing constraint bands behind the lines.
    plot_data = Any[]
    for case in cases
        run_data = (hasproperty(case, :resampled) && case.resampled !== nothing) ?
                    case.resampled : case.run_data
        label = case.label
        color = case.color
        sim_state = run_data["sim_state"]
        sim_time = run_data["sim_time"]
        if size(sim_state, 1) < 3
            @warn "plot_TRO_flightdata_figure: sim_state has fewer than 3 rows; skipping $(label)" size=size(sim_state)
            continue
        end
        t = vec(collect(sim_time))
        z = vec(collect(sim_state[3, :]))
        vz = _ddt(z, t)
        if reject_velocity_outliers
            mask = _outlier_mask(vz; window = outlier_window, k = outlier_k)
            n_removed = count(.!mask)
            @info "plot_TRO_flightdata_figure: velocity outlier rejection" label n_samples=length(vz) n_removed=n_removed window=outlier_window k=outlier_k
            t_v  = t[mask]
            vz_v = vz[mask]
        else
            t_v  = t
            vz_v = vz
        end
        err = compute_guidance_error_history(run_data)
        push!(plot_data, (case = case, label = label, color = color,
                          t_full = t, t_v = t_v, vz_v = vz_v, err = err))
    end

    # Compute the velocity-axis y-limits explicitly so the constraint boxes can
    # be drawn to the *exact* visible bottom (no gap between box and axis edge)
    # while still leaving a visual margin above/below the data. Includes the
    # velocity limits in the extents so the boundary lines are always in view.
    all_vz = isempty(plot_data) ? Float64[] : vcat([d.vz_v for d in plot_data]...)
    limit_vals = Float64[Float64(d.case.velocity_limit) for d in plot_data
                        if hasproperty(d.case, :velocity_limit) && d.case.velocity_limit !== nothing]
    y_vals = vcat(all_vz, limit_vals)
    if isempty(y_vals)
        y_lo_plot, y_hi_plot = -1.0, 1.0
    else
        y_lo_data, y_hi_data = extrema(y_vals)
        pad = 0.05 * max(y_hi_data - y_lo_data, eps())  # 5% margin (Makie default)
        y_lo_plot, y_hi_plot = y_lo_data - pad, y_hi_data + pad
    end
    ylims!(ax_vel, y_lo_plot, y_hi_plot)

    # Background: constraint bands and limit lines (drawn first so they sit
    # behind the velocity traces). The band floor sits exactly at the y-axis
    # bottom; the limit line spans only the case's own time range.
    for d in plot_data
        case = d.case
        if hasproperty(case, :velocity_limit) && case.velocity_limit !== nothing
            limit = Float64(case.velocity_limit)
            t_lo, t_hi = extrema(d.t_full)
            band!(ax_vel,
                  [t_lo, t_hi],
                  [y_lo_plot, y_lo_plot],
                  [limit, limit];
                  color = (d.color, constraint_alpha))
            lines!(ax_vel, [t_lo, t_hi], [limit, limit];
                   color = d.color, linewidth = 1.5, linestyle = :dot)
        end
    end

    # Foreground: velocity & tracking-error traces.
    for d in plot_data
        lines!(ax_vel, d.t_v,    d.vz_v; color = d.color, linewidth = 2.5, label = d.label)
        lines!(ax_err, d.t_full, d.err;  color = d.color, linewidth = 2.5)
    end

    axislegend(ax_vel; position = :rb, framevisible = true)

    if save_path !== nothing
        CairoMakie.save(save_path, f)
        println("Saved: ", save_path)
    end
    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    end
    return f
end
