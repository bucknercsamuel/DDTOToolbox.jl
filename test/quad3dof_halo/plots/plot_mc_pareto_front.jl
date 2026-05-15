using CairoMakie
using Colors
using LinearAlgebra
using Statistics
include("../../utils/plot_utils.jl")
include("plot_defaults.jl")

"""
Evaluate a 2D Gaussian KDE on a regular grid.

Returns `(x_grid, y_grid, density)` where `density[i,j]` corresponds to the
point `(x_grid[i], y_grid[j])`. Bandwidth defaults to Silverman's rule of thumb
applied per-axis. `gridpad` is the number of bandwidths of padding beyond the
data extent in each axis.
"""
function _kde_2d(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real};
                 gridsize::Int = 120, bandwidth = nothing, gridpad::Real = 3)
    n = length(xs)
    n == 0 && error("_kde_2d: empty input")
    σx = std(xs); σy = std(ys)
    hx, hy = if bandwidth === nothing
        # Silverman's rule of thumb (per-axis)
        (1.06 * (σx == 0 ? 1e-9 : σx) * n^(-1/5),
         1.06 * (σy == 0 ? 1e-9 : σy) * n^(-1/5))
    elseif bandwidth isa Tuple || bandwidth isa AbstractVector
        (Float64(bandwidth[1]), Float64(bandwidth[2]))
    else
        (Float64(bandwidth), Float64(bandwidth))
    end
    xmin = minimum(xs) - gridpad*hx; xmax = maximum(xs) + gridpad*hx
    ymin = minimum(ys) - gridpad*hy; ymax = maximum(ys) + gridpad*hy
    x_grid = collect(LinRange(xmin, xmax, gridsize))
    y_grid = collect(LinRange(ymin, ymax, gridsize))
    density = zeros(Float64, gridsize, gridsize)
    norm_const = 1.0 / (n * 2π * hx * hy)
    @inbounds for i in 1:gridsize, j in 1:gridsize
        gx = x_grid[i]; gy = y_grid[j]
        s = 0.0
        for k in 1:n
            dx = (gx - xs[k]) / hx
            dy = (gy - ys[k]) / hy
            s += exp(-0.5 * (dx*dx + dy*dy))
        end
        density[i,j] = s * norm_const
    end
    return x_grid, y_grid, density
end

"""
Compute a vector of density-level thresholds corresponding to each of the
requested HDR (highest-density region) coverage `percentiles` (0-100).

A point `(x,y)` is considered inside the `p`-percent HDR iff `density(x,y) ≥
level(p)`. The level is found by sorting the grid density values in
descending order and accumulating mass until the target fraction is reached.
"""
function _kde_levels_from_percentiles(x_grid, y_grid, density, percentiles)
    dx = x_grid[2] - x_grid[1]
    dy = y_grid[2] - y_grid[1]
    dA = dx * dy
    vals = sort(vec(density); rev = true)
    total = sum(vals) * dA
    levels = similar(percentiles, Float64)
    for (k, p) in enumerate(percentiles)
        target = (p/100) * total
        cum = 0.0
        level = vals[end]
        for v in vals
            cum += v * dA
            if cum >= target
                level = v
                break
            end
        end
        levels[k] = level
    end
    return levels
end

"""
Compute the 2D convex hull (counter-clockwise) of the points `(xs, ys)` using
Andrew's monotone chain algorithm. Returns a vector of `Point2f` forming the
closed hull (first point repeated at the end).
"""
function _convex_hull_2d(xs::AbstractVector{<:Real}, ys::AbstractVector{<:Real})
    pts = [(Float64(xs[i]), Float64(ys[i])) for i in eachindex(xs)]
    pts = unique(pts)
    length(pts) < 3 && return [Point2f(p...) for p in pts]
    sort!(pts; by = p -> (p[1], p[2]))

    cross2d(o, a, b) = (a[1]-o[1])*(b[2]-o[2]) - (a[2]-o[2])*(b[1]-o[1])

    # Lower hull
    lower = Vector{Tuple{Float64,Float64}}()
    for p in pts
        while length(lower) >= 2 && cross2d(lower[end-1], lower[end], p) <= 0
            pop!(lower)
        end
        push!(lower, p)
    end

    # Upper hull
    upper = Vector{Tuple{Float64,Float64}}()
    for p in reverse(pts)
        while length(upper) >= 2 && cross2d(upper[end-1], upper[end], p) <= 0
            pop!(upper)
        end
        push!(upper, p)
    end

    hull = vcat(lower[1:end-1], upper[1:end-1])
    push!(hull, hull[1]) # close the ring
    return [Point2f(p...) for p in hull]
end

"""
    plot_mc_pareto_front(data, label1, label2; region_type=:ellipse, n=3, ...)

Scatter plot of two MC metrics across the three solution types
(Graph-DDTO, Gr-1, Gr-∞) with a per-type uncertainty region drawn behind the
scatter points. Colors match the convention used in `plot_mc_statistics`
(dodgerblue / indianred / orange).

Arguments:
- `data`   : dict keyed by solution type ("Graph-DDTO", "Gr-1", "Gr-∞")
- `label1` : data key for the x-axis metric (e.g. "cum_thrust")
- `label2` : data key for the y-axis metric (e.g. "largest_radius_at_cutoff")

Keyword arguments:
- `region_type`         : `:ellipse` for an `n`-sigma confidence ellipse,
                          `:convex_hull` for the convex hull of the points,
                          or `:kde` for kernel density estimate contours at the
                          supplied `percentiles`
- `n`                   : number of standard deviations (only used for `:ellipse`)
- `percentiles`         : vector of HDR coverage percentiles in 0-100 (only used
                          for `:kde`). Example: `[50, 90]` draws two contours
                          enclosing 50% and 90% of the density mass per type.
- `kde_gridsize`        : grid resolution per axis for the KDE evaluation
- `kde_bandwidth`       : optional bandwidth override for the Gaussian KDE;
                          either a scalar (isotropic) or `(hx, hy)` tuple
- `interactive`         : show an interactive GLMakie window instead of saving
- `label`               : extra suffix appended to the saved figure filename
- `outlier_threshold_1` : if set, drop runs with `label1 > outlier_threshold_1`
                          (same convention as `plot_mc_statistics`)
- `outlier_threshold_2` : if set, drop runs with `label2 > outlier_threshold_2`
- `pareto_dir_1`        : `:increasing` or `:decreasing` — direction of
                          desirability for `label1`. If set, draws a small
                          arrow just above the x-axis pointing the desirable way.
- `pareto_dir_2`        : same as `pareto_dir_1` but for `label2` (drawn just
                          right of the y-axis).
- `arrow_pad`           : fraction of data range reserved as the arrow "lane"
                          on each arrow-hosting side. Smaller = tighter plot
                          and arrow sits closer to the data.
- `arrow_lane_pos`      : 0-1; where inside the lane the arrow sits. 0.0 = right
                          against the data, 1.0 = right against the outer spine.
                          Increase to press the arrow toward the axis bound.
"""
function plot_mc_pareto_front(
        data,
        label1,
        label2;
        region_type = :ellipse,
        n = 3,
        percentiles = [50, 90],
        kde_gridsize = 120,
        kde_bandwidth = nothing,
        xlabel = "",
        ylabel = "",
        interactive = true,
        label = "",
        outlier_threshold_1 = nothing,
        outlier_threshold_2 = nothing,
        pareto_dir_1 = nothing,
        pareto_dir_2 = nothing,
        arrow_pad = 0.13,
        arrow_lane_pos = 0.75,
    )
    region_type in (:ellipse, :convex_hull, :kde) ||
        error("plot_mc_pareto_front: region_type must be :ellipse, :convex_hull, or :kde, got $region_type")
    if region_type == :kde
        all(0 .< percentiles .< 100) ||
            error("plot_mc_pareto_front: percentiles must be in (0, 100); got $percentiles")
    end
    for (nm, pd) in (("pareto_dir_1", pareto_dir_1), ("pareto_dir_2", pareto_dir_2))
        pd === nothing || pd in (:increasing, :decreasing) ||
            error("plot_mc_pareto_front: $nm must be :increasing, :decreasing, or nothing; got $pd")
    end

    # Figure
    f = Figure(size=(600,500))
    defaults = Dict(
        :topspinevisible=>true,
        :rightspinevisible=>true,
        :xgridvisible=>false,
        :ygridvisible=>false,
    )
    ax = Axis(f[1,1], xlabel=xlabel, ylabel=ylabel; defaults...)

    # Color scheme (dark, light) — matches plot_mc_statistics
    colors = [
        (:dodgerblue4, :dodgerblue1),
        (:indianred4, :indianred1),
        (:orange4, :orange1),
    ]
    key_order = ["Graph-DDTO", "Gr-1", "Gr-∞"]

    # Collect per-type point clouds (filtering on feasibility, same as plot_mc_statistics).
    # Also drop:
    #   - runs with non-finite values (Inf / -Inf / NaN) on either axis, since these
    #     break covariance / eigendecomposition for the confidence ellipse
    #   - runs exceeding user-supplied `outlier_threshold_{1,2}` on either axis
    #     (same convention as plot_mc_statistics)
    point_clouds = Dict{String,Tuple{Vector{Float64},Vector{Float64}}}()
    for key in key_order
        value = data[key]
        idx_feas = findall(τ -> τ == 1, [value[k]["error_code"] for k in 1:length(value)])
        xs_raw = Float64[value[k][string(label1)] for k in idx_feas]
        ys_raw = Float64[value[k][string(label2)] for k in idx_feas]

        # Finite-value filter
        finite_mask = isfinite.(xs_raw) .& isfinite.(ys_raw)
        n_nonfinite = length(xs_raw) - count(finite_mask)
        if n_nonfinite > 0
            println("plot_mc_pareto_front: dropping $n_nonfinite non-finite run(s) for $key")
        end
        xs = xs_raw[finite_mask]
        ys = ys_raw[finite_mask]

        # Extreme-outlier filter (per-axis), mirroring plot_mc_statistics behavior
        outlier_mask = trues(length(xs))
        if !isnothing(outlier_threshold_1)
            outlier_mask .&= xs .<= outlier_threshold_1
        end
        if !isnothing(outlier_threshold_2)
            outlier_mask .&= ys .<= outlier_threshold_2
        end
        n_outliers = length(xs) - count(outlier_mask)
        if n_outliers > 0
            println("plot_mc_pareto_front: dropping $n_outliers extreme-outlier run(s) for $key")
        end

        point_clouds[key] = (xs[outlier_mask], ys[outlier_mask])
    end

    # Pass 1: draw uncertainty regions first so they sit behind the scatter.
    # Borders are removed (strokewidth=0) and fills use the brighter color.
    for (iter, key) in enumerate(key_order)
        xs, ys = point_clouds[key]
        length(xs) < 2 && continue

        if region_type == :kde
            x_grid, y_grid, density = _kde_2d(xs, ys;
                gridsize = kde_gridsize, bandwidth = kde_bandwidth)
            # Filled HDR bands at the requested percentile levels. We sort
            # percentiles descending so the largest-coverage (outermost) band
            # is drawn first and the smallest-coverage (innermost) band last,
            # giving the inner region a more opaque stacked fill.
            perc_sorted = sort(collect(percentiles); rev = true)
            levels_sorted = _kde_levels_from_percentiles(x_grid, y_grid, density, perc_sorted)
            dmax = maximum(density)
            for (p, lvl) in zip(perc_sorted, levels_sorted)
                α = 0.15 + 0.35 * (1 - p/100)
                # NOTE: do NOT pass extendlow/extendhigh — that would fill
                # cells outside [lvl, dmax] too (i.e. the whole bounding box).
                contourf!(ax, x_grid, y_grid, density;
                    levels = [lvl, dmax + eps(dmax)],
                    colormap = [(colors[iter][2], α)])
                # Border around the filled band
                contour!(ax, x_grid, y_grid, density;
                    levels = [lvl],
                    color = colors[iter][2],
                    linewidth = 1)
            end
        else
            ring = if region_type == :ellipse
                μ = [mean(xs), mean(ys)]
                Σ = cov(hcat(xs, ys))
                F = eigen(Symmetric(Σ))
                D = F.values
                V = F.vectors
                a = n * sqrt(max(D[1], 0))
                b = n * sqrt(max(D[2], 0))
                θs = LinRange(0, 2π, 200)
                [Point2f(μ .+ V * [a*cos(θ), b*sin(θ)]...) for θ in θs]
            else # :convex_hull
                _convex_hull_2d(xs, ys)
            end

            poly!(ax, ring;
                color = colors[iter][2],
                alpha = 0.25,
                strokewidth = 1)
        end
    end

    # Pass 2: scatter points on top, using the brighter color
    for (iter, key) in enumerate(key_order)
        xs, ys = point_clouds[key]
        isempty(xs) && continue
        scatter!(ax, xs, ys;
            color = colors[iter][2],
            markersize = 5)
    end

    # Pass 3: Pareto-direction arrows that span (nearly) the full axis, drawn
    # inside the axis box next to the corresponding spine. Extra padding is
    # added on the side where an arrow lives so data never overlaps the arrow.
    if !isnothing(pareto_dir_1) || !isnothing(pareto_dir_2)
        all_xs = vcat([point_clouds[k][1] for k in key_order]...)
        all_ys = vcat([point_clouds[k][2] for k in key_order]...)
        if !isempty(all_xs) && !isempty(all_ys)
            xmin, xmax = extrema(all_xs)
            ymin, ymax = extrema(all_ys)
            dx = xmax - xmin; dy = ymax - ymin
            dx = dx == 0 ? 1.0 : dx
            dy = dy == 0 ? 1.0 : dy

            # Asymmetric padding: larger on the side hosting an arrow so the
            # data never overlaps the arrow lane.
            # `arrow_pad` controls the arrow-lane width (fraction of data range);
            # shrink it to press the arrow in toward the outer spine.
            pad_arrow_side = arrow_pad
            pad_other_side = 0.04
            pad_bot = !isnothing(pareto_dir_1) ? pad_arrow_side : pad_other_side
            pad_top = pad_other_side
            pad_lft = !isnothing(pareto_dir_2) ? pad_arrow_side : pad_other_side
            pad_rgt = pad_other_side

            xlo = xmin - pad_lft*dx; xhi = xmax + pad_rgt*dx
            ylo = ymin - pad_bot*dy; yhi = ymax + pad_top*dy
            xlims!(ax, xlo, xhi)
            ylims!(ax, ylo, yhi)
            dxp = xhi - xlo; dyp = yhi - ylo

            # Draw arrows manually as a `lines!` shaft + triangle-marker tip.
            # This sidesteps CairoMakie's `Arrows2D` tip rendering which can
            # produce deformed / rectangular arrowheads in SVG output.
            arrow_color     = :black
            arrow_linewidth = 2
            arrow_tipsize   = 20  # triangle marker size (pixels)

            function _draw_arrow(x_tail, y_tail, x_tip, y_tip, marker)
                lines!(ax, [x_tail, x_tip], [y_tail, y_tip];
                    color = arrow_color, linewidth = arrow_linewidth)
                scatter!(ax, [x_tip], [y_tip];
                    marker = marker, markersize = arrow_tipsize,
                    color = arrow_color, strokewidth = 0)
            end

            # Crossover point of the two arrow lanes (even if one arrow is
            # absent, we can still define it so the other arrow's near end
            # lines up consistently).
            # `arrow_lane_pos` sets where the arrow sits inside its lane:
            # 0.0 → right against the data; 1.0 → right against the outer
            # spine. Larger values press the arrow toward the axis bound.
            y_pos = ymin - arrow_lane_pos * pad_bot * dy
            x_pos = xmin - arrow_lane_pos * pad_lft * dx

            # Margin used at both ends of each arrow. At the far end it's the
            # distance to the plot edge; at the near end it's the distance to
            # the arrow crossover point — by construction these are equal.
            end_margin = 0.04

            # X-axis arrow: horizontal, centered in the bottom pad lane. Its
            # near end is offset from the crossover `x_pos` by `end_margin*dxp`
            # (matching the `end_margin*dxp` gap at the far end), so the two
            # arrows never overlap near the origin.
            if !isnothing(pareto_dir_1)
                x_lo_a = x_pos + end_margin * dxp
                x_hi_a = xhi   - end_margin * dxp
                if pareto_dir_1 == :increasing
                    _draw_arrow(x_lo_a, y_pos, x_hi_a, y_pos, :rtriangle)
                else # :decreasing
                    _draw_arrow(x_hi_a, y_pos, x_lo_a, y_pos, :ltriangle)
                end
            end

            # Y-axis arrow: vertical, centered in the left pad lane. Same
            # symmetric-margin logic as the X-axis arrow.
            if !isnothing(pareto_dir_2)
                y_lo_a = y_pos + end_margin * dyp
                y_hi_a = yhi   - end_margin * dyp
                if pareto_dir_2 == :increasing
                    _draw_arrow(x_pos, y_lo_a, x_pos, y_hi_a, :utriangle)
                else # :decreasing
                    _draw_arrow(x_pos, y_hi_a, x_pos, y_lo_a, :dtriangle)
                end
            end
        end
    end

    # Build the legend manually so each entry matches the per-series region
    # color + alpha. This works around a Makie quirk where auto-generated
    # legend entries for `contourf!` (and other colormapped recipes) render as
    # a neutral gray `PolyElement` because the recipe represents multiple
    # levels, not a single color.
    legend_elements = Any[]
    legend_labels   = String[]
    for (iter, key) in enumerate(key_order)
        xs, _ = point_clouds[key]
        isempty(xs) && continue
        c = colors[iter][2]
        fill_alpha = region_type == :kde ? 0.40 : 0.25
        strokew    = region_type == :kde ? 1    : 1
        push!(legend_elements, [
            PolyElement(color = (c, fill_alpha), strokecolor = c, strokewidth = strokew),
            MarkerElement(marker = :circle, color = c, markersize = 6, strokewidth = 0),
        ])
        push!(legend_labels, key)
    end
    if !isempty(legend_elements)
        axislegend(ax, legend_elements, legend_labels; position = :rt)
    end

    if interactive
        GLMakie.activate!()
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    else
        fname = "pareto_$(label1)_vs_$(label2)"
        if !isempty(label)
            fname *= "_$(label)"
        end
        CairoMakie.activate!()
        CairoMakie.save(joinpath(fig_path, fname*fig_ext), f)
    end
end
