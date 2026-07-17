using CairoMakie
using GLMakie
using Colors
include("plot_defaults.jl")

# ----------------------------------------------------------------------------
# Shared helpers
# ----------------------------------------------------------------------------

"""
    _setup_terrain_map(map_data; kwargs...) -> (f, ax, hm)

Build the base figure/axis with the terrain heatmap rendered. Both
`plot_terrain_map` and `plot_terrain_per_algorithm` call into this helper to
keep the map-rendering logic (bounds derivation, downsampling, axis styling)
in exactly one place.

Indexing convention: `map_data["zlookup"]` is keyed by `(NED-y, NED-x)`
(i.e. first index = East, second index = North), matching the existing 3D
terrain code (`plot_paper_demo_traj_history.jl`, `data_proc_functions.jl`).
"""
function _setup_terrain_map(map_data;
        downsample = 500,
        terrain_alpha = 1.0,
        title = "",
    )
    zlookup = map_data["zlookup"]

    # NED bounds from the lookup keys.
    east_min  = typemax(Int); east_max  = typemin(Int)
    north_min = typemax(Int); north_max = typemin(Int)
    for k in keys(zlookup)
        e, n = k[1], k[2]
        e < east_min  && (east_min  = e); e > east_max  && (east_max  = e)
        n < north_min && (north_min = n); n > north_max && (north_max = n)
    end

    # Downsampled integer grids along each NED axis (integer steps avoid any
    # rounding ambiguity when querying the lookup dict).
    xs = collect(north_min:downsample:north_max)  # NED-X (North) along plot x-axis
    ys = collect(east_min:downsample:east_max)    # NED-Y (East)  along plot y-axis
    zs = [zlookup[e, n] for n in xs, e in ys]     # zs[i,j] = altitude at (xs[i], ys[j])

    f = Figure(size = (750, 700))
    ax = Axis(f[1, 1];
        xlabel = "NED-X / North [m]",
        ylabel = "NED-Y / East [m]",
        # title = title,
        aspect = DataAspect(),
        xgridvisible = false,
        ygridvisible = false,
        topspinevisible = true,
        rightspinevisible = true,
    )
    cmap = range(parse(Colorant, "seashell"), stop=parse(Colorant, "saddlebrown"), length=100)
    hm = heatmap!(ax, xs, ys, zs; colormap = cmap, alpha = terrain_alpha)
    return f, ax, hm
end

"""
    _finalize_terrain_plot(f, interactive, save_path, default_name) -> screen | f

Either pop the figure into a GLMakie window (interactive) or save via
CairoMakie. Shared between the terrain-overlay plots so they all behave the
same way w.r.t. the `interactive` / `save_path` arguments.
"""
function _finalize_terrain_plot(f, interactive, save_path, default_name)
    if interactive
        GLMakie.activate!()
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    else
        CairoMakie.activate!()
        out_path = save_path === nothing ?
            joinpath(fig_path, default_name * fig_ext) :
            save_path
        CairoMakie.save(out_path, f)
        return f
    end
end

# Filename-safe rendering of a spec name (e.g. "Gr-∞" -> "Gr_inf").
_spec_slug(spec::AbstractString) = replace(spec, "∞" => "inf", "-" => "_", " " => "_")

# ----------------------------------------------------------------------------
# Public plot functions
# ----------------------------------------------------------------------------

"""
    plot_terrain_map(map_data; kwargs...)

Plot a 2D top-down projection of the terrain map onto the lateral X-Y NED plane.
When `data` is supplied, the initial X-Y position of every unsafe run
(`safe_run == false`) is overlaid as a scatter marker color-coded by algorithm.

# Arguments
- `downsample::Int`: stride (in meters) used when sampling the lookup.
- `interactive::Bool`: GLMakie window vs. CairoMakie save.
- `mapid::String`: appended to the default save filename.
- `save_path`: explicit save path; overrides the default constructed from `mapid`.
- `colormap`: any Makie-recognized colormap.
- `terrain_alpha::Real`: opacity of the terrain heatmap layer in [0,1].
- `data`: MC result set (as built in `mc_analysis.jl`) keyed by algorithm/spec name.
- `unsafe_marker_size::Real`: marker size for the unsafe-run overlay.
"""
function plot_terrain_map(map_data;
        interactive = false,
        downsample = 4,
        save_path = nothing,
        mapid = "",
        colormap = :terrain,
        terrain_alpha = 1.0,
        data = nothing,
        unsafe_marker_size = 12,
    )
    f, ax, hm = _setup_terrain_map(map_data;
        downsample = downsample, terrain_alpha = terrain_alpha)
    Colorbar(f[1, 2], hm; label = "Elevation (NED-Z) [m]")

    # Overlay unsafe MC runs' initial X-Y positions, color-coded by algorithm.
    if data !== nothing
        spec_colors = Dict(
            "Graph-DDTO" => :dodgerblue3,
            "Gr-1"       => :indianred3,
            "Gr-∞"       => :orange3,
        )
        spec_order = ["Graph-DDTO", "Gr-1", "Gr-∞"]
        for spec in spec_order
            haskey(data, spec) || continue
            runs = data[spec]
            unsafe_idxs = findall(k -> !convert(Bool, runs[k]["safe_run"]), 1:length(runs))
            isempty(unsafe_idxs) && continue
            xs_init = [runs[k]["sim_state"][1, 1] for k in unsafe_idxs]
            ys_init = [runs[k]["sim_state"][2, 1] for k in unsafe_idxs]
            color = get(spec_colors, spec, :black)
            scatter!(ax, xs_init, ys_init;
                color = color,
                strokecolor = :black,
                strokewidth = 1,
                markersize = unsafe_marker_size,
                label = "$(spec) ($(length(unsafe_idxs)) unsafe)",
            )
        end
        axislegend(ax; position = :rt, framevisible = true)
    end

    return _finalize_terrain_plot(f, interactive, save_path, "terrain_map_$(mapid)")
end

"""
    plot_terrain_per_algorithm(map_data, runs; kwargs...)

Same base terrain projection as `plot_terrain_map`, but for a single algorithm.
Every run in `runs` is scattered at its initial X-Y NED position, colored green
if `safe_run == true` and red otherwise. Typically called once per algorithm
(Graph-DDTO, Gr-1, Gr-∞) to produce a side-by-side view of where each one
succeeds and fails.

# Arguments
- `runs::Vector`: list of run dicts (one algorithm's slice of `data_non67`).
- `spec::String`: algorithm name; used as the axis title and in the default filename.
- All other kwargs mirror `plot_terrain_map`.
"""
function plot_terrain_per_algorithm(map_data, runs;
        spec = "",
        interactive = false,
        downsample = 4,
        save_path = nothing,
        mapid = "",
        colormap = :terrain,
        terrain_alpha = 0.6,
        marker_size = 10,
    )
    f, ax, hm = _setup_terrain_map(map_data;
        downsample = downsample, terrain_alpha = terrain_alpha,
        title = spec)
    Colorbar(f[1, 2], hm; label = "Elevation (NED-Z) [m]")

    safe_mask   = Bool[convert(Bool, runs[k]["safe_run"]) for k in 1:length(runs)]
    safe_idxs   = findall(safe_mask)
    unsafe_idxs = findall(.!safe_mask)

    # Draw unsafe first so the green "safe" markers sit on top in dense regions.
    if !isempty(unsafe_idxs)
        xs_unsafe = [runs[k]["sim_state"][1, 1] for k in unsafe_idxs]
        ys_unsafe = [runs[k]["sim_state"][2, 1] for k in unsafe_idxs]
        scatter!(ax, xs_unsafe, ys_unsafe;
            color = :red,
            strokecolor = :black,
            strokewidth = 1,
            markersize = marker_size,
            label = "Unsafe",
        )
    end
    if !isempty(safe_idxs)
        xs_safe = [runs[k]["sim_state"][1, 1] for k in safe_idxs]
        ys_safe = [runs[k]["sim_state"][2, 1] for k in safe_idxs]
        scatter!(ax, xs_safe, ys_safe;
            color = :limegreen,
            strokecolor = :black,
            strokewidth = 1,
            markersize = marker_size,
            label = "Safe",
        )
    end
    axislegend(ax; position = :rt, framevisible = true)

    default_name = "terrain_runs_$(_spec_slug(spec))_$(mapid)"
    return _finalize_terrain_plot(f, interactive, save_path, default_name)
end
