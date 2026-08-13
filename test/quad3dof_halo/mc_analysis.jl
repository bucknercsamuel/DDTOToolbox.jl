using DDTOToolbox
using LinearAlgebra
using DataFrames
using Pandas
using Pickle
using Statistics
using Printf
using PrettyTables
include("plots/plot_defaults.jl")
include("plots/data_proc_functions.jl")
include("plots/plot_mc_statistics.jl")
include("plots/plot_mc_pareto_front.jl")
include("plots/plot_terrain_map.jl")

# Get map name from ID
map_id_to_name = Dict(
    "map1" => "msl_test_easy",
    "map2" => "msl_test_easy",
    "map3" => "dunes_test_hard",
)
alg_order = ["Gr-1", "Gr-∞", "Graph-DDTO"]
local_path = abspath(@__DIR__)

# `get_altitude_at_cutoff` reads this as a global (not an argument).
map_data = Dict()

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# IQR mean/std/outlier count on error_code == 1 runs, matching plot_mc_statistics.
function iqr_metric_stats(runs, label)
    idx_feas = findall(τ -> τ == 1, [runs[k]["error_code"] for k in 1:length(runs)])
    data = [runs[k][string(label)] for k in idx_feas]
    Q1, _, Q3 = quantile(data, [0.25, 0.5, 0.75])
    IQR = Q3 - Q1
    n_outliers = count(x -> (x < Q1 - 1.5 * IQR) || (x > Q3 + 1.5 * IQR), data)
    return (mean = mean(data), std = std(data), n_outliers = n_outliers)
end

function fmt_mean_std_out(s)
    @sprintf("%.2f ± %.2f (%d)", s.mean, s.std, s.n_outliers)
end

function print_mc_comparison_table(records)
    n = length(records)
    table = Matrix{String}(undef, n, 6)
    thrust_means  = Vector{Float64}(undef, n)
    radius_means  = Vector{Float64}(undef, n)
    divert_means  = Vector{Float64}(undef, n)
    success_pcts  = Vector{Float64}(undef, n)

    for (i, r) in enumerate(records)
        map_label = (i % 3 == 1) ? "Map #$(r.mapnum)" : ""
        table[i, 1] = map_label
        table[i, 2] = r.spec
        table[i, 3] = fmt_mean_std_out(r.thrust)
        table[i, 4] = fmt_mean_std_out(r.radius)
        table[i, 5] = fmt_mean_std_out(r.diverts)
        table[i, 6] = @sprintf("%.2f", r.success_pct)
        thrust_means[i] = r.thrust.mean
        radius_means[i] = r.radius.mean
        divert_means[i] = r.diverts.mean
        success_pcts[i] = r.success_pct
    end

    function best_in_map(i, values, sense)
        gstart = 3 * div(i - 1, 3) + 1
        rows = gstart:(gstart + 2)
        idx = sense == :min ? argmin(values[rows]) : argmax(values[rows])
        return i == rows[idx]
    end

    hl = TextHighlighter((d, i, j) -> begin
        j == 3 ? best_in_map(i, thrust_means, :min) :
        j == 4 ? best_in_map(i, radius_means, :max) :
        j == 5 ? best_in_map(i, divert_means, :min) :
        j == 6 ? best_in_map(i, success_pcts, :max) :
        false
    end, crayon"bold")

    println()
    pretty_table(
        table;
        column_labels = [
            "Map",
            "Algorithm",
            "Cumulative Thrust [N·s]",
            "Safe Radius [m]",
            "# Diverts",
            "% Success",
        ],
        highlighters = [hl],
        alignment = [:l, :l, :c, :c, :c, :c],
        fit_table_in_display_horizontally = false,
        display_size = (-1, -1),
    )
end

# Make function that iteratively calls plot_mc_statistics for a collection of labels
function plot_mc_statistics_collection(data, labels, saturations; interactive=false, mapid="")
    for (label, saturation) in zip(labels, saturations)
        plot_mc_statistics(data, label; saturation=saturation, interactive=interactive, mapid=mapid)
    end
end

# Plot a per-iteration scalar label (e.g. "altitude_at_cutoff", "cum_energy") for each
# algorithm/spec. X-axis is the MC iteration index, Y-axis is the value of `label`.
function plot_mc_per_iteration(data, label; interactive=false, mapid="", ylabel=label, save=true)
    f = Figure(size=(700, 400))
    ax = Axis(f[1, 1], xlabel="MC iteration", ylabel=ylabel,
              xgridvisible=false, ygridvisible=false,
              topspinevisible=true, rightspinevisible=true)

    spec_colors = Dict(
        "Graph-DDTO" => :dodgerblue3,
        "Gr-1"       => :indianred3,
        "Gr-∞"       => :orange3,
    )

    for spec in sort(collect(keys(data)))
        runs = data[spec]
        ys = [runs[k][label] for k in 1:length(runs)]
        xs = collect(1:length(ys))
        color = get(spec_colors, spec, :gray)
        scatter!(ax, xs, ys; color=color, markersize=8, label=spec)
    end
    axislegend(ax; position=:rt)

    if interactive
        GLMakie.activate!()
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    elseif save
        CairoMakie.activate!()
        CairoMakie.save(joinpath(fig_path, "mc_per_iter_$(label)_$(mapid)"*fig_ext), f)
    end
    return f
end

# ---------------------------------------------------------------------------
# Per-map processing (current configuration, applied to every map)
# ---------------------------------------------------------------------------

table_records = []
interactive = false
labels_mc = ["cum_thrust", "induced_energy", "mechanical_energy", "ATE", "num_recomputations", "radius_at_cutoff"]
saturations_mc = [450, Inf, Inf, Inf, Inf, Inf]

for (mapnum, mapid) in enumerate(["map1", "map2", "map3"])
    println("="^72)
    println("Map $mapid")
    println("="^72)

    mapname = map_id_to_name[mapid]

    # Specify relevant paths (hardcoded for now)
    path_mc  = "/data/$(mapid)_testFinal/"
    path_mc = abspath(@__DIR__)*abspath(path_mc)
    map_rel_path = "map_lookups\\maps\\$(mapname)\\lookup_table.pkl"

    # Parse mc data
    data = Dict()
    nfiles = 0
    for (_,_,files) in walkdir(path_mc)
        for file in files
            nfiles += 1
            file_ = replace(file, ".pkl" => "")
            contents = split(file_,("_"))
            spec = contents[2]
            spec = replace(spec, "gr" => "Gr-")
            spec = replace(spec, "Inf" => "∞")
            spec = replace(spec, "ddto" => "Graph-DDTO")
            if ~haskey(data,spec)
                data[spec] = []
            end
            data_ = read_pickle(joinpath(path_mc,file))
            append!(data[spec], [data_])
        end
    end
    println("Loaded $nfiles files from $path_mc")
    for spec in alg_order
        haskey(data, spec) && println("  $spec: $(length(data[spec])) runs")
    end

    # Parse map data
    println("Loading map data...")
    map_data["zlookup"] = read_pickle(joinpath(local_path, map_rel_path))
    println("Map data loaded successfully")

    # Additional data processing
    invalid_runs = Dict()
    for (spec, data_) in data
        invalid_runs[spec] = []
        for (idx, data__) in enumerate(data_)
            # # Fill in sim gaps
            # println("Filling in sim gaps...")
            # fill_sim_gaps!(data__)
            # println("Sim gaps filled successfully")

            # Validate the run
            if !validate_run(data__, map_data)
                push!(invalid_runs[spec], idx)
                data[spec][idx]["error_code"] = 67 # invalid run caught during post-analysis error code
            end

            # Cumulative thrust
            label = "cum_thrust"
            if ~haskey(data__,label) # only compute if not already computed
                data[spec][idx][label] = compute_cum_thrust(data[spec][idx])
            end

            # Induced energy
            label = "induced_energy"
            if ~haskey(data__,label) # only compute if not already computed
                data[spec][idx][label] = compute_induced_energy(data[spec][idx])
            end

            # Mechanical energy
            label = "mechanical_energy"
            if ~haskey(data__,label) # only compute if not already computed
                data[spec][idx][label] = compute_mechanical_energy(data[spec][idx])
            end

            # Average trajectory error (ATE)
            label = "ATE"
            if ~haskey(data__,label) # only compute if not already computed
                # data[spec][idx][label] = compute_ate(data[spec][idx])
                data[spec][idx][label] = 0. # not using for now
            end

            # Num recomputations
            label = "num_recomputations"
            if ~haskey(data__,label) # only compute if not already computed
                data[spec][idx][label] = length(data__["guid_update_times"]) - 1 # first update time is the initial time, so we don't count it
            end

            # Largest radius at cutoff time
            label = "radius_at_cutoff"
            if ~haskey(data__,label) # only compute if not already computed
                data[spec][idx][label] = compute_radius_at_cutoff(data[spec][idx])
            end

            # Cutoff altitude
            label = "altitude_at_cutoff"
            if ~haskey(data__,label) # only compute if not already computed
                data[spec][idx][label] = get_altitude_at_cutoff(data[spec][idx])
            end

            # Safe run
            label = "safe_run"
            if ~haskey(data__,label) # only compute if not already computed
                data[spec][idx][label] = compute_safety_of_run(data[spec][idx])
            end
        end
    end

    # Invalidate all runs on the superset of invalid runs
    invalid_runs_union = unique(vcat([invalid_runs[spec] for spec in keys(invalid_runs)]...))
    for (spec, data_) in data
        for (idx, data__) in enumerate(data_)
            if idx in invalid_runs_union
                data[spec][idx]["error_code"] = 67 # invalid run caught during post-analysis error code
            end
        end
    end
    println("Invalid-run union size: $(length(invalid_runs_union))")

    # Build filtered data subset that excludes runs flagged with error_code 67 (i.e.
    # runs caught as invalid during post-analysis). Downstream safety reporting and
    # plots operate on this subset so they all share a consistent denominator.
    data_non67 = Dict()
    for spec in keys(data)
        idx_non67_runs = findall(x -> x != 67, [data[spec][k]["error_code"] for k in 1:length(data[spec])])
        data_non67[spec] = data[spec][idx_non67_runs]
    end

    # For each spec, print out the percentage of runs that are safe using key safe_run
    for spec in alg_order
        haskey(data_non67, spec) || continue
        runs = data_non67[spec]
        num_runs = length(runs)
        num_valid_runs = count(r -> r["safe_run"] == true, runs)
        println("$(spec): $(num_valid_runs)/$(num_runs) ($(num_valid_runs/num_runs*100)%)")
    end

    # Collect table rows (paper order: Gr-1, Gr-∞, Graph-DDTO)
    for spec in alg_order
        push!(table_records, (
            mapid = mapid,
            mapnum = mapnum,
            spec = spec,
            thrust = iqr_metric_stats(data[spec], "cum_thrust"),
            radius = iqr_metric_stats(data[spec], "radius_at_cutoff"),
            diverts = iqr_metric_stats(data[spec], "num_recomputations"),
            success_pct = begin
                runs = data_non67[spec]
                count(r -> r["safe_run"] == true, runs) / length(runs) * 100
            end,
        ))
    end

    # Plot results
    with_theme(theme2d) do
        screens = [
            # plot_mc_statistics_collection(data, labels_mc, saturations_mc; interactive=interactive, mapid=mapid),
            # plot_mc_per_iteration(data, "altitude_at_cutoff"; interactive=interactive, mapid=mapid,
            #     ylabel="Altitude at cutoff [m]"),
            # plot_mc_statistics(data, "cum_thrust"; saturation=450, interactive=interactive, mapid=mapid),
            # plot_mc_pareto_front(data,
            #     "cum_thrust", "radius_at_cutoff";
            #     xlabel="Cumulative thrust [N]",
            #     ylabel="Cutoff Safety Radius [m]",
            #     n=3, interactive=interactive, label=mapid,
            #     region_type=:kde,
            #     percentiles=[90],
            #     outlier_threshold_1 = 450,
            #     pareto_dir_1 = :decreasing,
            #     pareto_dir_2 = :increasing
            # ),
            # plot_terrain_map(map_data; interactive=interactive, mapid=mapid,
            #     terrain_alpha=0.6, data=data_non67),
            [plot_terrain_per_algorithm(map_data, data_non67[spec];
                spec=spec, interactive=interactive, mapid=mapid, terrain_alpha=1, downsample=10)
                for spec in ["Graph-DDTO", "Gr-1", "Gr-∞"] if haskey(data_non67, spec)]...,
        ]
        if interactive
            hold_interactive(screens)
        end
    end
end

print_mc_comparison_table(table_records)
