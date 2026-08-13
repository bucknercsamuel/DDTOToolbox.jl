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
const E_CAP = 2500.0
const R_MIN = 1.0
alg_order = ["Gr-1", "Gr-∞", "Graph-DDTO"]
local_path = abspath(@__DIR__)

# `get_altitude_at_cutoff` reads this as a global (not an argument).
map_data = Dict()

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# IQR mean/std/outlier count on a filtered run list.
function iqr_metric_stats(runs, label)
    data = [r[string(label)] for r in runs if isfinite(r[string(label)])]
    isempty(data) && return (mean=NaN, std=NaN, n_outliers=0, n=0)
    Q1, _, Q3 = quantile(data, [0.25, 0.5, 0.75])
    IQR = Q3 - Q1
    n_outliers = count(x -> (x < Q1 - 1.5 * IQR) || (x > Q3 + 1.5 * IQR), data)
    return (mean=mean(data), std=std(data), n_outliers=n_outliers, n=length(data))
end

function fmt_mean_std_out(s)
    (!isfinite(s.mean) || get(s, :n, 1) == 0) && return "n/a"
    @sprintf("%.2f ± %.2f (%d)", s.mean, s.std, s.n_outliers)
end

function print_mc_comparison_table(records)
    n = length(records)
    table = Matrix{String}(undef, n, 6)
    thrust_means  = Vector{Float64}(undef, n)
    energy_means  = Vector{Float64}(undef, n)
    radius_means  = Vector{Float64}(undef, n)
    success_pcts  = Vector{Float64}(undef, n)

    for (i, r) in enumerate(records)
        map_label = (i % 3 == 1) ? "Map #$(r.mapnum)" : ""
        table[i, 1] = map_label
        table[i, 2] = r.spec
        table[i, 3] = fmt_mean_std_out(r.thrust)
        table[i, 4] = fmt_mean_std_out(r.energy)
        table[i, 5] = fmt_mean_std_out(r.radius)
        table[i, 6] = @sprintf("%.2f (%d/%d)", r.success_pct, r.n_succ, r.n_kept)
        thrust_means[i] = r.thrust.mean
        energy_means[i] = r.energy.mean
        radius_means[i] = r.radius.mean
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
        j == 4 ? best_in_map(i, energy_means, :min) :
        j == 5 ? best_in_map(i, radius_means, :max) :
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
            "Induced energy [J]",
            "Safe Radius [m]",
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
            endswith(file, ".pkl") || continue
            nfiles += 1
            file_ = replace(file, ".pkl" => "")
            contents = split(file_,("_"))
            spec = contents[2]
            spec = replace(spec, "gr" => "Gr-")
            spec = replace(spec, "Inf" => "∞")
            spec = replace(spec, "ddto" => "Graph-DDTO")
            iter = parse(Int, replace(contents[1], "iter" => ""))
            if ~haskey(data,spec)
                data[spec] = []
            end
            data_ = read_pickle(joinpath(path_mc,file))
            data_["mc_iter"] = iter
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

    # Event scoring on raw logs. No gap-fill. AGL/commit are not gates;
    # missing unique radius is imputed to 1 m. Operational = landed + R≥1 + energy cap.
    unusable_iters = Dict(s => Set{Int}() for s in alg_order)
    for (spec, data_) in data
        for (idx, data__) in enumerate(data_)
            score_run!(data__, map_data)
            if !data__["usable"]
                push!(unusable_iters[spec], data__["mc_iter"])
            end
            logged_ok = data__["committed"] && isfinite(data__["radius_at_cutoff"])
            data__["radius_imputed"] = !logged_ok
            if !logged_ok
                data__["radius_at_cutoff"] = 1.0
            end
            data__["cum_energy"] = data__["induced_energy"]
            data__["safe_run"] = data__["landed"] &&
                isfinite(data__["radius_at_cutoff"]) && data__["radius_at_cutoff"] >= R_MIN &&
                data__["induced_energy"] <= E_CAP
        end
    end

    # Trial-level drop: only unusable *logs* (corrupt / start off-map), keyed by iter.
    unusable_union = union((unusable_iters[s] for s in alg_order)...)
    println("Unusable-iter union size: $(length(unusable_union))")

    common = intersect((Set(r["mc_iter"] for r in data[s]) for s in alg_order)...)
    common = setdiff(common, unusable_union)
    println("Paired usable triples: $(length(common))")

    data_paired = Dict()
    for spec in alg_order
        data_paired[spec] = filter(r -> r["usable"] && (r["mc_iter"] in common), data[spec])
    end

    for spec in alg_order
        runs = data_paired[spec]
        num_runs = length(runs)
        num_valid_runs = count(r -> r["safe_run"] == true, runs)
        println("$(spec): $(num_valid_runs)/$(num_runs) ($(num_valid_runs/num_runs*100)%)")
    end

    for spec in alg_order
        kept = data_paired[spec]
        cost = filter(r -> r["safe_run"] == true, kept)
        n_succ = length(cost)
        push!(table_records, (
            mapid = mapid,
            mapnum = mapnum,
            spec = spec,
            thrust = iqr_metric_stats(cost, "cum_thrust"),
            energy = iqr_metric_stats(cost, "cum_energy"),
            radius = iqr_metric_stats(cost, "radius_at_cutoff"),
            n_succ = n_succ,
            n_kept = length(kept),
            success_pct = length(kept) == 0 ? NaN : n_succ / length(kept) * 100,
        ))
    end

    data_non67 = data_paired

    # Plot results: energy violin saturated at 7500 J, shared y-limits across maps.
    with_theme(theme2d) do
        plot_mc_statistics(data_paired, "cum_energy"; saturation=7500.0, ylims=(0.0, 7500.0), interactive=interactive, mapid=mapid)
    end
end

print_mc_comparison_table(table_records)
