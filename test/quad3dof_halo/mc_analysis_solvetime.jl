using LinearAlgebra
using Pandas
using Statistics
using Printf
using PrettyTables
include("plots/data_proc_functions.jl")

# Same configuration as mc_analysis.jl: pair by iteration ID across all three
# algorithms, drop the unusable-log union, no gap-fill. Solve times are every
# guidance recompute on those paired usable flights, pooled over all maps.
map_id_to_name = Dict(
    "map1" => "msl_test_easy",
    "map2" => "msl_test_easy",
    "map3" => "dunes_test_hard",
)
alg_order = ["Gr-1", "Gr-∞", "Graph-DDTO"]
local_path = abspath(@__DIR__)
map_data = Dict()

function run_solve_times(data)
    t = if haskey(data, "solve_times")
        data["solve_times"]
    elseif haskey(data, "guid_update_solve_time")
        data["guid_update_solve_time"]
    else
        return Float64[]
    end
    v = Float64.(vec(collect(t)))
    return filter(x -> isfinite(x) && x >= 0, v)
end

function solvetime_stats(times)
    isempty(times) && return (n=0, n_flights=0, mean=NaN, std=NaN, p95=NaN, worst=NaN)
    return (
        n = length(times),
        mean = mean(times),
        std = std(times),
        p95 = quantile(times, 0.95),
        worst = maximum(times),
    )
end

function fmt_sec(x)
    isfinite(x) || return "n/a"
    return @sprintf("%.3f", x)
end

pooled = Dict(s => Float64[] for s in alg_order)
n_flights = Dict(s => 0 for s in alg_order)

for mapid in ["map1", "map2", "map3"]
    println("="^72)
    println("Map $mapid")
    println("="^72)

    mapname = map_id_to_name[mapid]
    path_mc = joinpath(@__DIR__, "data", "$(mapid)_testFinal")
    map_rel_path = joinpath("map_lookups", "maps", mapname, "lookup_table.pkl")

    data = Dict()
    nfiles = 0
    for (_, _, files) in walkdir(path_mc)
        for file in files
            endswith(file, ".pkl") || continue
            nfiles += 1
            contents = split(replace(file, ".pkl" => ""), "_")
            spec = replace(replace(replace(contents[2], "gr" => "Gr-"), "Inf" => "∞"), "ddto" => "Graph-DDTO")
            iter = parse(Int, replace(contents[1], "iter" => ""))
            haskey(data, spec) || (data[spec] = [])
            data_ = read_pickle(joinpath(path_mc, file))
            data_["mc_iter"] = iter
            append!(data[spec], [data_])
        end
    end
    println("Loaded $nfiles files from $path_mc")

    println("Loading map data...")
    map_data["zlookup"] = read_pickle(joinpath(local_path, map_rel_path))
    println("Map data loaded successfully")

    unusable_iters = Dict(s => Set{Int}() for s in alg_order)
    for (spec, runs) in data
        for run in runs
            run["usable"] = usable_log(run, map_data)
            if !run["usable"]
                push!(unusable_iters[spec], run["mc_iter"])
            end
        end
    end

    unusable_union = union((unusable_iters[s] for s in alg_order)...)
    common = intersect((Set(r["mc_iter"] for r in data[s]) for s in alg_order)...)
    common = setdiff(common, unusable_union)
    println("Unusable-iter union size: $(length(unusable_union))")
    println("Paired usable triples: $(length(common))")

    for spec in alg_order
        kept = filter(r -> r["usable"] && (r["mc_iter"] in common), data[spec])
        n_flights[spec] += length(kept)
        n_solves = 0
        for run in kept
            ts = run_solve_times(run)
            append!(pooled[spec], ts)
            n_solves += length(ts)
        end
        println("  $spec: $(length(kept)) flights, $n_solves solves")
    end
end

println()
println("Guidance solve times [s], pooled over map1–map3")
println("Sample: paired usable triples, all recomputes (not success-filtered)")

n = length(alg_order)
table = Matrix{String}(undef, n, 7)
means = Vector{Float64}(undef, n)
p95s = Vector{Float64}(undef, n)
worsts = Vector{Float64}(undef, n)
for (i, spec) in enumerate(alg_order)
    s = solvetime_stats(pooled[spec])
    table[i, 1] = spec
    table[i, 2] = string(n_flights[spec])
    table[i, 3] = string(s.n)
    table[i, 4] = fmt_sec(s.mean)
    table[i, 5] = fmt_sec(s.std)
    table[i, 6] = fmt_sec(s.p95)
    table[i, 7] = fmt_sec(s.worst)
    means[i] = s.mean
    p95s[i] = s.p95
    worsts[i] = s.worst
end

hl = TextHighlighter((d, i, j) -> begin
    j == 4 ? i == argmin(means) :
    j == 6 ? i == argmin(p95s) :
    j == 7 ? i == argmin(worsts) :
    false
end, crayon"bold")

pretty_table(
    table;
    column_labels = [
        "Algorithm",
        "# Flights",
        "# Solves",
        "Mean [s]",
        "Std [s]",
        "95th pct [s]",
        "Worst [s]",
    ],
    highlighters = [hl],
    alignment = [:l, :c, :c, :c, :c, :c, :c],
    fit_table_in_display_horizontally = false,
    display_size = (-1, -1),
)
