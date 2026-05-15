using DDTOToolbox
using LinearAlgebra
using DataFrames
using Pandas
using Pickle
include("plots/plot_defaults.jl")
include("plots/data_proc_functions.jl")
include("plots/plot_mc_statistics.jl")
include("plots/plot_mc_pareto_front.jl")

# Specify map ID
mapid = "map3"

# Get map name from ID
map_id_to_name = Dict(
    "map1" => "msl_test_easy",
    "map2" => "msl_test_easy",
    "map3" => "dunes_test_hard",
)
mapname = map_id_to_name[mapid]

# Specify relevant paths (hardcoded for now)
path_mc  = "/data/$(mapid)_testFinal/"
path_mc = abspath(@__DIR__)*abspath(path_mc)
local_path = abspath(@__DIR__)
map_rel_path = "map_lookups\\maps\\$(mapname)\\lookup_table.pkl"

# Parse mc data
data = Dict()
for (_,_,files) in walkdir(path_mc)
    for file in files
        println("Processing file: $file")
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

# Parse map data
if !@isdefined map_data
    println("Loading map data...")
    map_data = Dict()
    map_data["zlookup"] = read_pickle(joinpath(local_path, map_rel_path))
    println("Map data loaded successfully")
end

# Additional data processing
invalid_runs = Dict()
for (spec, data_) in data
    invalid_runs[spec] = []
    for (idx, data__) in enumerate(data_)
        println("Processing ($spec,$idx)")

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

# For each spec, print out the percentgae of runs that are safe using key safe_run without error code 67
for spec in keys(data)
    idx_non67_runs = findall(x -> x != 67, [data[spec][k]["error_code"] for k in 1:length(data[spec])])
    num_runs = length(idx_non67_runs)
    num_valid_runs = length(findall(x -> x == true, [data[spec][k]["safe_run"] for k in idx_non67_runs]))
    println("$(spec): $(num_valid_runs)/$(num_runs) ($(num_valid_runs/num_runs*100)%)")
end

# Plot results
interactive = false
labels_mc = ["cum_thrust", "induced_energy", "mechanical_energy", "ATE", "num_recomputations", "radius_at_cutoff"]
saturations_mc = [450, Inf, Inf, Inf, Inf, Inf]
with_theme(theme2d) do
    screens = [
        plot_mc_statistics_collection(data, labels_mc, saturations_mc; interactive=interactive, mapid=mapid),
        # plot_mc_per_iteration(data, "altitude_at_cutoff"; interactive=interactive, mapid=mapid,
        #     ylabel="Altitude at cutoff [m]"),
        # plot_mc_statistics(data, "cum_thrust"; saturation=450, interactive=interactive, mapid=mapid),
        plot_mc_pareto_front(data, 
            "cum_thrust", "radius_at_cutoff";
            xlabel="Cumulative thrust [N]",
            ylabel="Cutoff Safety Radius [m]",
            n=3, interactive=interactive, label=mapid,
            region_type=:kde,
            percentiles=[90],
            outlier_threshold_1 = 450,
            pareto_dir_1 = :decreasing,
            pareto_dir_2 = :increasing
        ),
    ]
    if interactive
        hold_interactive(screens)
    end
end