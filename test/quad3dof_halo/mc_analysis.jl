using DDTOToolbox
using LinearAlgebra
using DataFrames
using Pandas
using Pickle
include("plots/plot_defaults.jl")
include("plots/plot_mc_statistics.jl")
include("plots/plot_mc_pareto_front.jl")

# Specify relevant paths (hardcoded for now)
mapid = "map3"
path_mc  = "/data/$(mapid)_testFinal/"
path_mc = abspath(@__DIR__)*abspath(path_mc)

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

# Additional data processing
function compute_ate(data; track_ahead_alt=1., downsample=10)
    # Compute ATE
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
function compute_cum_energy(data)
    cum_energy = 0.
    for k = 1:length(data["sim_time"])-1
        delta_time = data["sim_time"][k+1] - data["sim_time"][k]
        cur_control = data["sim_control"][:,k]
        cum_energy += norm(cur_control)^(3/2) * delta_time
    end
    return cum_energy
end
function compute_radius_at_cutoff(data; window_size=10)
    # Find the largest radius at the cutoff time
    cutoff_time = data["guid_update_times"][end]
    cutoff_idx = findlast(x -> x <= cutoff_time, data["sim_time"])
    largest_radius = maximum(data["sim_targs_radii"][:,cutoff_idx-window_size:cutoff_idx])
    return largest_radius
end
function compute_safety_at_cutoff(data; safe_radius=1.)
    # Find the largest radius at the cutoff time
    cutoff_time = data["guid_update_times"][end]
    cutoff_idx = findlast(x -> x <= cutoff_time, data["sim_time"])
    largest_radius = maximum(data["sim_targs_radii"][:,cutoff_idx])
    if largest_radius >= safe_radius
        return 1
    else
        return 0
    end
end

for (spec, data_) in data
    for (idx, data__) in enumerate(data_)
        # Average trajectory error (ATE)
        label = "ATE"
        if ~haskey(data__,label) # only compute if not already computed
            # compute ATE by taking every point in trajectory and computing the distance to the closest point in the reference trajectory at that time
            data[spec][idx][label] = compute_ate(data[spec][idx])
            println("Computed ATE for ($spec,$idx): $(data[spec][idx][label])")
        end
        # Num recomputations
        label = "num_recomputations"
        if ~haskey(data__,label) # only compute if not already computed
            data[spec][idx][label] = length(data__["guid_update_times"]) - 1 # first update time is the initial time, so we don't count it
            println("Computed num_recomputations for ($spec,$idx): $(data[spec][idx][label])")
        end
        # Cumulative energy
        label = "cum_energy"
        if ~haskey(data__,label) # only compute if not already computed
            data[spec][idx][label] = compute_cum_energy(data[spec][idx])
            println("Computed cum_energy for ($spec,$idx): $(data[spec][idx][label])")
        end
        # Largest radius at cutoff time
        label = "radius_at_cutoff"
        if ~haskey(data__,label) # only compute if not already computed
            data[spec][idx][label] = compute_radius_at_cutoff(data[spec][idx])
            println("Computed radius_at_cutoff for ($spec,$idx): $(data[spec][idx][label])")
        end
        # # Safety at cutoff time
        # label = "safety_at_cutoff"
        # if ~haskey(data__,label) # only compute if not already computed
        #     data[spec][idx][label] = compute_safety_at_cutoff(data[spec][idx])
        #     println("Computed safety_at_cutoff for ($spec,$idx): $(data[spec][idx][label])")
        # end
    end
end

# Make function that iteratively calls plot_mc_statistics for a collection of labels
function plot_mc_statistics_collection(data, labels, saturations; interactive=true, mapid="")
    for (label, saturation) in zip(labels, saturations)
        plot_mc_statistics(data, label; saturation=saturation, interactive=interactive, mapid=mapid)
    end
end

# Plot results
interactive = false
labels_mc = ["cum_thrust", "cum_energy", "ATE", "num_recomputations", "radius_at_cutoff"]
saturations_mc = [450, 1500, Inf, Inf, Inf]
with_theme(theme2d) do
    screens = [
        plot_mc_statistics_collection(data, labels_mc, saturations_mc; interactive=interactive, mapid=mapid),
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