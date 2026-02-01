using DDTOToolbox
using LinearAlgebra
using DataFrames
using Pandas
include("plots.jl")

# Specify relevant paths (hardcoded for now)
mapid = "map1"
path_mc  = "/data/$(mapid)_testFinal/"
path_mc = abspath(@__DIR__)*abspath(path_mc)

# Parse mc data
data = Dict()
for (_,_,files) in walkdir(path_mc)
    for file in files
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

# # Additional data processing
# function compute_ate(data; track_ahead_alt=1., downsample=10)
#     # Compute ATE
#     ate = 0.
#     iter = 0
#     for k = 1:downsample:length(data["sim_time"])
#         cur_time = data["sim_time"][k]
#         cur_guid_idx = findlast(x -> x <= cur_time, data["guid_update_times"])
#         guid_defer_idx = data["guid_prefer_vecs"][cur_guid_idx][end]
#         cur_guid = data["guid_ddto_trajs_sims"][cur_guid_idx,guid_defer_idx]
#         cur_state = data["sim_state"][1:3,k]
#         cur_state_shifted = reshape(cur_state + [0,0,-track_ahead_alt],3,1)
#         dists = [norm(cur_state_shifted - cur_guid[1:3,kk]) for kk = 1:size(cur_guid,2)]
#         track_idx = argmin(dists)
#         ate += dists[track_idx]
#         iter += 1
#     end
#     return ate/iter
# end

# for (spec, data_) in data
#     for (idx, data__) in enumerate(data_)
#         # Average trajectory error (ATE)
#         label = "ATE"
#         if ~haskey(data__,label) # only compute if not already computed
#             # compute ATE by taking every point in trajectory and computing the distance to the closest point in the reference trajectory at that time
#             data[spec][idx][label] = compute_ate(data[spec][idx])
#             println("Computed ATE for ($spec,$idx): $(data[spec][idx][label])")
#         end
#         # Num recomputations
#         label = "num_recomputations"
#         if ~haskey(data__,label) # only compute if not already computed
#             data[spec][idx][label] = length(data__["guid_update_times"])
#             println("Computed num_recomputations for ($spec,$idx): $(data[spec][idx][label])")
#         end
#     end
# end

# Plot results
with_theme(theme2d; fontsize=fontsize) do
    screens = [
        plot_mc_statistics(data, interactive=false, label=mapid)
    ]
    hold_interactive(screens)
end