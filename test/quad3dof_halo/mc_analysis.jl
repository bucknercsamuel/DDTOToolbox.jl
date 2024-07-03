using DDTOSCP
using DataFrames
using LinearAlgebra
using Pandas
include("plots.jl")

# Specify relevant paths (hardcoded for now)
path_mc  = "C:/users/chris/HALO_DDTOSCP/AirSim/mc_results/" 
path_map = "C:/users/chris/HALO_DDTOSCP/AirSim/maps/msl_test_easy/"

# Parse mc data
data = Dict()
for (_,_,files) in walkdir(path_mc)
    for file in files
        file_ = replace(file, ".pkl" => "")
        contents = split(file_,("_"))
        spec = contents[2]
        spec = replace(spec, "greedy" => "Gr-")
        spec = replace(spec, "ddto" => "DDTO")
        if ~haskey(data,spec)
            data[spec] = []
        end
        df = read_pickle(path_mc*file)
        append!(data[spec], [df])
    end
end

# Obtain map data
# map_data = Dict()
# map_data["xlims"] = [-904.178, 986.822]
# map_data["ylims"] = [-859.034, 905.966]
# map_data["zlookup"] = read_pickle(path_map*"lookup_table.pkl")
# map_data["imgpath"] = path_map*"cam.png"

# Plot results
with_theme(theme3d; fontsize=fontsize) do
    screens = [
        plot_mc_trajs(data["DDTO"], map_data),
        # plot_mc_trajs(data["Gr-0.1"], map_data),
        # plot_mc_trajs(data["Gr-100"], map_data),
        plot_mc_statistics(data)
    ]
    hold_interactive(screens)
end