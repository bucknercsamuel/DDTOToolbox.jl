using DDTOSCP
using LinearAlgebra
using DataFrames
using Pandas
# using PythonCall
include("plots.jl")

# Specify relevant paths (hardcoded for now)
path_mc  = "/home/samuelbuckner/HALO_ROS/halo_mc_results/test_20250210_205353" 

# Parse mc data
data = Dict()
for (_,_,files) in walkdir(path_mc)
    for file in files
        file_ = replace(file, ".pkl" => "")
        contents = split(file_,("_"))
        spec = contents[2]
        spec = replace(spec, "gr" => "Gr-")
        spec = replace(spec, "Inf" => "∞")
        spec = replace(spec, "ddto" => "DDTO")
        if ~haskey(data,spec)
            data[spec] = []
        end
        data_ = read_pickle(joinpath(path_mc,file))
        append!(data[spec], [data_])
    end
end

# Plot results
with_theme(theme2d; fontsize=fontsize) do
    screens = [
        plot_mc_statistics(data)
    ]
    hold_interactive(screens)
end
# plot_mc_statistics(data, interactive=false)