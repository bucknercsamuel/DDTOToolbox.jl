using DDTOToolbox
using LinearAlgebra
using DataFrames
using Pandas
using Pickle
using PyCall
const pypickle = pyimport("pickle")
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

# Output path: same as path_mc but with "_proc" appended to the leaf folder name
path_proc = abspath(@__DIR__)*abspath("/data/$(mapid)_testFinal_proc/")

# Parse mc data (mirrors mc_analysis.jl, but also tracks original file names so we
# can write each processed entry back out under the same name)
data = Dict()
filenames = Dict()
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
            filenames[spec] = String[]
        end
        data_ = read_pickle(joinpath(path_mc,file))
        append!(data[spec], [data_])
        push!(filenames[spec], file)
    end
end

# Fill sim gaps and save each file to the _proc folder, mirroring the original layout
mkpath(path_proc)
for (spec, data_) in data
    for (idx, data__) in enumerate(data_)
        file = filenames[spec][idx]
        println("Filling sim gaps for ($spec, $idx) -> $file")
        fill_sim_gaps!(data__)

        # Save to mirrored path. Pickle.jl can't serialize multi-dim Julia arrays,
        # so we go through Python's pickle via PyCall: PyCall auto-converts Julia
        # arrays to numpy arrays, matching what read_pickle originally returned.
        out_path = joinpath(path_proc, file)
        f = pybuiltin("open")(out_path, "wb")
        try
            pypickle.dump(data__, f)
        finally
            f.close()
        end
        println("Saved: $out_path")
    end
end

println("Post-processing complete. Processed data written to: $path_proc")
