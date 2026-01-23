using DDTOToolbox
using LinearAlgebra
using Pandas
using DataFrames
include("plots.jl")

# # Load relevant data
# local_path = "C:\\Users\\sbuck\\Documents\\ACL\\Code\\DDTOSCP.jl\\test\\quad3dof_halo"
# run_data = read_pickle(joinpath(local_path,"data\\paper_demo_data.pkl"))

# # Parse map data
# map_data = Dict()
# map_data["zlookup"] = read_pickle(joinpath(local_path,"map_lookups\\maps\\dunes_test_hard\\lookup_table.pkl"))
# # map_data["imgpath"] = path_map*"cam.png"

# Plot results
with_theme(theme3d; fontsize=fontsize) do
    screens = [
        plot_paper_demo_traj_history(run_data, map_data; interactive=true),
    ]
    hold_interactive(screens)
end
;