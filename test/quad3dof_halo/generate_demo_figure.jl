using DDTOToolbox
using LinearAlgebra
using DataFrames
using Pandas
using Pickle
include("plots.jl")
include("anim_halo_maneuver.jl")

# Paths
local_path = abspath(@__DIR__)
map_rel_path = "map_lookups\\maps\\dunes_test_hard\\lookup_table.pkl"
# demo_rel_path = "data\\paper_demo_data.pkl"
demo_rel_path = "data\\ddto2.pkl"

# Load data
println("Loading demo data...")
run_data = read_pickle(joinpath(local_path, demo_rel_path))
println("Demo data loaded successfully")

# HACKY: need to resolve logging bug in HALO_ROS
# Fill sim gaps at guidance updates (solve duration + interpolated state/control) so animation doesn't teleport
fill_sim_gaps!(run_data)

# Parse map data
if !@isdefined map_data
    println("Loading map data...")
    map_data = Dict()
    map_data["zlookup"] = read_pickle(joinpath(local_path, map_rel_path))
    # map_data["imgpath"] = path_map*"cam.png"
    println("Map data loaded successfully")
end

# Plot results (set save_path to a path string to record video with 10% progress prints)
with_theme(theme3d) do
    # result = animate_paper_demo_traj_history(run_data, map_data; fps=30, playback_speed=5.0, loop=true, show_time_label=false, save_path=joinpath(local_path, "figures", "demo_figure.mp4"))
    result = animate_paper_demo_traj_history(run_data, map_data; fps=30, playback_speed=1.0, loop=true, show_time_label=false, camera_rotation_rate=.2)
    # screens = [
    #     plot_paper_demo_traj_history(run_data, map_data; interactive=true),
    # ]
    # hold_interactive(screens)
end
;