using DDTOToolbox
using LinearAlgebra
using DataFrames
using Pandas
using Pickle
include("plots/data_proc_functions.jl")
include("plots/plot_defaults.jl")
include("plots/plot_paper_demo_traj_history.jl")
include("anim_halo_maneuver.jl")

# Paths
local_path = abspath(@__DIR__)
map_rel_path = "map_lookups\\maps\\dunes_test_hard\\lookup_table.pkl"
# demo_rel_path = "data\\paper_demo_data.pkl"
run_name = "ddto2"
# run_name = "grOne2"
# run_name = "grInf1"
demo_rel_path = "data\\$(run_name).pkl"

# map_rel_path = "map_lookups\\maps\\msl_test_easy\\lookup_table.pkl"
# run_name = "iter91_ddto"
# demo_rel_path = "data\\map1_testFinal\\$(run_name).pkl"

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

# # Get some extra data for printout
# cum_thrust = compute_cum_thrust(run_data)
# cum_momentum = compute_cum_momentum(run_data)
# cum_mechanical_power = compute_cum_mechanical_power(run_data)
# radius_at_cutoff = compute_radius_at_cutoff(run_data)
# altitude_at_cutoff = get_altitude_at_cutoff(run_data)[1]
# num_computations = length(run_data["guid_update_times"]) - 1

# # Print out results of the run
# println("Run name: $run_name")
# println("Cumulative thrust: $cum_thrust")
# println("Cumulative momentum: $cum_momentum")
# println("Cumulative mechanical power: $cum_mechanical_power")
# println("Radius at cutoff: $radius_at_cutoff")
# println("Altitude at cutoff: $altitude_at_cutoff")
# println("Number of computations: $num_computations")

# Plot results (set save_path to a path string to record video with 10% progress prints)
with_theme(theme3d) do
    # result = animate_paper_demo_traj_history(run_data, map_data; fps=30, playback_speed=5.0, loop=true, show_time_label=false)
    # result = animate_paper_demo_traj_history(run_data, map_data;
    #     fps=30,
    #     playback_speed=5.0,
    #     loop=false,
    #     show_time_label=false,
    #     # camera_rotation_rate=.1,
    #     show_guidance_error=true,
    #     azel=(pi/4,pi/8),
    #     map_downsample=1,
    #     save_path=joinpath(local_path, "figures", "$(run_name).mp4"),
    # )
    # plot_paper_demo_traj_history(run_data, map_data; 
    #     interactive=false, 
    #     azel=(3*pi/4,pi/6),
    #     save_path=joinpath(local_path, "figures", "$(run_name).png"),
    # )
    plot_paper_demo_traj_bundle(run_data, map_data;
        bundle_idx = 2,
        interactive = false,
        azel = (3pi/4, pi/6),
        defer_color = :darkgray,
        branch_colors = [:red, :blue, :forestgreen, :orange],
        save_path = joinpath(local_path, "figures", "$(run_name)_bundle1.png"),
    )
end
;