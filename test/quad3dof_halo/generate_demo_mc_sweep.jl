using DDTOToolbox
using LinearAlgebra
using DataFrames
using Pandas
using Pickle
include("plots/plot_defaults.jl")
include("anim_halo_maneuver.jl")

# Paths
local_path = abspath(@__DIR__)
map_rel_path = "map_lookups\\maps\\dunes_test_hard\\lookup_table.pkl"

# Parse map data
if !@isdefined map_data
    println("Loading map data...")
    map_data = Dict()
    map_data["zlookup"] = read_pickle(joinpath(local_path, map_rel_path))
    println("Map data loaded successfully")
end

idx_mc = [1,2,3,4,5,6,7,8]
for idx in idx_mc
    run_name = "iter$(idx)_ddto"
    demo_rel_path = "data\\map3_testFinal\\$(run_name).pkl"

    # Load run data
    println("Loading demo data...")
    run_data = read_pickle(joinpath(local_path, demo_rel_path))
    println("Demo data loaded successfully")

    # HACKY: need to resolve logging bug in HALO_ROS
    # Fill sim gaps at guidance updates (solve duration + interpolated state/control) so animation doesn't teleport
    fill_sim_gaps!(run_data)

    # Plot results (set save_path to a path string to record video with 10% progress prints)
    with_theme(theme3d) do
        result = animate_paper_demo_traj_history(run_data, map_data; 
            fps=30,
            playback_speed=5.0,
            loop=false,
            show_time_label=false,
            show_guidance_error=false,
            azel=(pi/4,pi/8),
            map_downsample=1,
            save_path=joinpath(local_path, "figures", "$(run_name).mp4"),
        )
    end
end
;