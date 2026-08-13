using DDTOToolbox
using LinearAlgebra
using DataFrames
using Pandas
using Pickle
using Printf
include("plots/data_proc_functions.jl")
include("plots/plot_defaults.jl")
include("plots/plot_paper_demo_traj_history.jl")
include("anim_halo_maneuver.jl")

# Paths
local_path = abspath(@__DIR__)
map_rel_path = "map_lookups\\maps\\dunes_test_hard\\lookup_table.pkl"
run_names = ["ddto2", "grOne2", "grInf1"]

# map_rel_path = "map_lookups\\maps\\msl_test_easy\\lookup_table.pkl"
# run_name = "iter91_ddto"
# demo_rel_path = "data\\map1_testFinal\\$(run_name).pkl"

# Parse map data
if !@isdefined map_data
    println("Loading map data...")
    map_data = Dict()
    map_data["zlookup"] = read_pickle(joinpath(local_path, map_rel_path))
    println("Map data loaded successfully")
end

summary = []
for run_name in run_names
    demo_rel_path = "data\\$(run_name).pkl"

    println("Loading demo data ($run_name)...")
    run_data = read_pickle(joinpath(local_path, demo_rel_path))
    println("Demo data loaded successfully")

    # Integrals on raw logged samples (no fill), same convention as T123.
    score_run!(run_data, map_data)
    t = run_data["sim_time"]
    flight_time = t[end] - t[1]
    cum_thrust = run_data["cum_thrust"]
    induced_energy = run_data["induced_energy"]
    mechanical_energy = run_data["mechanical_energy"]
    radius_at_cutoff = run_data["radius_at_cutoff"]
    altitude_at_cutoff = run_data["altitude_at_cutoff"]
    num_computations = run_data["num_recomputations"]

    println("Run name: $run_name")
    @printf("  Flight time: %.2f s\n", flight_time)
    @printf("  Cumulative thrust: %.2f N·s\n", cum_thrust)
    @printf("  Cumulative induced energy: %.2f J\n", induced_energy)
    @printf("  Mechanical energy: %.2f J\n", mechanical_energy)
    @printf("  Radius at cutoff: %.2f m\n", radius_at_cutoff)
    @printf("  Altitude at cutoff: %.2f m\n", altitude_at_cutoff)
    println("  Number of recomputations: $num_computations")
    println("  Landed: $(run_data["landed"])  committed: $(run_data["committed"])")

    push!(summary, (
        run = run_name,
        flight_time = flight_time,
        cum_thrust = cum_thrust,
        induced_energy = induced_energy,
        mechanical_energy = mechanical_energy,
        radius_at_cutoff = radius_at_cutoff,
        altitude_at_cutoff = altitude_at_cutoff,
        num_recomputations = num_computations,
        landed = run_data["landed"],
        committed = run_data["committed"],
    ))

    # HACKY: need to resolve logging bug in HALO_ROS
    # Fill sim gaps at guidance updates (solve duration + interpolated state/control) so animation doesn't teleport
    fill_sim_gaps!(run_data)

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
end

println("\nDemo run comparison (raw logged integrals):")
println("  Run       t [s]   Thrust [N·s]   Energy [J]   Mech [J]   R [m]   AGL [m]   #recomp")
for r in summary
    @printf("  %-8s %6.2f   %12.2f   %10.2f   %8.2f   %5.2f   %7.2f   %7d\n",
        r.run, r.flight_time, r.cum_thrust, r.induced_energy, r.mechanical_energy,
        r.radius_at_cutoff, r.altitude_at_cutoff, r.num_recomputations)
end
