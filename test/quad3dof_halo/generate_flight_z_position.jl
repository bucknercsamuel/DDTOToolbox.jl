using DDTOToolbox
using LinearAlgebra
using GLMakie
using CairoMakie
using Colors
include("plots/plot_defaults.jl")
include("plots/plot_flight_tracking_error.jl")
include("load_flight_data.jl")

local_path = abspath(@__DIR__)
data_dir = joinpath(local_path, "data", "flight_tests")

# Inspection plot: raw data with index labels on each point. Zoom in (mouse
# wheel + drag) on the GLMakie window to read indices near any jump/excursion,
# then copy those indices into the EXCLUDE_* lists at the top of
# `generate_flight_tracking_error.jl`.
println("Loading sim_results.pkl...")
sim_data  = load_flight_pickle(joinpath(data_dir, "sim_results.pkl"))
println("Loading slow_descent_results.pkl...")
slow_data = load_flight_pickle(joinpath(data_dir, "slow_descent_results.pkl"))
println("Loading fast_descent_results.pkl...")
fast_data = load_flight_pickle(joinpath(data_dir, "fast_descent_results.pkl"))
println("All data loaded.")

cases = [
    (run_data = sim_data,  label = "Simulated",    color = colorant"gray"),
    (run_data = slow_data, label = "Slow Descent", color = colorant"orange"),
    (run_data = fast_data, label = "Fast Descent", color = colorant"skyblue"),
]

# Set to true to overlay each scatter point with its 1-based sample index, so
# you can read off the indices of erroneous points and copy them into the
# EXCLUDE_* lists in `generate_flight_tracking_error.jl`. Zoom in (mouse wheel
# + drag) on the GLMakie window to read labels in crowded regions.
SHOW_INDICES = true

with_theme(theme2d) do
    plot_flight_z_position(
        cases;
        interactive = true,
        show_indices = SHOW_INDICES,
        index_fontsize = 8,
        # save_path = joinpath(local_path, "figures", "flight_z_position.svg"),
    )
end
;
