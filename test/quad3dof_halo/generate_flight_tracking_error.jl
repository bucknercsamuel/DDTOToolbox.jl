using DDTOToolbox
using LinearAlgebra
using GLMakie
using CairoMakie
using Colors
include("plots/plot_defaults.jl")
include("plots/plot_flight_tracking_error.jl")
include("plots/plot_paper_demo_traj_history.jl")
include("load_flight_data.jl")

# Minimal stand-in for the dunes/MSL terrain lookup: any integer (i, j) query
# returns the same z value, so `plot_paper_demo_traj_history`'s
# `[zlookup[round(x), round(y)] for x in xs, y in ys]` builds a flat plane.
struct FlatGround
    z::Float64
end
Base.getindex(g::FlatGround, ::Any, ::Any) = g.z

local_path = abspath(@__DIR__)
data_dir = joinpath(local_path, "data", "flight_tests")

# Bad sample indices to remove, per case (1-based into the original sim_time).
# Identify them by running `generate_flight_z_position.jl` (raw data with
# per-point index labels) and reading off the offending indices.
EXCLUDE_SIM  = Int[1,2]
EXCLUDE_SLOW = Int[1,2,3,63,376,420,479,611,878,879,880,969,920,971,972,973,974,975,1099,1131,1291]
EXCLUDE_FAST = Int[1,2,3,4,5,241,256,257,258,384,385,386,397,610]
# EXCLUDE_SLOW = Int[]
# EXCLUDE_FAST = Int[]

println("Loading sim_results.pkl...")
sim_data  = load_flight_pickle(joinpath(data_dir, "sim_results.pkl"))
println("Loading slow_descent_results.pkl...")
slow_data = load_flight_pickle(joinpath(data_dir, "slow_descent_results.pkl"))
println("Loading fast_descent_results.pkl...")
fast_data = load_flight_pickle(joinpath(data_dir, "fast_descent_results.pkl"))
println("All data loaded.")

println("Postprocessing...")
sim_data  = postprocess_flight_data(sim_data;  exclude_indices = EXCLUDE_SIM)
slow_data = postprocess_flight_data(slow_data; exclude_indices = EXCLUDE_SLOW)
fast_data = postprocess_flight_data(fast_data; exclude_indices = EXCLUDE_FAST)

# Simulated run also has logging gaps; close them by shifting time forward.
# Any gap > SIM_GAP_THRESHOLD seconds collapses to a 0-second step and all
# later timestamps shift back by the gap duration. Real flight runs are left
# alone — their sim_time is treated as ground truth modulo manual exclusions.
SIM_GAP_THRESHOLD = 0.025
sim_data = compress_time_gaps(sim_data; threshold = SIM_GAP_THRESHOLD)
println("Postprocessing complete.")

# Resample all three onto a uniform time grid. Anything downstream
# (tracking-error computation, plots) uses the resampled data; the
# cleaned-but-irregular data is kept as `*_clean` for the z-pos overlay
# comparison.
# INTERPOLATION:
#   :cubic   smooth, can overshoot across large gaps
#   :linear  kinked at samples, no overshoot
#   :hybrid  ALPHA*cubic + (1-ALPHA)*linear; ALPHA=1 ≡ :cubic, ALPHA=0 ≡ :linear
# ALPHA is ignored unless INTERPOLATION = :hybrid.
RESAMPLE_DT = 0.01
INTERPOLATION = :hybrid
ALPHA = 0.4
println("Resampling onto uniform $(RESAMPLE_DT)s grid via :$(INTERPOLATION)" *
        (INTERPOLATION === :hybrid ? " (alpha=$ALPHA)" : "") * "...")
sim_clean,  sim_resampled  = sim_data,  resample_flight_data(sim_data;  dt = RESAMPLE_DT, interpolation = :linear)
slow_clean, slow_resampled = slow_data, resample_flight_data(slow_data; dt = RESAMPLE_DT, interpolation = :hybrid, alpha = ALPHA)
fast_clean, fast_resampled = fast_data, resample_flight_data(fast_data; dt = RESAMPLE_DT, interpolation = :hybrid, alpha = ALPHA)
println("Resampling complete.")

# Velocity limits per scenario (m/s, with sign matching the data convention —
# negative = descent). Shown as transparent constraint bands on the velocity
# subplot of the TRO figure.
VELOCITY_LIMIT_SIM  = -0.2
VELOCITY_LIMIT_SLOW = -0.1
VELOCITY_LIMIT_FAST = -0.2

# For the tracking-error plot we want the analysis done on the uniform grid.
cases_resampled = [
    (run_data = sim_resampled,  label = "Simulated",    color = colorant"gray",    velocity_limit = VELOCITY_LIMIT_SIM),
    (run_data = slow_resampled, label = "Slow Descent", color = colorant"orange",  velocity_limit = VELOCITY_LIMIT_SLOW),
    (run_data = fast_resampled, label = "Fast Descent", color = colorant"skyblue", velocity_limit = VELOCITY_LIMIT_FAST),
]

# For the z-pos comparison plot, each case carries both the cleaned irregular
# data (as run_data, shown as dots+line) and the resampled uniform data (as
# resampled, shown as a dashed line overlay).
cases_z = [
    (run_data = sim_clean,  resampled = sim_resampled,  label = "Simulated",    color = colorant"gray"),
    (run_data = slow_clean, resampled = slow_resampled, label = "Slow Descent", color = colorant"orange"),
    (run_data = fast_clean, resampled = fast_resampled, label = "Fast Descent", color = colorant"skyblue"),
]

with_theme(theme2d) do
    # plot_flight_tracking_error(
    #     cases_resampled;
    #     interactive = true,
    #     align_start = false,  # already zeroed by postprocess_flight_data
    #     # save_path = joinpath(local_path, "figures", "flight_tracking_error.svg"),
    # )
    # plot_flight_z_position(
    #     cases_z;
    #     interactive = true,
    #     show_indices = true,
    #     show_resampled = true,
    #     title = "Z position: cleaned (dots+line) vs resampled (dashed)",
    #     # save_path = joinpath(local_path, "figures", "flight_z_position.svg"),
    # )
    # plot_TRO_flightdata_figure(
    #     cases_resampled;
    #     interactive = true,
    #     save_path = joinpath(local_path, "figures", "TRO_flightdata.svg"),
    # )
end

# 3D trajectory-history figure per case, drawn on a flat ground plane. The
# ground is placed a small buffer below each case's lowest sim_state z so the
# vehicle markers never visually intersect the surface.
GROUND_BUFFER = 0.5  # meters below the case's minimum z
with_theme(theme3d) do
    for case in cases_resampled
        z_min = minimum(case.run_data["sim_state"][3, :])
        map_data = Dict("zlookup" => FlatGround(z_min - GROUND_BUFFER))
        fname = replace(lowercase(case.label), " " => "_") * "_traj_history.png"
        plot_paper_demo_traj_history(case.run_data, map_data;
            interactive = false,
            azel = (3pi/4, pi/6),
            show_drone_timelapse = true,
            save_path = joinpath(local_path, "figures", fname),
        )
    end
end
;
