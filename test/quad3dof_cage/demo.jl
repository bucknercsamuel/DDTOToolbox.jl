using DDTOToolbox
using LinearAlgebra
using Random
include("scenarios.jl")
include("plots/plot_trajs.jl")

# Set the random seed
Random.seed!(123)

# Choose scenario
params = scenario_obstacles_hard()
# params = scenario_obstacles_easy()
# params = scenario_no_obstacles()

# Solve
scp_sol, scp_sim, ddtoscp_sol, ddtoscp_sim = solve(params)

# Plot results
screens = []
interactive = false
with_theme(theme2d; fontsize=fontsize) do
    push!(screens, plot_trajs([scp_sol],     [scp_sim],     params; interactive=interactive, ddto=false))
    push!(screens, plot_trajs([ddtoscp_sol], [ddtoscp_sim], params; interactive=interactive))
end
if interactive
    hold_interactive(screens)
end
;