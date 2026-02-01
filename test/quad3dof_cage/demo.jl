using DDTOToolbox
using LinearAlgebra
using Random
include("scenarios.jl")
include("plots/plot_trajs.jl")

# Set the random seed
Random.seed!(123)

# Test DDTO-LEX
lex = true
# lex = false

params = scenario_obstacles_hard(lex)
# params = scenario_obstacles_easy()
# params = scenario_no_obstacles()
# params = scenario_obstacles_hard_random_targets(lex=lex, n_targets=4, min_distance_from_obstacle=0.01)

# Solve
if lex
    scp_sol, scp_sim, ddtoscp_sol, ddtoscp_sim = solve_lex(params)
else
    scp_sol, scp_sim, ddtoscp_sol, ddtoscp_sim = solve(params)
end

# Plot results
screens = []
interactive = false
with_theme(theme2d; fontsize=fontsize) do
    # push!(screens, plot_trajs([scp_sol],     [scp_sim],     params; interactive=interactive, ddto=false))
    push!(screens, plot_trajs([ddtoscp_sol], [ddtoscp_sim], params; interactive=interactive))
end
if interactive
    hold_interactive(screens)
end
;