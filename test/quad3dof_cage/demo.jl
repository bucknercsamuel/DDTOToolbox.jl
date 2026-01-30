using DDTOToolbox
using LinearAlgebra
include("scenarios.jl")
include("plots.jl")

params = scenario_obstacles_hard()
# params = scenario_obstacles_easy()
# params = scenario_no_obstacles()

# Solve
scp_sol, scp_sim, ddtoscp_sol, ddtoscp_sim = solve(params)
# scp_sol, scp_sim, ddtoscp_sol, ddtoscp_sim = solve_lex(params)

# Plot results
screens = []
interactive = true
with_theme(theme2d; fontsize=fontsize) do
    push!(screens, plot_trajs([scp_sol],     [scp_sim],     params; interactive=interactive, ddto=false))
    push!(screens, plot_trajs([ddtoscp_sol], [ddtoscp_sim], params; interactive=interactive))
end
hold_interactive(screens)
;