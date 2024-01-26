using DDTOSCP
using LinearAlgebra

include("scenarios.jl")
include("plots.jl")

# params = scenario_obstacles_hard()
# params = scenario_obstacles_easy()
params = scenario_no_obstacles()

scp_solutions, scp_simulations, ddtoscp_solutions, ddtoscp_simulations = solve(params)
build_plots(params, scp_solutions, scp_simulations, ddtoscp_solutions, ddtoscp_simulations)
;