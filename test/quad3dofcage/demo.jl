using DDTOSCP
using LinearAlgebra
include("scenarios.jl")
include("plots.jl")

params = scenario_obstacles_hard()
# params = scenario_obstacles_easy()
# params = scenario_no_obstacles()

scp_solution, scp_simulation, ddtoscp_solution, ddtoscp_simulation = solve(params)
build_plots([scp_solution], [scp_simulation], [ddtoscp_solution], [ddtoscp_simulation], params)

;