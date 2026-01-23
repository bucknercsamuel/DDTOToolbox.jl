using DDTOToolbox
using LinearAlgebra
include("scenarios.jl")
include("plots.jl")

# params = scenario_obstacles_hard()

# params.a.ctcs_enabled = false
# params.a.w_trust = 1e2
# _,_, sols_ctcs_off, sims_ctcs_off = solve(params)
# params.a.ctcs_enabled = true
# params.a.w_trust = 1e3
# _,_, sols_ctcs_on , sims_ctcs_on  = solve(params)

build_plots(sols_ctcs_off, sims_ctcs_off, sols_ctcs_on, sims_ctcs_on, params; interactive=true)
;