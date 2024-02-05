using DDTOSCP
using LinearAlgebra
include("plots.jl")

params = DIntegrator2DoFParams()
params_cvx = copy(params)
params_scp = copy(params)
cvx_solution, cvx_simulation, ddtocvx_solution, ddtocvx_simulation = solve_cvx(params_cvx)
scp_solution, scp_simulation, ddtoscp_solution, ddtoscp_simulation = solve(params_scp)

# build_plots_compare_cvx_scp(ddtocvx_solution, ddtocvx_simulation, ddtoscp_solution, ddtoscp_simulation, params_cvx, params_scp)
build_plots_single([scp_solution], [scp_simulation], params_scp; ddto=false)
;