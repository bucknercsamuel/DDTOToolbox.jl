using DDTOToolbox
using LinearAlgebra
using Printf
include("plots.jl")

params_cvx = DIntegrator2DoFParams()
sol_cvx, sim_cvx, sol_ddtocvx, sim_ddtocvx = solve_cvx(params_cvx)

params_scp = DIntegrator2DoFParams()
sol_scp, sim_scp, sol_ddtoscp, sim_ddtoscp, _ = solve(params_scp)

build_plots_compare_cvx_scp(
    sol_ddtocvx, 
    sim_ddtocvx, 
    sol_ddtoscp,
    sim_ddtoscp,
    params_cvx,
    params_scp)