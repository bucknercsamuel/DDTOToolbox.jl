using DDTOSCP
using LinearAlgebra
using Statistics
using Random
using Printf
include("plots.jl")

MersenneTwister(0)
params = DIntegrator2DoFParams(autogen_targs=true, autogen_targ_count=50)

# Modify params
params.w_obj_sing = 1e0
params.w_obj_ddto = 1e0
params.w_ctrl = 1e4
params.w_trust = 1e2
params.scp_iters = 10

# Solve
_,_,sol,sim = solve(params)

# Plot a solution
build_plots_single([sol],[sim],params)