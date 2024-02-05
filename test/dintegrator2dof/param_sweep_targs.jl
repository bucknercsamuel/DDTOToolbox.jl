using DDTOSCP
using LinearAlgebra
using Statistics
using Random
using Printf
include("plots.jl")

MersenneTwister(0)
params = DIntegrator2DoFParams()

function generate_random_targets(N::Int, radius, vertex)
    rf_targs = Matrix(undef, 2, N)
    for j = 1:N
        r_targ = radius * rand(Float64)
        θ_targ = 2 * pi * rand(Float64)
        rf_targs[:,j] = vertex + r_targ*cos(θ_targ)*[1;0] + r_targ*sin(θ_targ)*[0;1]
    end
    return rf_targs
end

# Baseline params
ϵ = 0.3
params.w_obj = 1e0
params.w_ctrl = 1e3
params.w_trust = 1e1

# Generate n targets
N = 100
rf_targs = generate_random_targets(N, 20, [30;30])
zf_targs = vcat(rf_targs, zeros(2,N))

# Update params
params.n_targs = N
params.zf_targs = zf_targs
params.λ_targs = collect(1:N)
params.T_targs = collect(1:N)
params.τ_targs = zeros(N)
params.α_targs = ones(N)
params.ϵ_targs = ϵ*ones(N)

# Solve
_,_,sol,sim = solve(params)

# Plot a solution
build_plots_single([sol],[sim],params)