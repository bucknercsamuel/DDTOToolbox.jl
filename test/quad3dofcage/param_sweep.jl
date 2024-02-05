using DDTOSCP
using LinearAlgebra

include("scenarios.jl")
include("plots.jl")

# params = scenario_obstacles_hard()
# params = scenario_obstacles_easy()
params = scenario_no_obstacles()

logrange(x1, x2, n) = (10^y for y in range(log10(x1), log10(x2), length=n))

N = 30
# param_range = collect(range(1e-3,1e0,N))
param_range = collect(logrange(1e-3,1e3,N))
# param_range = collect(logrange(1e-2,1e0,N))
sols = []
sims = []
for j = 1:N
    # params.ϵ_targs = fill(param_range[j], params.n_targs)
    params.α_targs = fill(1, params.n_targs)
    params.α_targs[params.λ_targs[1]] = param_range[j]
    _,_,sol,sim = solve(params)
    push!(sols, sol)
    push!(sims, sim)
end

build_plots([], [], sols, sims, params)