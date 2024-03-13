using DDTOSCP
using LinearAlgebra
using Statistics
using Random
using Printf
include("plots.jl")

# MersenneTwister(0)

# # Test params
# N_targ_perms = 10
# N_iters_per_targ_perm = 10
# N_min = 10
# N_max = 1000
# logrange(x1, x2, n) = (10^y for y in range(log10(x1), log10(x2), length=n))
# targ_count_range = collect(Int.(round.(logrange(N_min,N_max,N_targ_perms))))

# # Core loop
# solve_times_cvx = zeros(N_targ_perms, N_iters_per_targ_perm)
# solve_times_scp = zeros(N_targ_perms, N_iters_per_targ_perm)
# for j = 1:N_targ_perms
#     n = targ_count_range[j] 
#     for k = 1:N_iters_per_targ_perm

#         # DDTO-Cvx
#         try
#             params_cvx = DIntegrator2DoFParams(autogen_targs=true, autogen_targ_count=n)
#             solve_times_cvx[j,k] = @elapsed solve_cvx(params_cvx; simulate_solutions=false, process_the_solutions=false)
#         catch
#             solve_times_cvx[j,k] = Inf
#         end

#         # DDTO-SCP
#         try
#             params_scp = DIntegrator2DoFParams(autogen_targs=true, autogen_targ_count=n)
#             solve_times_scp[j,k] = @elapsed solve(params_scp; simulate_solutions=false, process_the_solutions=false)
#         catch
#             solve_times_scp[j,k] = Inf
#         end
#     end
# end

timing_comparison(targ_count_range, solve_times_cvx, solve_times_scp)

# # Plot a solution
# build_plots_single([sol],[sim],params_cvx)