using DDTOSCP
using LinearAlgebra
using Statistics
using Random
using Printf
include("plots.jl")

MersenneTwister(0)

# SCP params
w_obj_sing = 1e-1
w_obj_ddto = 1e-2
w_ctrl = 1e4
w_trust = 1e2
scp_iters = 50

# Test params
N_targ_perms = 5
N_iters_per_targ_perm = 5
N_min = 10
N_max = 200
logrange(x1, x2, n) = (10^y for y in range(log10(x1), log10(x2), length=n))
targ_count_range = collect(Int.(round.(logrange(N_min,N_max,N_targ_perms))))

# Core loop
solve_times_cvx = zeros(N_targ_perms, N_iters_per_targ_perm)
solve_times_scp = zeros(N_targ_perms, N_iters_per_targ_perm)
for j = 1:N_targ_perms
    n = targ_count_range[j] 
    for k = 1:N_iters_per_targ_perm
        # Generate DDTO-Cvx params
        params_cvx = DIntegrator2DoFParams(autogen_targs=true, autogen_targ_count=n, scp=false)

        # Time DDTO-Cvx solve
        try
            solve_times_cvx[j,k] = @elapsed solve_cvx(params_cvx; simulate_solutions=false, process_the_solutions=false)
        catch
            solve_times_cvx[j,k] = Inf
        end

        # Generate DDTO-SCP params
        params_scp = DIntegrator2DoFParams(autogen_targs=true, autogen_targ_count=n)
        params_scp.w_obj_sing = w_obj_sing
        params_scp.w_obj_ddto = w_obj_ddto
        params_scp.w_ctrl = w_ctrl
        params_scp.w_trust = w_trust
        params_scp.scp_iters = scp_iters

        # Time DDTO-SCP solve
        try
            solve_times_scp[j,k] = @elapsed solve(params_scp; simulate_solutions=false, process_the_solutions=false)
        catch
            solve_times_scp[j,k] = Inf
        end
    end
end

timing_comparison(targ_count_range, solve_times_cvx, solve_times_scp)

# # Plot a solution
# build_plots_single([sol],[sim],params_cvx)