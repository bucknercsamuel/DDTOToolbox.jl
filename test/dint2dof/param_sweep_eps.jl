using DDTOSCP
using LinearAlgebra
using Statistics
using Printf
include("plots.jl")

params = DIntegrator2DoFParams()

logrange(x1, x2, n) = (10^y for y in range(log10(x1), log10(x2), length=n))

N = 20
param_range = collect(logrange(1e-3,1,N))
defer_times_ddtocvx = Matrix(undef, params.a.n_targs-1, N)
defer_times_ddtoscp = Matrix(undef, params.a.n_targs-1, N)
τ_lu(j) = params.a.τ_targs[findfirst(i->i==j, params.a.λ_targs)] # obtain the deferrability index of the j-th target (solution)
for k = 1:N
    # Apply epsilon change
    params.a.ϵ_targs = fill(param_range[k], params.a.n_targs)

    # Solve DDTO-Cvx, record deferral times
    _,_,sol,sim = solve_cvx(params)
    for j∈1:params.a.n_targs
        if j != params.a.λ_targs[end]
            defer_times_ddtocvx[j,k] = sol.targs[j].t[τ_lu(j)]
        end
    end

    # Solve DDTO-SCP, record deferral times
    _,_,sol,sim = solve(params)
    for j∈1:params.a.n_targs
        if j != params.a.λ_targs[end]
            defer_times_ddtoscp[j,k] = sol.targs[j].t[τ_lu(j)]
        end
    end
end

mean_cvx = mean(defer_times_ddtocvx, dims=2)
mean_scp = mean(defer_times_ddtoscp, dims=2)
std_cvx = std(defer_times_ddtocvx, dims=2)
std_scp = std(defer_times_ddtoscp, dims=2)

println("=== Results ===")
println("DDTO-Cvx:")
@printf("   > Target 1 | mean: %.2f s, std: %.2f s\n", mean_cvx[1], std_cvx[1])
@printf("   > Target 2 | mean: %.2f s, std: %.2f s\n", mean_cvx[2], std_cvx[2])
@printf("   > Target 3 | mean: %.2f s, std: %.2f s\n", mean_cvx[3], std_cvx[3])
println("DDTO-SCP:")
@printf("   > Target 1 | mean: %.2f s, std: %.2f s\n", mean_scp[1], std_scp[1])
@printf("   > Target 2 | mean: %.2f s, std: %.2f s\n", mean_scp[2], std_scp[2])
@printf("   > Target 3 | mean: %.2f s, std: %.2f s\n", mean_scp[3], std_scp[3])



# build_plots([], [], sols, sims, params)