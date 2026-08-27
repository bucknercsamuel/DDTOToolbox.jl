using DDTOToolbox
using LinearAlgebra
using Printf
using Random
include("scenarios.jl")
include("plots/plot_trajs.jl")

function ctcs_branch_totals(sim)
    return [sim.targs[j].x[end, end] - sim.targs[j].x[end, 1] for j in eachindex(sim.targs)]
end

function print_ctcs_verification(label, sim, ϵ)
    ξ = ctcs_branch_totals(sim)
    println("\nCTCS nonlinear-sim path violation — ", label)
    @printf("  tolerance ϵ_ctcs = %.1e\n", ϵ)
    for (j, ξj) in enumerate(ξ)
        flag = ξj <= ϵ ? "≤ ϵ" : "> ϵ"
        @printf("  target %d: %9.3e  (%s)\n", j, ξj, flag)
    end
    ξmax = maximum(ξ)
    @printf("  max:      %9.3e  (%s)\n", ξmax, ξmax <= ϵ ? "PASS" : "FAIL")
end

# Set the random seed
Random.seed!(123)

# Choose scenario
lex = false
params = scenario_obstacles_hard(lex=lex)
# params = scenario_no_obstacles(lex=lex)

# Solve
if lex
    scp_sol, scp_sim, ddtoscp_sol, ddtoscp_sim = solve_lex(params)
else
    scp_sol, scp_sim, ddtoscp_sol, ddtoscp_sim = solve(params)
end

print_ctcs_verification("decoupled SCP", scp_sim, params.a.ϵ_ctcs)
print_ctcs_verification(lex ? "DDTO-lex" : "Graph-DDTO", ddtoscp_sim, params.a.ϵ_ctcs)

# Plot results
screens = []
interactive = true
fontsize = 12
with_theme(theme2d; fontsize=fontsize) do
    push!(screens, plot_trajs([scp_sol],     [scp_sim],     params; interactive=interactive, ddto=false))
    push!(screens, plot_trajs([ddtoscp_sol], [ddtoscp_sim], params; interactive=interactive))
end
if interactive
    hold_interactive(screens)
end
;