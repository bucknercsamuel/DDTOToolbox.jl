using DDTOToolbox
using LinearAlgebra
using Printf
include("scenarios.jl")
include("plots/plot_ddto_compare.jl")

# =============================================================================
# Nominal (non-random) comparison: Lex-DDTO vs Graph-DDTO vs DDTO-QCvx
# =============================================================================

# :hard -> scenario_obstacles_hard; :none -> scenario_no_obstacles
obstacle_mode = :none

interactive = false

function load_nominal_params(; lex::Bool)
    if obstacle_mode == :hard
        return scenario_obstacles_hard(lex)
    elseif obstacle_mode == :none
        return scenario_no_obstacles(lex)
    else
        error("Unknown obstacle_mode=$obstacle_mode (use :hard or :none)")
    end
end

function unpack_solve(out)
    # Compatible with 7- or 8-element returns (n_iters optional)
    sol_scp, sim_scp, sol_ddto, sim_ddto = out[1], out[2], out[3], out[4]
    converged, solver_time, deferral = out[5], out[6], out[7]
    n_iters = length(out) >= 8 ? out[8] : NaN
    return sol_scp, sim_scp, sol_ddto, sim_ddto, converged, solver_time, deferral, n_iters
end

function print_method_summary(name, params, sol_ddto, converged, solver_time, deferral, n_iters)
    println("\n========== $name ==========")
    println("  converged    = $converged")
    println("  solver_time  = $(round(solver_time; digits=3)) s")
    println("  n_iters      = $n_iters")
    println("  warmstart    = $(params.a.warmstart_method)")
    println("  N            = $(params.a.N)")
    println("  τ_targs      = $(params.a.τ_targs)")
    println("  Σ deferral   = $(round(sum(deferral); digits=4)) s")
    for j = 1:params.a.n_targs
        @printf("  Target %d: defer=%.4f s, cost=%.4f\n",
            j, deferral[j], sol_ddto.targs[j].cost)
    end
end

# ----- Solve all three on the same nominal geometry -----
params_graph = load_nominal_params(lex=false)
params_lex   = load_nominal_params(lex=true)
params_qcvx  = load_nominal_params(lex=false)

println("Obstacle mode: $obstacle_mode")
println("Targets (Graph/QCvx N=$(params_graph.a.N), Lex N=$(params_lex.a.N)):")
for j = 1:params_graph.a.n_targs
    println("  j=$j  rf=$(params_graph.a.zf_targs[1:3,j])  λ-order pos=$(findfirst(==(j), params_graph.a.λ_targs))")
end

println("\n>>> Solving Lex-DDTO...")
lex_out = solve_lex(params_lex)
_, _, lex_sol, lex_sim, lex_conv, lex_time, lex_defer, lex_iters = unpack_solve(lex_out)
print_method_summary("Lex-DDTO", params_lex, lex_sol, lex_conv, lex_time, lex_defer, lex_iters)

println("\n>>> Solving Graph-DDTO (SCP)...")
scp_out = solve(params_graph)
_, _, scp_sol, scp_sim, scp_conv, scp_time, scp_defer, scp_iters = unpack_solve(scp_out)
print_method_summary("Graph-DDTO", params_graph, scp_sol, scp_conv, scp_time, scp_defer, scp_iters)

println("\n>>> Solving DDTO-QCvx...")
qcvx_time = @elapsed begin
    _, _, qcvx_sol, qcvx_sim = solve_cvx(params_qcvx; simulate_solutions=true, process_the_solutions=true)
end
n = params_qcvx.a.n_targs
qcvx_defer = zeros(n)
for j = 1:n
    τ = max(Int(params_qcvx.a.τ_targs[j]), 1)
    qcvx_defer[j] = qcvx_sol.targs[j].t[τ]
end
qcvx_conv = all(isfinite(qcvx_sol.targs[j].cost) for j = 1:n)
print_method_summary("DDTO-QCvx", params_qcvx, qcvx_sol, qcvx_conv, qcvx_time, qcvx_defer, NaN)

# ----- Plot -----
screens = []
titles = [
    "Lex-DDTO\nΣ defer=$(round(sum(lex_defer); digits=2)) s, conv=$(lex_conv)",
    "Graph-DDTO\nΣ defer=$(round(sum(scp_defer); digits=2)) s, conv=$(scp_conv)",
    "DDTO-QCvx\nΣ defer=$(round(sum(qcvx_defer); digits=2)) s, conv=$(qcvx_conv)",
]
with_theme(theme_latexfonts()) do
    push!(screens, plot_ddto_compare(
        [lex_sol, scp_sol, qcvx_sol],
        [lex_sim, scp_sim, qcvx_sim],
        [params_lex, params_graph, params_qcvx],
        titles;
        interactive=interactive,
        fig_name="ddto_compare_$(obstacle_mode)",
    ))
end
if interactive
    hold_interactive(screens)
end
;
