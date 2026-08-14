using DDTOToolbox
using LinearAlgebra
using Random
using Statistics
using JLD2
include("scenarios.jl")
include("plots/plot_mc_compare.jl")

# =============================================================================
# Monte Carlo: Lex-DDTO vs Graph-DDTO (+ DDTO-QCvx when obstacle_mode == :none)
# Saves after every trial so results survive a killed terminal.
# =============================================================================

Random.seed!(123)

# ----- User knobs -----
# :hard -> hard-obstacle arena (Lex + Graph only)
# :none -> free space (Lex + Graph + QCvx)
obstacle_mode = :hard
# obstacle_mode = :none

n_sweep = [2, 3, 4]
n_trials = 10

# Fixed wall-clock ToF [s] shared by Lex / Graph / QCvx (nothing -> scenario default)
# Defaults: :none -> 10 s, :hard -> 12 s
ToF = nothing

# Set true to skip solves and only plot/print from an existing results file
analyze_existing = false

# ----- Paths / methods -----
include_qcvx = obstacle_mode == :none
methods = include_qcvx ? ["lex", "scp", "qcvx"] : ["lex", "scp"]

data_dir = joinpath(@__DIR__, "data")
mkpath(data_dir)
fname = joinpath(data_dir, "montecarlo_compare_$(obstacle_mode)_n$(join(n_sweep, "-"))_t$(n_trials).jld2")

function make_scenario(n_targets::Int; lex::Bool)
    kwargs = isnothing(ToF) ? NamedTuple() : (; ToF=ToF)
    if obstacle_mode == :hard
        return scenario_with_random_targets(scenario_obstacles_hard;
            n_targets=n_targets, min_distance_from_obstacle=0.01, lex=lex, kwargs...)
    elseif obstacle_mode == :none
        return scenario_with_random_targets(scenario_no_obstacles;
            n_targets=n_targets, min_distance_from_obstacle=0.01, lex=lex, kwargs...)
    else
        error("Unknown obstacle_mode=$obstacle_mode (use :hard or :none)")
    end
end

"""Run DDTO-QCvx and return (converged, solver_time, deferral_times, n_iters)."""
function run_qcvx(params)
    # QCvx has no PTR outer loop; n_iters left as NaN for the iterations plot.
    solver_time = @elapsed begin
        _, _, ddto_sol, _ = solve_cvx(params; simulate_solutions=true, process_the_solutions=true)
    end
    n = params.a.n_targs
    deferral_times = zeros(n)
    for j = 1:n
        τ = Int(params.a.τ_targs[j])
        deferral_times[j] = ddto_sol.targs[j].t[τ]
    end
    converged = all(isfinite(ddto_sol.targs[j].cost) for j = 1:n)
    return converged, solver_time, deferral_times, NaN
end

function empty_results(n_sweep, n_trials, methods)
    results = Dict{String,Any}(
        "obstacle_mode" => String(obstacle_mode),
        "n_sweep" => collect(n_sweep),
        "n_trials" => n_trials,
        "methods" => collect(methods),
        "convergence" => Dict{String,Dict{Int,Vector{Float64}}}(),
        "solver_time" => Dict{String,Dict{Int,Vector{Float64}}}(),
        "deferral_time" => Dict{String,Dict{Int,Matrix{Float64}}}(),
        "iters" => Dict{String,Dict{Int,Vector{Float64}}}(),
        "completed" => Dict{Int,Vector{Bool}}(),
    )
    for m in methods
        results["convergence"][m] = Dict{Int,Vector{Float64}}()
        results["solver_time"][m] = Dict{Int,Vector{Float64}}()
        results["deferral_time"][m] = Dict{Int,Matrix{Float64}}()
        results["iters"][m] = Dict{Int,Vector{Float64}}()
        for n in n_sweep
            results["convergence"][m][n] = fill(NaN, n_trials)
            results["solver_time"][m][n] = fill(NaN, n_trials)
            results["deferral_time"][m][n] = fill(NaN, n_trials, n)
            results["iters"][m][n] = fill(NaN, n_trials)
        end
    end
    for n in n_sweep
        results["completed"][n] = falses(n_trials)
    end
    return results
end

function save_results!(results_dict, fname)
    jldsave(fname; results_dict)
    println("Saved checkpoint -> $fname")
end

function print_mc_summary(results_dict)
    methods = results_dict["methods"]
    n_sweep = results_dict["n_sweep"]
    n_trials = results_dict["n_trials"]
    println("\n========== Monte Carlo summary (mode=$(results_dict["obstacle_mode"])) ==========")
    for n in n_sweep
        println("--- n_targets = $n ---")
        for m in methods
            conv = results_dict["convergence"][m][n]
            times = results_dict["solver_time"][m][n]
            iters = results_dict["iters"][m][n]
            defer = results_dict["deferral_time"][m][n]
            done = findall(isfinite, conv)
            isempty(done) && (println("  $m: no completed trials"); continue)
            conv_pct = 100 * mean(conv[done])
            conv_idx = findall(i -> i in done && conv[i] > 0, 1:n_trials)
            time_str = "$(mean(times[done])) ± $(std(times[done])) s"
            if isempty(conv_idx)
                defer_str = "n/a (0 converged)"
                iters_str = "n/a (0 converged)"
            else
                defer_sums = [sum(defer[t, :]) for t in conv_idx]
                iters_use = filter(isfinite, iters[conv_idx])
                defer_str = "$(mean(defer_sums)) ± $(std(defer_sums)) s (n=$(length(conv_idx)))"
                iters_str = isempty(iters_use) ? "n/a" :
                    "$(mean(iters_use)) ± $(std(iters_use)) (n=$(length(iters_use)))"
            end
            println("  $m: conv=$(conv_pct)%, time=$time_str, iters=$iters_str, Σ deferral=$defer_str")
        end
    end
end

# ----- Load or initialize -----
if analyze_existing || isfile(fname)
    if !isfile(fname)
        error("analyze_existing=true but results file not found: $fname")
    end
    results_dict = load(fname)["results_dict"]
    println("Loaded results from $fname")
    # Back-compat: older files used per-method top-level keys
    if !haskey(results_dict, "methods")
        error("Results file is from an older format; delete or rename it and re-run.")
    end
else
    results_dict = empty_results(n_sweep, n_trials, methods)
end

# ----- Solve loop (skips completed trials; checkpoints after each trial) -----
if !analyze_existing
    for n_targets in n_sweep
        for trial = 1:n_trials
            if results_dict["completed"][n_targets][trial]
                println("Skipping completed trial $trial / $n_trials (n=$n_targets)")
                continue
            end
            println("========== Trial $trial / $n_trials for $n_targets targets (mode=$obstacle_mode) ==========")

            # Shared random target set for all methods in this trial
            params_base = make_scenario(n_targets; lex=false)
            params_lex = make_scenario(n_targets; lex=true)
            params_lex.a.zf_targs = copy(params_base.a.zf_targs)
            params_lex.a.uf_targs = copy(params_base.a.uf_targs)
            params_lex.a.λ_targs = copy(params_base.a.λ_targs)
            params_lex.a.J_targs = copy(params_base.a.J_targs)
            params_lex.a.α_targs = copy(params_base.a.α_targs)
            params_lex.a.ϵ_targs = fill(params_lex.a.ϵ_targs[1], n_targets)

            # Lex-DDTO
            println("========== Lex-DDTO ==========")
            lex_out = solve_lex(params_lex)
            converged, solver_time, deferral_times = lex_out[5], lex_out[6], lex_out[7]
            n_iters = length(lex_out) >= 8 ? lex_out[8] : NaN
            results_dict["convergence"]["lex"][n_targets][trial] = Float64(converged)
            results_dict["solver_time"]["lex"][n_targets][trial] = solver_time
            results_dict["deferral_time"]["lex"][n_targets][trial, :] = deferral_times
            results_dict["iters"]["lex"][n_targets][trial] = Float64(n_iters)

            # Graph-DDTO (SCP)
            println("========== Graph-DDTO (SCP) ==========")
            params = deepcopy(params_base)
            scp_out = solve(params)
            converged, solver_time, deferral_times = scp_out[5], scp_out[6], scp_out[7]
            n_iters = length(scp_out) >= 8 ? scp_out[8] : NaN
            results_dict["convergence"]["scp"][n_targets][trial] = Float64(converged)
            results_dict["solver_time"]["scp"][n_targets][trial] = solver_time
            results_dict["deferral_time"]["scp"][n_targets][trial, :] = deferral_times
            results_dict["iters"]["scp"][n_targets][trial] = Float64(n_iters)

            # DDTO-QCvx (no-obstacle mode only)
            if include_qcvx
                println("========== DDTO-QCvx ==========")
                params = deepcopy(params_base)
                converged, solver_time, deferral_times, n_iters = run_qcvx(params)
                results_dict["convergence"]["qcvx"][n_targets][trial] = Float64(converged)
                results_dict["solver_time"]["qcvx"][n_targets][trial] = solver_time
                results_dict["deferral_time"]["qcvx"][n_targets][trial, :] = deferral_times
                results_dict["iters"]["qcvx"][n_targets][trial] = Float64(n_iters)
            end

            results_dict["completed"][n_targets][trial] = true
            save_results!(results_dict, fname)
        end
    end
end

# ----- Analysis -----
print_mc_summary(results_dict)

screens = []
interactive = false
with_theme(theme_latexfonts()) do
    push!(screens, plot_mc_compare(results_dict; interactive=interactive))
end
if interactive
    hold_interactive(screens)
end
;
