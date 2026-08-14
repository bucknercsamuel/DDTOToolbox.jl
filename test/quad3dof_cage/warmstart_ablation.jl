using DDTOToolbox
using LinearAlgebra
using Random
using Statistics
using JLD2
include("scenarios.jl")

function print_warmstart_ablation_summary(results_dict)
    methods = results_dict["warmstart_methods"]
    n_trials = results_dict["n_trials"]
    convergence_container = results_dict["convergence_container"]
    solver_time_container = results_dict["solver_time_container"]
    deferral_time_container = results_dict["deferral_time_container"]
    iters_container = results_dict["iters_container"]

    println("\n========== Warmstart ablation summary ==========")
    for method in methods
        conv = convergence_container[method]
        converged_idx = findall(c -> c > 0, conv)
        times = solver_time_container[method]
        iters = iters_container[method]

        if isempty(converged_idx)
            deferral_str = "n/a (0 converged)"
            iters_str = "n/a (0 converged)"
        else
            defer_sums = [sum(deferral_time_container[method][t, :]) for t in converged_idx]
            iters_use = iters[converged_idx]
            n_conv = length(converged_idx)
            defer_std = n_conv > 1 ? std(defer_sums) : 0.0
            iters_std = n_conv > 1 ? std(iters_use) : 0.0
            deferral_str = "$(mean(defer_sums)) ± $defer_std s (n=$n_conv)"
            iters_str = "$(mean(iters_use)) ± $iters_std (n=$n_conv)"
        end

        time_std = n_trials > 1 ? std(times) : 0.0
        println("Method \"$method\": conv=$(100 * mean(conv))%, " *
                "solver_time=$(mean(times)) ± $time_std s, " *
                "iters=$iters_str, " *
                "Σ deferral=$deferral_str")
    end
end

# =============================================================================
# Warmstart ablation: DDTO-SCP only, fixed n=4, methods = linear / single / ddto
# Same scenario config as montecarlo_compare.jl (random targets, hard obstacles).
# =============================================================================

# Set true to skip solves and only plot from an existing results file in data/
analyze_existing = false

# Ablation parameters
n_targets = 4
n_trials = 50
warmstart_methods = ["linear", "single", "ddto"]
Random.seed!(123)

data_dir = joinpath(@__DIR__, "data")
mkpath(data_dir)
fname = joinpath(data_dir, "warmstart_ablation_n$(n_targets)_t$(n_trials).jld2")

if analyze_existing
    if !isfile(fname)
        error("analyze_existing=true but results file not found: $fname")
    end
    results_dict = load(fname)["results_dict"]
    println("Loaded existing results from $fname")
else
    # Containers keyed by warmstart method
    convergence_container = Dict{String,Vector{Float64}}()
    solver_time_container = Dict{String,Vector{Float64}}()
    deferral_time_container = Dict{String,Matrix{Float64}}()
    iters_container = Dict{String,Vector{Float64}}()
    for method in warmstart_methods
        convergence_container[method] = zeros(n_trials)
        solver_time_container[method] = zeros(n_trials)
        deferral_time_container[method] = zeros(n_trials, n_targets)
        iters_container[method] = zeros(n_trials)
    end

    for trial = 1:n_trials
        println("========== Trial $trial / $n_trials (n=$n_targets) ==========")

        # Draw one random target set per trial so methods are paired fairly
        # params_base = scenario_obstacles_hard_random_targets(
        #     lex=false,
        #     n_targets=n_targets,
        #     min_distance_from_obstacle=0.01,
        # )
        params_base = scenario_with_random_targets(scenario_obstacles_hard;
            lex=false,
            n_targets=n_targets,
            min_distance_from_obstacle=0.01,
        )

        for method in warmstart_methods
            println("========== DDTO-SCP warmstart=\"$method\" ==========")
            params = deepcopy(params_base)
            params.a.warmstart_method = method
            _, _, _, _, converged, solver_time, deferral_times, n_iters = solve(params)
            convergence_container[method][trial] = Float64(converged)
            solver_time_container[method][trial] = solver_time
            deferral_time_container[method][trial, :] = deferral_times
            iters_container[method][trial] = Float64(n_iters)
        end
    end

    results_dict = Dict(
        "n_targets" => n_targets,
        "n_trials" => n_trials,
        "warmstart_methods" => warmstart_methods,
        "convergence_container" => convergence_container,
        "solver_time_container" => solver_time_container,
        "deferral_time_container" => deferral_time_container,
        "iters_container" => iters_container,
    )
    jldsave(fname; results_dict)
    println("Saved results to $fname")
end

# Analysis
print_warmstart_ablation_summary(results_dict)
;
