using DDTOToolbox
using LinearAlgebra
using Random
using JLD2
include("scenarios.jl")
include("plots/plot_warmstart_ablation.jl")

# =============================================================================
# Warmstart ablation: DDTO-SCP only, fixed n=4, methods = linear / single / ddto
# Same scenario config as montecarlo_compare.jl (random targets, hard obstacles).
# =============================================================================

# Set true to skip solves and only plot from an existing results file in data/
analyze_existing = false

# Ablation parameters
n_targets = 4
n_trials = 10
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
    for method in warmstart_methods
        convergence_container[method] = zeros(n_trials)
        solver_time_container[method] = zeros(n_trials)
        deferral_time_container[method] = zeros(n_trials, n_targets)
    end

    for trial = 1:n_trials
        println("========== Trial $trial / $n_trials (n=$n_targets) ==========")

        # Draw one random target set per trial so methods are paired fairly
        params_base = scenario_obstacles_hard_random_targets(
            lex=false,
            n_targets=n_targets,
            min_distance_from_obstacle=0.01,
        )

        for method in warmstart_methods
            println("========== DDTO-SCP warmstart=\"$method\" ==========")
            params = deepcopy(params_base)
            params.a.warmstart_method = method
            _, _, _, _, converged, solver_time, deferral_times = solve(params)
            convergence_container[method][trial] = Float64(converged)
            solver_time_container[method][trial] = solver_time
            deferral_time_container[method][trial, :] = deferral_times
        end
    end

    results_dict = Dict(
        "n_targets" => n_targets,
        "n_trials" => n_trials,
        "warmstart_methods" => warmstart_methods,
        "convergence_container" => convergence_container,
        "solver_time_container" => solver_time_container,
        "deferral_time_container" => deferral_time_container,
    )
    jldsave(fname; results_dict)
    println("Saved results to $fname")
end

# Analysis / plots
screens = []
interactive = false
with_theme(theme_latexfonts()) do
    push!(screens, plot_warmstart_ablation(results_dict; interactive=interactive))
end
if interactive
    hold_interactive(screens)
end
;
