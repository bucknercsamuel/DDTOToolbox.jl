using DDTOToolbox
using LinearAlgebra
using Random
using Statistics
using Printf
using PrettyTables
using JLD2
include("scenarios.jl")
include("plots/plot_mc_compare.jl")

# User-facing knobs (same fair-compare setup as compare_methods.jl)
const TF_MIN = 5.0       # [s] minimum physical time-of-flight
const TF_MAX = 20.0      # [s] maximum physical time-of-flight
const J_UB = 8.0         # normalized total-thrust upper bound (∫||T|| dt / ρ_max)
const N_NODES = 10
const SCP_ITERS = 100
const N_TARGS = [2, 3, 4, 5, 6, 7, 8]
const N_TRIALS = 10
const SEED = 123
const RESULTS_FILE = "quad3dof_cage/montecarlo_compare_results.jld2"
const SCHEMA = 2

# :scratch  resolve everything from scratch
# :resume   pick up missing (case, method, n_targets, trial) cells
# :plot     analysis / figures only
const RUN_MODE = :resume

const PLOT_INTERACTIVE = false

const CASES = [
    (name = "no_obstacles", obstacles = false, methods = ["qcvx", "lex", "scp"]),
    (name = "obstacles",    obstacles = true,  methods = ["lex", "scp"]),
]

function qcvx_deferral_times(params)
    Δt = params.a.Δt_cvx
    return [(max(Int(round(params.a.τ_targs[j])), 1) - 1) * Δt for j in 1:params.a.n_targs]
end

function run_method(method, params)
    if method == "qcvx"
        elapsed = @elapsed begin
            solve_cvx(params; simulate_solutions=false, process_the_solutions=false)
        end
        t_defer = qcvx_deferral_times(params)
        return (converged = true, elapsed = elapsed, t_defer = t_defer, n_iters = NaN)
    elseif method == "lex"
        _, _, converged, elapsed, t_defer = solve_lex(params; simulate_solutions=false, process_the_solutions=false)
        return (converged = converged, elapsed = elapsed, t_defer = collect(t_defer), n_iters = Float64(params.a.last_iters))
    elseif method == "scp"
        _, _, converged, elapsed, t_defer = solve(params; simulate_solutions=false, process_the_solutions=false)
        return (converged = converged, elapsed = elapsed, t_defer = collect(t_defer), n_iters = Float64(params.a.last_iters))
    else
        error("Unknown method $method")
    end
end

function empty_bucket()
    return Dict{String,Any}(
        "elapsed" => Float64[],
        "defer_obj" => Float64[],
        "converged" => Bool[],
        "n_iters" => Float64[],
    )
end

function trial_seed(obstacles::Bool, n_targets::Int, trial::Int)
    return SEED + 1000 * Int(obstacles) + 100 * n_targets + trial
end

function new_results()
    cases = Dict{String,Any}()
    for case in CASES
        cases[case.name] = Dict{String,Any}(m => Dict{Int,Any}() for m in case.methods)
    end
    return Dict{String,Any}(
        "schema" => SCHEMA,
        "meta" => Dict{String,Any}(
            "tf_min" => TF_MIN,
            "tf_max" => TF_MAX,
            "j_ub" => J_UB,
            "n_nodes" => N_NODES,
            "scp_iters" => SCP_ITERS,
            "seed" => SEED,
        ),
        "cases" => cases,
    )
end

function normalize_results!(results)
    cases = results["cases"]
    for (cname, methods) in cases
        for (mname, by_n) in methods
            normalized = Dict{Int,Any}()
            for (n, bucket) in by_n
                normalized[Int(n)] = bucket
            end
            methods[mname] = normalized
        end
    end
    return results
end

function load_results(path)
    if !isfile(path)
        return new_results()
    end
    raw = load(path)
    results = haskey(raw, "results") ? raw["results"] : raw
    if get(results, "schema", 0) != SCHEMA
        @warn "Existing results file has an incompatible schema; starting a new container" path=path
        return new_results()
    end
    return normalize_results!(results)
end

function save_results(path, results)
    mkpath(dirname(path))
    jldsave(path; results)
    return nothing
end

function ensure_bucket!(results, case_name, method, n_targets)
    by_n = results["cases"][case_name][method]
    if !haskey(by_n, n_targets)
        by_n[n_targets] = empty_bucket()
    end
    return by_n[n_targets]
end

function n_done(bucket)
    return length(bucket["elapsed"])
end

function append_trial!(bucket, rec, defer_obj)
    push!(bucket["elapsed"], rec.elapsed)
    push!(bucket["defer_obj"], defer_obj)
    push!(bucket["converged"], rec.converged)
    push!(bucket["n_iters"], rec.n_iters)
    return nothing
end

results = RUN_MODE == :scratch ? new_results() : load_results(RESULTS_FILE)

if RUN_MODE != :plot
    for case in CASES
        if !haskey(results["cases"], case.name)
            results["cases"][case.name] = Dict{String,Any}(m => Dict{Int,Any}() for m in case.methods)
        end
        for method in case.methods
            if !haskey(results["cases"][case.name], method)
                results["cases"][case.name][method] = Dict{Int,Any}()
            end
        end
        println("="^72)
        println("Case: $(case.name)   ToF∈[$TF_MIN, $TF_MAX] s   J_ub=$J_UB   mode=$RUN_MODE")
        println("="^72)
        for n_targets in N_TARGS
            for trial in 1:N_TRIALS
                pending = String[]
                for method in case.methods
                    bucket = ensure_bucket!(results, case.name, method, n_targets)
                    n_done(bucket) < trial && push!(pending, method)
                end
                isempty(pending) && continue

                Random.seed!(trial_seed(case.obstacles, n_targets, trial))
                params_base = scenario_fair_compare(
                    n_targets = n_targets,
                    obstacles = case.obstacles,
                    tf_min = TF_MIN,
                    tf_max = TF_MAX,
                    J_ub = J_UB,
                    N = N_NODES,
                    scp_iters = SCP_ITERS,
                )
                println("\n----- $(case.name)  n=$n_targets  trial $trial / $N_TRIALS -----")
                for method in pending
                    println("\n>>> $(case.name) / $method / n=$n_targets / trial $trial")
                    params = deepcopy(params_base)
                    rec = try
                        run_method(method, params)
                    catch err
                        @warn "Solve failed" case=case.name method=method n_targets=n_targets trial=trial exception=err
                        (converged = false, elapsed = NaN, t_defer = fill(NaN, n_targets), n_iters = NaN)
                    end
                    J_def = deferrability_objective(params_base.a.α_targs, rec.t_defer)
                    append_trial!(ensure_bucket!(results, case.name, method, n_targets), rec, J_def)
                    @printf("    elapsed = %.2f s   deferrability = %.2f s   converged = %s   iters = %s\n",
                        rec.elapsed, J_def, rec.converged,
                        isfinite(rec.n_iters) ? string(Int(round(rec.n_iters))) : "n/a")
                    save_results(RESULTS_FILE, results)
                end
            end
        end
    end
    save_results(RESULTS_FILE, results)
    println("\nSaved results → $RESULTS_FILE")
else
    println("Plot-only mode; loaded $RESULTS_FILE")
end

screens = []
with_theme(theme2d) do
    for case in CASES
        haskey(results["cases"], case.name) || continue
        screen = plot_mc_compare(results; case_name=case.name, interactive=PLOT_INTERACTIVE)
        screen !== nothing && push!(screens, screen)
    end
end
if PLOT_INTERACTIVE && !isempty(screens)
    hold_interactive(screens)
end
