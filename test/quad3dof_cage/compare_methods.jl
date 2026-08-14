using DDTOToolbox
using LinearAlgebra
using Random
using Statistics
using Printf
using PrettyTables
include("scenarios.jl")
include("plots/plot_trajs.jl")

# User-facing knobs (free final time + fixed thrust bound, no ε-suboptimality)
const TF_MIN = 5.0       # [s] minimum physical time-of-flight
const TF_MAX = 20.0      # [s] maximum physical time-of-flight
const J_UB = 8.0         # normalized total-thrust upper bound (∫||T|| dt / ρ_max)
const N_NODES = 10
const SCP_ITERS = 100
const N_TARGETS = 4
const N_TRIALS = 1
const SEED = 123
const PLOT_FIRST_TRIAL = true
const PLOT_INTERACTIVE = true

const CASES = [
    (name = "no_obstacles", obstacles = false, methods = ["qcvx", "lex", "scp"]),
    # (name = "obstacles",    obstacles = true,  methods = ["lex", "scp"]),
]

function qcvx_deferral_times(params)
    Δt = params.a.Δt_cvx
    return [(max(Int(round(params.a.τ_targs[j])), 1) - 1) * Δt for j in 1:params.a.n_targs]
end

function run_method(method, params; keep_plots=false)
    sim = keep_plots
    proc = keep_plots
    if method == "qcvx"
        elapsed = @elapsed begin
            out = solve_cvx(params; simulate_solutions=sim, process_the_solutions=proc)
        end
        t_defer = qcvx_deferral_times(params)
        sol = sim ? out[3] : nothing
        simu = sim ? out[4] : nothing
        return (converged = true, elapsed = elapsed, t_defer = t_defer, sol = sol, sim = simu)
    elseif method == "lex"
        out = solve_lex(params; simulate_solutions=sim, process_the_solutions=proc)
        if sim
            _, _, sol, simu, converged, elapsed, t_defer = out
        else
            _, _, converged, elapsed, t_defer = out
            sol = nothing
            simu = nothing
        end
        return (converged = converged, elapsed = elapsed, t_defer = collect(t_defer), sol = sol, sim = simu)
    elseif method == "scp"
        out = solve(params; simulate_solutions=sim, process_the_solutions=proc)
        if sim
            _, _, sol, simu, converged, elapsed, t_defer = out
        else
            _, _, converged, elapsed, t_defer = out
            sol = nothing
            simu = nothing
        end
        return (converged = converged, elapsed = elapsed, t_defer = collect(t_defer), sol = sol, sim = simu)
    else
        error("Unknown method $method")
    end
end

function fmt_msd(xs)
    data = filter(isfinite, xs)
    isempty(data) && return "n/a"
    @sprintf("%.2f ± %.2f", mean(data), std(data))
end

results = Dict()
for case in CASES
    results[case.name] = Dict(m => (
        elapsed = fill(NaN, N_TRIALS),
        defer_obj = fill(NaN, N_TRIALS),
        converged = fill(false, N_TRIALS),
    ) for m in case.methods)
end

plot_screens = []
for case in CASES
    println("="^72)
    println("Case: $(case.name)   n_targets=$N_TARGETS   ToF∈[$TF_MIN, $TF_MAX] s   J_ub=$J_UB")
    println("="^72)
    for trial in 1:N_TRIALS
        Random.seed!(SEED + 1000 * Int(case.obstacles) + trial)
        params_base = scenario_fair_compare(
            n_targets = N_TARGETS,
            obstacles = case.obstacles,
            tf_min = TF_MIN,
            tf_max = TF_MAX,
            J_ub = J_UB,
            N = N_NODES,
            scp_iters = SCP_ITERS,
        )
        println("\n----- Trial $trial / $N_TRIALS -----")
        for method in case.methods
            println("\n>>> $(case.name) / $method / trial $trial")
            params = deepcopy(params_base)
            keep_plots = PLOT_FIRST_TRIAL && trial == 1
            rec = try
                run_method(method, params; keep_plots=keep_plots)
            catch err
                @warn "Solve failed" case=case.name method=method trial=trial exception=err
                (converged = false, elapsed = NaN, t_defer = fill(NaN, N_TARGETS), sol = nothing, sim = nothing)
            end
            J_def = deferrability_objective(params_base.a.α_targs, rec.t_defer)
            results[case.name][method].elapsed[trial] = rec.elapsed
            results[case.name][method].defer_obj[trial] = J_def
            results[case.name][method].converged[trial] = rec.converged
            @printf("    elapsed = %.2f s   deferrability = %.2f s   converged = %s\n",
                rec.elapsed, J_def, rec.converged)
            if keep_plots && rec.sol !== nothing && rec.sim !== nothing
                save_name = "compare_$(case.name)_$(method)_trial1"
                println("    plotting first-trial DDTO solution → $save_name")
                with_theme(theme2d) do
                    screen = plot_trajs([rec.sol], [rec.sim], params;
                        interactive = PLOT_INTERACTIVE,
                        ddto = true,
                        obstacles = case.obstacles,
                        save_name = save_name,
                    )
                    PLOT_INTERACTIVE && push!(plot_screens, screen)
                end
            end
        end
    end
end

println("\n" * "="^72)
println("Summary  (mean ± std over $N_TRIALS trials; ToF∈[$TF_MIN, $TF_MAX] s, J_ub=$J_UB, N=$N_TARGETS)")
println("="^72)
for case in CASES
    n = length(case.methods)
    table = Matrix{String}(undef, n, 4)
    for (i, method) in enumerate(case.methods)
        r = results[case.name][method]
        n_ok = count(r.converged)
        table[i, 1] = method
        table[i, 2] = fmt_msd(r.elapsed)
        table[i, 3] = fmt_msd(r.defer_obj)
        table[i, 4] = @sprintf("%d/%d", n_ok, N_TRIALS)
    end
    println("\n$(case.name):")
    pretty_table(
        table;
        column_labels = ["Method", "Solve time [s]", "Deferrability obj [s]", "Converged"],
        alignment = [:l, :c, :c, :c],
    )
end

if PLOT_FIRST_TRIAL && PLOT_INTERACTIVE && !isempty(plot_screens)
    hold_interactive(plot_screens)
end
