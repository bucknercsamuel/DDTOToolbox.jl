using DDTOToolbox
using LinearAlgebra
using Random
using JLD2
include("scenarios.jl")
include("plots/plot_mc_compare.jl")

# Set the random seed
Random.seed!(123)

# Monte-carlo parameters
# n_sweep = [2,4,6,8,10,12,15,20]
# n_trials = 10
n_sweep = [2,3,4]
n_trials = 5
# n_sweep = [2]
# n_trials = 1

# Monte-Carlo main loop
fname = "quad3dof_cage/montecarlo_compare_results.jld2"
if isfile(fname)
    results_dict = load(fname)
    results_dict = results_dict["results_dict"] # due to JLD2 loading implementation
else
    # Create container to record convergence and solver time for all trials at each target count level for both DDTO-LEX and DDTO-SCP
    convergence_container_lex = Dict()
    solver_time_container_lex = Dict()
    deferral_time_container_lex = Dict()
    convergence_container_scp = Dict()
    solver_time_container_scp = Dict()
    deferral_time_container_scp = Dict()
    for n_targets in n_sweep
        convergence_container_lex[n_targets] = zeros(n_trials)
        solver_time_container_lex[n_targets] = zeros(n_trials)
        deferral_time_container_lex[n_targets] = zeros(n_trials, n_targets)
        convergence_container_scp[n_targets] = zeros(n_trials)
        solver_time_container_scp[n_targets] = zeros(n_trials)
        deferral_time_container_scp[n_targets] = zeros(n_trials, n_targets)
    end

    # Run Monte Carlo
    for n_targets in n_sweep
        for trial = 1:n_trials
            println("========== Trial $trial for $n_targets targets ==========")
            
            println("========== DDTO-LEX Trial ==========")
            params = scenario_obstacles_hard_random_targets(lex=true, n_targets=n_targets, min_distance_from_obstacle=0.01)
            _,_,_,_, converged, solver_time, deferral_times = solve_lex(params)
            convergence_container_lex[n_targets][trial] = converged
            solver_time_container_lex[n_targets][trial] = solver_time
            deferral_time_container_lex[n_targets][trial,:] = deferral_times

            println("========== DDTO-SCP Trial ==========")
            params = scenario_obstacles_hard_random_targets(lex=false, n_targets=n_targets, min_distance_from_obstacle=0.01)
            _,_,_,_, converged, solver_time, deferral_times = solve(params)
            convergence_container_scp[n_targets][trial] = converged
            solver_time_container_scp[n_targets][trial] = solver_time
            deferral_time_container_scp[n_targets][trial,:] = deferral_times
        end
    end

    # Save results to JLD2 file
    results_dict = Dict()
    results_dict["convergence_container_lex"] = convergence_container_lex
    results_dict["solver_time_container_lex"] = solver_time_container_lex
    results_dict["deferral_time_container_lex"] = deferral_time_container_lex
    results_dict["convergence_container_scp"] = convergence_container_scp
    results_dict["solver_time_container_scp"] = solver_time_container_scp
    results_dict["deferral_time_container_scp"] = deferral_time_container_scp
    jldsave(fname; results_dict)
end

# Plot results
screens = []
interactive = false
with_theme(theme2d; fontsize=fontsize) do
    push!(screens, plot_mc_compare(results_dict; interactive=interactive))
end
hold_interactive(screens)
;