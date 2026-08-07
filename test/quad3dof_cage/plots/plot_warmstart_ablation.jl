using CairoMakie
using Colors
using Statistics
include("../../utils/plot_utils.jl")
include("plot_defaults.jl")

"""
    plot_warmstart_ablation(results_dict; interactive=true)

Bar + scatter summary of warmstart ablation metrics:
convergence rate, solver time, and sum of deferral times.
"""
function plot_warmstart_ablation(results_dict; interactive=true)
    axis_defaults = Dict(
        :topspinevisible => true,
        :rightspinevisible => true,
        :xgridvisible => false,
        :ygridvisible => false,
    )
    ax_label_size = 15

    methods = results_dict["warmstart_methods"]
    n_methods = length(methods)
    n_trials = results_dict["n_trials"]
    convergence_container = results_dict["convergence_container"]
    solver_time_container = results_dict["solver_time_container"]
    deferral_time_container = results_dict["deferral_time_container"]

    # Aggregate per method
    conv_pct = Float64[]
    time_mean = Float64[]
    time_std = Float64[]
    defer_mean = Float64[]
    defer_std = Float64[]
    time_trials = Vector{Float64}[]
    defer_trials = Vector{Float64}[]

    for method in methods
        conv = convergence_container[method]
        push!(conv_pct, 100 * mean(conv))

        times = solver_time_container[method]
        push!(time_mean, mean(times))
        push!(time_std, std(times))
        push!(time_trials, times)

        defer_sums = [sum(deferral_time_container[method][t, :]) for t = 1:n_trials]
        # Report deferral objective only on converged trials when available
        converged_idx = findall(c -> c > 0, conv)
        defer_use = isempty(converged_idx) ? defer_sums : defer_sums[converged_idx]
        push!(defer_mean, mean(defer_use))
        push!(defer_std, length(defer_use) > 1 ? std(defer_use) : 0.0)
        push!(defer_trials, defer_sums)

        println("Method \"$method\": conv=$(conv_pct[end])%, " *
                "solver_time=$(time_mean[end]) ± $(time_std[end]) s, " *
                "Σ deferral=$(defer_mean[end]) ± $(defer_std[end]) s")
    end

    xs = 1:n_methods
    colors = [:steelblue, :darkorange, :seagreen][1:n_methods]
    f = Figure(size=(1400, 400))

    # Convergence
    ax = Axis(f[1, 1],
        xlabel="Warmstart Method",
        ylabel="% Converged";
        axis_defaults...)
    barplot!(ax, xs, conv_pct; color=colors, width=0.6)
    ax.xticks = (xs, methods)
    ax.yticks = 0:25:100
    ylims!(ax, 0, 105)

    # Solver time
    ax = Axis(f[1, 2],
        xlabel="Warmstart Method",
        ylabel="Solver Time [s]";
        axis_defaults...)
    barplot!(ax, xs, time_mean; color=colors, width=0.6)
    errorbars!(ax, xs, time_mean, time_std; color=:black, whiskerwidth=10)
    for (i, trials) in enumerate(time_trials)
        scatter!(ax, fill(xs[i], length(trials)), trials;
            color=:black, markersize=6, alpha=0.45)
    end
    ax.xticks = (xs, methods)

    # Deferral objective
    ax = Axis(f[1, 3],
        xlabel="Warmstart Method",
        ylabel=L"$\Sigma$ Deferral Times [s]";
        axis_defaults...)
    barplot!(ax, xs, defer_mean; color=colors, width=0.6)
    errorbars!(ax, xs, defer_mean, defer_std; color=:black, whiskerwidth=10)
    for (i, trials) in enumerate(defer_trials)
        # Mark only converged trials as filled; non-converged as open
        conv = convergence_container[methods[i]]
        for t = 1:n_trials
            mcolor = conv[t] > 0 ? :black : :gray
            marker = conv[t] > 0 ? :circle : :xcross
            scatter!(ax, [xs[i]], [trials[t]];
                color=mcolor, marker=marker, markersize=7, alpha=0.55)
        end
    end
    ax.xticks = (xs, methods)

    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    else
        out_dir = joinpath(@__DIR__, "..", "figures")
        mkpath(out_dir)
        out = joinpath(out_dir, "warmstart_ablation" * fig_ext)
        CairoMakie.save(out, f)
        println("Saved figure to $out")
        return nothing
    end
end
