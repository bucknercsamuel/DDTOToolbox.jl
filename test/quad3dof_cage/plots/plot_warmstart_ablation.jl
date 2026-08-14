using CairoMakie
using Colors
using Statistics
include("../../utils/plot_utils.jl")
include("plot_defaults.jl")

"""
    plot_warmstart_ablation(results_dict; interactive=true)

Bar + scatter summary of warmstart ablation metrics:
convergence rate, solver time, iterations to convergence, and sum of deferral times.
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
    iters_container = results_dict["iters_container"]

    # Aggregate per method
    conv_pct = Float64[]
    time_mean = Float64[]
    time_std = Float64[]
    defer_mean = Float64[]
    defer_std = Float64[]
    iters_mean = Float64[]
    iters_std = Float64[]
    time_trials = Vector{Float64}[]
    defer_trials = Vector{Float64}[]
    iters_trials = Vector{Float64}[]

    for method in methods
        conv = convergence_container[method]
        push!(conv_pct, 100 * mean(conv))

        times = solver_time_container[method]
        push!(time_mean, mean(times))
        push!(time_std, std(times))
        push!(time_trials, times)

        iters = iters_container[method]
        converged_idx = findall(c -> c > 0, conv)
        iters_use = isempty(converged_idx) ? iters : iters[converged_idx]
        push!(iters_mean, mean(iters_use))
        push!(iters_std, length(iters_use) > 1 ? std(iters_use) : 0.0)
        push!(iters_trials, iters)

        # Σ deferral is only meaningful on converged trials
        defer_sums = [sum(deferral_time_container[method][t, :]) for t = 1:n_trials]
        if isempty(converged_idx)
            push!(defer_mean, NaN)
            push!(defer_std, NaN)
            push!(defer_trials, Float64[])
            deferral_str = "n/a (0 converged)"
        else
            defer_use = defer_sums[converged_idx]
            push!(defer_mean, mean(defer_use))
            push!(defer_std, length(defer_use) > 1 ? std(defer_use) : 0.0)
            push!(defer_trials, defer_use)
            deferral_str = "$(defer_mean[end]) ± $(defer_std[end]) s (n=$(length(defer_use)))"
        end

        println("Method \"$method\": conv=$(conv_pct[end])%, " *
                "solver_time=$(time_mean[end]) ± $(time_std[end]) s, " *
                "iters=$(iters_mean[end]) ± $(iters_std[end]), " *
                "Σ deferral=$deferral_str")
    end

    xs = 1:n_methods
    colors = [:steelblue, :darkorange, :seagreen][1:n_methods]
    f = Figure(size=(1800, 400))

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

    # Iterations to convergence
    ax = Axis(f[1, 3],
        xlabel="Warmstart Method",
        ylabel="# Iterations to Convergence";
        axis_defaults...)
    barplot!(ax, xs, iters_mean; color=colors, width=0.6)
    errorbars!(ax, xs, iters_mean, iters_std; color=:black, whiskerwidth=10)
    for (i, trials) in enumerate(iters_trials)
        conv = convergence_container[methods[i]]
        for t = 1:n_trials
            mcolor = conv[t] > 0 ? :black : :gray
            marker = conv[t] > 0 ? :circle : :xcross
            scatter!(ax, [xs[i]], [trials[t]];
                color=mcolor, marker=marker, markersize=7, alpha=0.55)
        end
    end
    ax.xticks = (xs, methods)

    # Deferral objective (converged trials only)
    ax = Axis(f[1, 4],
        xlabel="Warmstart Method",
        ylabel=L"$\Sigma$ Deferral Times [s]";
        axis_defaults...)
    defer_finite = findall(isfinite, defer_mean)
    if !isempty(defer_finite)
        barplot!(ax, xs[defer_finite], defer_mean[defer_finite];
            color=colors[defer_finite], width=0.6)
        errorbars!(ax, xs[defer_finite], defer_mean[defer_finite], defer_std[defer_finite];
            color=:black, whiskerwidth=10)
    end
    for (i, trials) in enumerate(defer_trials)
        isempty(trials) && continue
        scatter!(ax, fill(xs[i], length(trials)), trials;
            color=:black, markersize=7, alpha=0.55)
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
