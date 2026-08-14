using CairoMakie
using Colors
using Statistics
include("../../utils/plot_utils.jl")
include("plot_defaults.jl")

"""
    plot_mc_compare(results_dict; interactive=true)

Monte Carlo statistical comparison plots:
Σ deferral, solver time, % converged, and # iterations to convergence.
"""
function plot_mc_compare(results_dict; interactive=true)
    axis_defaults = Dict(
        :topspinevisible => true,
        :rightspinevisible => true,
        :xgridvisible => false,
        :ygridvisible => false,
    )
    ax_label_size = 15

    methods = results_dict["methods"]
    n_sweep = sort(Int.(results_dict["n_sweep"]))
    n_trials = results_dict["n_trials"]
    x_range = n_sweep

    method_style = Dict(
        "lex" => (label="Lex-DDTO", color=:red),
        "scp" => (label="Graph-DDTO", color=:blue),
        "qcvx" => (label="DDTO-QCvx", color=:green),
    )

    # Aggregate Σ deferral per trial
    deferral_obj = Dict{String,Dict{Int,Vector{Float64}}}()
    conv_pct = Dict{String,Dict{Int,Vector{Float64}}}()
    for m in methods
        deferral_obj[m] = Dict{Int,Vector{Float64}}()
        conv_pct[m] = Dict{Int,Vector{Float64}}()
        for n in n_sweep
            deferral_obj[m][n] = [sum(results_dict["deferral_time"][m][n][t, :]) for t = 1:n_trials]
            conv_pct[m][n] = results_dict["convergence"][m][n] .* 100
        end
    end

    function plot_mean_and_funnel(ax, data, label, colors; funnel=true, convergences=nothing, saturate_zero=false)
        means = Float64[]
        stds = Float64[]
        for n in n_sweep
            data_trials = Float64[]
            for i = 1:n_trials
                proceed = isnothing(convergences) ? true : (convergences[n][i] > 0.0)
                val = data[n][i]
                if proceed && isfinite(val)
                    push!(data_trials, val)
                end
            end
            if !isempty(data_trials)
                push!(means, mean(data_trials))
                push!(stds, length(data_trials) > 1 ? std(data_trials) : 0.0)
            else
                push!(means, NaN)
                push!(stds, NaN)
            end
        end
        means_upper = means .+ stds
        means_lower = means .- stds
        if saturate_zero
            means = [max(m, 1e-10) for m in means]
            means_upper = [max(m, 1e-10) for m in means_upper]
            means_lower = [max(m, 1e-10) for m in means_lower]
        end
        finite = findall(isfinite, means)
        isempty(finite) && return
        if funnel
            band!(ax, x_range[finite], means_lower[finite], means_upper[finite];
                color=colors, alpha=0.2)
        end
        lines!(ax, x_range[finite], means[finite];
            color=colors, linewidth=2, label=label)
        ax.xticks = (x_range, string.(x_range))
        for (k, n) in enumerate(n_sweep)
            println("  $label @ n=$n: $(means[k])")
        end
    end

    f = Figure(size=(1800, 400))

    # Σ deferral (converged only)
    ax = Axis(f[1, 1], xlabel="Number of Targets", ylabel=L"$\Sigma$ Deferral Times [s]"; axis_defaults...)
    println("Σ deferral means:")
    for m in methods
        st = method_style[m]
        plot_mean_and_funnel(ax, deferral_obj[m], st.label, st.color; convergences=conv_pct[m])
    end
    axislegend(ax, position=:lt, labelsize=ax_label_size)

    # Solver time (all completed trials)
    ax = Axis(f[1, 2], xlabel="Number of Targets", ylabel="Solver Time [s]"; axis_defaults...)
    println("Solver time means:")
    for m in methods
        st = method_style[m]
        plot_mean_and_funnel(ax, results_dict["solver_time"][m], st.label, st.color)
    end
    axislegend(ax, position=:lt, labelsize=ax_label_size)

    # % converged
    ax = Axis(f[1, 3], xlabel="Number of Targets", ylabel="% Converged"; axis_defaults...)
    println("Convergence means:")
    for m in methods
        st = method_style[m]
        plot_mean_and_funnel(ax, conv_pct[m], st.label, st.color; funnel=false)
    end
    axislegend(ax, position=:lb, labelsize=ax_label_size)
    ylims!(ax, 0, 105)

    # Iterations to convergence (PTR methods; QCvx is NaN)
    ax = Axis(f[1, 4], xlabel="Number of Targets", ylabel="# Iterations to Convergence"; axis_defaults...)
    println("Iteration means (converged PTR trials):")
    for m in methods
        st = method_style[m]
        plot_mean_and_funnel(ax, results_dict["iters"][m], st.label, st.color; convergences=conv_pct[m])
    end
    axislegend(ax, position=:lt, labelsize=ax_label_size)

    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    else
        out_dir = joinpath(@__DIR__, "..", "figures")
        mkpath(out_dir)
        mode = get(results_dict, "obstacle_mode", "unknown")
        out = joinpath(out_dir, "mc_compare_$(mode)" * fig_ext)
        CairoMakie.save(out, f)
        println("Saved figure to $out")
        return nothing
    end
end
