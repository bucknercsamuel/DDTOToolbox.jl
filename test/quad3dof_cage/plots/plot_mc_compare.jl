using CairoMakie
using Colors
using Statistics
using Printf
using PrettyTables
include("../../utils/plot_utils.jl")
include("plot_defaults.jl")

const METHOD_STYLE = Dict(
    "qcvx" => (label = "QCVX-DDTO", color = :seagreen),
    "lex"  => (label = "Lex-DDTO",  color = :red),
    "scp"  => (label = "Graph-DDTO", color = :dodgerblue),
)

function _int_keys(by_n)
    return sort(Int[Int(k) for k in keys(by_n)])
end

function _case_n_range(case_data)
    ns = Int[]
    for by_n in values(case_data)
        append!(ns, _int_keys(by_n))
    end
    return sort(unique(ns))
end

function _mean_std(values; only_converged=false, converged=nothing)
    data = Float64[]
    for i in eachindex(values)
        v = values[i]
        if only_converged
            (converged === nothing || !converged[i]) && continue
        end
        isfinite(v) && push!(data, v)
    end
    isempty(data) && return (NaN, NaN)
    length(data) == 1 && return (data[1], 0.0)
    return (mean(data), std(data))
end

function _series(case_data, method, n_range; field, only_converged=false)
    means = Float64[]
    stds = Float64[]
    haskey(case_data, method) || return (fill(NaN, length(n_range)), fill(NaN, length(n_range)))
    by_n = case_data[method]
    for n in n_range
        if !haskey(by_n, n) || isempty(by_n[n][field])
            push!(means, NaN)
            push!(stds, NaN)
            continue
        end
        bucket = by_n[n]
        μ, σ = _mean_std(bucket[field]; only_converged=only_converged, converged=bucket["converged"])
        push!(means, μ)
        push!(stds, σ)
    end
    return (means, stds)
end

function _conv_pct(case_data, method, n_range)
    pcts = Float64[]
    haskey(case_data, method) || return fill(NaN, length(n_range))
    by_n = case_data[method]
    for n in n_range
        if !haskey(by_n, n) || isempty(by_n[n]["converged"])
            push!(pcts, NaN)
            continue
        end
        push!(pcts, 100 * mean(by_n[n]["converged"]))
    end
    return pcts
end

function plot_mean_and_funnel!(ax, x, means, stds, color, label; funnel=true, saturate_zero=false)
    y = copy(means)
    lo = means .- stds
    hi = means .+ stds
    if saturate_zero
        y = [isfinite(v) ? max(v, 1e-10) : v for v in y]
        lo = [isfinite(v) ? max(v, 1e-10) : v for v in lo]
        hi = [isfinite(v) ? max(v, 1e-10) : v for v in hi]
    end
    finite = findall(isfinite, y)
    isempty(finite) && return
    xf = x[finite]
    if funnel
        band!(ax, xf, lo[finite], hi[finite]; color=color, alpha=0.2)
    end
    lines!(ax, xf, y[finite]; color=color, linewidth=2, label=label)
    scatter!(ax, xf, y[finite]; color=color, markersize=8)
end

function print_mc_summary(results; case_name)
    case_data = results["cases"][case_name]
    n_range = _case_n_range(case_data)
    isempty(n_range) && return
    methods = [m for m in ("qcvx", "lex", "scp") if haskey(case_data, m)]
    println("\n", case_name, "  (mean ± std; iters logged, not plotted)")
    header = ["n", "Method", "Defer [s]", "Time [s]", "% Conv", "Iters"]
    rows = Matrix{String}(undef, length(n_range) * length(methods), length(header))
    r = 0
    for n in n_range
        for method in methods
            r += 1
            style = get(METHOD_STYLE, method, (label=method, color=:black))
            defer_μ, defer_σ = _series(case_data, method, [n]; field="defer_obj", only_converged=true)
            time_μ, time_σ = _series(case_data, method, [n]; field="elapsed")
            iter_μ, iter_σ = _series(case_data, method, [n]; field="n_iters", only_converged=true)
            pct = _conv_pct(case_data, method, [n])[1]
            rows[r, 1] = string(n)
            rows[r, 2] = style.label
            rows[r, 3] = isfinite(defer_μ[1]) ? @sprintf("%.2f ± %.2f", defer_μ[1], defer_σ[1]) : "n/a"
            rows[r, 4] = isfinite(time_μ[1]) ? @sprintf("%.2f ± %.2f", time_μ[1], time_σ[1]) : "n/a"
            rows[r, 5] = isfinite(pct) ? @sprintf("%.1f", pct) : "n/a"
            rows[r, 6] = isfinite(iter_μ[1]) ? @sprintf("%.1f ± %.1f", iter_μ[1], iter_σ[1]) : "n/a"
        end
    end
    pretty_table(rows; column_labels=header, alignment=[:c, :l, :c, :c, :c, :c])
end

"""
    plot_mc_compare(results; case_name, interactive=true)

One three-panel statistical figure for a single CASE: deferrability (mean ± std
on converged trials), solve time (log-y, all trials), and percent converged.
"""
function plot_mc_compare(results; case_name, interactive=true)
    haskey(results, "cases") || error("results dict is missing the 'cases' schema")
    haskey(results["cases"], case_name) || error("no data for case $case_name")
    case_data = results["cases"][case_name]
    n_range = _case_n_range(case_data)
    if isempty(n_range)
        @warn "No Monte Carlo cells to plot" case_name
        return nothing
    end

    print_mc_summary(results; case_name=case_name)

    axis_defaults = Dict(
        :topspinevisible => true,
        :rightspinevisible => true,
        :xgridvisible => false,
        :ygridvisible => false,
    )
    ax_label_size = 15
    f = Figure(size=(1400, 400))
    methods = [m for m in ("qcvx", "lex", "scp") if haskey(case_data, m)]

    ax_def = Axis(f[1, 1], xlabel="Number of Targets", ylabel=L"Deferrability obj. $[s]$"; axis_defaults...)
    ax_time = Axis(f[1, 2], xlabel="Number of Targets", ylabel="Solver Time [s]", yscale=log10; axis_defaults...)
    ax_conv = Axis(f[1, 3], xlabel="Number of Targets", ylabel="% Converged"; axis_defaults...)

    for method in methods
        style = METHOD_STYLE[method]
        defer_μ, defer_σ = _series(case_data, method, n_range; field="defer_obj", only_converged=true)
        time_μ, time_σ = _series(case_data, method, n_range; field="elapsed")
        pct = _conv_pct(case_data, method, n_range)
        plot_mean_and_funnel!(ax_def, n_range, defer_μ, defer_σ, style.color, style.label)
        plot_mean_and_funnel!(ax_time, n_range, time_μ, time_σ, style.color, style.label; saturate_zero=true)
        plot_mean_and_funnel!(ax_conv, n_range, pct, zero(pct), style.color, style.label; funnel=false)
    end

    ax_def.xticks = (n_range, string.(n_range))
    ax_time.xticks = (n_range, string.(n_range))
    ax_conv.xticks = (n_range, string.(n_range))
    ylims!(ax_conv, 0, 105)
    axislegend(ax_def, position=:lt, labelsize=ax_label_size)
    axislegend(ax_time, position=:lt, labelsize=ax_label_size)
    axislegend(ax_conv, position=:lb, labelsize=ax_label_size)

    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    else
        mkpath(fig_path)
        out = joinpath(fig_path, "mc_compare_$(case_name)" * fig_ext)
        CairoMakie.save(out, f)
        println("Saved figure → $out")
        return nothing
    end
end
