using CairoMakie
using Colors
using Statistics
using Printf
using PrettyTables
include("../../utils/plot_utils.jl")
include("plot_defaults.jl")

"""
    print_stats_table(data_name, stats)

Print a formatted summary table for the given MC statistics, bolding the best
(lowest, since each metric is `↓`) value in each numeric column.
"""
function print_stats_table(data_name, stats)
    order = ["Gr-1", "Gr-∞", "Graph-DDTO"]
    rows  = [stats[findfirst(s -> s.name == n, stats)] for n in order]

    fmt_mean(s) = s.name == "Graph-DDTO" ?
        @sprintf("%.2f ± %.2f", s.mean, s.std) :
        @sprintf("%.2f ± %.2f  (%+.2f%%)", s.mean, s.std, s.percent_increase)

    data = hcat(
        [s.name               for s in rows],
        [fmt_mean(s)          for s in rows],
        [@sprintf("%.2f", s.median) for s in rows],
        [@sprintf("%.2f", s.min)    for s in rows],
        [@sprintf("%.2f", s.max)    for s in rows],
        [s.num_outliers       for s in rows],
    )

    # Best (minimum) row per numeric column (2..6 in the displayed table)
    best = Dict(
        2 => argmin([s.mean         for s in rows]),
        3 => argmin([s.median       for s in rows]),
        4 => argmin([s.min          for s in rows]),
        5 => argmin([s.max          for s in rows]),
        6 => argmin([s.num_outliers for s in rows]),
    )
    hl = TextHighlighter((d, i, j) -> get(best, j, 0) == i, crayon"bold")

    println("\n$(data_name) statistics:")
    pretty_table(
        data;
        column_labels = ["Algorithm", "↓ Mean ± Std.", "↓ Median", "↓ Min", "↓ Max", "↓ # Outliers"],
        highlighters  = [hl],
    )
end

function plot_mc_statistics(solution_set, label; saturation=Inf, interactive=true, groupings::Vector = [], mapid="")

    # Build figure
    f = Figure(size=(500,400))
    defaults = Dict(
        :topspinevisible=>true, 
        :rightspinevisible=>true,
        :xgridvisible=>false,
        :ygridvisible=>false,
        # # :xminorticksvisible=>true,
        # :yminorticksvisible=>true,
    )
    
    # Custom colors for each solution type (dark, light)
    colors = [
        (:dodgerblue4, :dodgerblue1),
        (:indianred4, :indianred1),
        (:orange4, :orange1),
    ]

    function add_box_plot_entry(ax, idx, Q1, md, Q3, outliers; width=.5, color_dark=:red, color_light=:pink, saturate_zero=false, alpha_fill=0.5, linewidth_scale=.5)
        w = width
        IQR = Q3 - Q1
        mn = saturate_zero ? max(Q1 - 1.5*IQR,1e-3) : Q1 - 1.5*IQR
        mx = Q3 + 1.5*IQR

        # Make IQR box filled color_dark with white circle in center to represent median
        w_box = w/8
        band!(ax, [idx-w_box/2, idx+w_box/2], [Q1, Q1], [Q3, Q3]; color=color_dark, alpha=1)
        scatter!(ax, [idx], [md]; color=:white, markersize=10)

        # Add line component below the IQR to min
        lines!(ax, [idx, idx], [Q1, mn]; color=color_dark, linewidth=2*linewidth_scale)
        lines!(ax, [idx-w_box/4, idx+w_box/4], [mn, mn]; color=color_dark, linewidth=2*linewidth_scale)

        # Add line component above the IQR to max
        lines!(ax, [idx, idx], [Q3, mx]; color=color_dark, linewidth=2*linewidth_scale)
        lines!(ax, [idx-w_box/4, idx+w_box/4], [mx, mx]; color=color_dark, linewidth=2*linewidth_scale)

        # # fill
        # band!(ax, [idx-w/2, idx+w/2], [md, md], [Q3, Q3]; color=color_light, alpha=alpha_fill)
        # band!(ax, [idx-w/2, idx+w/2], [Q1, Q1], [md, md]; color=color_light, alpha=alpha_fill)

        # # horizontal lines
        # lines!(ax, [idx-w/4, idx+w/4], [mn, mn]; color=color_dark, linewidth=2*linewidth_scale)
        # lines!(ax, [idx-w/2, idx+w/2], [Q1, Q1]; color=color_dark, linewidth=2*linewidth_scale)
        # lines!(ax, [idx-w/2, idx+w/2], [md, md]; color=color_dark, linewidth=4*linewidth_scale)
        # lines!(ax, [idx-w/2, idx+w/2], [Q3, Q3]; color=color_dark, linewidth=2*linewidth_scale)
        # lines!(ax, [idx-w/4, idx+w/4], [mx, mx]; color=color_dark, linewidth=2*linewidth_scale)

        # # vertical lines
        # lines!(ax, [idx-w/2, idx-w/2], [Q1, Q3]; color=color_dark, linewidth=2*linewidth_scale)
        # lines!(ax, [idx+w/2, idx+w/2], [Q1, Q3]; color=color_dark, linewidth=2*linewidth_scale)
        # lines!(ax, [idx, idx], [mn, Q1]; color=color_dark, linewidth=2*linewidth_scale)
        # lines!(ax, [idx, idx], [Q3, mx]; color=color_dark, linewidth=2*linewidth_scale)

        # Outliers
        if length(outliers) > 0
            scatter!(ax, fill(idx, length(outliers)), outliers; color=color_dark, markersize=10*linewidth_scale)
        end
    end

    function add_plot_entries(ax, solution_set, data_name; colors=[], saturate_zero=false, groupings=[], width_factor=.4, show_violin=true, outlier_threshold=Inf, plot_saturation=Inf)
        if length(groupings) == 0
            groupings = [(i,) for i in 1:length(solution_set)]
        end
        box_pos = 1.
        box_poses = []
        key_order = ["Graph-DDTO", "Gr-1", "Gr-∞"]
        mean_graphSCvx = nothing
        data_plot = nothing
        stats_table = []
        for (iter,key) in enumerate(key_order)
            # Raw data statistical analysis (uses IQR band to identify outliers)
            value = solution_set[key]
            append!(box_poses, box_pos)
            idx_feas = findall(τ->τ==1, [value[k]["error_code"] for k∈1:length(value)])
            data = [value[k][string(data_name)] for k∈idx_feas]            
            quant_data(p) = quantile(data, p)
            Q1,median,Q3 = map(quant_data, [.25,.5,.75])
            IQR = Q3-Q1

            # Identify extreme outliers and disregard these from plots
            @assert outlier_threshold > Q3+1.5*IQR # outlier threshold must be greater than upper bound of outlier detection metric
            idx_outliers = findall(x->(x<Q1-1.5*IQR).|(x>Q3+1.5*IQR), data)
            idx_extreme_outliers = findall(x->x>outlier_threshold, data)
            data_outliers = data[idx_outliers]
            data_plot = copy(data)
            if length(idx_extreme_outliers) > 0
                data_plot = data_plot[setdiff(1:length(data_plot), idx_extreme_outliers)]
            end
            data_plot_outliers = data_plot[findall(x->(x<Q1-1.5*IQR).|(x>Q3+1.5*IQR), data_plot)]
            data_plot_inliers = setdiff(data_plot, data_plot_outliers)

            # Add violin plot (with inliers only)
            if show_violin
                violin!(ax, fill(iter,length(data_plot_inliers)), data_plot_inliers; color=colors[iter][2], scale=:width, width=width_factor*length(solution_set))
            end

            # Add box plot overlay (with outliers)
            add_box_plot_entry(ax, box_pos, Q1, median, Q3, data_plot_outliers; width=width_factor*length(solution_set), color_dark=colors[iter][1], color_light=colors[iter][2], saturate_zero=saturate_zero, alpha_fill=0.5)
            grouping_idx = findfirst(g->iter in g, groupings)
            group_idx = findfirst(g-> iter in g, groupings[grouping_idx])
            if group_idx < length(groupings[grouping_idx])
                box_pos += 1.05*width_factor*length(solution_set)
            else
                box_pos += 1.
            end

            # Cache GraphSCvx data for later comparison
            if key == "Graph-DDTO"
                mean_graphSCvx = mean(data)
            end
            percent_increase = (mean(data) - mean_graphSCvx) / mean_graphSCvx * 100

            # Collect stats for tabular output after the loop
            push!(stats_table, (
                name             = key,
                mean             = mean(data),
                std              = std(data),
                median           = median,
                min              = minimum(data),
                max              = maximum(data),
                num_outliers     = length(data_outliers),
                percent_increase = percent_increase,
            ))
        end

        # Print stats table (bold = best, i.e. lowest, in each column)
        print_stats_table(data_name, stats_table)

        # Customize ticks
        labels = key_order
        label_pointers = Dict([(box_poses[k], labels[k]) for k in 1:length(labels)])
        ax.xticks = box_poses
        ax.xtickformat = values -> [label_pointers[value] for value in values]
        # hidedecorations!(ax, label=false, ticklabels=false, ticks=false, minorticks=false)
    end

    if label == "cum_thrust"
        ylabel = "Cumulative thrust [N]"
    elseif label == "cum_energy"
        ylabel = "Cumulative energy [N^(3/2)]"
    elseif label == "ATE"
        ylabel = "Average trajectory error [m]"
    elseif label == "num_recomputations"
        ylabel = "Number of recomputations"
    elseif label == "radius_at_cutoff"
        ylabel = "Largest radius at cutoff [m]"
    elseif label == "safety_at_cutoff"
        ylabel = "Safety at cutoff [m]"
    end

    ax = Axis(f[1,1], ylabel=ylabel; defaults...)
    add_plot_entries(ax, solution_set, label; colors=colors, groupings=groupings, outlier_threshold=saturation)
    
    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    else
        CairoMakie.save(joinpath(fig_path, "mc_$(label)_$(mapid)"*fig_ext), f)
    end
end
