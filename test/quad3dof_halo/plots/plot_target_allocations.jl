using CairoMakie
using Colors
include("../../utils/plot_utils.jl")
include("plot_defaults.jl")

function plot_target_allocations(params, results; interactive=true, f=nothing, ax_idx=[1,1])

    if isnothing(f)
        f = Figure(size=(800,500))
    end
    defaults = Dict(:xgridvisible=>false, :ygridvisible=>false, :xautolimitmargin=>(0,0), :yautolimitmargin=>(0,0), :topspinevisible=>true, :rightspinevisible=>true)

    # Obtain target allocations
    target_allocations = get_target_allocations(results)
    R_targs_max = maximum(results["targpool_radii"])
    # colors = distinguishable_colors(length(results["targpool_ID"]), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)

    function map_bounds_saturate(x, bounds_in, bounds_out)
        # Saturate
        for k = 1:length(x)
            x[k] = max(min(x[k], bounds_in[2]), bounds_in[1])
        end
        # Map
        x = (x .- bounds_in[1]) ./ (bounds_in[2] - bounds_in[1]) * (bounds_out[2] - bounds_out[1]) .+ bounds_out[1]
    end

    function plot_target_history(ax, results, target_allocations, pool_idx; pad_in = 0.1, rad_min = 0, rad_max = 1, color=RGB(0,0,0))
        target_id = results["targpool_ID"][pool_idx]
        time_history = results["sim_time"]
        radii_history = results["targpool_radii"][pool_idx,:]
        plot_bounds = [pool_idx - .5 + pad_in, pool_idx + .5 - pad_in]
        color_dark = dark_color(color, fraction=0.3)
        color_light = bright_color(color, fraction=0.6)

        # Plot unallocated band
        band!(ax, [time_history[1], time_history[end]], [plot_bounds[1], plot_bounds[1]], [plot_bounds[2], plot_bounds[2]]; color=:gray95)

        # Plot allocation bands
        if target_id in keys(target_allocations)
            for span in target_allocations[target_id]
                id0,idf = span
                band!(ax, [time_history[id0], time_history[idf]], [plot_bounds[1], plot_bounds[1]], [plot_bounds[2], plot_bounds[2]]; color=color_light)
            end
        end

        # Plot radii history
        lines!(ax, time_history, map_bounds_saturate(radii_history, [rad_min, rad_max], plot_bounds); color=color_dark)
    end

    ax = Axis(f[ax_idx...], xlabel="Time [s]", ylabel="Target ID"; defaults...)
    for k = 1:length(results["targpool_ID"])
        id = results["targpool_ID"][k]
        plot_target_history(ax, results, target_allocations, k, rad_min=params.R_targs_min, rad_max=R_targs_max, color=target_colors[id])
    end
    ax.yticks = results["targpool_ID"]
    
    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    end
end
