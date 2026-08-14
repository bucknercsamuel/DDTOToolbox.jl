using CairoMakie
using Colors
using InvertedIndices
include("../../utils/plot_utils.jl")
include("plot_defaults.jl")

"""
    plot_ddto_compare(solutions, simulations, params_list, titles; interactive=false)

Side-by-side 2D trajectory panels for Lex / Graph / QCvx DDTO solutions.
"""
function plot_ddto_compare(
        solutions,
        simulations,
        params_list,
        titles;
        interactive=false,
        fig_name="ddto_compare",
    )
    axis_defaults = Dict(
        :topspinevisible => true,
        :rightspinevisible => true,
        :xgridvisible => false,
        :ygridvisible => false,
        :aspect => DataAspect(),
    )
    style2D_dt = Dict(:color => :gray, :marker => :circle, :markersize => 8,
                      :strokecolor => :black, :strokewidth => 2)
    style2D_ct = Dict(:color => :black, :linewidth => 3)

    n = length(solutions)
    f = Figure(size=(420 * n, 520))
    J = [1, 2]
    color_map = cgrad(:rainbow, params_list[1].a.n_targs, categorical=true)
    color_branch = k -> color_map[k]

    for m = 1:n
        ax = Axis(f[1, m],
            xlabel=L"$x$-position [m]",
            ylabel=L"$y$-position [m]",
            title=titles[m];
            axis_defaults...)
        params = params_list[m]

        for o = 1:params.n_obstacles
            draw2d_circle(ax, params.p_obstacles[1:2, o], params.R_obstacles[o]; color=:red)
        end

        plot2D_bundle(ax,
            [solutions[m].targs[j].r[J[1], :] for j = 1:params.a.n_targs],
            [simulations[m].targs[j].r[J[1], :] for j = 1:params.a.n_targs],
            [solutions[m].targs[j].r[J[2], :] for j = 1:params.a.n_targs],
            [simulations[m].targs[j].r[J[2], :] for j = 1:params.a.n_targs],
            params,
            style2D_ct,
            style2D_dt;
            color_branch=color_branch,
            show_sol_nodes=true,
            show_defer_nodes=true,
            show_ddto_split=true,
        )
    end

    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    else
        out_dir = joinpath(@__DIR__, "..", "figures")
        mkpath(out_dir)
        out = joinpath(out_dir, fig_name * fig_ext)
        CairoMakie.save(out, f)
        println("Saved figure to $out")
        return nothing
    end
end
