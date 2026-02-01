# using CairoMakie
using GLMakie
using Colors
using InvertedIndices
include("../../../utils/plot_utils.jl")

# Generic styling
style2D_dt = Dict(:color=>:gray, :marker=>:circle, :markersize=>7.5, :strokecolor=>:black, :strokewidth=>3) 
style2D_ct = Dict(:color=>:black, :linewidth=>3)

# Themes
theme2d = merge(theme_minimal(), theme_latexfonts())
fontsize = 20

function plot_time_dilation(
        solutions,
        simulations,
        params;
        interactive = true,
        ddto = true
    )
    # Axis settings
    axis_defaults = Dict(
        :xautolimitmargin=>(0,0), 
        :topspinevisible=>true, 
        :rightspinevisible=>true,
        :xgridvisible=>true,
        :ygridvisible=>true)

    # Setup
    f = Figure(size=(1600,1000))

    # Color conditions
    color_branch = n -> 0
    if length(solutions) == 1
        # if params.a.n_targs <= 4
        #     colors = Colors.JULIA_LOGO_COLORS
        # else
        #     colors = range(HSV(0,1,1), stop=HSV(-360,1,1), length=params.a.n_targs)
        # end
        color_map_bundles = cgrad(:rainbow, params.a.n_targs, categorical=true)
        color_branch = n -> color_map_bundles[n] # contains all colors
    else
        color_map_bundles = cgrad(:rainbow, length(solutions), categorical=true)
    end

    # Flag conditions
    n_solutions = length(solutions)
    show_ddto_split = n_solutions > 1 || !ddto ? false : true
    # show_ddto_split = true
    show_sol_nodes = n_solutions > 1 ? false : true
    # show_sol_nodes = true
    show_defer_nodes = ddto ? true : false

    # Plot time
    ax = Axis(f[1,1], xlabel=L"$\tau$", ylabel=L"$t(\tau)$ [s]"; axis_defaults...)
    for k = 1:n_solutions
        if length(solutions) > 1
            color_branch = n -> color_map_bundles[k]
        end
        plot2D_bundle(ax,
            [solutions[k].targs[j].τ for j∈1:params.a.n_targs],
            [simulations[k].targs[j].τ for j∈1:params.a.n_targs],
            [solutions[k].targs[j].t for j∈1:params.a.n_targs],
            [simulations[k].targs[j].t for j∈1:params.a.n_targs],
            params,
            style2D_ct,
            style2D_dt;
            color_branch = color_branch,
            show_sol_nodes = show_sol_nodes,
            show_defer_nodes = show_defer_nodes,
            show_ddto_split = show_ddto_split
        )
    end

    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    end
end