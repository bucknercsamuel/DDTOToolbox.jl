using CairoMakie
# using GLMakie
using Colors
using InvertedIndices
include("../../utils/plot_utils.jl")
include("plot_defaults.jl")

# Generic styling
style2D_dt = Dict(:color=>:gray, :marker=>:circle, :markersize=>10, :strokecolor=>:black, :strokewidth=>3) 
style2D_ct = Dict(:color=>:black, :linewidth=>4)

function plot_trajs(
        solutions,
        simulations,
        params;
        interactive = true,
        ddto = true, 
        obstacles = true,
        projection_indices = [1,2] # x-y 2D projection
    )
    # Axis settings
    axis_defaults = Dict(
        # :xautolimitmargin=>(0,0), 
        :topspinevisible=>true, 
        :rightspinevisible=>true,
        :xgridvisible=>false,
        :ygridvisible=>false,
        :aspect=>DataAspect())

    # Setup
    f = Figure(size=(520,800))
    ax = Axis(f[1,1], xlabel=L"$x$-position [m]", ylabel=L"$y$-position [m]"; axis_defaults...)
    J = projection_indices

    # Color conditions
    color_branch = n -> 0
    if length(solutions) == 1
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

    # Plot obstacles
    if obstacles
        flag_labeled = false
        for o = 1:params.n_obstacles
            if !flag_labeled
                label = "Obstacles"
                flag_labeled = true
            else
                label = ""
            end
            draw2d_circle(ax, params.p_obstacles[1:2,o], params.R_obstacles[o]; color="red")
            # plt.text(params.p_obstacles[1,o], params.p_obstacles[2,o], string(o), color="black", fontsize=12, ha="center", va="center", fontweight="extra bold")
        end
    end

    # Plot trajectories
    for k = 1:n_solutions
        if length(solutions) > 1
            color_branch = n -> color_map_bundles[k]
        end
        plot2D_bundle(ax,
            [solutions[k].targs[j].r[J[1],:] for j∈1:params.a.n_targs],
            [simulations[k].targs[j].r[J[1],:] for j∈1:params.a.n_targs],
            [solutions[k].targs[j].r[J[2],:] for j∈1:params.a.n_targs],
            [simulations[k].targs[j].r[J[2],:] for j∈1:params.a.n_targs],
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
    else
        CairoMakie.save(joinpath(fig_path, "trajs"*fig_ext), f)
    end
end