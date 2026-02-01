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

function plot_traj_compare(
        solutions,
        simulations,
        params;
        interactive = true,
        ddto = true, 
        obstacles = true,
        projection_indices = [1,2], # x-y 2D projection
        titles = nothing
    )
    # Axis settings
    axis_defaults = Dict(
        # :xautolimitmargin=>(0,0), 
        :topspinevisible=>true, 
        :rightspinevisible=>true,
        :xgridvisible=>true,
        :ygridvisible=>true,
        :aspect=>DataAspect())

    # Setup
    f = Figure(size=(1200,800))
    J = projection_indices
    n_solutions = length(solutions)
    if isnothing(titles)
        titles = fill("", n_solutions)
    end

    # Color conditions
    # color_map_bundles = range(HSV(colorant"blue"), stop=HSV(colorant"orange"), length=params.a.n_targs)
    color_map_bundles = cgrad(:rainbow, params[1].a.n_targs, categorical=true)
    color_branch = n -> color_map_bundles[n] # contains all colors

    # Flag conditions
    show_ddto_split = !ddto ? false : true
    show_sol_nodes = true
    show_defer_nodes = ddto ? true : false

    τ_lu(k,j) = params[k].τ_targs[findfirst(i->i==j, params[k].λ_targs)] # obtain the deferrability index of the j-th target (solution)
    for m = 1:n_solutions
        ax = Axis(f[1,m], xlabel=L"$x$-position [m]", ylabel=L"$y$-position [m]", title=titles[m]; axis_defaults...)

        # Plot obstacles
        if obstacles
            flag_labeled = false
            for o = 1:params[1].n_obstacles
                if !flag_labeled
                    label = "Obstacles"
                    flag_labeled = true
                else
                    label = ""
                end
                draw2d_circle(ax, params[1].p_obstacles[1:2,o], params[1].R_obstacles[o]; color=:red)
            end
        end

        # Plot trajectories
        for k ∈ [collect(1:n_solutions)[1:end .!= m]; m]
            show_defer_times = k == m ? true : false
            alpha = k == m ? 1 : 0.1
            plot2D_bundle(ax,
                [solutions[k].targs[j].r[J[1],:] for j∈1:params[k].a.n_targs],
                [simulations[k].targs[j].r[J[1],:] for j∈1:params[k].a.n_targs],
                [solutions[k].targs[j].r[J[2],:] for j∈1:params[k].a.n_targs],
                [simulations[k].targs[j].r[J[2],:] for j∈1:params[k].a.n_targs],
                params[k],
                style2D_ct,
                style2D_dt;
                color_branch = color_branch,
                show_sol_nodes = show_sol_nodes,
                show_defer_nodes = show_defer_nodes,
                show_ddto_split = show_ddto_split,
                show_defer_times = show_defer_times,
                defer_times = [solutions[k].targs[j].t[τ_lu(k,j)] for j∈1:params[k].a.n_targs],
                alpha = alpha
            )
        end
    end

    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        save("quad3dofcage/figures/plot_compare.png", f; px_per_unit = 4)
        return screen
    else
        save("quad3dofcage/figures/plot_compare.png", f; px_per_unit = 4)
    end
end