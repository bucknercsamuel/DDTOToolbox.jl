# using CairoMakie
using GLMakie
using Colors
using InvertedIndices
include("../utils/plot_utils.jl")

# Generic styling
style2D_dt = Dict(:color=>:gray, :marker=>:circle, :markersize=>15, :strokecolor=>:black, :strokewidth=>3) 
style2D_ct = Dict(:color=>:black, :linewidth=>3)

# Themes
theme2d = merge(theme_minimal(), theme_latexfonts())
fontsize = 20

function build_plots_single(sols, sims, params; ddto=true, interactive=true)
    screens = []
    with_theme(theme2d; fontsize=fontsize) do
        push!(screens, plot_trajs(sols, sims, params; interactive=interactive, ddto=ddto))
    end

    if interactive
        println("\nPress any key when finished using plots...")
        readline() # Wait for user to finish plotting
        [GLMakie.destroy!(screen) for screen in screens]
    end
end

function build_plots_compare_cvx_scp(ddtocvx_sols, ddtocvx_sims, ddtoscp_sols, ddtoscp_sims, params_cvx, params_scp; ddto=true, interactive=true)
    screens = []
    labels = ["DDTO-Cvx", "DDTO-SCP"]
    with_theme(theme2d; fontsize=fontsize) do
        push!(screens, plot_compare([ddtocvx_sols, ddtoscp_sols], [ddtocvx_sims, ddtoscp_sims], [params_cvx, params_scp]; interactive=interactive, titles=labels, ddto=ddto))
    end

    if interactive
        println("\nPress any key when finished using plots...")
        readline() # Wait for user to finish plotting
        [GLMakie.destroy!(screen) for screen in screens]
    end
end

function plot_trajs(
        solutions,
        simulations,
        params;
        interactive = true,
        ddto = true, 
        projection_indices = [1,2] # x-y 2D projection
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
    ax = Axis(f[1,1], xlabel=L"$x$-position [m]", ylabel=L"$y$-position [m]"; axis_defaults...)
    J = projection_indices

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
        save("dintegrator2dof/figures/plot_trajs.png", f; px_per_unit = 4)
        return screen
    else
        save("dintegrator2dof/figures/plot_trajs.svg", f)
    end
end

function plot_compare(
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

    # Plot trajectories
    τ_lu(k,j) = params[k].τ_targs[findfirst(i->i==j, params[k].λ_targs)] # obtain the deferrability index of the j-th target (solution)
    for m = 1:n_solutions
        ax = Axis(f[1,m], xlabel=L"$x$-position [m]", ylabel=L"$y$-position [m]", title=titles[m]; axis_defaults...)
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
        save("dintegrator2dof/figures/plot_compare.png", f; px_per_unit = 4)
        # save("dintegrator2dof/figures/compare_ddto_weight1.png", f; px_per_unit = 4)
        # save("dintegrator2dof/figures/compare_ddto_weight2.png", f; px_per_unit = 4)
        # save("dintegrator2dof/figures/compare_ddto_weight3.png", f; px_per_unit = 4)
        # save("dintegrator2dof/figures/compare_ddto_equalweights.png", f; px_per_unit = 4)
        return screen
    else
        save("dintegrator2dof/figures/plot_compare.svg", f)
    end
end

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
    f = Figure(size=(1200,800))

    # Color conditions
    color_branch = n -> 0
    if length(solutions) == 1
        # if params.a.n_targs <= 4
        #     colors = Colors.JULIA_LOGO_COLORS
        # else
        #     colors = range(HSV(0,1,1), stop=HSV(-360,1,1), length=params.a.n_targs)
        # end
        # color_map_bundles = cgrad(:thermal, params.a.n_targs, categorical=true)
        color_map_bundles = range(HSV(colorant"blue"), stop=HSV(colorant"orange"), length=params.a.n_targs)
        color_branch = n -> color_map_bundles[n] # contains all colors
    else
        color_map_bundles = cgrad(:thermal, length(solutions), categorical=true)
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

function plot_accel_norm(
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
    f = Figure(size=(1200,800))

    # Color conditions
    color_branch = n -> 0
    if length(solutions) == 1
        # if params.a.n_targs <= 4
        #     colors = Colors.JULIA_LOGO_COLORS
        # else
        #     colors = range(HSV(0,1,1), stop=HSV(-360,1,1), length=params.a.n_targs)
        # end
        # color_map_bundles = cgrad(:thermal, params.a.n_targs, categorical=true)
        color_map_bundles = range(HSV(colorant"blue"), stop=HSV(colorant"orange"), length=params.a.n_targs)
        color_branch = n -> color_map_bundles[n] # contains all colors
    else
        color_map_bundles = cgrad(:thermal, length(solutions), categorical=true)
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
        if params.a.disc == 0
            norms_sol = [[[norm(solutions[k].targs[j].a[:,l]) for l=1:size(solutions[k].targs[j].a)[2]]...; norm(solutions[k].targs[j].a[:,end])] for j∈1:params.a.n_targs]
        elseif params.a.disc == 1
            norms_sol = [[norm(solutions[k].targs[j].a[:,l]) for l=1:size(solutions[k].targs[j].a)[2]] for j∈1:params.a.n_targs]
        end
        norms_sim = [[norm(simulations[k].targs[j].a[:,l]) for l=1:size(simulations[k].targs[j].a)[2]] for j∈1:params.a.n_targs]
        plot2D_bundle(ax,
            [solutions[k].targs[j].t for j∈1:params.a.n_targs],
            [simulations[k].targs[j].t for j∈1:params.a.n_targs],
            norms_sol,
            norms_sim,
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

function timing_comparison(
        x_range,
        solve_times_cvx,
        solve_times_scp
    )
    with_theme(theme2d; fontsize=fontsize) do
        # Axis settings
        axis_defaults = Dict(
            :topspinevisible=>true, 
            :rightspinevisible=>true,
            :xgridvisible=>true,
            :ygridvisible=>true,
            :xminorticksvisible=>true,
            :yminorticksvisible=>true,
            :xminorticks=>IntervalsBetween(10),
            :yminorticks=>IntervalsBetween(10),
            # :xscale=>log10,
            # :yscale=>log10
        )

        # Setup
        f = Figure(size=(800,600))
        ax = Axis(f[1,1], xlabel=L"$n$", ylabel="Solve time [s]"; axis_defaults...)

        # Get statistics on data
        means = Vector{Vector{Float64}}(undef,2)
        stds  = Vector{Vector{Float64}}(undef,2)
        solve_times = Vector{Matrix{Float64}}(undef,2)
        solve_times[1] = solve_times_cvx
        solve_times[2] = solve_times_scp
        for j = 1:2
            means[j] = Vector(undef,size(solve_times[j])[1])
            stds[j] = Vector(undef,size(solve_times[j])[1])
            for k = 1:length(x_range)
                vec = solve_times[j][k,:]
                vec_noinfs = []
                n_vec = length(vec)
                for m = 1:n_vec
                    if !isinf(vec[m])
                        append!(vec_noinfs, vec[m])
                    end
                end
                if length(vec_noinfs) == 0
                    error("can't have all samples be inf..")
                end
                means[j][k] = mean(vec_noinfs)
                stds[j][k] = std(vec_noinfs)
            end
        end

        # Styling
        style_scatter = Dict(:markersize=>15, :strokewidth=>0) 
        style_line    = Dict(:linewidth=>2)
        style_fill    = Dict(:alpha=>.25)

        # Plot loop
        colors = [:red,:blue]
        labels = ["DDTO-Cvx", "DDTO-SCP"]
        for k = 1:length(means)
            # Means plot
            lines!(ax,
                x_range,
                means[k];
                style_line..., :color=>colors[k], label=labels[k])
            scatter!(ax,
                x_range,
                means[k];
                style_scatter..., :color=>colors[k], :strokecolor=>colors[k])

            # Stdev plot
            std_upper = means[k] .+ stds[k]
            std_lower = means[k] -+ stds[k]
            std_lower = [max(std_lower[j],means[k][j]-1e-1) for j=1:length(means[k])]
            band!(ax,
                x_range,
                std_lower,
                std_upper;
                style_fill..., :color=>(colors[k], style_fill[:alpha])
            )
        end
        Legend(f[1,2], ax, framevisible=false)

        screen = GLMakie.Screen()
        display(screen, f)
        save("dintegrator2dof/figures/solve_time_compare.png", f; px_per_unit = 4)
    end
end