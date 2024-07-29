# using CairoMakie
using GLMakie
using Colors
using InvertedIndices
using GeometryBasics
using LinearAlgebra
include("../utils/plot_utils.jl")

# Generic styling
style2D_dt = Dict(:color=>:gray, :marker=>:circle, :markersize=>15, :strokecolor=>:black, :strokewidth=>3) 
style2D_ct = Dict(:color=>:black, :linewidth=>3)
style2D_ct_ddto = Dict(:color=>:black, :linewidth=>1)
style3D_dt = Dict(:color=>:gray, :marker=>:circle, :markersize=>15, :strokecolor=>:black, :strokewidth=>3)
style3D_ct = Dict(:color=>:black, :linewidth=>3)
style3D_ground_base = Dict(:color=>bright_color(:orange), :transparency=>false, :alpha=>1)
style3D_ground_base_frame = Dict(:color=>bright_color(:saddlebrown))

# Themes
theme2d = merge(theme_minimal(), theme_latexfonts())
theme3d = theme_latexfonts()
fontsize = 20

function ddto_color_scheme(n_ddto_sols, n_targs)
    base_colors = ["red", "gold", "blue", "green", "purple", "pink", "brown", "cyan", "orange", "yellow"]
    idx_colors = k -> mod(k-1,length(base_colors))+1
    color_map_bundles = []
    for k = 1:n_ddto_sols
        color = base_colors[idx_colors(k)]
        color1 = parse(Colorant, color*"1")
        color2 = parse(Colorant, color*"4")
        if n_targs > 1
            cmap = range(color1, color2, length=n_targs)
        else
            cmap = color1
        end
        append!(color_map_bundles, [cmap])
    end
    return base_colors, color_map_bundles, idx_colors
end

# function ddto_color_scheme(n_ddto_sols, n_targs)
#     base_colors = range(parse(Colorant, "red1"), parse(Colorant, "red4"), length=n_ddto_sols)
#     color_map_bundles = base_colors
#     idx_colors = k -> k

#     return base_colors, color_map_bundles, idx_colors
# end

function plot_3d_trajs(
        results;
        interactive = true,
        azel=(pi/4,pi/4)
    )
    # Setup
    f = Figure(size=(800,800))
    # Create axis
    ax = Axis3(
        f[1,1],
        xlabel = "East [m]",
        ylabel = "North [m]",
        zlabel = "Up [m]",
        aspect = :equal,
        azimuth = azel[1],
        elevation = azel[2],
        xgridvisible = false,
        ygridvisible = false,
        zgridvisible = false)

    # Results parsing
    ddto_params = results["guid_update_ddto_params"]
    ddto_bundles_sol = results["guid_update_ddto_bundles"]
    ddto_bundles_sim = results["guid_update_ddto_bundles_sims"]
    sim_solution = results["sim_state"]
    n_ddto_sols = length(ddto_bundles_sol)

    # Zoom onto right area
    rmin = [min(vcat([[params.a.z0[k,:]; params.a.zf_targs[k,:]] for params in ddto_params]...)...) for k∈1:3]
    rmax = [max(vcat([[params.a.z0[k,:]; params.a.zf_targs[k,:]] for params in ddto_params]...)...) for k∈1:3]
    xLims,yLims,zLims = get_equal_3d_lims(rmin, rmax)
    xlims!(ax, xLims...)
    ylims!(ax, yLims...)
    zlims!(ax, zLims...)

    # Plot the ground
    ϵ = (zLims[2] - zLims[1]) * 0.001 # Epsilon in altitude for objects that are stacked
    ΔL(L) = L[2] - L[1]
    box_lower = [xLims[1], yLims[1], zLims[1]]
    box_upper = [ΔL(xLims), ΔL(yLims), -zLims[1]-ϵ]
    groundBaseDef = Rect3f(box_lower, box_upper)
    groundBaseMesh = GeometryBasics.mesh(groundBaseDef)
    mesh!(groundBaseMesh; style3D_ground_base...)
    boxframe_3D(ax, box_lower, box_upper; style=style3D_ground_base_frame)
    
    # Plot DDTO trajectories
    base_colors, color_map_bundles, idx_colors = ddto_color_scheme(n_ddto_sols, ddto_params[1].a.n_targs)
    proj_idxs = [1,2,3]
    for k = 1:n_ddto_sols
        params = ddto_params[k]
        color_branch = params.a.n_targs > 1 ? j -> color_map_bundles[k][j] : j -> color_map_bundles[k]
        plot_bundle(ax,
            [[ddto_bundles_sol[k].targs[j].r[c,:] for j∈1:params.a.n_targs] for c∈proj_idxs],
            [[ddto_bundles_sim[k].targs[j].r[c,:] for j∈1:params.a.n_targs] for c∈proj_idxs],
            params,
            style3D_ct,
            style3D_dt;
            color_branch = color_branch,
            color_trunk = base_colors[idx_colors(k)],
            show_sol_nodes = false,
            show_defer_nodes = true,
            show_ddto_split = true,
            alpha=0.5
        )
    end

    # Plot solution
    positions = sim_solution[1:3,:]
    lines!(ax,
           positions[1,:], positions[2,:], positions[3,:];
           style3D_ct..., :alpha=>1, :color=>:black)


    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    end
end

function plot_states(
        results;
        interactive = true,
        integrated_sim = true
    )
    f = Figure(size=(1600,1000))

    # Results parsing
    ddto_params = results["guid_update_ddto_params"]
    ddto_bundles_sol = results["guid_update_ddto_bundles"]
    ddto_bundles_sim = results["guid_update_ddto_bundles_sims"]
    sim_time = results["sim_time"]
    sim_state = results["sim_state"]
    update_times = results["guid_update_time"]
    n_ddto_sols = length(ddto_bundles_sol) - 1 # don't include last solution which is just the guidance lock
    params = ddto_params[1]

    # Default variables
    proj_idxs = [1,2,3]
    base_colors, color_map_bundles, idx_colors = ddto_color_scheme(n_ddto_sols, params.a.n_targs)
    axis_defaults_2d = Dict(
        :xautolimitmargin=>(0,0), 
        :topspinevisible=>true, 
        :rightspinevisible=>true,
        :xgridvisible=>true,
        :ygridvisible=>true)

    # Positions axes
    labels = ["Pos-East [m]", "Pos-North [m]", "Pos-Up [m]"]
    for (k,c) in enumerate(proj_idxs)
        ax = Axis(f[k,1], ylabel=labels[k]; axis_defaults_2d...)
        lines!(ax, sim_time, sim_state[c,:];
            style2D_ct..., :alpha=>1, :color=>:black)
        for k = 1:n_ddto_sols
            params = ddto_params[k]
            color_branch = params.a.n_targs > 1 ? j -> color_map_bundles[k][j] : j -> color_map_bundles[k]
            plot_bundle(ax,
                [[ddto_bundles_sol[k].targs[j].t .+ update_times[k] for j∈1:params.a.n_targs], [ddto_bundles_sol[k].targs[j].r[c,:] for j∈1:params.a.n_targs]],
                [[ddto_bundles_sim[k].targs[j].t .+ update_times[k] for j∈1:params.a.n_targs], [ddto_bundles_sim[k].targs[j].r[c,:] for j∈1:params.a.n_targs]],
                params,
                style2D_ct_ddto,
                style2D_dt;
                color_branch = color_branch,
                color_trunk = base_colors[idx_colors(k)],
                show_sol_nodes = false,
                show_defer_nodes = true,
                show_ddto_split = true,
                alpha=0.5
            )
        end
    end
    # ax_labels = ax

    # Velocities axes
    labels = ["Vel-East [m]", "Vel-North [m]", "Vel-Up [m]"]
    for (k,c) in enumerate(proj_idxs)
        ax = Axis(f[k,2], ylabel=labels[k]; axis_defaults_2d...)
        lines!(ax, sim_time, sim_state[c+3,:];
            style2D_ct..., :alpha=>1, :color=>:black)
        for k = 1:n_ddto_sols
            params = ddto_params[k]
            color_branch = params.a.n_targs > 1 ? j -> color_map_bundles[k][j] : j -> color_map_bundles[k]
            plot_bundle(ax,
                [[ddto_bundles_sol[k].targs[j].t .+ update_times[k] for j∈1:params.a.n_targs], [ddto_bundles_sol[k].targs[j].v[c,:] for j∈1:params.a.n_targs]],
                [[ddto_bundles_sim[k].targs[j].t .+ update_times[k] for j∈1:params.a.n_targs], [ddto_bundles_sim[k].targs[j].v[c,:] for j∈1:params.a.n_targs]],
                params,
                style2D_ct_ddto,
                style2D_dt;
                color_branch = color_branch,
                color_trunk = base_colors[idx_colors(k)],
                show_sol_nodes = false,
                show_defer_nodes = true,
                show_ddto_split = true,
                alpha=0.5
            )
        end
    end

    # Thrust norm axis
    ax = Axis(f[4,1:2], ylabel="Thrust Norm [N]"; axis_defaults_2d...)
    if integrated_sim
        data = sim_state[7,:]
    else
        data = [norm(results["sim_control"][1:3,l]) for l=1:length(sim_time)]
    end
    lines!(ax, sim_time, data;
        style2D_ct..., :alpha=>1, :color=>:black)
    for k = 1:n_ddto_sols
        params = ddto_params[k]
        color_branch = params.a.n_targs > 1 ? j -> color_map_bundles[k][j] : j -> color_map_bundles[k]
        plot_bundle(ax,
            [[ddto_bundles_sol[k].targs[j].t .+ update_times[k] for j∈1:params.a.n_targs], [[norm(ddto_bundles_sol[k].targs[j].T[:,l]) for l∈1:length(ddto_bundles_sol[k].targs[j].t)] for j∈1:params.a.n_targs]],
            [[ddto_bundles_sim[k].targs[j].t .+ update_times[k] for j∈1:params.a.n_targs], [[norm(ddto_bundles_sim[k].targs[j].T[:,l]) for l∈1:length(ddto_bundles_sim[k].targs[j].t)] for j∈1:params.a.n_targs]],
            params,
            style2D_ct_ddto,
            style2D_dt;
            color_branch = color_branch,
            color_trunk = base_colors[idx_colors(k)],
            show_sol_nodes = false,
            show_defer_nodes = true,
            show_ddto_split = true,
            alpha=0.5
        )
    end

    # ax = Legend(f[1:4,3], ax_labels)

    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    end
end

function plot_greedy_compare(
        results_all;
        interactive = true,
        azel=(pi/4,pi/4)
    )
    # Setup
    f = Figure(size=(800,800))
    # Create axis
    ax = Axis3(
        f[1,1],
        xlabel = "East [m]",
        ylabel = "North [m]",
        zlabel = "Up [m]",
        aspect = :equal,
        azimuth = azel[1],
        elevation = azel[2],
        xgridvisible = false,
        ygridvisible = false,
        zgridvisible = false)

    # Zoom onto right area
    params_all = []
    [[append!(params_all, [params]) for params in results["guid_update_ddto_params"]] for results in results_all]
    rmin = [min(vcat([[params.a.z0[k,:]; params.a.zf_targs[k,:]] for params in params_all]...)...) for k∈1:3]
    rmax = [max(vcat([[params.a.z0[k,:]; params.a.zf_targs[k,:]] for params in params_all]...)...) for k∈1:3]
    xLims,yLims,zLims = get_equal_3d_lims(rmin, rmax)
    xlims!(ax, xLims...)
    ylims!(ax, yLims...)
    zlims!(ax, zLims...)

    # Plot the ground
    ϵ = (zLims[2] - zLims[1]) * 0.001 # Epsilon in altitude for objects that are stacked
    ΔL(L) = L[2] - L[1]
    box_lower = [xLims[1], yLims[1], zLims[1]]
    box_upper = [ΔL(xLims), ΔL(yLims), -zLims[1]-ϵ]
    groundBaseDef = Rect3f(box_lower, box_upper)
    groundBaseMesh = GeometryBasics.mesh(groundBaseDef)
    mesh!(groundBaseMesh; style3D_ground_base...)
    boxframe_3D(ax, box_lower, box_upper; style=style3D_ground_base_frame)

    # DDTO Color conditions
    color_ddto = parse(Colorant, "blue")
    color1_greedy = parse(Colorant, "red1")
    color2_greedy = parse(Colorant, "red4")
    color_greedy = range(color1_greedy, color2_greedy, length=length(results_all)-1)
    colors = [color_ddto, color_greedy...]

    # Plot solutions
    for (k,results) in enumerate(results_all)
        sim_solution = results["sim_state"]
        positions = sim_solution[1:3,:]
        alpha = 1
        lines!(ax,
            positions[1,:], positions[2,:], positions[3,:];
            style3D_ct..., :alpha=>alpha, :color=>colors[k])
        update_times = results["guid_update_time"]
        update_idxs = [findfirst(τ->τ>=update_times[k], results["sim_time"]) for k=1:length(update_times)]
        final_radius = max(results["targs_radii"][:,end]...)
        scatter!(ax,
            positions[1,update_idxs], positions[2,update_idxs], positions[3,update_idxs];
            style3D_dt..., :alpha=>alpha, :color=>bright_color(colors[k]), :strokecolor=>(colors[k],alpha),
            :markersize=>7)
        draw_circle_3d(ax, positions[:,end], final_radius, pointing_direction=e_z, color=bright_color(colors[k]))
        end

    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    end
end

function plot_mc_trajs(
        data,
        map_data;
        interactive = true,
        hidef_ground = true,
        azel=(pi/4,pi/4)
    )
    # Setup
    f = Figure(size=(800,800))
    # Create axis
    ax = Axis3(
        f[1,1],
        xlabel = "East [m]",
        ylabel = "North [m]",
        zlabel = "Up [m]",
        # aspect = :equal,
        # aspect = DataAspect(),
        azimuth = azel[1],
        elevation = azel[2],
        xgridvisible = false,
        ygridvisible = false,
        zgridvisible = false)
    hidespines!(ax)

    # Zoom onto right area
    (xmin,xmax) = map_data["xlims"]
    (ymin,ymax) = map_data["ylims"]
    zmin = min(map_data["zlookup"]...)[2]
    zmax = max(map_data["zlookup"]...)[2] + 300.
    # xLims,yLims,zLims = get_equal_3d_lims([xmin,ymin,zmin], [xmax,ymax,zmax])
    # xlims!(ax, xLims...)
    # ylims!(ax, yLims...)
    # zlims!(ax, zLims...)
    xlims!(ax, (xmin,xmax))
    ylims!(ax, (ymin,ymax))
    zlims!(ax, (zmin,zmax))

    # Plot the ground
    if hidef_ground
        pad = 50.
        xs = LinRange(xmin+pad, xmax-pad, 1000)
        ys = LinRange(ymin+pad, ymax-pad, 1000)
        zs = [map_data["zlookup"][Int(round(x)),Int(round(y))] for x in xs, y in ys]
        surface!(ax,
            xs, ys, zs,
            colormap = :darkterrain)
    else
        ϵ = (zLims[2] - zLims[1]) * 0.001 # Epsilon in altitude for objects that are stacked
        ΔL(L) = L[2] - L[1]
        box_lower = [xLims[1], yLims[1], zLims[1]]
        box_upper = [ΔL(xLims), ΔL(yLims), -zLims[1]-ϵ]
        groundBaseDef = Rect3f(box_lower, box_upper)
        groundBaseMesh = GeometryBasics.mesh(groundBaseDef)
        mesh!(groundBaseMesh; style3D_ground_base...)
        boxframe_3D(ax, box_lower, box_upper; style=style3D_ground_base_frame)
    end

    # Color settings
    color_success = :limegreen
    color_failure = bright_color(:red, fraction=0.6)
    alpha_sim = 1
    alpha_ddto = 0.5
    
    # Plot all MC solutions
    style3D_ddto = copy(style3D_ct)
    style3D_ddto[:linewidth] = 1
    for (k,results) in enumerate(data)
        if results["error_code"] == 1
            sim_solution = results["sim_state"]
            positions = sim_solution[1:3,:]
            color = convert(Bool, results["safe_landing"]) ? color_success : color_failure
            lines!(ax,
                positions[1,:], positions[2,:], positions[3,:];
                style3D_ct..., :alpha=>alpha_sim, :color=>color)
            draw_circle_3d(ax, results["terminal_pos"], results["final_radius_truth"], pointing_direction=e_z, color=color)
        
            # Plot DDTO trajectories
            virt_params = Quad3DoFHaloParams() # virtual params for plotting only
            ddto_sols = results["guid_ddto_trajs_sols"]
            ddto_sims = results["guid_ddto_trajs_sims"]
            n_ddto_sols = size(ddto_sols,1)
            proj_idxs = [2,1,3]
            dark_frac = 0
            for k = 1:n_ddto_sols
                virt_params.a.τ_targs = results["guid_defer_vecs"][k]
                virt_params.a.λ_targs = results["guid_prefer_vecs"][k]
                virt_params.a.n_targs = virt_params.n_targs_max
                color_ = dark_color(color, fraction=dark_frac)
                plot_bundle(ax,
                    [[ddto_sols[k,j][c,:] for j∈1:virt_params.n_targs_max] for c∈proj_idxs],
                    [[ddto_sims[k,j][c,:] for j∈1:virt_params.n_targs_max] for c∈proj_idxs],
                    virt_params,
                    style3D_ddto,
                    style3D_dt;
                    color_branch = _ -> color_,
                    color_trunk = color_,
                    show_sol_nodes = false,
                    show_defer_nodes = true,
                    show_ddto_split = true,
                    alpha=alpha_ddto
                )
                dark_frac = min(dark_frac + 0.1, 0.5)
            end
        end
    end

    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    end
end

function plot_mc_statistics(solution_set, interactive=true)

    # Build figure
    f = Figure(size=(1400,1000))
    # defaults = Dict(:xgridvisible=>false, :ygridvisible=>false)
    defaults = Dict()
    n_mc = length(solution_set[collect(keys(solution_set))[1]])
    display(n_mc)

    # Custom colors for each solution type (dark, light)
    colors = [
        (:blue4, :blue1),
        (:red4, :red1),
        (:orange4, :orange1),
    ]

    function add_box_plot_entry(ax, idx, Q1, md, Q3, outliers; width=.5, color_dark=:red, color_light=:pink, saturate_zero=false)
        w = width
        IQR = Q3 - Q1
        mn = saturate_zero ? max(Q1 - 1.5*IQR,1e-3) : Q1 - 1.5*IQR
        mx = Q3 + 1.5*IQR

        # fill
        band!(ax, [idx-w/2, idx+w/2], [md, md], [Q3, Q3]; color=color_light)
        band!(ax, [idx-w/2, idx+w/2], [Q1, Q1], [md, Q3]; color=color_light)

        # horizontal lines
        lines!(ax, [idx-w/4, idx+w/4], [mn, mn]; color=color_dark, linewidth=2)
        lines!(ax, [idx-w/2, idx+w/2], [Q1, Q1]; color=color_dark, linewidth=2)
        lines!(ax, [idx-w/2, idx+w/2], [md, md]; color=color_dark, linewidth=4)
        lines!(ax, [idx-w/2, idx+w/2], [Q3, Q3]; color=color_dark, linewidth=2)
        lines!(ax, [idx-w/4, idx+w/4], [mx, mx]; color=color_dark, linewidth=2)

        # vertical lines
        lines!(ax, [idx-w/2, idx-w/2], [Q1, Q3]; color=color_dark, linewidth=2)
        lines!(ax, [idx+w/2, idx+w/2], [Q1, Q3]; color=color_dark, linewidth=2)
        lines!(ax, [idx, idx], [mn, Q1]; color=color_dark, linewidth=2)
        lines!(ax, [idx, idx], [Q3, mx]; color=color_dark, linewidth=2)

        # Outliers
        if length(outliers) > 0
            scatter!(ax, fill(idx, length(outliers)), outliers; color=color_dark)
        end
    end

    function add_box_plot_entries(ax, solution_set, data_name; colors=[], saturate_zero=false)
        for (iter,(_, value)) in enumerate(solution_set)
            idx_feas = findall(τ->τ==1, [value[k]["error_code"] for k∈1:n_mc])
            data = [value[k][string(data_name)] for k∈idx_feas]
            quant_data(p) = quantile(data, p)
            Q1,median,Q3 = map(quant_data, [.25,.5,.75])
            IQR = Q3-Q1
            outliers = findall(x->(x<Q1-1.5*IQR).|(x>Q3+1.5*IQR), data)
            add_box_plot_entry(ax, iter, Q1, median, Q3, data[outliers]; width=length(solution_set)/5, color_dark=colors[iter][1], color_light=colors[iter][2], saturate_zero=saturate_zero)
        end
    end

    axes = []
    ax = Axis(f[1,1], title="Maneuver Time", ylabel="[s]"; defaults...)
    add_box_plot_entries(ax, solution_set, "final_time"; colors=colors, saturate_zero=true)
    push!(axes, ax)

    ax = Axis(f[1,2], title="Cumulative Thrust", ylabel="[N]"; defaults...)
    add_box_plot_entries(ax, solution_set, "cum_thrust"; colors=colors, saturate_zero=true)
    push!(axes, ax)

    ax = Axis(f[1,3], title="Average Thrust", ylabel="[N]"; defaults...)
    add_box_plot_entries(ax, solution_set, "avg_thrust"; colors=colors, saturate_zero=true)
    push!(axes, ax)

    ax = Axis(f[1,4], title="Final Site Radius", ylabel="[m]"; defaults...)
    add_box_plot_entries(ax, solution_set, "final_radius_truth"; colors=colors, saturate_zero=true)
    push!(axes, ax)
    
    # Customize ticks
    labels = collect(keys(solution_set))
    n_mc = length(solution_set[labels[1]])
    for ax in axes
        ax.xticks = 1:length(labels)
        ax.xtickformat = k -> labels[Int.(k)]
        hidedecorations!(ax, label=false, ticklabels=false, ticks=false, minorticks=false)
    end

    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    end
end