# using CairoMakie
using GLMakie
using Colors
using InvertedIndices
using GeometryBasics
include("../utils/plot_utils.jl")

# Generic styling
style2D_dt = Dict(:color=>:gray, :marker=>:circle, :markersize=>15, :strokecolor=>:black, :strokewidth=>3) 
style2D_ct = Dict(:color=>:black, :linewidth=>3)
style3D_dt = Dict(:color=>:gray, :marker=>:circle, :markersize=>15, :strokecolor=>:black, :strokewidth=>3)
style3D_ct = Dict(:color=>:black, :linewidth=>3)
style3D_ground_base = Dict(:color=>bright_color(:orange), :transparency=>false, :alpha=>1)
style3D_ground_base_frame = Dict(:color=>bright_color(:saddlebrown))

# Themes
theme2d = merge(theme_minimal(), theme_latexfonts())
theme3d = theme_latexfonts()
fontsize = 20

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

    # DDTO Color conditions
    base_colors = ["red", "gold", "blue", "green", "purple", "pink", "brown", "cyan", "orange", "yellow"]
    idx_colors = k -> mod(k-1,length(base_colors))+1
    color_map_bundles = []
    for k = 1:n_ddto_sols
        color = base_colors[idx_colors(k)]
        color1 = parse(Colorant, color*"1")
        color2 = parse(Colorant, color*"4")
        if length(ddto_bundles_sol[k].targs) > 1
            cmap = range(color1, color2, length=length(ddto_bundles_sol[k].targs))
        else
            cmap = color1
        end
        append!(color_map_bundles, [cmap])
    end
    
    # Plot DDTO trajectories
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