using CairoMakie
# using GLMakie
using Colors
using InvertedIndices
using GeometryBasics
using LinearAlgebra
using Statistics
include("../utils/plot_utils.jl")

# Generic styling
markersize = 15
style2D_dt = Dict(:color=>:gray, :marker=>:circle, :markersize=>markersize, :strokecolor=>:black, :strokewidth=>3) 
style2D_ct = Dict(:color=>:black, :linewidth=>3)
style2D_ct_ddto = Dict(:color=>:black, :linewidth=>3)
style3D_dt = Dict(:color=>:gray, :marker=>:circle, :markersize=>markersize, :strokecolor=>:black, :strokewidth=>3)
style3D_ct = Dict(:color=>:black, :linewidth=>5, :overdraw=>true)
style3D_ct_ddto = Dict(:color=>:black, :linewidth=>5)
# style3D_ground_base = Dict(:color=>bright_color(:orange), :transparency=>false, :alpha=>1)
# style3D_ground_base_frame = Dict(:color=>bright_color(:saddlebrown))
style3D_ground_base = Dict(:color=>bright_color(:gray95), :transparency=>false, :alpha=>1)
style3D_ground_base_frame = Dict(:color=>bright_color(:gray90))

# Themes
theme2d = merge(theme_minimal(), theme_latexfonts())
theme3d = theme_latexfonts()
fontsize = 20

fig_path = "quad3dof_halo\\figures"
fig_ext = ".png"

CairoMakie.activate!(type="svg")
fig_ext = ".svg"

function generate_custom_colors(max_targs)
    # target_colors = range(colorant"magenta", stop=colorant"cyan", length=max_targs)
    target_colors = range(HSV(45,1,1), stop=HSV(-360,1,1), length=max_targs)

    return target_colors
end

b_ref = rand(3)
b_ref = b_ref / norm(b_ref)
function plot_drone(
        ax,
        position,
        thrust_direction;
        scale=1
    )
    # Sizing parameters for vehicle
    b_x = thrust_direction
    b_y = cross(b_x, b_ref)
    b_y = b_y / norm(b_y)
    b_z = cross(b_x, b_y)
    body_radius = 0.15
    body_height = 0.05
    arm_length = 0.4
    arm_radius = 0.02
    prop_radius = 0.2
    prop_height = 0.05
    cmap_frame(N) = colormap("Grays",N)
    cmap_arm(N) = colormap("Grays", N)
    cmap_prop(N) = colormap("Reds", N)
    style_prop = Dict(:alpha=>.3)

    # Create a drone object with cylinders
    # inputs: ax, vertex, pointing_direction, radius; length, cmap
    draw_cylinder_3d(ax, position, b_x, scale*body_radius; length=scale*body_height, cmap=cmap_frame)
    arm_dirs = [b_y, b_z, -b_y, -b_z]
    for arm_dir in arm_dirs
        draw_cylinder_3d(ax, position, arm_dir, scale*arm_radius; length=scale*(arm_length-prop_radius), cmap=cmap_arm)
        draw_cylinder_3d(ax, position + scale*arm_length*arm_dir + scale*arm_radius*b_x, b_x, scale*prop_radius; length=scale*prop_height, cmap=cmap_prop, style=style_prop)
    end
end

function plot_3d_trajs(
        results;
        interactive=true,
        azel=(pi/4,pi/4),
        f=nothing,
        ax_idx=[1,1],
        show_drones=false,
        drone_frame_count=10,
        save_fig=true
    )
    # Setup
    if isnothing(f)
        f = Figure(size=(800,800))
    end
    # Create axis
    ax = Axis3(
        f[ax_idx...],
        xlabel = "East [m]",
        ylabel = "North [m]",
        zlabel = "Up [m]",
        aspect = :equal,
        azimuth = azel[1],
        elevation = azel[2],
        xgridvisible = false,
        ygridvisible = false,
        zgridvisible = false,
        # xypanelcolor = :gray95,
        # yzpanelcolor = :gray95,
        # xzpanelcolor = :gray95
        )

    # Results parsing
    ddto_params = results["guid_update_ddto_params"]
    ddto_bundles_sol = results["guid_update_ddto_bundles"]
    ddto_bundles_sim = results["guid_update_ddto_bundles_sims"]
    sim_solution = results["sim_state"]
    sim_control = results["sim_control"]
    n_ddto_sols = length(ddto_bundles_sol)

    # Zoom onto right area
    rmin = [min(vcat([[params.a.z0[k,:]; params.a.zf_targs[k,:]] for params in ddto_params]...)...) for k∈1:3]
    rmax = [max(vcat([[params.a.z0[k,:]; params.a.zf_targs[k,:]] for params in ddto_params]...)...) for k∈1:3]
    rmin = [min(vcat(rmin[k], results["targpool_positions"][k,:])...) for k∈1:3]
    rmax = [max(vcat(rmax[k], results["targpool_positions"][k,:])...) for k∈1:3]
    xLims,yLims,zLims = get_equal_3d_lims(rmin, rmax)
    xlims!(ax, xLims...)
    ylims!(ax, yLims...)
    zlims!(ax, zLims...)

    # Plot the ground
    ϵ = (zLims[2] - zLims[1]) * 0.001 # Epsilon in altitude for objects that are stacked
    fac_exp = 0.01
    ΔL(L) = L[2] - L[1]
    box_adjust = [ΔL(xLims), ΔL(yLims), 0] * fac_exp
    box_lower = [xLims[1], yLims[1], zLims[1]] - box_adjust
    box_upper = [ΔL(xLims), ΔL(yLims), -zLims[1]-ϵ] + 2*box_adjust
    groundBaseDef = Rect3f(box_lower, box_upper)
    groundBaseMesh = GeometryBasics.mesh(groundBaseDef)
    mesh!(groundBaseMesh; style3D_ground_base...)
    boxframe_3D(ax, box_lower, box_upper; style=style3D_ground_base_frame)
    
    # Plot circular patches on the ground for each target
    radius = 0.025 * ΔL(xLims)
    for k = 1:length(results["targpool_ID"])
        id = results["targpool_ID"][k]
        draw_circle_3d(ax, results["targpool_positions"][:,k], radius, pointing_direction=e_z, color=bright_color(target_colors[id], fraction=0.5))
    end

    # Plot the cylindrical obstacles from the ground to the initial altitude
    params = ddto_params[1]
    for k = 1:params.n_obstacles
        style = Dict(:transparency=>true, :alpha=>0.2)
        draw_cylinder_3d(ax, params.p_obstacles[:,k], e_z, params.R_obstacles[k], length=params.a.z0[3], cmap=N->colormap("Reds",N), style=style)
    end

    # Plot DDTO trajectories
    proj_idxs = [1,2,3]
    darken_frac_start = 0
    darken_frac_end = .3
    for k = 1:n_ddto_sols
        params = ddto_params[k]
        darken_frac = darken_frac_start+(darken_frac_end-darken_frac_start)*(k-1)/(n_ddto_sols-1)
        color_branch = j -> dark_color(target_colors[params.a.ID_targs[j]], fraction=darken_frac)
        plot_bundle(ax,
            [[ddto_bundles_sol[k].targs[j].r[c,:] for j∈1:params.a.n_targs] for c∈proj_idxs],
            [[ddto_bundles_sim[k].targs[j].r[c,:] for j∈1:params.a.n_targs] for c∈proj_idxs],
            params,
            style3D_ct_ddto,
            style3D_dt;
            color_branch = color_branch,
            show_sol_nodes = false,
            show_defer_nodes = false,
            show_ddto_split = false,
            alpha=0.5
        )
    end

    # Plot the trajectory that was actually taken
    positions = sim_solution[1:3,:]
    lines!(ax,
           positions[1,:], positions[2,:], positions[3,:];
           style3D_ct..., :alpha=>1, :color=>:black)

    # Plot the vehicle body timelapse
    if show_drones
        idx_frames = round.(Int, LinRange(1, size(sim_solution,2), drone_frame_count))
        for idx in idx_frames
            position = sim_solution[1:3,idx]
            thrust = sim_control[1:3,idx]
            thrust_direction = thrust / norm(thrust)
            plot_drone(ax, position, thrust_direction; scale=15)
        end
    end

    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    else
        if save_fig
            save(joinpath(fig_path, "3d_trajs"*fig_ext), f)
        end
    end
end

function plot_2d_trajs_XY(
        results;
        interactive=true,
        f=nothing,
        ax_idx=[1,1],
        save_fig=true
    )
    # Setup
    if isnothing(f)
        f = Figure(size=(800,800))
    end
    # Create axis
    ax = Axis(
        f[ax_idx...],
        xlabel = "East [m]",
        ylabel = "North [m]",
        aspect = AxisAspect(1),
        # xgridvisible = false,
        # ygridvisible = false,
        )

    # Results parsing
    ddto_params = results["guid_update_ddto_params"]
    ddto_bundles_sol = results["guid_update_ddto_bundles"]
    ddto_bundles_sim = results["guid_update_ddto_bundles_sims"]
    sim_solution = results["sim_state"]
    n_ddto_sols = length(ddto_bundles_sol)

    # Zoom onto right area
    rmin = [min(vcat([[params.a.z0[k,:]; params.a.zf_targs[k,:]] for params in ddto_params]...)...) for k∈1:3]
    rmax = [max(vcat([[params.a.z0[k,:]; params.a.zf_targs[k,:]] for params in ddto_params]...)...) for k∈1:3]
    rmin = [min(vcat(rmin[k], results["targpool_positions"][k,:])...) for k∈1:3]
    rmax = [max(vcat(rmax[k], results["targpool_positions"][k,:])...) for k∈1:3]
    rmin[3] = 0
    rmax[3] = 0
    xLims,yLims,_ = get_equal_3d_lims(rmin, rmax)
    xlims!(ax, xLims...)
    ylims!(ax, yLims...)

    # Plot circular patches on the ground for each target
    ΔL(L) = L[2] - L[1]
    radius = 0.025 * ΔL(xLims)
    for k = 1:length(results["targpool_ID"])
        id = results["targpool_ID"][k]
        draw_circle_3d(ax, results["targpool_positions"][:,k], radius, pointing_direction=e_z, color=bright_color(target_colors[id], fraction=0.5))
    end

    # Plot the cylindrical obstacles from the ground to the initial altitude (projected)
    params = ddto_params[1]
    for k = 1:params.n_obstacles
        style = Dict(:transparency=>true, :alpha=>0.4)
        draw_circle_3d(ax, params.p_obstacles[:,k], params.R_obstacles[k], pointing_direction=e_z, color=:red, style=style)
    end

    # Plot DDTO trajectories
    proj_idxs = [1,2]
    darken_frac_start = 0
    darken_frac_end = .3
    for k = 1:n_ddto_sols
        params = ddto_params[k]
        darken_frac = darken_frac_start+(darken_frac_end-darken_frac_start)*(k-1)/(n_ddto_sols-1)
        color_branch = j -> dark_color(target_colors[params.a.ID_targs[j]], fraction=darken_frac)
        plot_bundle(ax,
            [[ddto_bundles_sol[k].targs[j].r[c,:] for j∈1:params.a.n_targs] for c∈proj_idxs],
            [[ddto_bundles_sim[k].targs[j].r[c,:] for j∈1:params.a.n_targs] for c∈proj_idxs],
            params,
            style3D_ct_ddto,
            style3D_dt;
            color_branch = color_branch,
            show_sol_nodes = false,
            show_defer_nodes = false,
            show_ddto_split = false,
            alpha=0.5
        )
    end

    # Plot the trajectory that was actually taken
    positions = sim_solution[1:2,:]
    lines!(ax,
        positions[1,:], positions[2,:];
        style3D_ct..., :alpha=>1, :color=>:black)


    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    else
        if save_fig
            save(joinpath(fig_path, "2d_trajs_XY"*fig_ext), f)
        end
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
    axis_defaults_2d = Dict(
        :xautolimitmargin=>(0,0), 
        :topspinevisible=>true, 
        :rightspinevisible=>true,
        :xgridvisible=>true,
        :ygridvisible=>true,
        )
    color_branch = j -> target_colors[params.a.ID_targs[j]]

    # Positions axes
    labels = ["Pos-East [m]", "Pos-North [m]", "Pos-Up [m]"]
    for (k,c) in enumerate(proj_idxs)
        ax = Axis(f[k,1], ylabel=labels[k]; axis_defaults_2d...)
        for k = 1:n_ddto_sols
            params = ddto_params[k]
            plot_bundle(ax,
                [[ddto_bundles_sol[k].targs[j].t .+ update_times[k] for j∈1:params.a.n_targs], [ddto_bundles_sol[k].targs[j].r[c,:] for j∈1:params.a.n_targs]],
                [[ddto_bundles_sim[k].targs[j].t .+ update_times[k] for j∈1:params.a.n_targs], [ddto_bundles_sim[k].targs[j].r[c,:] for j∈1:params.a.n_targs]],
                params,
                style2D_ct_ddto,
                style2D_dt;
                color_branch = color_branch,
                show_sol_nodes = false,
                show_defer_nodes = false,
                show_ddto_split = false,
                alpha=0.5
            )
        end
        lines!(ax, sim_time, sim_state[c,:];
            style2D_ct..., :alpha=>1, :color=>:black)
    end
    # ax_labels = ax

    # Velocities axes
    labels = ["Vel-East [m]", "Vel-North [m]", "Vel-Up [m]"]
    for (k,c) in enumerate(proj_idxs)
        ax = Axis(f[k,2], ylabel=labels[k]; axis_defaults_2d...)
        for k = 1:n_ddto_sols
            params = ddto_params[k]
            plot_bundle(ax,
                [[ddto_bundles_sol[k].targs[j].t .+ update_times[k] for j∈1:params.a.n_targs], [ddto_bundles_sol[k].targs[j].v[c,:] for j∈1:params.a.n_targs]],
                [[ddto_bundles_sim[k].targs[j].t .+ update_times[k] for j∈1:params.a.n_targs], [ddto_bundles_sim[k].targs[j].v[c,:] for j∈1:params.a.n_targs]],
                params,
                style2D_ct_ddto,
                style2D_dt;
                color_branch = color_branch,
                show_sol_nodes = false,
                show_defer_nodes = false,
                show_ddto_split = false,
                alpha=0.5
            )
        end
        lines!(ax, sim_time, sim_state[c+3,:];
            style2D_ct..., :alpha=>1, :color=>:black)
    end

    # Thrust norm axis
    ax = Axis(f[4,1:2], ylabel="Thrust Norm [N]"; axis_defaults_2d...)
    if integrated_sim
        data = sim_state[7,:]
    else
        data = [norm(results["sim_control"][1:3,l]) for l=1:length(sim_time)]
    end
    for k = 1:n_ddto_sols
        params = ddto_params[k]
        plot_bundle(ax,
            [[ddto_bundles_sol[k].targs[j].t .+ update_times[k] for j∈1:params.a.n_targs], [[norm(ddto_bundles_sol[k].targs[j].T[:,l]) for l∈1:length(ddto_bundles_sol[k].targs[j].t)] for j∈1:params.a.n_targs]],
            [[ddto_bundles_sim[k].targs[j].t .+ update_times[k] for j∈1:params.a.n_targs], [[norm(ddto_bundles_sim[k].targs[j].T[:,l]) for l∈1:length(ddto_bundles_sim[k].targs[j].t)] for j∈1:params.a.n_targs]],
            params,
            style2D_ct_ddto,
            style2D_dt;
            color_branch = color_branch,
            show_sol_nodes = false,
            show_defer_nodes = false,
            show_ddto_split = false,
            alpha=0.5
        )
    end
    lines!(ax, sim_time, data;
        style2D_ct..., :alpha=>1, :color=>:black)

    # ax = Legend(f[1:4,3], ax_labels)

    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    else
        save(joinpath(fig_path, "states"*fig_ext), f)
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
    else
        save(joinpath(fig_path, "mc_trajs"*fig_ext), f)
    end
end

function plot_mc_statistics(solution_set; interactive=true, groupings::Vector = [], label="")

    # Build figure
    f = Figure(size=(500,400))
    defaults = Dict(
        :topspinevisible=>true, 
        :rightspinevisible=>true,
        # :xgridvisible=>true,
        # :ygridvisible=>true,
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

    function add_plot_entries(ax, solution_set, data_name; colors=[], saturate_zero=false, groupings=[], width_factor=.4, show_violin=true, outlier_threshold=nothing)
        if length(groupings) == 0
            groupings = [(i,) for i in 1:length(solution_set)]
        end
        box_pos = 1.
        box_poses = []
        key_order = ["GraphSCvx", "Gr-1", "Gr-∞"]
        mean_graphSCvx = nothing
        for (iter,key) in enumerate(key_order)
            # Process data
            value = solution_set[key]
            append!(box_poses, box_pos)
            idx_feas = findall(τ->τ==1, [value[k]["error_code"] for k∈1:length(value)])
            data = [value[k][string(data_name)] for k∈idx_feas]            
            quant_data(p) = quantile(data, p)
            Q1,median,Q3 = map(quant_data, [.25,.5,.75])
            IQR = Q3-Q1
            outliers = findall(x->(x<Q1-1.5*IQR).|(x>Q3+1.5*IQR), data)
            if !isnothing(outlier_threshold)
                extreme_outliers = findall(x->x>outlier_threshold, data)
                if length(extreme_outliers) > 0
                    data = data[setdiff(1:length(data), extreme_outliers)]
                    outliers = findall(x->(x<Q1-1.5*IQR).|(x>Q3+1.5*IQR), data)
                end
            end

            # Add violin plot (with inliers only)
            data_inliers = data[setdiff(1:length(data), outliers)]
            if show_violin
                violin!(ax, fill(iter,length(data_inliers)), data_inliers; color=colors[iter][2], scale=:width, width=width_factor*length(solution_set))
            end
            # Add box plot overlay (with outliers)
            add_box_plot_entry(ax, box_pos, Q1, median, Q3, data[outliers]; width=width_factor*length(solution_set), color_dark=colors[iter][1], color_light=colors[iter][2], saturate_zero=saturate_zero, alpha_fill=0.5)
            grouping_idx = findfirst(g->iter in g, groupings)
            group_idx = findfirst(g-> iter in g, groupings[grouping_idx])
            if group_idx < length(groupings[grouping_idx])
                box_pos += 1.05*width_factor*length(solution_set)
            else
                box_pos += 1.
            end
            # Cache GraphSCvx data for later comparison
            if key == "GraphSCvx"
                mean_graphSCvx = mean(data)
            end
            percent_increase = (mean(data) - mean_graphSCvx) / mean_graphSCvx * 100
            # Print out stats:
            round2(x) = round(x, digits=2)
            println("$(data_name)::$(key)")
            println("Q1: $(round2(Q1)), Median: $(round2(median)), Q3: $(round2(Q3)), IQR: $(round2(IQR))")
            println("Num outliers: $(length(outliers))")
            println("Mean data: $(round2(mean(data))) ± $(round2(std(data))) (%$(round2(percent_increase)))")
            println("Bounds data: [$(round2(minimum(data))), $(round2(maximum(data)))]")
            println()
        end
        # Customize ticks
        labels = key_order
        label_pointers = Dict([(box_poses[k], labels[k]) for k in 1:length(labels)])
        ax.xticks = box_poses
        ax.xtickformat = values -> [label_pointers[value] for value in values]
        # hidedecorations!(ax, label=false, ticklabels=false, ticks=false, minorticks=false)
    end

    # Instantiation
    types = ""
    axes = []
    idx = 1

    ax = Axis(f[1,idx], ylabel="Cumulative thrust [N]"; defaults...)
    add_plot_entries(ax, solution_set, "cum_thrust"; colors=colors, saturate_zero=true, groupings=groupings, outlier_threshold=500)
    push!(axes, ax)
    idx += 1
    types *= "cum_thrust_"

    # ax = Axis(f[1,idx], title="Cumulative Jerk", ylabel="[N]"; defaults...)
    # add_plot_entries(ax, solution_set, "cum_jerk"; colors=colors, saturate_zero=true, groupings=groupings, outlier_threshold=1000)
    # push!(axes, ax)
    # idx += 1
    # types *= "cum_jerk_"

    # ax = Axis(f[1,idx], ylabel="Average solve time [s]"; defaults...)
    # add_plot_entries(ax, solution_set, "avg_solve_time"; colors=colors, saturate_zero=true, groupings=groupings)
    # push!(axes, ax)
    # idx += 1
    # types *= "avg_solve_time_"

    # ax = Axis(f[1,idx], title="ATE", ylabel="[m]"; defaults...)
    # add_plot_entries(ax, solution_set, "ATE"; colors=colors, saturate_zero=true, groupings=groupings)
    # push!(axes, ax)
    # idx += 1
    # types *= "ATE_"

    # ax = Axis(f[1,idx], title="Num Recomputations", ylabel="[-]"; defaults...)
    # add_plot_entries(ax, solution_set, "num_recomputations"; colors=colors, saturate_zero=true, groupings=groupings)
    # push!(axes, ax)
    # idx += 1
    # types *= "num_recomputations_"

    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    else
        CairoMakie.save(joinpath(fig_path, "mc_$(types)$(label)"*fig_ext), f)
    end
end

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

function paper_plot_trajallocation(
    params,
    results;
    interactive = true,
    azel=(pi/4,pi/4)
    )

    with_theme(theme3d; fontsize=fontsize) do
        f = Figure(size=(1500,700))
        plot_3d_trajs(results; interactive=false, f=f, ax_idx=[1,1], azel=azel)
        plot_target_allocations(params, results; interactive=false, f=f, ax_idx=[1,2])

        if interactive
            screen = GLMakie.Screen()
            display(screen, f)
            return screen
        else
            save(joinpath(fig_path, "trajallocation"*fig_ext), f)
        end
    end
end

function paper_plot_traj_with_top_projection(
    params,
    results;
    interactive = true,
    azel=(pi/4,pi/6)
    )

    with_theme(theme3d; fontsize=fontsize) do
        f = Figure(size=(1500,700))
        plot_3d_trajs(results; interactive=false, f=f, ax_idx=[1,1], azel=azel)
        plot_2d_trajs_XY(results; interactive=false, f=f, ax_idx=[1,2])

        if interactive
            screen = GLMakie.Screen()
            display(screen, f)
            return screen
        else
            save(joinpath(fig_path, "traj_with_top_projection"*fig_ext), f)
        end
    end
end

function plot_paper_demo_traj_history(
        run_data,
        map_data;
        interactive=true,
        azel=(pi/4,pi/4)
    )
    # Setup
    f = Figure(size=(800,800))
    ax = Axis3(
        f[1,1],
        xlabel = "",
        ylabel = "",
        zlabel = "",
        aspect = :equal,
        azimuth = azel[1],
        elevation = azel[2],
        xgridvisible = false,
        ygridvisible = false,
        zgridvisible = false,
        xticklabelsvisible = false, 
        yticklabelsvisible = false,
        zticklabelsvisible = false,
        xgridcolor = :transparent,
        ygridcolor = :transparent,
        zgridcolor = :transparent,
        xticksvisible = false,
        yticksvisible = false,
        zticksvisible = false,
        # xypanelcolor = :gray95,
        # yzpanelcolor = :gray95,
        # xzpanelcolor = :gray95
        )
    hidespines!(ax)

    # Results parsing
    sim_state = run_data["sim_state"]
    # sim_control = run_data["sim_control"]
    guid_ddto_trajs = run_data["guid_ddto_trajs_sims"]
    guid_defer_vecs = run_data["guid_defer_vecs"]
    n_ddto_sols = length(guid_ddto_trajs)

    # # Clean part of sim state where the vehicle resets to above its position (should not be plotted)
    # idx_reset = findfirst(dh->dh>0, diff(sim_state[3,:]))
    # sim_state = sim_state[:,1:idx_reset]

    # Results parsing
    rmin = [min(sim_state[k,:]...) for k∈1:3]
    rmax = [max(sim_state[k,:]...) for k∈1:3]
    xLims,yLims,zLims = get_equal_3d_lims(rmin, rmax)
    # xlims!(ax, xLims...)
    # ylims!(ax, yLims...)
    # zlims!(ax, zLims...)

    # Load terrain data (requires conversions between NED and ENU)
    xlims_terrain = yLims # limit swap to convert from ENU to NED
    ylims_terrain = xLims # limit swap to convert from ENU to NED
    xs = LinRange(xlims_terrain[1], xlims_terrain[2], Int(round(xlims_terrain[2]-xlims_terrain[1])))
    ys = LinRange(ylims_terrain[1], ylims_terrain[2], Int(round(ylims_terrain[2]-ylims_terrain[1])))
    center = [(xlims_terrain[1]+xlims_terrain[2])/2, (ylims_terrain[1]+ylims_terrain[2])/2]
    zs = [map_data["zlookup"][Int(round(x)),Int(round(y))] for x in xs, y in ys]
    xs,ys,zs = ys,xs,transpose(zs)
    
    # Plot the terrain
    surface!(ax,
        xs, ys, zs,
        colormap = range(parse(Colorant, "orangered4"), stop=parse(Colorant, "darkorange2"), length=100),
        alpha=0.8)

    # Plot the recorded sim state
    positions = sim_state[1:3,:]
    lines!(ax,
        positions[1,:], positions[2,:], positions[3,:];
        style3D_ct..., :alpha=>.5, :color=>:black)

    # Plot the vehicle body timelapse
    num_frames = 10
    idx_frames = round.(Int, LinRange(1, size(sim_state,2), num_frames))
    for idx in idx_frames
        position = sim_state[1:3,idx]
        thrust = run_data["sim_control"][1:3,idx]
        thrust_direction = thrust / norm(thrust)
        plot_drone(ax, position, thrust_direction; scale=15)
    end

    # # Plot DDTO trajectories
    # proj_idxs = [1,2,3]
    # darken_frac_start = 0
    # darken_frac_end = .3
    # for k = 1:n_ddto_sols
    #     params = ddto_params[k]
    #     darken_frac = darken_frac_start+(darken_frac_end-darken_frac_start)*(k-1)/(n_ddto_sols-1)
    #     color_branch = j -> dark_color(target_colors[params.a.ID_targs[j]], fraction=darken_frac)
    #     plot_bundle(ax,
    #         [[ddto_bundles_sol[k].targs[j].r[c,:] for j∈1:params.a.n_targs] for c∈proj_idxs],
    #         [[ddto_bundles_sim[k].targs[j].r[c,:] for j∈1:params.a.n_targs] for c∈proj_idxs],
    #         params,
    #         style3D_ct_ddto,
    #         style3D_dt;
    #         color_branch = color_branch,
    #         show_sol_nodes = false,
    #         show_defer_nodes = false,
    #         show_ddto_split = false,
    #         alpha=0.5
    #     )
    # end

    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    else
        save(joinpath(fig_path, "paper_demo_traj"*fig_ext), f)
    end
end

function paper_plot_greedy_compare(
    results_all;
    interactive = true,
    azel=(pi/4,pi/4)
    )

    with_theme(theme3d; fontsize=fontsize) do
        f = Figure(size=(1200,700))
        plot_3d_trajs(results_all[1]; interactive=false, f=f, ax_idx=[1,1], azel=azel, save_fig=false)
        plot_3d_trajs(results_all[2]; interactive=false, f=f, ax_idx=[1,2], azel=azel, save_fig=false)
        plot_3d_trajs(results_all[3]; interactive=false, f=f, ax_idx=[1,3], azel=azel, save_fig=false)
        plot_2d_trajs_XY(results_all[1]; interactive=false, f=f, ax_idx=[2,1], save_fig=false)
        plot_2d_trajs_XY(results_all[2]; interactive=false, f=f, ax_idx=[2,2], save_fig=false)
        plot_2d_trajs_XY(results_all[3]; interactive=false, f=f, ax_idx=[2,3], save_fig=false)

        if interactive
            screen = GLMakie.Screen()
            display(screen, f)
            return screen
        else
            save(joinpath(fig_path, "greedy_compare"*fig_ext), f)
        end
    end
end