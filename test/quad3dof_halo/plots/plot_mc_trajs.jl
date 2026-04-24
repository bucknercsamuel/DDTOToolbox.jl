using CairoMakie
using Colors
using GeometryBasics
using LinearAlgebra
include("../../utils/plot_utils.jl")
include("plot_defaults.jl")

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
