using CairoMakie
using Colors
using LinearAlgebra
include("../../utils/plot_utils.jl")
include("plot_defaults.jl")

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
