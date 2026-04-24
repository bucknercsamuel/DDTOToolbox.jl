using CairoMakie
using Colors
using LinearAlgebra
include("../../utils/plot_utils.jl")
include("plot_defaults.jl")
include("plot_drone.jl")

function plot_paper_demo_traj_history(
        run_data,
        map_data;
        interactive=true,
        azel=(pi/4,pi/4),
        save_path=nothing
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
    sim_time = run_data["sim_time"]
    sim_state = run_data["sim_state"]
    # sim_control = run_data["sim_control"]
    guid_ddto_trajs = run_data["guid_ddto_trajs_sims"]
    guid_defer_vecs = run_data["guid_defer_vecs"]
    guid_update_times = run_data["guid_update_times"]
    n_ddto_sols = length(guid_update_times)

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

    # Plot DDTO trajectory bundles, color-coded by temporal progress along the
    # trajectory (guid_update_times[k] mapped onto [sim_time[1], sim_time[end]]
    # via the viridis colormap).
    sim_t0, sim_tf = sim_time[1], sim_time[end]
    Δsim = sim_tf - sim_t0
    bundle_cmap = cgrad(:viridis)
    is_matrix_layout = ndims(guid_ddto_trajs) == 2
    for k = 1:n_ddto_sols
        frac = Δsim > 0 ? clamp((guid_update_times[k] - sim_t0) / Δsim, 0.0, 1.0) : 0.0
        color_k = bundle_cmap[frac]

        # Number of targets per bundle is data-driven (no reliance on params)
        n_targs = is_matrix_layout ? size(guid_ddto_trajs, 2) : length(guid_ddto_trajs[k])

        for j = 1:n_targs
            traj = collect(is_matrix_layout ? guid_ddto_trajs[k, j] : guid_ddto_trajs[k][j])
            (size(traj, 1) >= 3 && size(traj, 2) >= 1) || continue
            lines!(ax,
                traj[1,:], traj[2,:], traj[3,:];
                style3D_ct_ddto..., :alpha=>0.9, :color=>color_k)
        end
    end

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

    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    else
        if !isnothing(save_path)
            CairoMakie.save(save_path, f, px_per_unit=4)
        else
            CairoMakie.save(joinpath(fig_path, "paper_demo_traj"*fig_ext), f, px_per_unit=4)
        end
    end
end
