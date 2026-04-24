using CairoMakie
using Colors
using LinearAlgebra
include("../../utils/plot_utils.jl")
include("plot_defaults.jl")

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
