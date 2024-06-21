# using CairoMakie
using GLMakie
using Colors
using InvertedIndices
include("../utils/plot_utils.jl")

# Generic styling
style2D_dt = Dict(:color=>:gray, :marker=>:circle, :markersize=>15, :strokecolor=>:black, :strokewidth=>3) 
style2D_ct = Dict(:color=>:black, :linewidth=>3)
style3D_dt = Dict(:color=>:gray, :marker=>:circle, :markersize=>15, :strokecolor=>:black, :strokewidth=>3)
style3D_ct = Dict(:color=>:black, :linewidth=>3)

# Themes
theme2d = merge(theme_minimal(), theme_latexfonts())
theme3d = theme_latexfonts()
fontsize = 20

function build_plots(results; interactive=true)
    screens = []
    with_theme(theme3d; fontsize=fontsize) do
        push!(screens, plot_3d_trajs(results))
    end

    if interactive
        println("\nPress any key when finished using plots...")
        readline() # Wait for user to finish plotting
        [GLMakie.destroy!(screen) for screen in screens]
    end
end

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

    # DDTO Color conditions
    base_colors = ["red", "gold", "blue", "green", "purple", "pink", "brown", "cyan", "orange", "yellow"]
    color_map_bundles = []
    for k = 1:n_ddto_sols
        color = base_colors[k]
        color1 = parse(Colorant, color*"1")
        color2 = parse(Colorant, color*"4")
        cmap = range(color1, color2, length=length(ddto_bundles_sol[k].targs))
        append!(color_map_bundles, [cmap])
    end
    
    # Plot DDTO trajectories
    proj_idxs = [1,2,3]
    for k = 1:n_ddto_sols
        color_branch = j -> color_map_bundles[k][j]
        params = ddto_params[k]
        plot_bundle(ax,
            [[ddto_bundles_sol[k].targs[j].r[c,:] for j∈1:params.a.n_targs] for c∈proj_idxs],
            [[ddto_bundles_sim[k].targs[j].r[c,:] for j∈1:params.a.n_targs] for c∈proj_idxs],
            params,
            style3D_ct,
            style3D_dt;
            color_branch = color_branch,
            color_trunk = base_colors[k],
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