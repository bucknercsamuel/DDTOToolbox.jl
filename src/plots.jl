#= DDTO for Landing -- Plotting Functions.

Author: Samuel Buckner (UW-ACL)
=#

# ..:: Plotting parameters ::..

# Setup
PyPlot.svg(true)
fig_path = "figures"

# Plot styling dictionaries
style_ct = Dict(:color=>"black",:linewidth=>1.5)
style_dt = Dict(:color=>"limegreen",:linestyle=>"none",:marker=>"o",:markersize=>5,:markeredgecolor=>"none")
style_relax = Dict(:color=>"red",:linestyle=>"none",:marker=>"o",:markersize=>3,:markeredgecolor=>"black")
style_hover = Dict(:color=>"orange",:linewidth=>1.5)
style_ground = Dict(:facecolor=>"gray",:edgecolor=>"black",:alpha=>0.7)
style_constraint = Dict(:color=>"red",:linestyle=>"--",:linewidth=>2)
style_constraint_fill = Dict(:edgecolor=>"none",:facecolor=>"black",:alpha=>0.1)

# Indexing offset based on interpolation method, 1 for ZOH
ctrl_offset = 1

# 3D plot camera viewing angle parameters
traj3d_elev = 90
traj3d_azim = 90

# Set trunk/target colors
# color_trunk = "black"
# colormap_targs = (n) -> collect(range(colorant"blue", stop=colorant"limegreen", length=n))

# Target colors
targ_colors = [
    ["blue", "cyan"],
    ["maroon", "red"],
    ["green", "limegreen"],
    ["purple", "magenta"],
]


# ..:: Convenience Functions ::..

function modify_styling_dict(styling_dict::Dict, key::String, value::Any)::Dict
    # Modify a parameter of an input styling dict
    styling_dict_mod = copy(styling_dict)
    styling_dict_mod[Symbol(key)] = value
    return styling_dict_mod
end

function set_fonts()::Nothing
    # Set the figure fonts.
    fig_smaller_sz = 13
    fig_small_sz = 14
    fig_med_sz = 15
    fig_big_sz = 17
    plt.rc("text", usetex=true)
    plt.rc("font", size=fig_small_sz, family="serif")
    plt.rc("axes", titlesize=fig_small_sz)
    plt.rc("axes", labelsize=fig_small_sz)
    plt.rc("xtick", labelsize=fig_smaller_sz)
    plt.rc("ytick", labelsize=fig_smaller_sz)
    plt.rc("grid", alpha=0.25)
    plt.rc("legend", fontsize=fig_smaller_sz)
    plt.rc("figure", titlesize=fig_big_sz)
    plt.rc("figure", dpi=300) 
    plt.rc("figure", figsize = [8,6])
    plt.rc("animation", html="html5")
    return nothing
end

function set_axes_equal(ax)
    #= Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    =#

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[2] - x_limits[1])
    x_middle = mean(x_limits)
    y_range = abs(y_limits[2] - y_limits[1])
    y_middle = mean(y_limits)
    z_range = abs(z_limits[2] - z_limits[1])
    z_middle = mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*maximum([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
end

# ..:: Plotting Functions ::..

function plot_parametric_optimal_trajectories(lander::Lander, sols_optimal::Array{Solution})

    # Create figure
    fun_name = nameof(var"#self#")
    fig = plt.figure(facecolor="white", figsize=[8,8])
    plt.clf()

    # Create and format subplot
    fig = plt.figure(facecolor="white")
    ax = plt.gca()
    plt.axis("equal")
    plt.grid(true)

    # Trajectory plots
    for j = 1:lander.n_targs
        dark_color = targ_colors[j][1]
        light_color = targ_colors[j][2]
        ax.plot(sols_optimal[j].r[1,:], sols_optimal[j].r[2,:], color=dark_color, markersize=5, marker="o", markerfacecolor=light_color, markeredgecolor=dark_color, label="Target "*string(j))
    end

    # Obstacles
    flag_labeled = false
    for o = 1:lander.n_obstacles
        if !flag_labeled
            label = "Obstacles"
            flag_labeled = true
        else
            label = ""
        end
        patch = plt.Circle((lander.p_obstacles[1,o], lander.p_obstacles[2,o]), lander.R_obstacles[o], facecolor="lightcoral", edgecolor="indianred", label=label)
        ax.add_patch(patch)
        plt.text(lander.p_obstacles[1,o], lander.p_obstacles[2,o], string(o), color="black", fontsize=12, ha="center", va="center", fontweight="extra bold")
    end

    # Set cage boundaries and plot limits
    fill_color = "gray"
    bound_color = "darkgray"
    alpha = 0.35
    xlims = [-100, 100]
    ylims = [-100, 100]
    ax.fill_between([xlims[1], lander.x_arena_lims[1]], 0, 1, facecolor=fill_color, edgecolor="none", alpha=alpha, transform=ax.get_xaxis_transform(), label="Cage bounds")
    ax.fill_between([lander.x_arena_lims[2], xlims[2]], 0, 1, facecolor=fill_color, edgecolor="none", alpha=alpha, transform=ax.get_xaxis_transform())
    ax.fill_between([lander.x_arena_lims[1], lander.x_arena_lims[2]], ylims[1], lander.y_arena_lims[1], facecolor=fill_color, edgecolor="none", alpha=alpha)
    ax.fill_between([lander.x_arena_lims[1], lander.x_arena_lims[2]], lander.y_arena_lims[2], ylims[2], facecolor=fill_color, edgecolor="none", alpha=alpha)
    ax.plot([lander.x_arena_lims[1], lander.x_arena_lims[1]], [lander.y_arena_lims[1], lander.y_arena_lims[2]], color=bound_color)
    ax.plot([lander.x_arena_lims[2], lander.x_arena_lims[2]], [lander.y_arena_lims[1], lander.y_arena_lims[2]], color=bound_color)
    ax.plot([lander.x_arena_lims[1], lander.x_arena_lims[2]], [lander.y_arena_lims[1], lander.y_arena_lims[1]], color=bound_color)
    ax.plot([lander.x_arena_lims[1], lander.x_arena_lims[2]], [lander.y_arena_lims[2], lander.y_arena_lims[2]], color=bound_color)

    pad = 1
    ax.set_xlim((lander.x_arena_lims[1] - pad, lander.x_arena_lims[2] + pad))
    ax.set_ylim((lander.y_arena_lims[1] - pad, lander.y_arena_lims[2] + pad))

    # Extra formatting
    plt.xlabel("X [m]")
    plt.ylabel("Y [m]")
    plt.legend(loc="upper right")

    # Save and show figure
    # fig.savefig("$fig_path/$fun_name.pdf", bbox_inches="tight")
    ;

end

function plot_parametric_ddto_trajectories(lander::Lander, sols_ddto::Array{BranchSolution})

    # Create figure
    fun_name = nameof(var"#self#")
    fig = plt.figure(facecolor="white", figsize=[8,8])
    plt.clf()

    # Create and format subplot
    fig = plt.figure(facecolor="white")
    ax = plt.gca()
    plt.axis("equal")
    plt.grid(true)

    # Trajectory plots
    for j = 1:lander.n_targs
        dark_color = targ_colors[j][1]
        light_color = targ_colors[j][2]
        ax.plot(sols_ddto[j].sol.r[1,:], sols_ddto[j].sol.r[2,:], color=dark_color, markersize=5, marker="o", markerfacecolor=light_color, markeredgecolor=dark_color, label="Target "*string(j))
    end

    # Obstacles
    flag_labeled = false
    for o = 1:lander.n_obstacles
        if !flag_labeled
            label = "Obstacles"
            flag_labeled = true
        else
            label = ""
        end
        patch = plt.Circle((lander.p_obstacles[1,o], lander.p_obstacles[2,o]), lander.R_obstacles[o], facecolor="lightcoral", edgecolor="indianred", label=label)
        ax.add_patch(patch)
        plt.text(lander.p_obstacles[1,o], lander.p_obstacles[2,o], string(o), color="black", fontsize=12, ha="center", va="center", fontweight="extra bold")
    end

    # Set cage boundaries and plot limits
    fill_color = "gray"
    bound_color = "darkgray"
    alpha = 0.35
    pad = 1
    xlims = [-100, 100]
    ylims = [-100, 100]
    ax.fill_between([xlims[1], lander.x_arena_lims[1]], 0, 1, facecolor=fill_color, edgecolor="none", alpha=alpha, transform=ax.get_xaxis_transform(), label="Cage bounds")
    ax.fill_between([lander.x_arena_lims[2], xlims[2]], 0, 1, facecolor=fill_color, edgecolor="none", alpha=alpha, transform=ax.get_xaxis_transform())
    ax.fill_between([lander.x_arena_lims[1], lander.x_arena_lims[2]], ylims[1], lander.y_arena_lims[1], facecolor=fill_color, edgecolor="none", alpha=alpha)
    ax.fill_between([lander.x_arena_lims[1], lander.x_arena_lims[2]], lander.y_arena_lims[2], ylims[2], facecolor=fill_color, edgecolor="none", alpha=alpha)
    ax.plot([lander.x_arena_lims[1], lander.x_arena_lims[1]], [lander.y_arena_lims[1], lander.y_arena_lims[2]], color=bound_color)
    ax.plot([lander.x_arena_lims[2], lander.x_arena_lims[2]], [lander.y_arena_lims[1], lander.y_arena_lims[2]], color=bound_color)
    ax.plot([lander.x_arena_lims[1], lander.x_arena_lims[2]], [lander.y_arena_lims[1], lander.y_arena_lims[1]], color=bound_color)
    ax.plot([lander.x_arena_lims[1], lander.x_arena_lims[2]], [lander.y_arena_lims[2], lander.y_arena_lims[2]], color=bound_color)
    ax.set_xlim((lander.x_arena_lims[1] - pad, lander.x_arena_lims[2] + pad))
    ax.set_ylim((lander.y_arena_lims[1] - pad, lander.y_arena_lims[2] + pad))


    # Extra formatting
    plt.xlabel("X [m]")
    plt.ylabel("Y [m]")
    plt.legend(loc="upper right")

    # Save and show figure
    # fig.savefig("$fig_path/$fun_name.pdf", bbox_inches="tight")
    ;

end

function plot_states(lander::Lander, sols_ddto::Array{BranchSolution})

    labels = ["pos-X [m]","pos-Y [m]","pos-Z [m]","vel-X [m/s]","vel-Y [m/s]","vel-Z [m/s]","acc-X [m/s2]","acc-Y [m/s2]","acc-Z [m/s2]"]
    data = Vector{Matrix}()
    for j in 1:lander.n_targs
        data_targ = vcat(
            sols_ddto[j].sol.r,
            sols_ddto[j].sol.v,
            hcat(sols_ddto[j].sol.T / lander.mass, [0;0;0]),
        )
        push!(data, data_targ)
    end

    # Create figure
    fun_name = nameof(var"#self#")
    fig, axs = plt.subplots(3, 3, facecolor="white", constrained_layout=true, figsize=[8,6])

    for (ind,ax) in enumerate(axs)
        for j in 1:lander.n_targs
            dark_color = targ_colors[j][1]
            light_color = targ_colors[j][2]
            ax.plot(sols_ddto[j].sol.t, data[j][ind,:], color=dark_color, markersize=5, marker="o", markerfacecolor=light_color, markeredgecolor=dark_color)
        end
        ax.grid(true)
        ax.set_xlabel("Time [s]")
        ax.set_ylabel(labels[ind])
    end
    ;
end