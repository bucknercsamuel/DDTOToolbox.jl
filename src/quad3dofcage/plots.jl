#= DDTO for Landing -- Plotting Functions.

Author: Samuel Buckner (UW-ACL)
=#

# ..:: Plotting parameters ::..

# Setup
PyPlot.svg(true)
fig_path = "figures"

# Target colors
targ_colors = [
    ["blue", "cyan"],
    ["maroon", "red"],
    ["green", "limegreen"],
    ["purple", "magenta"],
]


# ..:: Convenience Functions ::..

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

# ..:: Plotting Functions ::..

function plot_parametric_trajectories(
    params::Params, 
    solutions::Array{ProcessedBranchSolution}, 
    simulations::Array{ProcessedBranchSolution};
    display_cage::Bool=false,
    display_obstacles::Bool=false, 
    fname::String="default_name")

    # Create figure
    fig = plt.figure(facecolor="white", figsize=[8,8])
    plt.clf()

    # Create and format subplot
    fig = plt.figure(facecolor="white")
    ax = plt.gca()
    plt.axis("equal")
    plt.grid(true)

    # Trajectory plots
    for j = 1:params.n_targs
        dark_color = targ_colors[j][1]
        light_color = targ_colors[j][2]
        ax.plot(simulations[j].sol.r[1,:], simulations[j].sol.r[2,:], color=dark_color)
        ax.plot(solutions[j].sol.r[1,:], solutions[j].sol.r[2,:], color="none", markersize=5, marker="o", markerfacecolor=light_color, markeredgecolor=dark_color, label="Target "*string(j))
    end

    # Obstacles
    if display_obstacles
        flag_labeled = false
        for o = 1:params.n_obstacles
            if !flag_labeled
                label = "Obstacles"
                flag_labeled = true
            else
                label = ""
            end
            patch = plt.Circle((params.p_obstacles[1,o], params.p_obstacles[2,o]), params.R_obstacles[o], facecolor="lightcoral", edgecolor="indianred", label=label)
            ax.add_patch(patch)
            plt.text(params.p_obstacles[1,o], params.p_obstacles[2,o], string(o), color="black", fontsize=12, ha="center", va="center", fontweight="extra bold")
        end
    end

    # Set cage boundaries and plot limits
    if display_cage
        fill_color = "gray"
        bound_color = "darkgray"
        alpha = 0.35
        pad = 1
        xlims = [-100, 100]
        ylims = [-100, 100]
        ax.fill_between([xlims[1], params.x_arena_lims[1]], 0, 1, facecolor=fill_color, edgecolor="none", alpha=alpha, transform=ax.get_xaxis_transform(), label="Cage bounds")
        ax.fill_between([params.x_arena_lims[2], xlims[2]], 0, 1, facecolor=fill_color, edgecolor="none", alpha=alpha, transform=ax.get_xaxis_transform())
        ax.fill_between([params.x_arena_lims[1], params.x_arena_lims[2]], ylims[1], params.y_arena_lims[1], facecolor=fill_color, edgecolor="none", alpha=alpha)
        ax.fill_between([params.x_arena_lims[1], params.x_arena_lims[2]], params.y_arena_lims[2], ylims[2], facecolor=fill_color, edgecolor="none", alpha=alpha)
        ax.plot([params.x_arena_lims[1], params.x_arena_lims[1]], [params.y_arena_lims[1], params.y_arena_lims[2]], color=bound_color)
        ax.plot([params.x_arena_lims[2], params.x_arena_lims[2]], [params.y_arena_lims[1], params.y_arena_lims[2]], color=bound_color)
        ax.plot([params.x_arena_lims[1], params.x_arena_lims[2]], [params.y_arena_lims[1], params.y_arena_lims[1]], color=bound_color)
        ax.plot([params.x_arena_lims[1], params.x_arena_lims[2]], [params.y_arena_lims[2], params.y_arena_lims[2]], color=bound_color)
        ax.set_xlim((params.x_arena_lims[1] - pad, params.x_arena_lims[2] + pad))
        ax.set_ylim((params.y_arena_lims[1] - pad, params.y_arena_lims[2] + pad))
    end

    # Extra formatting
    plt.xlabel("X [m]")
    plt.ylabel("Y [m]")
    plt.legend(loc="upper right")

    # Save and show figure
    fig.savefig("$fig_path/$fname.pdf", bbox_inches="tight")
    ;
end

function plot_time_dilation(
    params::Params, 
    solutions::Array{ProcessedBranchSolution}, 
    simulations::Array{ProcessedBranchSolution};
    fname::String="default_name")

    # Create figure
    fig = plt.figure(facecolor="white", figsize=[8,8])
    plt.clf()

    # Create and format subplot
    fig = plt.figure(facecolor="white")
    ax = plt.gca()
    plt.grid(true)

    for j = 1:params.n_targs
        # Obtain uniformly-sampled τ
        τ_sim = CVector(range(0, stop=1, length=length(simulations[j].sol.t)))
        τ_sol = CVector(range(0, stop=1, length=length(solutions[j].sol.t)))

        # Core plots
        dark_color = targ_colors[j][1]
        light_color = targ_colors[j][2]
        ax.plot(τ_sim, simulations[j].sol.t, color=dark_color)
        ax.plot(τ_sol, solutions[j].sol.t, color="none", markersize=5, marker="o", markerfacecolor=light_color, markeredgecolor=dark_color, label="Target "*string(j))
    end

    # Extra formatting
    plt.xlabel(L"\tau")
    plt.ylabel("t [s]")
    plt.legend(loc="upper left")

    # Save and show figure
    fig.savefig("$fig_path/$fname.pdf", bbox_inches="tight")
    ;
end

function plot_thrust_magnitude(
    params::Params, 
    solutions::Array{ProcessedBranchSolution}, 
    simulations::Array{ProcessedBranchSolution};
    fname::String="default_name")

    # Create figure
    fig = plt.figure(facecolor="white", figsize=[8,8])
    plt.clf()

    # Create and format subplot
    fig = plt.figure(facecolor="white")
    ax = plt.gca()
    plt.grid(true)

    for j = 1:params.n_targs
        # Obtain uniformly-sampled τ
        τ_sim = CVector(range(0, stop=1, length=length(simulations[j].sol.t)))
        τ_sol = CVector(range(0, stop=1, length=length(solutions[j].sol.t)))

        # Core plots
        dark_color = targ_colors[j][1]
        light_color = targ_colors[j][2]
        ax.plot(τ_sim, simulations[j].sol.T_nrm, color=dark_color)
        ax.plot(τ_sol, vcat(solutions[j].sol.T_nrm, solutions[j].sol.T_nrm[end]), color="none", markersize=5, marker="o", markerfacecolor=light_color, markeredgecolor=dark_color, label="Target "*string(j))
    end

    # Extra formatting
    plt.xlabel(L"\tau")
    plt.ylabel(L"\|T\|_2~[N]")
    plt.legend(loc="upper left")

    # Save and show figure
    fig.savefig("$fig_path/$fname.pdf", bbox_inches="tight")
    ;
end

function plot_3vec(
    params::Params, 
    solutions::Array{ProcessedBranchSolution}, 
    simulations::Array{ProcessedBranchSolution},
    vec_name::String="r";
    fname::String="default_name")

    # Create figure
    fig = plt.figure(facecolor="white", figsize=[8,8])
    plt.clf()

    # Create and format subplot
    fig, axs = plt.subplots(1, 3, facecolor="white", constrained_layout=true, figsize=[12,4])
    plt.grid(true)

    # Core plots
    for j = 1:params.n_targs
        # Obtain uniformly-sampled τ
        τ_sim = CVector(range(0, stop=1, length=length(simulations[j].sol.t)))
        τ_sol = CVector(range(0, stop=1, length=length(solutions[j].sol.t)))
        dark_color = targ_colors[j][1]
        light_color = targ_colors[j][2]

        for (ind,ax) in enumerate(axs)
            if vec_name == "r"
                y_sim = simulations[j].sol.r[ind,:]
                y_sol = solutions[j].sol.r[ind,:]
            elseif vec_name == "v"
                y_sim = simulations[j].sol.v[ind,:]
                y_sol = solutions[j].sol.v[ind,:]
            elseif vec_name == "T"
                y_sim = simulations[j].sol.T[ind,:]
                y_sol = solutions[j].sol.T[ind,:]
            end
            ax.plot(τ_sim, y_sim, color=dark_color)
            ax.plot(τ_sol, y_sol, color="none", markersize=5, marker="o", markerfacecolor=light_color, markeredgecolor=dark_color, label="Target "*string(j))
            ax.set_xlabel(L"\tau")
            ax.set_ylabel(vec_name*"["*string(ind)*"]")
        end
    end


    # Save and show figure
    fig.savefig("$fig_path/$fname.pdf", bbox_inches="tight")
    ;
end

# function plot_states(params::Params, solutions::Array{BranchSolution})

#     labels = ["pos-X [m]","pos-Y [m]","pos-Z [m]","vel-X [m/s]","vel-Y [m/s]","vel-Z [m/s]","acc-X [m/s2]","acc-Y [m/s2]","acc-Z [m/s2]"]
#     data = Vector{Matrix}()
#     for j in 1:params.n_targs
#         data_targ = vcat(
#             solutions[j].sol.r,
#             solutions[j].sol.v,
#             hcat(solutions[j].sol.T / params.mass, [0;0;0]),
#         )
#         push!(data, data_targ)
#     end

#     # Create figure
#     fun_name = nameof(var"#self#")
#     fig, axs = plt.subplots(3, 3, facecolor="white", constrained_layout=true, figsize=[8,6])

#     for (ind,ax) in enumerate(axs)
#         for j in 1:params.n_targs
#             dark_color = targ_colors[j][1]
#             light_color = targ_colors[j][2]
#             ax.plot(solutions[j].sol.t, data[j][ind,:], color=dark_color, markersize=5, marker="o", markerfacecolor=light_color, markeredgecolor=dark_color)
#         end
#         ax.grid(true)
#         ax.set_xlabel("Time [s]")
#         ax.set_ylabel(labels[ind])
#     end
#     ;
# end