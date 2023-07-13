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
    defer_solution::ProcessedSolution=EmptyProcessedSolution(),
    defer_simulation::ProcessedSolution=EmptyProcessedSolution(),
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

    # Deferrable segment plot
    if !isempty(defer_solution.t) && !isempty(defer_simulation.t)
        ax.plot(defer_simulation.r[1,:], defer_simulation.r[2,:], color="black")
        ax.plot(defer_solution.r[1,:], defer_solution.r[2,:], color="none", markersize=5, marker="o", markerfacecolor="gray", markeredgecolor="black", label="Deferred")
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
    defer_solution::ProcessedSolution=EmptyProcessedSolution(),
    defer_simulation::ProcessedSolution=EmptyProcessedSolution(),
    fname::String="default_name")

    # Create figure
    fig = plt.figure(facecolor="white", figsize=[8,8])
    plt.clf()

    # Create and format subplot
    fig = plt.figure(facecolor="white")
    ax = plt.gca()
    plt.grid(true)

    # Trajectory plots
    for j = 1:params.n_targs
        N_sim = length(simulations[j].sol.t)
        N_sol = length(solutions[j].sol.t)
        N_sim_ctrl = length(simulations[j].sol.s)
        N_sol_ctrl = length(solutions[j].sol.s)
        τ_sim = CVector(range(0, stop=1, length=N_sim))
        τ_sol = CVector(range(0, stop=1, length=N_sol))
        dτ_sim = diff(τ_sim)
        dτ_sol = diff(τ_sol)
        dτ_sim = vcat(dτ_sim, dτ_sim[end])
        dτ_sol = vcat(dτ_sol, dτ_sol[end])
        if params.disc == 0
            t_sim = cumsum([simulations[j].sol.s[k] * dτ_sim[k] for k = 1:N_sim-1])
            t_sol = cumsum([solutions[j].sol.s[k] * dτ_sol[k] for k = 1:N_sol-1])
        elseif params.disc == 1
            t_sim = cumsum([(1/2) * (simulations[j].sol.s[k] + simulations[j].sol.s[k+1]) * dτ_sim[k] for k = 1:N_sim-1])
            t_sol = cumsum([(1/2) * (solutions[j].sol.s[k] + solutions[j].sol.s[k+1]) * dτ_sol[k] for k = 1:N_sol-1])
        end
        t_sim = vcat(0, t_sim)
        t_sol = vcat(0, t_sol)

        # Core plots
        dark_color = targ_colors[j][1]
        light_color = targ_colors[j][2]
        ax.plot(τ_sim, t_sim, color=dark_color)
        ax.plot(τ_sol, t_sol, color="none", markersize=5, marker="o", markerfacecolor=light_color, markeredgecolor=dark_color, label="Target "*string(j))
    end

    # Deferrable segment plot
    if !isempty(defer_solution.t) && !isempty(defer_simulation.t)
        N_sim = length(defer_simulation.t)
        N_sol = length(defer_solution.t)
        N_sim_ctrl = length(defer_simulation.s)
        N_sol_ctrl = length(defer_solution.s)
        ratio_sim = (N_sim-1) / (length(simulations[1].sol.t)-1)
        ratio_sol = (N_sol-1) / (length(solutions[1].sol.t)-1)
        τ_sim = CVector(range(0, stop=ratio_sim, length=N_sim))
        τ_sol = CVector(range(0, stop=ratio_sol, length=N_sol))
        dτ_sim = diff(τ_sim)
        dτ_sol = diff(τ_sol)
        dτ_sim = vcat(dτ_sim, dτ_sim[end])
        dτ_sol = vcat(dτ_sol, dτ_sol[end])
        if params.disc == 0
            t_sim = cumsum([defer_simulation.s[k] * dτ_sim[k] for k = 1:N_sim-1])
            t_sol = cumsum([defer_solution.s[k] * dτ_sol[k] for k = 1:N_sol-1])
        elseif params.disc == 1
            t_sim = cumsum([(1/2) * (defer_simulation.s[k] + defer_simulation.s[k+1]) * dτ_sim[k] for k = 1:N_sim-1])
            t_sol = cumsum([(1/2) * (defer_solution.s[k] + defer_solution.s[k+1]) * dτ_sol[k] for k = 1:N_sol-1])
        end
        t_sim = vcat(0, t_sim)
        t_sol = vcat(0, t_sol)

        ax.plot(τ_sim, t_sim, color="black")
        ax.plot(τ_sol, t_sol, color="none", markersize=5, marker="o", markerfacecolor="gray", markeredgecolor="black", label="Deferred")
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
    defer_solution::ProcessedSolution=EmptyProcessedSolution(),
    defer_simulation::ProcessedSolution=EmptyProcessedSolution(),
    fname::String="default_name")

    # Create figure
    fig = plt.figure(facecolor="white", figsize=[8,8])
    plt.clf()

    # Create and format subplot
    fig = plt.figure(facecolor="white")
    ax = plt.gca()
    plt.grid(true)

    # Trajectory plots
    for j = 1:params.n_targs
        # Obtain uniformly-sampled τ
        τ_sim = CVector(range(0, stop=1, length=length(simulations[j].sol.t)))
        τ_sol = CVector(range(0, stop=1, length=length(solutions[j].sol.t)))

        # Obtain thrust
        T_sim = simulations[j].sol.T_nrm
        T_sol = solutions[j].sol.T_nrm
        if params.disc == 0
            T_sol = vcat(T_sol, T_sol[end])
        end

        # Core plots
        dark_color = targ_colors[j][1]
        light_color = targ_colors[j][2]
        ax.plot(τ_sim, T_sim, color=dark_color)
        ax.plot(τ_sol, T_sol, color="none", markersize=5, marker="o", markerfacecolor=light_color, markeredgecolor=dark_color, label="Target "*string(j))
    end

    # Deferrable segment plot
    if !isempty(defer_solution.t) && !isempty(defer_simulation.t)
        # Obtain uniformly-sampled τ
        N_sim = length(defer_simulation.t)
        N_sol = length(defer_solution.t)
        ratio_sim = (N_sim-1) / (length(simulations[1].sol.t)-1)
        ratio_sol = (N_sol-1) / (length(solutions[1].sol.t)-1)
        τ_sim = CVector(range(0, stop=ratio_sim, length=N_sim))
        τ_sol = CVector(range(0, stop=ratio_sol, length=N_sol))

        # Obtain thrust
        T_sim = defer_simulation.T_nrm
        T_sol = defer_solution.T_nrm
        if params.disc == 0
            T_sol = vcat(T_sol, T_sol[end])
        end

        ax.plot(τ_sim, T_sim, color="black")
        ax.plot(τ_sol, T_sol, color="none", markersize=5, marker="o", markerfacecolor="gray", markeredgecolor="black", label="Deferred")
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
    defer_solution::ProcessedSolution=EmptyProcessedSolution(),
    defer_simulation::ProcessedSolution=EmptyProcessedSolution(),
    fname::String="default_name")

    # Create figure
    fig = plt.figure(facecolor="white", figsize=[8,8])
    plt.clf()

    # Create and format subplot
    num_comps = 2
    fig, axs = plt.subplots(1, num_comps, facecolor="white", constrained_layout=true, figsize=[4*num_comps,4])

    # Trajectory plots
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
            ax.grid(true)
        end
    end

    # Deferrable segment plot
    if !isempty(defer_solution.t) && !isempty(defer_simulation.t)
        # Obtain uniformly-sampled τ
        N_sim = length(defer_simulation.t)
        N_sol = length(defer_solution.t)
        ratio_sim = (N_sim-1) / (length(simulations[1].sol.t)-1)
        ratio_sol = (N_sol-1) / (length(solutions[1].sol.t)-1)
        τ_sim = CVector(range(0, stop=ratio_sim, length=N_sim))
        τ_sol = CVector(range(0, stop=ratio_sol, length=N_sol))

        for (ind,ax) in enumerate(axs)
            if vec_name == "r"
                y_sim = defer_simulation.r[ind,:]
                y_sol = defer_solution.r[ind,:]
            elseif vec_name == "v"
                y_sim = defer_simulation.v[ind,:]
                y_sol = defer_solution.v[ind,:]
            elseif vec_name == "T"
                y_sim = defer_simulation.T[ind,:]
                y_sol = defer_solution.T[ind,:]
            end
            ax.plot(τ_sim, y_sim, color="black")
            ax.plot(τ_sol, y_sol, color="none", markersize=5, marker="o", markerfacecolor="gray", markeredgecolor="black", label="Deferred")
            ax.set_xlabel(L"\tau")
            ax.set_ylabel(vec_name*"["*string(ind)*"]")
            ax.grid(true)
        end
    end

    # Save and show figure
    fig.savefig("$fig_path/$fname.pdf", bbox_inches="tight")
    ;
end