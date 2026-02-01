using CairoMakie
# using GLMakie
using Colors
using InvertedIndices
include("../../utils/plot_utils.jl")
include("plot_defaults.jl")

function plot_mc_compare(
        results_dict;
        interactive = true
    )
    # Show tick marks on x and y axis
    # Axis settings
    axis_defaults = Dict(
        # :xautolimitmargin=>(0,0), 
        :topspinevisible=>true, 
        :rightspinevisible=>true,
        :xgridvisible=>false,
        :ygridvisible=>false,
    )
    ax_label_size = 15

    # Figure setup
    f = Figure(size=(1400,400))

    # Unpack results
    convergence_container_lex = results_dict["convergence_container_lex"]
    solver_time_container_lex = results_dict["solver_time_container_lex"]
    deferral_time_container_lex = results_dict["deferral_time_container_lex"]
    convergence_container_scp = results_dict["convergence_container_scp"]
    solver_time_container_scp = results_dict["solver_time_container_scp"]
    deferral_time_container_scp = results_dict["deferral_time_container_scp"]
    num_targ_levels = length(convergence_container_lex)
    x_range = collect(keys(convergence_container_lex))
    x_range = [Int(j) for j in x_range]
    x_range = sort(x_range)
    targ_levels = x_range
    n_trials = length(convergence_container_lex[targ_levels[1]])

    # Obtain objective proxy of deferral time
    deferral_obj_container_lex = Dict()
    deferral_obj_container_scp = Dict()
    for j in targ_levels
        deferral_obj_container_lex[j] = zeros(n_trials)
        deferral_obj_container_scp[j] = zeros(n_trials)
        for trial = 1:n_trials
            deferral_obj_container_lex[j][trial] = sum(deferral_time_container_lex[j][trial,:])
            deferral_obj_container_scp[j][trial] = sum(deferral_time_container_scp[j][trial,:])
        end
    end

    # Convert convergence values to percentages
    for j in targ_levels
        convergence_container_lex[j] = [convergence_container_lex[j][i] * 100 for i in 1:n_trials]
        convergence_container_scp[j] = [convergence_container_scp[j][i] * 100 for i in 1:n_trials]
    end

    # Define function to compute mean across trials and use a 1-sigma funnel to show variability
    function plot_mean_and_funnel(ax, data, label, colors; funnel=true)
        means = [mean(data[j]) for j in targ_levels]
        stds = [std(data[j]) for j in targ_levels]
        means_upper = means .+ stds
        means_lower = means .- stds
        if funnel
            band!(ax,
                x_range,
                means_lower,
                means_upper;
                color=colors,
                alpha=0.2)
        end
        lines!(ax,
            x_range,
            means;
            color=colors,
            linewidth=2,
            label=label)

        # update axis ticks to only include target levels
        ax.xticks = (x_range, string.(x_range))
    end

    scp_label = "Graph-DDTO"
    lex_label = "Lex-DDTO"

    # Plot objective
    ax = Axis(f[1,1], xlabel="Number of Targets", ylabel=L"$\Sigma$ Deferral Times [s]"; axis_defaults...)
    plot_mean_and_funnel(ax, deferral_obj_container_lex, lex_label, :red)
    plot_mean_and_funnel(ax, deferral_obj_container_scp, scp_label, :blue)
    axislegend(ax, position=:lt, labelsize=ax_label_size)

    # Plot solver time
    ax = Axis(f[1,2], xlabel="Number of Targets", ylabel="Solver Time [s]"; axis_defaults...)
    plot_mean_and_funnel(ax, solver_time_container_lex, lex_label, :red)
    plot_mean_and_funnel(ax, solver_time_container_scp, scp_label, :blue)
    axislegend(ax, position=:lt, labelsize=ax_label_size)

    # Plot convergence
    ax = Axis(f[1,3], xlabel="Number of Targets", ylabel="% Converged"; axis_defaults...)
    plot_mean_and_funnel(ax, convergence_container_lex, lex_label, :red; funnel=false)
    plot_mean_and_funnel(ax, convergence_container_scp, scp_label, :blue; funnel=false)
    axislegend(ax, position=:lb, labelsize=ax_label_size)

    # Legend
    # Legend(f[2,1:3], ax, framevisible=true)

    if interactive
        screen = GLMakie.Screen()
        display(screen, f)
        return screen
    else
        CairoMakie.save(joinpath(fig_path, "mc_compare"*fig_ext), f)
    end
end