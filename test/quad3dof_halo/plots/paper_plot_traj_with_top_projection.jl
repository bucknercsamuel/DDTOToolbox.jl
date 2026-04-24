using CairoMakie
include("plot_defaults.jl")
include("plot_3d_trajs.jl")
include("plot_2d_trajs_XY.jl")

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
