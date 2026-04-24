using CairoMakie
include("plot_defaults.jl")
include("plot_3d_trajs.jl")
include("plot_2d_trajs_XY.jl")

function paper_plot_greedy_compare(
    results_all;
    interactive = true,
    azel=(pi/4,pi/4)
    )

    with_theme(theme3d; fontsize=fontsize) do
        f = Figure(size=(1200,700))
        plot_3d_trajs(results_all[1]; interactive=false, f=f, ax_idx=[1,1], azel=azel, save_fig=false)
        plot_3d_trajs(results_all[2]; interactive=false, f=f, ax_idx=[1,2], azel=azel, save_fig=false)
        plot_3d_trajs(results_all[3]; interactive=false, f=f, ax_idx=[1,3], azel=azel, save_fig=false)
        plot_2d_trajs_XY(results_all[1]; interactive=false, f=f, ax_idx=[2,1], save_fig=false)
        plot_2d_trajs_XY(results_all[2]; interactive=false, f=f, ax_idx=[2,2], save_fig=false)
        plot_2d_trajs_XY(results_all[3]; interactive=false, f=f, ax_idx=[2,3], save_fig=false)

        if interactive
            screen = GLMakie.Screen()
            display(screen, f)
            return screen
        else
            save(joinpath(fig_path, "greedy_compare"*fig_ext), f)
        end
    end
end
