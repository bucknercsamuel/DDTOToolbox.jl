using CairoMakie
include("plot_defaults.jl")
include("plot_3d_trajs.jl")
include("plot_target_allocations.jl")

function paper_plot_trajallocation(
    params,
    results;
    interactive = true,
    azel=(pi/4,pi/4)
    )

    with_theme(theme3d; fontsize=fontsize) do
        f = Figure(size=(1500,700))
        plot_3d_trajs(results; interactive=false, f=f, ax_idx=[1,1], azel=azel)
        plot_target_allocations(params, results; interactive=false, f=f, ax_idx=[1,2])

        if interactive
            screen = GLMakie.Screen()
            display(screen, f)
            return screen
        else
            save(joinpath(fig_path, "trajallocation"*fig_ext), f)
        end
    end
end
