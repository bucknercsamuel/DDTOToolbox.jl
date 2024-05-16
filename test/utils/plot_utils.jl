vectorize_components(V) = [[v] for v in V]

# ..:: General plotting functions ::..

"""
Generic function for 2D DDTO-formatted trajectory bundles
"""
function plot2D_bundle(ax,
        x_sols, 
        x_sims, 
        y_sols, 
        y_sims,
        params,
        style_sim,
        style_sol; 
        color_trunk = :black,
        color_branch = N->:gray, # color range across data set (function called on iterate)
        show_sol_nodes = true,
        show_defer_nodes = true,
        show_ddto_split = true, # determines if we should split data into a trunk and branches
        show_defer_times = false,
        defer_node_marker = :diamond,
        alpha = 1,
        defer_times = []
    )
    # Helper functions
    τ_split_sol_lookup(j) = params.a.τ_targs[findfirst(i->i==j, params.a.λ_targs)] # obtain the deferrability index of the j-th target (solution)
    τ_split_sim_lookup(j) = max((τ_split_sol_lookup(j)-1)*Int(round((length(x_sims[j])/(length(x_sols[j])-1))))+1,1) |> round |> Int

    # Extract DDTO-segmented solutions from traj bundles
    if show_ddto_split
        # Extract trunk from solution by finding the second-to-last deferred traj (this is the final "split point")
        τ_split_sol = τ_split_sol_lookup(params.a.λ_targs[end-1])
        τ_split_sim = τ_split_sim_lookup(params.a.λ_targs[end-1])
        x_sols_trunk = x_sols[params.a.λ_targs[end-1]][1:τ_split_sol]
        x_sims_trunk = x_sims[params.a.λ_targs[end-1]][1:τ_split_sim]
        y_sols_trunk = y_sols[params.a.λ_targs[end-1]][1:τ_split_sol]
        y_sims_trunk = y_sims[params.a.λ_targs[end-1]][1:τ_split_sim]
        
        # Extract branches from solution
        x_sols_branch = []
        x_sims_branch = []
        y_sols_branch = []
        y_sims_branch = []
        for j = 1:params.a.n_targs
            if j != params.a.λ_targs[end]
                idx = j
            else
                idx = params.a.λ_targs[end-1]
            end
            τ_split_sol = τ_split_sol_lookup(idx)
            τ_split_sim = τ_split_sim_lookup(idx)
            push!(x_sols_branch, x_sols[j][τ_split_sol:end])
            push!(x_sims_branch, x_sims[j][τ_split_sim:end])
            push!(y_sols_branch, y_sols[j][τ_split_sol:end])
            push!(y_sims_branch, y_sims[j][τ_split_sim:end])
        end
    else
        x_sols_branch = x_sols
        x_sims_branch = x_sims
        y_sols_branch = y_sols
        y_sims_branch = y_sims
    end

    # Plot simulated data
    for j = 1:params.a.n_targs
        lines!(ax,
            x_sims_branch[j],
            y_sims_branch[j];
            style_sim..., :alpha=>alpha, :color=>color_branch(j))
    end
    if show_ddto_split
        lines!(ax,
            x_sims_trunk,
            y_sims_trunk;
            style_sim..., :alpha=>alpha, :color=>color_trunk)
    end

    # Plot solution (optimization) data
    if show_sol_nodes
        for j = 1:params.a.n_targs
            scatter!(ax,
                x_sols_branch[j],
                y_sols_branch[j];
                style_sol..., :alpha=>alpha, :strokealpha=>alpha, :color=>bright_color(color_branch(j)), :strokecolor=>(color_branch(j),alpha))
        end
        if show_ddto_split
            if show_defer_nodes
                x_sols_trunk_ = x_sols_trunk[Not(params.a.τ_targs[1:end-1])]
                y_sols_trunk_ = y_sols_trunk[Not(params.a.τ_targs[1:end-1])]
            else
                x_sols_trunk = x_sols_trunk
                y_sols_trunk = y_sols_trunk
            end
            scatter!(ax,
                x_sols_trunk_,
                y_sols_trunk_;
                style_sol..., :alpha=>alpha, :color=>bright_color(color_trunk), :strokecolor=>(color_trunk,alpha))
        end
    end

    # Plot deferred nodes on trunk
    if show_defer_nodes
        for j = 1:params.a.n_targs
            τ_split = τ_split_sol_lookup(j)
            if j != params.a.λ_targs[end] # don't plot final deferrable node
                scatter!(ax,
                    x_sols[j][τ_split],
                    y_sols[j][τ_split];
                    style_sol..., :alpha=>alpha, :strokealpha=>alpha, :color=>bright_color(color_branch(j)), :strokecolor=>(color_branch(j),alpha), :marker=>defer_node_marker)
                if show_defer_times
                    if isempty(defer_times)
                        println("Warning: must specify defer times if using `show_defer_times` option.")
                        t_defer = 0
                    else
                        t_defer = defer_times[j]
                    end
                    defer_time_str = string(round(t_defer,digits=2)) * "s"
                    xlim = ax.xaxis.attributes.limits[]
                    ylim = ax.yaxis.attributes.limits[]
                    Δx = xlim[2]-xlim[1]
                    Δy = ylim[2]-ylim[1]
                    text!(ax,
                        defer_time_str,
                        position = tuple([x_sols[j][τ_split], y_sols[j][τ_split]]...) .+ (0.05*Δx,0),
                        align = (:left, :top),
                        # position = tuple([x_sols[j][τ_split], y_sols[j][τ_split]]...) .+ (0,.01*Δy),
                        # align = (:center, :bottom),
                        font = :bold,
                        color = (color_branch(j),alpha),
                        glowwidth = 5,
                        glowcolor = (:white,1)
                    )
                end
            end
        end
    end
end

"""
Overlays constraint bounds on an axis
"""
function overlay2D_constraint_bounds(ax,
        time::Vector,
        lowerbound::Vector,
        upperbound::Vector;
        lowerbound_lim::Float64=0.,
        upperbound_lim::Float64=1.,
        scaling::Float64=1.0, 
        padding::Float64=0.2,
        style_fill=Dict(),
        style_edge=Dict(),
        labels=false
    )
    # Scale constraint bounds
    lowerbound     *= scaling
    upperbound     *= scaling
    upperbound_lim *= scaling
    lowerbound_lim *= scaling

    # Obtain plot bounds
    Δy = upperbound_lim - lowerbound_lim
    lowerbound_lim -= Δy*padding
    upperbound_lim += Δy*padding
    lb_plot = fill(lowerbound_lim, length(time))
    ub_plot = fill(upperbound_lim, length(time))
    lb_plot[findall(x->x==Inf, lowerbound)] .= Inf
    ub_plot[findall(x->x==Inf, upperbound)] .= Inf

    # Fill in constraint
    some(type, arr) = ~all(type, arr) && any(type, arr)
    label1 = labels ? "Constraint Set" : nothing
    label2 = labels ? "Constraint Bound" : nothing
    if ~all(isinf, lowerbound)
        if ~some(isinf, lowerbound)
            band!(ax, Float64.(time), lb_plot, Float64.(lowerbound); style_fill..., :label=>label1)
            lines!(ax, Float64.(time), Float64.(lowerbound); style_edge..., :label=>label2)
        else
            inf_elems = findall(x->x==Inf, lowerbound)
            switch_points = findall(x->x>1, diff(inf_elems))
            for k = 1:length(switch_points)
                k0 = inf_elems[switch_points[k]]+1
                kf = inf_elems[switch_points[k]+1]-1
                lim_frac(k) = (lowerbound[k] - lowerbound_lim)/(upperbound_lim - lowerbound_lim)
                band!(ax, Float64.(time)[k0:kf], lb_plot[k0:kf], Float64.(lowerbound)[k0:kf]; style_fill..., :label=>label1)
                lines!(ax, Float64.(time)[k0:kf], Float64.(lowerbound)[k0:kf]; style_edge..., :label=>label2)
                vlines!(ax, [time[k0], time[kf]], ymin=[0,0], ymax=[lim_frac(k0), lim_frac(kf)]; style_edge...)
            end
        end
        labels = false
        label1 = labels ? "Constraint Set" : nothing
        label2 = labels ? "Constraint Bound" : nothing
    end
    if ~all(isinf, upperbound)
        if ~some(isinf, upperbound)
            band!(ax, Float64.(time), Float64.(upperbound), ub_plot; style_fill..., :label=>label1)
            lines!(ax, Float64.(time), Float64.(upperbound); style_edge..., :label=>label2)
        else
            inf_elems = findall(x->x==Inf, upperbound)
            switch_points = findall(x->abs(x)>1, diff(inf_elems))
            for k = 1:length(switch_points)
                k0 = inf_elems[switch_points[k]]+1
                kf = inf_elems[switch_points[k]+1]-1
                lim_frac(k) = (upperbound[k] - lowerbound_lim)/(upperbound_lim - lowerbound_lim)
                band!(ax, Float64.(time)[k0:kf], Float64.(upperbound)[k0:kf], ub_plot[k0:kf]; style_fill..., :label=>label1)
                lines!(ax, Float64.(time)[k0:kf], Float64.(upperbound)[k0:kf]; style_edge..., :label=>label2)
                vlines!(ax, [time[k0], time[kf]], ymin=[lim_frac(k0), lim_frac(kf)], ymax=[1,1]; style_edge...)
            end
        end
    end

    # Set plot limits
    xlims!(ax, (time[1], time[end]))
    ylims!(ax, (lowerbound_lim, upperbound_lim))
end

"""
Gets a brighter version of a color (interpolated by some fraction between white and the color)
"""
function bright_color(color; fraction=.5)
    if color isa Symbol
        colorant = parse(Colorant, color)
    else
        colorant = convert(RGB, color)
    end
    bright_color = weighted_color_mean(fraction, colorant"white", colorant)
    return bright_color
end

"""
Draws a circle in R2
"""
function draw2d_circle(ax, center, radius; color=:red, alpha=0.5, N=100)
    # circle body
    c = center
    r = radius
    lower = fill(Point2f(c), N)
    upper = [Point2f(c + r*[sin(ψ),cos(ψ)]) for ψ∈range(0,2pi,N)]
    band!(ax, lower, upper; color=(color,alpha))

    # circle edge
    arc!(ax, center, radius, 0, 2pi; color=color, alpha=alpha)
end