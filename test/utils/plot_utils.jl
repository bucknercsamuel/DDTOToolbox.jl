vectorize_components(V) = [[v] for v in V]

# ..:: General plotting functions ::..

"""
Generic function for DDTO-formatted trajectory bundles
"""
function plot_bundle(ax,
        data_sols, 
        data_sims, 
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
    τ_split_sol_lookup(j) = params.a.τ_targs[findfirst(i->i==j, params.a.λ_targs)]
    τ_split_sim_lookup(j) = max((τ_split_sol_lookup(j)-1)*params.a.N_sim + 1, 1) |> Int
    comps = 1:length(data_sols)
    
    # Extract DDTO-segmented solutions from traj bundles
    if show_ddto_split
        # Extract trunk from solution by finding the second-to-last deferred traj (this is the final "split point")
        idx_trunk = params.a.n_targs > 1 ? params.a.λ_targs[end-1] : 1
        τ_split_sol = τ_split_sol_lookup(idx_trunk)
        τ_split_sim = τ_split_sim_lookup(idx_trunk)
        data_sols_trunk = []
        data_sims_trunk = []
        for c∈comps
            append!(data_sols_trunk, [data_sols[c][idx_trunk][1:τ_split_sol]])
            append!(data_sims_trunk, [data_sims[c][idx_trunk][1:τ_split_sim]])
        end
        
        # Extract branches from solution
        data_sols_branch = []
        data_sims_branch = []
        for c∈comps
            data_sols_branch_c = []
            data_sims_branch_c = []
            for j = 1:params.a.n_targs
                if params.a.n_targs == 1
                    idx = 1
                elseif j != params.a.λ_targs[end]
                    idx = j
                else
                    idx = params.a.λ_targs[end-1]
                end
                τ_split_sol = τ_split_sol_lookup(idx)
                τ_split_sim = τ_split_sim_lookup(idx)
                append!(data_sols_branch_c, [data_sols[c][j][τ_split_sol:end]])
                append!(data_sims_branch_c, [data_sims[c][j][τ_split_sim:end]])
            end
            append!(data_sols_branch, [data_sols_branch_c])
            append!(data_sims_branch, [data_sims_branch_c])
        end
    else
        data_sols_branch = data_sols
        data_sims_branch = data_sims
    end

    # Plot simulated data
    for j = 1:params.a.n_targs
        lines!(ax,
            [data_sims_branch[c][j] for c∈comps]...;
            style_sim..., :alpha=>alpha, :color=>color_branch(j))
    end
    if show_ddto_split
        lines!(ax,
            data_sims_trunk...;
            style_sim..., :alpha=>alpha, :color=>color_trunk)
    end

    # Plot solution (optimization) data
    if show_sol_nodes
        for j = 1:params.a.n_targs
            scatter!(ax,
                [data_sols_branch[c][j] for c∈comps]...;
                style_sol..., :alpha=>alpha, :color=>bright_color(color_branch(j)), :strokecolor=>(color_branch(j),alpha))
        end
        if show_ddto_split
            if show_defer_nodes
                data_sols_trunk_ = [data_sols_trunk[c][Not(params.a.τ_targs[1:end-1])] for c∈comps]
            else
                data_sols_trunk_ = data_sols_trunk
            end
            scatter!(ax,
                data_sols_trunk_...;
                style_sol..., :alpha=>alpha, :color=>bright_color(color_trunk), :strokecolor=>(color_trunk,alpha))
        end
    end

    # Plot deferred nodes on trunk
    if show_defer_nodes
        for j = 1:params.a.n_targs
            τ_split = τ_split_sol_lookup(j)
            if j != params.a.λ_targs[end] # don't plot final deferrable node
                scatter!(ax,
                    [data_sols[c][j][τ_split] for c∈comps]...;
                    style_sol..., :alpha=>alpha, :color=>bright_color(color_branch(j)), :strokecolor=>(color_branch(j),alpha), :marker=>defer_node_marker)
                if show_defer_times
                    if isempty(defer_times)
                        println("Warning: must specify defer times if using `show_defer_times` option.")
                        t_defer = 0
                    else
                        t_defer = defer_times[j]
                    end
                    defer_time_str = string(round(t_defer,digits=2)) * "s"
                    if nc == 2
                        xlim = ax.xaxis.attributes.limits[]
                        ylim = ax.yaxis.attributes.limits[]
                        Δx = xlim[2]-xlim[1]
                        Δy = ylim[2]-ylim[1]
                        text!(ax,
                            defer_time_str,
                            position = tuple([data_sols[1][j][τ_split], data_sols[2][j][τ_split]]...) .+ (0.05*Δx,0),
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
end

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
    return plot_bundle(ax,
        [x_sols,y_sols], 
        [x_sims,y_sims], 
        params,
        style_sim,
        style_sol;
        color_trunk=color_trunk,
        color_branch=color_branch,
        show_sol_nodes=show_sol_nodes,
        show_defer_nodes=show_defer_nodes,
        show_ddto_split=show_ddto_split,
        show_defer_times=show_defer_times,
        defer_node_marker=defer_node_marker,
        alpha=alpha,
        defer_times=defer_times,
    )
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
Gets a darker version of a color (interpolated by some fraction between black and the color)
"""
function dark_color(color; fraction=.5)
    if color isa Symbol
        colorant = parse(Colorant, color)
    else
        colorant = convert(RGB, color)
    end
    bright_color = weighted_color_mean(fraction, colorant"black", colorant)
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

function boxframe_3D(ax, lower_corner, Δcorner; style=Dict())
    # Draws the outline of a box in R3

    l1,l2,l3 = lower_corner
    u1,u2,u3 = lower_corner + Δcorner
    p = [
        [l1, l2, l3],
        [l1, l2, u3],
        [l1, u2, l3],
        [u1, l2, l3],
        [u1, u2, l3],
        [l1, u2, u3],
        [u1, l2, u3],
        [u1, u2, u3]
    ]
    for p1 in p
        for p2 in p
            if count([p1[k] == p2[k] for k=1:3]) == 2
                lines!(ax, [[p1[k],p2[k]] for k=1:3]...; style...)
            end
        end
    end
end

function draw_cone_3d(ax, vertex, pointing_direction, half_angle; style=Dict(), number_circle_elems=100, length=1, cmap=Any, draw=true, project_z0=false)
   # Draws a cone in R3

   v = vertex
   n = normalize(pointing_direction)
   θ = half_angle
   L = length
   N = number_circle_elems

   # Obtain an (arbitrary) vector perpendicular to "n" 
   # (make sure it is not equivalent to "n" or this will fail!)
   ϵ = 1e-4
   rand_vec = [1,0,0]
   if n ≈ rand_vec
       rand_vec = normalize(rand_vec + [ϵ,ϵ,ϵ])
   end
   np = normalize(cross(n,rand_vec))

   # Obtain a DCM that rotates around "n" by some angle "ψ"
   R(ψ) = quat_to_dcm([cos(ψ/2), sin(ψ/2)*n...])

   # Used for band!
   lower = fill(Point3f(v), N)
   upper = [Point3f(v + L*n + R(ψ)*np*L*tan(θ)) for ψ∈range(0,2pi, length=N)]

   # TODO: hacky solution, may bug out under certain conditions, rework when needed
   if project_z0
       for k = 1:N
           l = lower[k]
           u = upper[k]
           dist = norm(u-l)
           dir = normalize(u-l)
           Δh = l[3] - u[3]
           u_new = Point3f(v + l[3] / Δh * dist * dir)
           upper[k] = u_new 
       end
   end

   map = vcat(cmap(Int(N/2)), reverse(cmap(Int(N/2))))
   col = repeat(map,outer=2)
   if draw
       band = band!(ax, lower, upper; style..., color=col, rasterize=true)
   else
       band = undef
   end

   return band, lower, upper
end

"""
Compute cylinder band geometry (lower and upper point rings) for band!.
Returns (lower::Vector{Point3f}, upper::Vector{Point3f}) so that Observables can drive updates.
"""
function cylinder_band_points(vertex, pointing_direction, radius; length=1, N=100)
   n = normalize(pointing_direction)
   v = vertex
   ρ = radius
   L = length
   rand_vec = [1,0,0]
   ϵ = 1e-4
   if n ≈ rand_vec
       rand_vec = normalize(rand_vec + [ϵ,ϵ,ϵ])
   end
   np = normalize(cross(n, rand_vec))
   R(ψ) = quat_to_dcm([cos(ψ/2), sin(ψ/2)*n...])
   lower = [Point3f(v +       R(ψ)*np*ρ) for ψ∈range(0,2pi, length=N)]
   upper = [Point3f(v + L*n + R(ψ)*np*ρ) for ψ∈range(0,2pi, length=N)]
   return lower, upper
end

"""
Compute circle (disk) band geometry for band!.
Returns (lower::Vector{Point3f}, upper::Vector{Point3f}).
"""
function circle_band_points(vertex, pointing_direction, radius; N=100)
   n = normalize(pointing_direction)
   v = vertex
   ρ = radius
   rand_vec = [1,0,0]
   ϵ = 1e-4
   if n ≈ rand_vec
       rand_vec = normalize(rand_vec + [ϵ,ϵ,ϵ])
   end
   np = normalize(cross(n, rand_vec))
   R(ψ) = quat_to_dcm([cos(ψ/2), sin(ψ/2)*n...])
   lower = [Point3f(v)             for ψ∈range(0,2pi, length=N)]
   upper = [Point3f(v + R(ψ)*np*ρ) for ψ∈range(0,2pi, length=N)]
   return lower, upper
end

function draw_cylinder_3d(ax, vertex, pointing_direction, radius; style=Dict(), number_circle_elems=100, length=1, cmap=Any, draw=true, draw_caps=true)
   # Draws a cylinder in R3
   v = vertex
   n = normalize(pointing_direction)
   ρ = radius
   L = length
   N = number_circle_elems

   lower, upper = cylinder_band_points(v, n, ρ; length=L, N=N)

   map = vcat(cmap(Int(N/2)), reverse(cmap(Int(N/2))))
   col = repeat(map,outer=2)
   if draw
       band = band!(ax, lower, upper; style..., color=col, rasterize=true)
   else
       band = undef
   end

    # Add caps to the cylinder
    if draw_caps
        draw_circle_3d(ax, v, ρ; pointing_direction=n, style=style, color=col, number_circle_elems=number_circle_elems, draw=draw)
        draw_circle_3d(ax, v + L*n, ρ; pointing_direction=n, style=style, color=col, number_circle_elems=number_circle_elems, draw=draw)
    end

   return band, lower, upper
end

function draw_circle_3d(ax, vertex, radius; pointing_direction=[0,0,1], style=Dict(), color=:yellow, number_circle_elems=100, draw=true)
   # Draws a 2D circle in R3 (defaults to XY plane with pointing_direction=[0,0,1])
   v = vertex
   n = normalize(pointing_direction)
   ρ = radius
   N = number_circle_elems
   lower, upper = circle_band_points(v, n, ρ; N=N)
   if draw
       band = band!(ax, lower, upper; style..., color=color, rasterize=true)
   else
       band = undef
   end
   return band, lower, upper
end

function get_equal_3d_lims(initial_position, final_position; pad=0.2)
    # Obtain equally-spaced 3D limits based on guidance boundary conditions
    #
    # :in initial_position: Guidance initial position
    # :in final_position: Guidance final position

    r0 = initial_position
    rf = final_position

    Δr = r0 - rf
    bounding_box_width = max(abs.(Δr)...)
    view_width = bounding_box_width * (1 + pad)

    centroid = (r0 + rf)/2
    Δ_lim = [-view_width/2, view_width/2]
    xlims = centroid[1] .+ Δ_lim
    ylims = centroid[2] .+ Δ_lim
    zlims = centroid[3] .+ Δ_lim

    return xlims, ylims, zlims
end

function get_target_allocations(results)
    # For each target ID, determine periods of active allocation
    #
    # :in results: Results from a simulation
    #
    # :out target_allocations: Target allocations

    target_allocations = Dict()
    active_targets = []
    for idx in 1:length(results["sim_time"])
        for idx2 in findall(x->x==1, results["targpool_allocated"][:,idx])
            id = results["targpool_ID"][idx2]

            # Add new entry to target_allocations if target hasn't been added yet
            if !(id in keys(target_allocations))
                target_allocations[id] = []
            end

            # Add target to active_targets if it hasn't been added yet
            if !(id in active_targets)
                push!(active_targets, id)
                push!(target_allocations[id], [idx])
            end
        end
        # Remove target from active_targets if it is no longer active
        for id in active_targets
            if !(id in results["targs_ID"][:,idx])
                push!(target_allocations[id][end], idx-1)
                deleteat!(active_targets, findfirst(x->x==id, active_targets))
            end
        end
    end

    # Remove remaining targets from activity once the final time is reached
    for id in active_targets
        push!(target_allocations[id][end], length(results["sim_time"]))
    end

    # Stitch segments together that differ by less than 3 time-steps (accounts for recomputation discontinuities)
    for id in keys(target_allocations)
        k = 1  
        for j = 1:length(target_allocations[id])-1
            if target_allocations[id][k+1][1] - target_allocations[id][k][2] <= 3
                target_allocations[id][k][2] = target_allocations[id][k+1][2]
                deleteat!(target_allocations[id], k+1)
            else
                k += 1
            end
        end
    end

    return target_allocations
end

function hold_interactive(screens)
    println("\nPress any key when finished using plots...")
    readline() # Wait for user to finish plotting
    [GLMakie.destroy!(screen) for screen in screens]
end

"""
Write me a function that, given a figure and 3D axis object with Makie.jl, 
produces a new axis that projects the 3D plot onto a 2D plane with a specified normal vector
and also displays the projection on the original 3D plot
"""
function add_2D_projection(fig, ax_3d, normal_vector, ax_idx_new)
    # Create a new axis for the 2D projection
    ax_proj = Axis(fig[ax_idx_new...])

    # Take all data from the 3D plot and project it onto the 2D plot if its a line object
    display(ax_3d)
    plots = ax_3d.children
    display(plots)


end