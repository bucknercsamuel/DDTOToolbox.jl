using CairoMakie
using Colors
using LinearAlgebra
include("../../utils/plot_utils.jl")
include("plot_defaults.jl")

function plot_drone(
        ax,
        position,
        thrust_direction;
        scale=1
    )
    # Sizing parameters for vehicle
    b_x = thrust_direction
    b_y = cross(b_x, b_ref)
    b_y = b_y / norm(b_y)
    b_z = cross(b_x, b_y)
    body_radius = 0.15
    body_height = 0.05
    arm_length = 0.4
    arm_radius = 0.02
    prop_radius = 0.2
    prop_height = 0.05
    cmap_frame(N) = colormap("Grays",N)
    cmap_arm(N) = colormap("Grays", N)
    cmap_prop(N) = colormap("Reds", N)
    style_prop = Dict(:alpha=>.3)

    # Create a drone object with cylinders
    # inputs: ax, vertex, pointing_direction, radius; length, cmap
    draw_cylinder_3d(ax, position, b_x, scale*body_radius; length=scale*body_height, cmap=cmap_frame)
    arm_dirs = [b_y, b_z, -b_y, -b_z]
    for arm_dir in arm_dirs
        draw_cylinder_3d(ax, position, arm_dir, scale*arm_radius; length=scale*(arm_length-prop_radius), cmap=cmap_arm)
        draw_cylinder_3d(ax, position + scale*arm_length*arm_dir + scale*arm_radius*b_x, b_x, scale*prop_radius; length=scale*prop_height, cmap=cmap_prop, style=style_prop)
    end
end

"""
Draw a drone whose position and thrust direction are given by Observables, so the
drone can be updated dynamically (e.g. for real-time animation).
"""
function plot_drone_observable(
        ax,
        position_obs,
        thrust_direction_obs;
        scale=1,
        number_circle_elems=100
    )
    body_radius = 0.15
    body_height = 0.05
    arm_length = 0.4
    arm_radius = 0.02
    prop_radius = 0.2
    prop_height = 0.05
    cmap_frame(N) = colormap("Grays", N)
    cmap_arm(N) = colormap("Grays", N)
    cmap_prop(N) = colormap("Reds", N)
    style_prop = Dict(:alpha=>0.3)

    # Thrust direction (body x) and ensure normalized
    b_x_obs = lift(thrust_direction_obs) do t
        n = normalize(t)
        iszero(norm(n)) ? [1.0, 0.0, 0.0] : n
    end

    # Body y from cross(b_x, b_ref)
    b_y_obs = lift(b_x_obs) do b_x
        by = cross(b_x, b_ref)
        nby = norm(by)
        iszero(nby) ? [0.0, 1.0, 0.0] : by / nby
    end

    b_z_obs = lift(b_x_obs, b_y_obs) do b_x, b_y
        cross(b_x, b_y)
    end

    # Body cylinder
    body_band_obs = lift(position_obs, b_x_obs) do pos, b_x
        cylinder_band_points(pos, b_x, scale*body_radius; length=scale*body_height, N=number_circle_elems)
    end
    lower_body = lift(first, body_band_obs)
    upper_body = lift(last, body_band_obs)
    col_body = repeat(vcat(cmap_frame(Int(number_circle_elems/2)), reverse(cmap_frame(Int(number_circle_elems/2)))), outer=2)
    band!(ax, lower_body, upper_body; color=col_body, rasterize=true)

    # Body caps
    cap1_band_obs = lift(position_obs, b_x_obs) do pos, b_x
        circle_band_points(pos, b_x, scale*body_radius; N=number_circle_elems)
    end
    band!(ax, lift(first, cap1_band_obs), lift(last, cap1_band_obs); color=col_body, rasterize=true)
    cap2_center = lift(position_obs, b_x_obs) do pos, b_x
        pos + scale*body_height * (norm(b_x) > 1e-10 ? normalize(b_x) : [1.0, 0.0, 0.0])
    end
    cap2_band_obs = lift(cap2_center, b_x_obs) do v, b_x
        circle_band_points(v, b_x, scale*body_radius; N=number_circle_elems)
    end
    band!(ax, lift(first, cap2_band_obs), lift(last, cap2_band_obs); color=col_body, rasterize=true)

    # Arms and props: b_y, b_z, -b_y, -b_z
    arm_dirs_obs = [
        b_y_obs,
        b_z_obs,
        lift(x -> -x, b_y_obs),
        lift(x -> -x, b_z_obs)
    ]
    col_arm = repeat(vcat(cmap_arm(Int(number_circle_elems/2)), reverse(cmap_arm(Int(number_circle_elems/2)))), outer=2)
    col_prop = repeat(vcat(cmap_prop(Int(number_circle_elems/2)), reverse(cmap_prop(Int(number_circle_elems/2)))), outer=2)

    for arm_dir_obs in arm_dirs_obs
        # Arm cylinder
        arm_band_obs = lift(position_obs, arm_dir_obs) do pos, ad
            cylinder_band_points(pos, ad, scale*arm_radius; length=scale*(arm_length - prop_radius), N=number_circle_elems)
        end
        band!(ax, lift(first, arm_band_obs), lift(last, arm_band_obs); color=col_arm, rasterize=true)
        # Arm cap at base
        arm_cap1_obs = lift(position_obs, arm_dir_obs) do pos, ad
            circle_band_points(pos, ad, scale*arm_radius; N=number_circle_elems)
        end
        band!(ax, lift(first, arm_cap1_obs), lift(last, arm_cap1_obs); color=col_arm, rasterize=true)
        # Arm cap at tip (before prop)
        arm_tip = lift(position_obs, arm_dir_obs) do pos, ad
            pos + scale*(arm_length - prop_radius) * (norm(ad) > 1e-10 ? normalize(ad) : ad)
        end
        arm_cap2_obs = lift(arm_tip, arm_dir_obs) do v, ad
            circle_band_points(v, ad, scale*arm_radius; N=number_circle_elems)
        end
        band!(ax, lift(first, arm_cap2_obs), lift(last, arm_cap2_obs); color=col_arm, rasterize=true)

        # Prop at tip
        prop_center_obs = lift(position_obs, arm_dir_obs, b_x_obs) do pos, ad, b_x
            nad = norm(ad) > 1e-10 ? normalize(ad) : ad
            pos + scale*arm_length*nad + scale*arm_radius*b_x
        end
        prop_band_obs = lift(prop_center_obs, b_x_obs) do v, b_x
            cylinder_band_points(v, b_x, scale*prop_radius; length=scale*prop_height, N=number_circle_elems)
        end
        band!(ax, lift(first, prop_band_obs), lift(last, prop_band_obs); style_prop..., color=col_prop, rasterize=true)
        prop_cap1_obs = lift(prop_center_obs, b_x_obs) do v, b_x
            circle_band_points(v, b_x, scale*prop_radius; N=number_circle_elems)
        end
        band!(ax, lift(first, prop_cap1_obs), lift(last, prop_cap1_obs); style_prop..., color=col_prop, rasterize=true)
        prop_cap2_center = lift(prop_center_obs, b_x_obs) do v, b_x
            v + scale*prop_height * (norm(b_x) > 1e-10 ? normalize(b_x) : [1.0, 0.0, 0.0])
        end
        prop_cap2_obs = lift(prop_cap2_center, b_x_obs) do v, b_x
            circle_band_points(v, b_x, scale*prop_radius; N=number_circle_elems)
        end
        band!(ax, lift(first, prop_cap2_obs), lift(last, prop_cap2_obs); style_prop..., color=col_prop, rasterize=true)
    end
    return nothing
end
