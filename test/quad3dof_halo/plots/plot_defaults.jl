using CairoMakie
using Colors
using LinearAlgebra
include("../../utils/plot_utils.jl")

# Generic styling
markersize = 15
style2D_dt = Dict(:color=>:gray, :marker=>:circle, :markersize=>markersize, :strokecolor=>:black, :strokewidth=>3)
style2D_ct = Dict(:color=>:black, :linewidth=>3)
style2D_ct_ddto = Dict(:color=>:black, :linewidth=>3)
style3D_dt = Dict(:color=>:gray, :marker=>:circle, :markersize=>markersize, :strokecolor=>:black, :strokewidth=>3)
style3D_ct = Dict(:color=>:black, :linewidth=>5, :overdraw=>true)
style3D_ct_ddto = Dict(:color=>:black, :linewidth=>5)
style3D_ground_base = Dict(:color=>bright_color(:gray95), :transparency=>false, :alpha=>1)
style3D_ground_base_frame = Dict(:color=>bright_color(:gray90))

# Themes
theme2d = theme_latexfonts()
theme2d.fontsize = 20
theme3d = theme_latexfonts()
theme3d.fontsize = 20

# Figure saving setup
fig_path = "quad3dof_halo\\figures"
fig_ext = ".svg"

function generate_custom_colors(max_targs)
    # target_colors = range(colorant"magenta", stop=colorant"cyan", length=max_targs)
    target_colors = range(HSV(45,1,1), stop=HSV(-360,1,1), length=max_targs)

    return target_colors
end

# Reference vector for drone body-frame construction (used by plot_drone helpers)
b_ref = rand(3)
b_ref = b_ref / norm(b_ref)
