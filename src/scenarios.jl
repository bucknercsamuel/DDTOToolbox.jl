function scenario_onr_demo!(quad::Lander)

    # High-level settings
    eps = 0.0035  # Accepted level of suboptimality
    obs_rad = 0.6 # [m] Radius of all cylindrical obstacles
    height = 1 # [m] Height of the maneuver

    filename = "input_demo_config/states.mat"
    input = matread(filename)

    n_targs = size(input["r_tar"])[1]
    n_obs   = size(input["r_obs"])[1]

    # >> Initial condition state <<
    quad.r0 = vec(input["r_drone"])
    quad.v0 = zeros(3)

    # >> Obstacle parameters <<
    quad.n_obstacles = n_obs
    quad.R_obstacles = fill(obs_rad,n_obs)
    quad.p_obstacles = input["r_obs"]
    quad.H_obstacles = repeat([I(3)],quad.n_obstacles)

    # >> Target conditions <<
    quad.n_targs = n_targs
    quad.rf_targs = input["r_tar"]
    quad.vf_targs = zeros(3,n_targs)
    quad.N_targs = fill(41,n_targs)
    quad.λ_targs = [3,1,2]
    quad.T_targs = 1:n_targs
    quad.ϵ_targs = fill(eps, n_targs)

    # >> Adjust all z-components to height << 
    quad.r0[3] = -height
    quad.p_obstacles[3,:] .= -height
    quad.rf_targs[3,:] .= -height

end

function scenario_toy1!(quad::Lander)

     # High-level settings
     eps = 0.2  # Accepted level of suboptimality
     obs_rad = 0.6 # [m] Radius of all cylindrical obstacles
     height = 1 # [m] Height of the maneuver

    # >> Obstacle parameters <<
    quad.n_obstacles = 14 # Number of obstacles
    quad.R_obstacles = fill(obs_rad, quad.n_obstacles) # Radii of all circular obstacles
    quad.p_obstacles = hcat( # Positions of circular obstacless
       -3*e_x + 1.5*e_y - height*e_z,
       -1*e_x + 1.5*e_y - height*e_z,
       +1*e_x + 1.5*e_y - height*e_z,
       +3*e_x + 1.5*e_y - height*e_z,
       -2*e_x + 0.5*e_y - height*e_z,
       +0*e_x + 0.5*e_y - height*e_z,
       +2*e_x + 0.5*e_y - height*e_z,
       -3*e_x - 0.5*e_y - height*e_z,
       -1*e_x - 0.5*e_y - height*e_z,
       +1*e_x - 0.5*e_y - height*e_z,
       +3*e_x - 0.5*e_y - height*e_z,
       -2*e_x - 1.5*e_y - height*e_z,
       +0*e_x - 1.5*e_y - height*e_z,
       +2*e_x - 1.5*e_y - height*e_z,
    )
    quad.H_obstacles = repeat([I(3)],quad.n_obstacles)

    # >> Dynamics <<
    Δt = 0.2

    # >> Initial condition state <<
    quad.r0 = -3*e_x + 0.5*e_y - height*e_z
    quad.v0 =  0*e_x + 0*e_y + 0*e_z

    # >> Target conditions <<
    quad.n_targs = 4
    quad.rf_targs = hcat(
        -1*e_x - 1.5*e_y - height*e_z,
        +3*e_x - 1.5*e_y - height*e_z,
        +3*e_x + 0.5*e_y - height*e_z,
        +0*e_x + 1.5*e_y - height*e_z,
    )
    quad.vf_targs = zeros(3,quad.n_targs)
    quad.N_targs = fill(51, quad.n_targs)
    quad.λ_targs = [3, 2, 4, 1]
    quad.T_targs = 1:quad.n_targs
    quad.ϵ_targs = fill(eps, quad.n_targs)

end