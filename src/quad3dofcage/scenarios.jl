function scenario_toy1!(quad::Params)

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
    quad.N_targs = fill(21, quad.n_targs)
    quad.λ_targs = [3, 2, 4, 1]
    quad.T_targs = 1:quad.n_targs
    quad.ϵ_targs = fill(eps, quad.n_targs)

end