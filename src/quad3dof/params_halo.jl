#= Parameter Structures and Functions.

Author: Samuel Buckner (UW-ACL)
=#

# ..:: Quadcopter Object ::..

"""
`Quad3DoFHaloParams` holds the quadcopter parameters.
"""
mutable struct Quad3DoFHaloParams{TF,TI} # TF: Type of Float, TI: Type of Int
    # >> Environmental parameters <<
    g::Vector{TF}                   # [m/s²] Acceleration due to gravity
    ρ::TF                           # [kg/m^3] Air density

    # >> Vehicle parameters <<
    C_d::TF                         # Linear drag coefficient
    S_A::TF                         # Average frontal surface area
    n_rotor::TI                     # Number of quadcopter rotors
    mass::TF                        # [kg] Mass of params
    ρ_min::TF                       # [N] Minimum thrust
    ρ_max::TF                       # [N] Maximum thrust
    drag_term_enabled::Bool         # Indicate if drag term should be enabled in dynamics

    # >> Constraint parameters <<
    ϵ_subopt::TF                    # Global suboptimality tolerance for all targets
    γ_gs::TF                        # [rad] Maximum approach angle
    γ_p::TF                         # [rad] Maximum pointing angle
    v_min_V::TF                     # [m/s] Minimum vertical velocity
    v_max_V::TF                     # [m/s] Maximum vertical velocity
    v_max_L::TF                     # [m/s] Maximum lateral velocity
    n_obstacles::TI                 # Number of obstacles
    R_obstacles::Vector{TF}         # [m] Radii of obstacles
    p_obstacles::Matrix{TF}         # [m] Positions of obstacles
    H_obstacles::Vector{Matrix{TF}} # Ellipse geometry

    # >> HALO-specific target parameters <<
    n_targs_min::TI                 # Minimum number of targets
    n_targs_max::TI                 # Maximum number of targets
    R_targs::Vector{TF}             # [m] Current bounding radii of all targets
    R_targs_min::TF                 # [m] Minimum safe bounding radius for a target
    p_targs::Dict                   # Hyperparameters for each target (pcd, prox_veh, prox_clust, µ_99)
    w_des::Vector{TF}               # Desirability score weights (pcd, prox_veh, prox_clust, µ_99, R_targs)

    # >> Algorithm parameters <<
    a::AlgorithmParams
    w_obj_decay_factor::TF          # Objective decay factor per PTR iteration
end

# ..:: Default Quad3DoFHaloParams Constructor ::..

function Quad3DoFHaloParams()::Quad3DoFHaloParams{CReal,Int}
    # >> Environmental parameters <<
    g = -9.807*e_z
    ρ = 1.225

    # Rotor Parameters (not stored)
    # (See line 21 of https://github.com/microsoft/AirSim/blob/master/AirLib/include/vehicles/multirotor/RotorParams.hpp)
    C_T = 0.109919
    RPM_max = 6396.667
    d_prop = 0.2286

    # >> Vehicle parameters <<
    C_D = 1.3/4.0
    S_A = .18*.11 # overhead rectangular area assuming vehicle's velocity is mostly aligned with body -Z, not including arms
    n_rotor = 4
    mass = 1
    ρ_min = 5                                             # [N] AirSim throttle lower bound default
    ρ_max = n_rotor * C_T * ρ * (RPM_max/60)^2 * d_prop^4 # [N] Max physical thrust of single engine
    drag_term_enabled = true

    # >> Constraint parameters <<
    ϵ_subopt = 1e-2
    γ_gs = 89 * DEG_2_RAD
    γ_p = 45 * DEG_2_RAD
    v_max_V = 1e-3
    v_min_V = -10
    v_max_L = 10

    # Obstacle and boundary parameters 
    # (defaults to empty, scenario-specific)
    n_obstacles = 0
    R_obstacles = CVector(undef,0)
    p_obstacles = CMatrix(undef,0,0)
    H_obstacles = Vector{CMatrix}(undef,0)

    # >> Algorithm parameters <<
    a = AlgorithmParams()
    a.nx = 7 # (position, velocity, thrust 2-norm)
    a.nu = 3 # (thrust)
    a.z0 = Inf * ones(a.nx) # empty initial state (to be populated with current state)
    a.u0 = Inf * ones(a.nu) # empty initial control (to be populated with current control)
    w_obj_decay_factor = 1.4

    # SCP parameters
    a.ctcs_enabled = true
    a.warmstart_method = "single" # types: (linear, single, ddto)
    a.w_obj_sing = 1.
    a.w_ctrl = 100.
    a.w_trust = 1.
    a.w_buff = a.w_ctrl
    a.ϵ_ctrl = 5e-3
    a.ϵ_buff = 5e-3
    a.ϵ_trust = 5e-3
    a.scp_iters = 50

    # Time dilation & discretization
    a.N = 20
    a.ToF_min = 10.
    a.ToF_max = 60.
    a.Δt_min = .5*a.ToF_min/(a.N-1)
    a.Δt_max =  2*a.ToF_max/(a.N-1)
    a.gss_cvx = true
    a.Δt_cvx = (a.Δt_min + a.Δt_max)/2.
    a.differentiator = "forwarddiff"

    # >> HALO-specific parameters <<
    n_targs_min = 2
    n_targs_max = 4
    R_targs_min = 1.
    R_targs = CVector(undef, a.n_targs)
    p_targs = Dict(
        "pcd" => CVector(undef, a.n_targs),         # Point cloud density
        "prox_veh" => CVector(undef, a.n_targs),    # Proximity of landing site to vehicle
        "prox_clust" => CVector(undef, a.n_targs),  # Proximity to other landing sites ("cluster proximity")
        "µ_99" => CVector(undef, a.n_targs),        # 99th percentile uncertainty
    )
    w_des = [0,0,1,0,0] # Weights for: [pcd, prox_veh, prox_clust, µ_99, R_targs]

    params = Quad3DoFHaloParams{CReal,Int}(
        g,
        ρ,
        C_D,
        S_A,
        n_rotor,
        mass,
        ρ_min,
        ρ_max,
        drag_term_enabled,
        ϵ_subopt,
        γ_gs,
        γ_p,
        v_min_V,
        v_max_V,
        v_max_L,
        n_obstacles,
        R_obstacles,
        p_obstacles,
        H_obstacles,
        n_targs_min,
        n_targs_max,
        R_targs,
        R_targs_min,
        p_targs,
        w_des,
        a,
        w_obj_decay_factor,
    )

    return params
end