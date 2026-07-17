#=
Parameter structures and defaults for 3-DOF quadcopter DDTO in an indoor
netted-cage arena (obstacles, cage bounds, constant-altitude flight).
=#

# ..:: Quadcopter Object ::..

"""
    Quad3DoFCageParams{TF,TI}

Vehicle, environment, constraint, and algorithm parameters for the indoor
netted-cage 3-DOF quadcopter DDTO scenario.
"""
mutable struct Quad3DoFCageParams{TF,TI}
    # >> Environmental parameters <<
    g::Vector{TF} # [m/s²] Acceleration due to gravity
    ρ::TF         # [kg/m^3] Air density

    # >> Vehicle parameters <<
    n_rotor::TI # Number of quadcopter rotors
    mass::TF    # [kg] Mass of params
    ρ_min::TF   # [N] Minimum thrust
    ρ_max::TF   # [N] Maximum thrust
    drag_term_enabled::Bool         # Indicate if drag term should be enabled in dynamics

    # >> Constraint parameters <<
    γ_p::TF                         # [rad] Maximum pointing angle
    v_max_V::TF                     # [m/s] Maximum vertical velocity
    v_max_L::TF                     # [m/s] Maximum lateral velocity
    h_constant::TF                  # [m] Fixed constant altitude
    n_obstacles::TI                 # Number of obstacles
    R_obstacles::Vector{TF}         # [m] Radii of obstacles
    p_obstacles::Matrix{TF}         # [m] Positions of obstacles
    H_obstacles::Vector{Matrix{TF}} # Ellipse geometry
    x_arena_lims::Vector{TF}        # [m] X limits for the indoor netted arena
    y_arena_lims::Vector{TF}        # [m] Y limits for the indoor netted arena
    z_arena_lims::Vector{TF}        # [m] Z limits for the indoor netted arena
    cage_bounds_enabled::Bool       # Determine if we should enable cage bound constraints 

    # >> Algorithm parameters <<
    a::AlgorithmParams
    w_obj_decay_factor::TF          # Objective decay factor per PTR iteration
end

# ..:: Default Quad3DoFCageParams Constructor ::..

"""
    Quad3DoFCageParams() -> Quad3DoFCageParams{CReal,Int}

Construct cage-scenario parameters with toolbox defaults (mass, thrust limits,
arena bounds, SCP settings).

# Arguments
- none

# Returns
- Default `Quad3DoFCageParams` instance ready for scenario customization.
"""
function Quad3DoFCageParams()::Quad3DoFCageParams{CReal,Int}
    # >> Environmental parameters <<
    g = -9.81*e_z
    ρ = 1.225

    # >> Vehicle parameters <<
    n_rotor = 4
    mass = 0.35
    ρ_min = 1.0
    ρ_max = 7.0
    # x_arena_lims = CVector([-4.5,+4.5])
    # y_arena_lims = CVector([-2.5,+2.5])
    x_arena_lims = CVector([-2,+2])
    y_arena_lims = CVector([-4,+4.5])
    z_arena_lims = CVector([-2,+0])
    drag_term_enabled = false

    # >> Constraint parameters <<
    γ_p = 45 * DEG_2_RAD
    v_max_V = 0.
    v_max_L = 2.
    h_constant = 0.

    # Obstacle and boundary parameters 
    # (defaults to empty, scenario-specific)
    n_obstacles = 0
    R_obstacles = CVector(undef,0)
    p_obstacles = CMatrix(undef,0,0)
    H_obstacles = Vector{CMatrix}(undef,0)
    cage_bounds_enabled = true

    # >> Algorithm parameters <<
    a = AlgorithmParams()
    a.nx = 7 # (position, velocity, thrust 2-norm)
    a.nu = 3 # (thrust)
    w_obj_decay_factor = 1.4

    # Set initial thrust input for all scenarios to be hover condition
    a.u0 = -g*mass

    params = Quad3DoFCageParams{CReal,Int}(
        g,
        ρ,
        n_rotor,
        mass,
        ρ_min,
        ρ_max,
        drag_term_enabled,
        γ_p,
        v_max_V,
        v_max_L,
        h_constant,
        n_obstacles,
        R_obstacles,
        p_obstacles,
        H_obstacles,
        x_arena_lims,
        y_arena_lims,
        z_arena_lims,
        cage_bounds_enabled,
        a,
        w_obj_decay_factor
    )

    return params
end

# ..:: Sample Scenario (needed for precompile purposes) ::..

"""
    Quad3DoFCageSampleScenario() -> Quad3DoFCageParams

Return a fully populated sample cage scenario (obstacles, IC/TCs, SCP settings)
suitable for demos and precompilation.

# Arguments
- none

# Returns
- Fully configured `Quad3DoFCageParams` with obstacles, targets, and SCP settings.
"""
function Quad3DoFCageSampleScenario()
    # Load default params first
    params = Quad3DoFCageParams()

    # High-level settings
    eps = .01  # Accepted level of suboptimality
    obs_rad = 0.6 # [m] Radius of all cylindrical obstacles
    height = 1 # [m] Height of the maneuver

    # >> Obstacle parameters <<
    params.n_obstacles = 14 # Number of obstacles
    params.R_obstacles = fill(obs_rad, params.n_obstacles) # Radii of all circular obstacles
    params.p_obstacles = hcat( # Positions of circular obstacless
       -3*e_y + 1.5*e_x - height*e_z,
       -1*e_y + 1.5*e_x - height*e_z,
       +1*e_y + 1.5*e_x - height*e_z,
       +3*e_y + 1.5*e_x - height*e_z,
       -2*e_y + 0.5*e_x - height*e_z,
       +0*e_y + 0.5*e_x - height*e_z,
       +2*e_y + 0.5*e_x - height*e_z,
       -3*e_y - 0.5*e_x - height*e_z,
       -1*e_y - 0.5*e_x - height*e_z,
       +1*e_y - 0.5*e_x - height*e_z,
       +3*e_y - 0.5*e_x - height*e_z,
       -2*e_y - 1.5*e_x - height*e_z,
       +0*e_y - 1.5*e_x - height*e_z,
       +2*e_y - 1.5*e_x - height*e_z,
    )
    params.H_obstacles = repeat([I(3)],params.n_obstacles)

    # >> Initial condition state <<
    r0 = -3*e_y + 0.5*e_x - height*e_z
    v0 =  0*e_y + 0*e_x + 0*e_z
    params.a.z0 = [r0;v0;0]
    params.h_constant = params.a.z0[3]

    # >> Target conditions <<
    params.a.n_targs = 4
    rf_targs = hcat(
        -1*e_y - 1.5*e_x - height*e_z,
        +3*e_y - 1.5*e_x - height*e_z,
        +3*e_y + 0.5*e_x - height*e_z,
        +0*e_y + 1.5*e_x - height*e_z,
    )
    vf_targs = zeros(3,params.a.n_targs)
    params.a.zf_targs = vcat(rf_targs,vf_targs,Inf*ones(1,params.a.n_targs)) # Inf: not constraining this state
    params.a.uf_targs = repeat(params.a.u0,1,params.a.n_targs) # repeat initial input cond
    params.a.λ_targs = [1,4,2,3]
    params.a.J_targs = 1:params.a.n_targs
    params.a.α_targs = ones(params.a.n_targs)
    params.a.ϵ_targs = fill(eps, params.a.n_targs)

    # >> SCP Params <<
    params.a.ctcs_enabled = true
    params.a.w_obj_sing = .001
    params.a.w_obj_ddto = params.a.w_obj_sing/params.a.n_targs
    params.a.w_ctrl = 50
    params.a.w_buff = params.a.w_ctrl
    params.a.w_trust = 2
    params.a.ϵ_ctrl = 1e-4
    params.a.ϵ_buff = 1e-4
    params.a.ϵ_trust = 1e-4
    params.a.scp_iters = 100

    # >> Time dilation & discretization <<
    params.a.N = 10
    params.a.Δt_min = 0.2
    params.a.Δt_max = 0.7
    params.a.ToF_max = 10.

    # >> Build custom scaling matrices <<
    custom_scaling!(params)

    return params
end