#= DDTO for double integrator landing -- Parameter Structures and Functions.

Author: Samuel Buckner (UW-ACL)
=#

# ..:: Quadcopter Object ::..

"""
`Quad3DoFCageParams` holds the quadcopter parameters.
"""
mutable struct Quad3DoFCageParams

    # >> Environmental parameters <<
    g::CVector     # [m/s²] Acceleration due to gravity
    ρ::CReal       # [kg/m^3] Air density

    # >> Vehicle parameters <<
    n_rotor::Int   # Number of quadcopter rotors
    mass::CReal    # [kg] Mass of params
    ρ_min::CReal   # [N] Minimum thrust
    ρ_max::CReal   # [N] Maximum thrust

    # >> Constraint parameters <<
    γ_p::CReal                   # [rad] Maximum pointing angle
    v_max_V::CReal               # [m/s] Maximum vertical velocity
    v_max_L::CReal               # [m/s] Maximum lateral velocity
    n_obstacles::Int             # Number of obstacles
    R_obstacles::CVector         # [m] Radii of obstacles
    p_obstacles::CMatrix         # [m] Positions of obstacles
    H_obstacles::Vector{CMatrix} # Ellipse geometry
    x_arena_lims::CVector        # [m] X limits for the indoor netted arena
    y_arena_lims::CVector        # [m] Y limits for the indoor netted arena
    z_arena_lims::CVector        # [m] Z limits for the indoor netted arena
    z0::CVector                  # [m] Initial state
    nx::Int                      # [-] Number of states
    nu::Int                      # [-] Number of controls

    # >> DDTO target conditions <<
    n_targs::Int          # Current number of targets
    zf_targs::CMatrix     # [m] Terminal state of each target
    λ_targs::Vector{Int}  # Order of target rejection
    T_targs::Vector{Int}  # Tag for each target
    τ_targs::Vector{Int}  # Deferrability index allocation (in order specified by λ_targs) -- set automatically in `solve_tree_ddto`
    α_targs::CVector      # Relative weight for deferrability of each target
    ϵ_targs::CVector      # Optimality tolerances

    # >> SCP Params <<
    w_obj::CReal          # Objective penalty weight
    w_ctrl::CReal         # Virtual control penalty weight
    w_buff::CReal         # Virtual buffer penalty weight
    w_trust::CReal        # Trust region penalty weight
    ϵ_ctrl::CReal         # Convergence threshold for virtual control penalty
    ϵ_buff::CReal         # Convergence threshold for virtual buffer penalty
    ϵ_trust::CReal        # Convergence threshold for trust region penalty
    scp_iters::Int        # Number of SCP subproblem iterations

    # >> Time dilation & discretization <<
    N::Int                # Number of nodes (for all targets)
    Δt_min::CReal         # [-] Minimum wall time step
    Δt_max::CReal         # [-] Maximum wall time step
    ToF_max::CReal        # [s] Maximum physical time-of-flight for all targets    
    disc::Int             # Discretization hold order (currently can either choose 0 or 1)

    # >> Affine scaling parameters <<
    Sx::CMatrix                  # Scaling transformation matrix for state "x"
    sx::CVector                  # Scaling affine vector for state "x"
    Su::CMatrix                  # Scaling transformation matrix for state "u"
    su::CVector                  # Scaling affine vector for state "u"
end

# ..:: Default Quad3DoFCageParams Constructor ::..

function Quad3DoFCageParams()::Quad3DoFCageParams
    # >> Environmental parameters <<
    g = -9.81*e_z
    ρ = 1.225

    # >> Vehicle parameters <<
    n_rotor = 4
    mass = 0.35
    ρ_min = 1.0
    ρ_max = 7.0
    x_arena_lims = CVector([-4.5,+4.5])
    y_arena_lims = CVector([-2.5,+2.5])
    z_arena_lims = CVector([-2,+0])

    # >> Constraint parameters <<
    γ_p = 45 * DEG_2_RAD
    v_max_V = 0.
    v_max_L = 2.
    nx = 7 # (position, velocity, thrust 2-norm)
    nu = 4 # (thrust, time dilation)

    # Obstacle and boundary parameters 
    # (defaults to empty, scenario-specific)
    n_obstacles = 0
    R_obstacles = CVector(undef,0)
    p_obstacles = CMatrix(undef,0,0)
    H_obstacles = Vector(undef,0)
    n_targs = 0
    z0 = CVector(undef,0)
    zf_targs = CMatrix(undef,0,0)
    λ_targs = Array{Int}(undef,0)
    T_targs = Array{Int}(undef,0)
    τ_targs = Array{Int}(undef,0)
    α_targs = CVector(undef,0)
    ϵ_targs = CVector(undef,0)
    
    # >> SCP Params <<
    w_obj = 1
    w_ctrl = 1e7
    w_buff = 1e-2
    w_trust = 1e3
    ϵ_ctrl = 1e-2
    ϵ_buff = 1e-2
    ϵ_trust = 1e-2
    scp_iters = 10

    # >> Time dilation & discretization <<
    N = 11
    Δt_min = 0.01
    Δt_max = 2.
    ToF_max = 10.
    disc = 1

    # >> Affine scaling parameters <<
    rmin = [x_arena_lims[1]; y_arena_lims[1]; z_arena_lims[1]]
    rmax = [x_arena_lims[2]; y_arena_lims[2]; z_arena_lims[2]]
    Δτ = 1/(N-1)
    Sx,sx = scaling_matrices([rmin; -v_max_L*ones(3); 0], [rmax; v_max_L*ones(3); ToF_max*ρ_max])
    Su,su = scaling_matrices([-ρ_max*ones(3); 0], [ρ_max*ones(3); Δt_max/Δτ])

    # >> Make params object <<
    params = Quad3DoFCageParams(
        g,
        ρ,
        n_rotor,
        mass,
        ρ_min,
        ρ_max,
        γ_p,
        v_max_V,
        v_max_L,
        n_obstacles,
        R_obstacles,
        p_obstacles,
        H_obstacles,
        x_arena_lims,
        y_arena_lims,
        z_arena_lims,
        z0,
        nx,
        nu,
        n_targs,
        zf_targs,
        λ_targs,
        T_targs,
        τ_targs,
        α_targs,
        ϵ_targs,
        w_obj,
        w_ctrl,
        w_buff,
        w_trust,
        ϵ_ctrl,
        ϵ_buff,
        ϵ_trust,
        scp_iters,
        N,
        Δt_min,
        Δt_max,
        ToF_max,
        disc,
        Sx,
        sx,
        Su,
        su
    )

    return params
end

# ..:: Sample Scenario (needed for precompile purposes) ::..
function Quad3DoFCageSampleScenario()
    # Load default params first
    params = Quad3DoFCageParams()

    # High-level settings
    eps = .1  # Accepted level of suboptimality
    obs_rad = 0.6 # [m] Radius of all cylindrical obstacles
    height = 1 # [m] Height of the maneuver

    # >> Obstacle parameters <<
    params.n_obstacles = 1 # Number of obstacles
    params.R_obstacles = fill(obs_rad, params.n_obstacles) # Radii of all circular obstacles
    params.p_obstacles = hcat( # Positions of circular obstacless
       +2*e_x + 0.5*e_y - height*e_z,
       -2*e_x + 0.5*e_y - height*e_z,
       +0*e_x - 0.5*e_y - height*e_z,
    )
    params.H_obstacles = repeat([I(3)],params.n_obstacles)

    # >> Initial condition state <<
    r0 = -3*e_x + 0.5*e_y - height*e_z
    v0 =  0*e_x + 0*e_y + 0*e_z
    params.z0 = [r0;v0;0]

    # >> Target conditions <<
    N = 10 # number of nodes for each targ
    params.n_targs = 4
    rf_targs = hcat(
        -1*e_x - 1.5*e_y - height*e_z,
        +3*e_x - 1.5*e_y - height*e_z,
        +3*e_x + 0.5*e_y - height*e_z,
        +0*e_x + 1.5*e_y - height*e_z,
    )
    vf_targs = zeros(3,params.n_targs)
    params.zf_targs = vcat(rf_targs,vf_targs,Inf*ones(1,params.n_targs)) # Inf: not constraining this state
    params.λ_targs = [3, 2, 4, 1]
    params.T_targs = 1:params.n_targs
    params.α_targs = [1,1,1000,1]
    params.ϵ_targs = fill(eps, params.n_targs)

    # >> SCP Params <<
    params.w_obj = 1e0
    params.w_ctrl = 1e5
    params.w_buff = 1e4
    params.w_trust = 1e3
    params.ϵ_ctrl = 1e-2
    params.ϵ_buff = 1e-2
    params.ϵ_trust = 1e-2
    params.scp_iters = 15

    # >> Time dilation & discretization <<
    params.N = N
    params.Δt_min = 0.001
    params.Δt_max = 1
    params.ToF_max = 10

    return params
end

# ..:: Post-processed Solution Structure ::..
mutable struct Quad3DoFCageSolution
    τ::CVector       # [s] Dilated Time Vector
    t::CVector       # [s] Wall-clock Time vector
    r::CMatrix       # [m] Position trajectory
    v::CMatrix       # [m/s] Velocity trajectory
    T::CMatrix       # [m/s^2] Thrust vector
    s::CVector       # Time dilation factor (if using free-final-time)
    T_nrm::CVector   # [N] Thrust magnitude
    ∫T::CVector      # [N] Cumulative thrust
    γ::CVector       # [rad] Pointing angle
    cost::CReal      # Optimization's optimal cost
end

mutable struct Quad3DoFCageDDTOSolution
    targs::Vector{Quad3DoFCageSolution}  # Contains the `Quad3DoFCageSolution` for each target
end

# ..:: Constructors for empty `*Solution` structs ::..

function EmptyQuad3DoFCageSolution()::Quad3DoFCageSolution

    τ = CVector(undef,0)
    t = CVector(undef,0)
    r = CMatrix(undef,0,0)
    v = CMatrix(undef,0,0)
    T = CMatrix(undef,0,0)
    s = CVector(undef,0)
    T_nrm = CVector(undef,0)
    ∫T = CVector(undef,0)
    γ = CVector(undef,0)
    cost = Inf

    return Quad3DoFCageSolution(τ,t,r,v,T,s,T_nrm,∫T,γ,cost)
end

function EmptyQuad3DoFCageDDTOSolution(n_targs)::Quad3DoFCageDDTOSolution
    targs = Vector{Quad3DoFCageSolution}(undef, n_targs)
    for j = 1:n_targs
        targs[j] = EmptyQuad3DoFCageSolution()
    end
    return Quad3DoFCageDDTOSolution(targs)
end

# ..:: Function to convert raw `Solution` data for each branch to a `Quad3DoFCageSolution` ::..

function process_solutions(solution::DDTOSolution, params::Quad3DoFCageParams)::Quad3DoFCageDDTOSolution
    solution_proc = EmptyQuad3DoFCageDDTOSolution(params.n_targs)
    for k = 1:params.n_targs
        # Obtain raw data from solution
        cost = solution.targs[k].cost
        τ = solution.targs[k].t
        x = solution.targs[k].x
        u = solution.targs[k].u
        N = size(x,2)
        Δτ = 1 / (N-1)
        if ~isempty(u)
            t = time_dilation_control_to_wall_clock_time(u[end,:], τ, params.disc)
        else
            t = 0
        end

        # Post-processing
        r = x[1:3,:]
        v = x[4:6,:]
        ∫T = x[7,:]
        T = u[1:3,:]
        s = u[4,:]
        T_nrm = CVector([norm(T[:,i],2) for i=1:length(T[1,:])])
        γ = CVector([acos(dot(T[:,k],e_z)/norm(T[:,k],2)) for k=1:length(T[1,:])])
        solution_proc.targs[k] = Quad3DoFCageSolution(τ,t,r,v,T,s,T_nrm,∫T,γ,cost)
    end

    return solution_proc
end