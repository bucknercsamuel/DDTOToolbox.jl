#= DDTO for double integrator landing -- Parameter Structures and Functions.

Author: Samuel Buckner (UW-ACL)
=#

# ..:: Quadcopter Object ::..

"""
`Params` holds the quadcopter parameters.
"""
mutable struct Params

    # >> Environmental parameters <<
    g::CVector     # [m/s²] Acceleration due to gravity
    ρ::CReal       # [kg/m^3] Air density

    # >> Params parameters <<
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

    # >> Target conditions <<
    n_targs::Int          # Current number of targets
    zf_targs::CMatrix     # [m] Terminal state of each target
    λ_targs::Vector{Int}  # Order of target rejection
    T_targs::Vector{Int}  # Tag for each target
    ϵ_targs::CVector      # Optimality tolerances

    # >> Dynamics <<
    n::Int         # Number of states
    m::Int         # Number of inputs
    A_c::CMatrix   # Continuous-time dynamics A matrix
    B_c::CMatrix   # Continuous-time dynamics B matrix
    p_c::CVector   # Continuous-time dynamics p vector

    # >> SCP Params <<
    w_ctrl::CReal         # Virtual control penalty weight
    w_buff::CReal         # Virtual buffer penalty weight
    w_trust::CReal        # Trust region penalty weight
    ϵ_ctrl::CReal         # Convergence threshold for virtual control penalty
    ϵ_buff::CReal         # Convergence threshold for virtual buffer penalty
    ϵ_trust::CReal        # Convergence threshold for trust region penalty
    scp_iters::Int        # Number of SCP subproblem iterations

    # >> Time Discretization <<
    free_final_time::Bool # Choose wether to use fixed- or free-final-time
    disc::Int             # Discretization hold order (currently can either choose 0 or 1)

    # Fixed-final-time
    Δt::CReal             # Discretization time-step
    N_targs::Vector{Int}  # Horizon length of each target

    # Free-final-time
    N_fft::Int            # Number of nodes (for all targets)
    τ::CVector            # [s] Normalized time grid
    Δτ::CVector           # [s] Vector diff on τ
    Δt_min::CReal         # [s] Minimum physical time step
    Δt_max::CReal         # [s] Maximum physical time step
    s_min::CReal          # [s] Minimum time mapping derivative value
    s_max::CReal          # [s] Maximum time mapping derivative value
    ToF_max::CReal        # [s] Maximum physical time-of-flight for all targets    

    # >> Other <<
    τ_max::Int     # Artificial maximum deferrability
end

"""
    Params()

Constructor for the quadcopter parameters.

# Returns
- `params`: the quadcopter definition.
"""
function Params()::Params

    # >> Environmental parameters <<
    g = -9.81*e_z
    ρ = 1.225

    # >> Params parameters <<
    n_rotor = 4
    mass = 0.35
    ρ_min = 1.0
    ρ_max = 7.0
    x_arena_lims = CVector([-4.5,+4.5])
    y_arena_lims = CVector([-2.5,+2.5])
    z_arena_lims = CVector([-2,+0])

    # >> Constraint parameters <<
    γ_p = 45 * DEG_2_RAD
    v_max_V = 0
    v_max_L = 5

    # >> Dynamics <<
    A_c = CMatrix([
        zeros(3,3) I(3);
        zeros(3,3) zeros(3,3)
    ])
    B_c = CMatrix([
        zeros(3,3) zeros(3);
        I(3)/mass  zeros(3)
    ])
    p_c = CVector(vcat(zeros(3),g))
    n,m = size(B_c)

    # Default scenario parameters
    # >> Obstacle parameters <<
    R_obstacles = [
        0.3,
        0.5,
        0.6,
        0.2,
        0.2
    ] # Radii of all circular obstacles
    n_obstacles = length(R_obstacles) # Number of obstacles
    p_obstacles = hcat( # Positions of circular obstacless
       -1.25*e_x + 0.5*e_y - 1*e_z,
        0*e_x    + 0*e_y   - 1*e_z,
        1*e_x    + 1*e_y   - 1*e_z,
        1*e_x    - 0.7*e_y - 1*e_z,
        1.8*e_x  + 0.3*e_y - 1*e_z,
    )
    H_obstacles = repeat([I(3)],n_obstacles)

    # >> Initial condition state <<
    r0 = -3*e_x + 1*e_y - 1*e_z
    v0 =  0*e_x + 0*e_y + 0*e_z
    z0 = [r0;v0]

    # >> Target conditions <<
    n_targs = 3
    rf_targs = hcat(
         0.0*e_x - 1.5*e_y - 1*e_z, # Target 1
         2.0*e_x - 1.0*e_y - 1*e_z, # Target 2
         3.0*e_x + 1*e_y   - 1*e_z, # Target 3
    ) 
    vf_targs = hcat(
        0*e_x + 0*e_y + 0*e_z, # Target 1
        0*e_x + 0*e_y + 0*e_z, # Target 2
        0*e_x + 0*e_y + 0*e_z, # Target 3
    )
    zf_targs = vcat(rf_targs,vf_targs)
    λ_targs = [1, 2, 3]
    T_targs = 1:n_targs
    ϵ_targs = CVector([0.2, 0.2, 0.2])

    # >> SCP Params <<
    w_ctrl = 1e7
    w_buff = 1e-2
    w_trust = 1e3
    ϵ_ctrl = 1e-2
    ϵ_buff = 1e-2
    ϵ_trust = 1e-2
    scp_iters = 10

    # >> Time Discretization <<
    free_final_time = true
    disc = 0

    # Fixed-final-time
    Δt = 0.5
    N_targs = [21, 21, 21]

    # Free-final-time
    N_fft = 11
    τ = CVector(range(0, stop=1, length=N_fft))
    Δτ = diff(τ)
    Δt_min = 0.01
    Δt_max = 2
    s_min = 0.01
    s_max = 10
    ToF_max = 10

    # >> Other <<
    τ_max = 1e10

    # >> Make params object <<
    params = Params(
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
        n_targs,
        zf_targs,
        λ_targs,
        T_targs,
        ϵ_targs,
        n,
        m,
        A_c,
        B_c,
        p_c,
        w_ctrl,
        w_buff,
        w_trust,
        ϵ_ctrl,
        ϵ_buff,
        ϵ_trust,
        scp_iters,
        free_final_time,
        disc,
        Δt,
        N_targs,
        N_fft,
        τ,
        Δτ,
        Δt_min,
        Δt_max,
        s_min,
        s_max,
        ToF_max,
        τ_max,
    )

    return params
end

# ..:: Post-processed Solution Structure ::..
"""
`ProcessedSolution` stores post-processed data from the `Solution` struct specific to this problem
`ProcessedBranchSolution` forms this into an equivalent `BranchSolution`-style object
"""

mutable struct ProcessedSolution
    t::CVector       # [s] Time vector
    r::CMatrix       # [m] Position trajectory
    v::CMatrix       # [m/s] Velocity trajectory
    T::CMatrix       # [m/s^2] Thrust vector
    Γ::CVector       # [m/s^2] Slack thrust magnitude
    T_nrm::CVector   # [N] Thrust magnitude
    γ::CVector       # [rad] Pointing angle
    cost::CReal      # Optimization's optimal cost
end

mutable struct ProcessedBranchSolution
    sol::ProcessedSolution  # Contains the `ProcessedSolution` for the branch
    cost_dd::CReal          # Cost for deferred decision
    idx_dd::Int             # Deferred decision branch point index
end

# ..:: Function to convert raw `Solution` data for each branch to a `ProcessedSolution` ::..

function process_solutions(branchsolutions::Vector{BranchSolution})::Vector{ProcessedBranchSolution}
    processed_branchsolutions = Vector{ProcessedBranchSolution}(undef, length(branchsolutions))
    
    for k = 1:length(branchsolutions)
        # Obtain raw data from solution
        cost = branchsolutions[k].sol.cost
        t = branchsolutions[k].sol.t
        x = branchsolutions[k].sol.x
        u = branchsolutions[k].sol.u
        cost_dd = branchsolutions[k].cost_dd
        idx_dd = branchsolutions[k].idx_dd

        # Post-processing
        r = x[1:3,:]
        v = x[4:6,:]
        T = u[1:3,:]
        Γ = u[4,:]
        T_nrm = CVector([norm(T[:,i],2) for i=1:length(Γ)])
        γ = CVector([acos(dot(T[:,k],e_z)/norm(T[:,k],2)) for k=1:length(Γ)])

        processed_solution = ProcessedSolution(t,r,v,T,Γ,T_nrm,γ,cost)
        processed_branchsolutions[k] = ProcessedBranchSolution(processed_solution, cost_dd, idx_dd)
    end

    return processed_branchsolutions
end