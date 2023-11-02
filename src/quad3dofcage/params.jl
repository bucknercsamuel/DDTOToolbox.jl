#= DDTO for double integrator landing -- Parameter Structures and Functions.

Author: Samuel Buckner (UW-ACL)
=#

# ..:: Quadcopter Object ::..

"""
`Params` holds the quadcopter parameters.
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

    # >> Target conditions <<
    n_targs::Int          # Current number of targets
    zf_targs::CMatrix     # [m] Terminal state of each target
    λ_targs::Vector{Int}  # Order of target rejection
    T_targs::Vector{Int}  # Tag for each target
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
    τ::CVector            # [s] Normalized time grid
    Δτ::CVector           # [s] Vector diff on τ
    Δt_min::CReal         # [-] Minimum wall time step
    Δt_max::CReal         # [-] Maximum wall time step
    s_min::CReal          # [-] Minimum time dilation factor
    s_max::CReal          # [-] Maximum time dilation factor
    ToF_max::CReal        # [s] Maximum physical time-of-flight for all targets    
    disc::Int             # Discretization hold order (currently can either choose 0 or 1)

    # >> Other <<
    τ_max::Int     # Artificial maximum deferrability
end

# ..:: Post-processed Solution Structure ::..
mutable struct Quad3DoFCageSolution
    t::CVector       # [s] Time vector
    r::CMatrix       # [m] Position trajectory
    v::CMatrix       # [m/s] Velocity trajectory
    T::CMatrix       # [m/s^2] Thrust vector
    s::CVector       # Time dilation factor (if using free-final-time)
    T_nrm::CVector   # [N] Thrust magnitude
    γ::CVector       # [rad] Pointing angle
    cost::CReal      # Optimization's optimal cost
end

mutable struct Quad3DoFCageBranchSolution
    sol::Quad3DoFCageSolution  # Contains the `ProcessedSolution` for the branch
    cost_dd::CReal          # Cost for deferred decision
    idx_dd::Int             # Deferred decision branch point index
end

# ..:: Constructors for empty `*Solution` structs ::..

function EmptyQuad3DoFCageSolution()::Quad3DoFCageSolution

    t = CVector(undef,0)
    r = CMatrix(undef,0,0)
    v = CMatrix(undef,0,0)
    T = CMatrix(undef,0,0)
    s = CVector(undef,0)
    T_nrm = CVector(undef,0)
    γ = CVector(undef,0)
    cost = Inf

    return ProcessedSolution(t,r,v,T,s,T_nrm,γ,cost)
end

# ..:: Function to convert raw `Solution` data for each branch to a `Quad3DoFCageSolution` ::..

function process_solutions(branchsolutions::Vector{BranchSolution}, params::Quad3DoFCageParams)::Vector{Quad3DoFCageBranchSolution}
    processed_branchsolutions = Vector{Quad3DoFCageBranchSolution}(undef, length(branchsolutions))
    
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
        s = u[4,:]
        T_nrm = CVector([norm(T[:,i],2) for i=1:length(T[1,:])])
        γ = CVector([acos(dot(T[:,k],e_z)/norm(T[:,k],2)) for k=1:length(T[1,:])])

        processed_solution = Quad3DoFCageSolution(t,r,v,T,s,T_nrm,γ,cost)
        processed_branchsolutions[k] = Quad3DoFCageBranchSolution(processed_solution, cost_dd, idx_dd)
    end

    return processed_branchsolutions
end