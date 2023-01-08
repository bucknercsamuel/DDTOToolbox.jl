#= DDTO for double integrator landing -- Parameter Structures and Functions.

Author: Samuel Buckner (UW-ACL)
=#

using LinearAlgebra

# ..:: Quadcopter Object ::..

"""
`Lander` holds the quadcopter parameters.
"""
mutable struct Lander

    # >> Environmental parameters <<
    g::CVector     # [m/s²] Acceleration due to gravity
    ρ::CReal       # [kg/m^3] Air density

    # >> Vehicle parameters <<
    n_rotor::Int   # Number of quadcopter rotors
    mass::CReal    # [kg] Mass of vehicle
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

    # >> Initial condition state <<
    r0::CVector    # [m] Initial position
    v0::CVector    # [m] Initial velocity

    # >> Target conditions <<
    n_targs::Int          # Current number of targets
    rf_targs::CMatrix     # [m] Position of each target
    vf_targs::CMatrix     # [m] Velocity of each target
    N_targs::Vector{Int}  # Horizon length of each target
    λ_targs::Vector{Int}  # Order of target rejection
    T_targs::Vector{Int}  # Tag for each target
    ϵ_targs::CVector      # Optimality tolerances

    # >> Dynamics <<
    Δt::CReal      # Discretization time-step
    n::Int         # Number of states
    m::Int         # Number of inputs
    A_c::CMatrix   # Continuous-time dynamics A matrix
    B_c::CMatrix   # Continuous-time dynamics B matrix
    p_c::CVector   # Continuous-time dynamics p vector

    # >> SCP Params <<
    w_buff::CReal  # Virtual buffer weight
    w_trust::CReal # Trust region weight
    w_r0::CReal    # Initial position relaxation weight
    w_rf::CReal    # Final position relaxation weight
    sub_iters::Int # Number of SCP subproblem iterations
    ϵ_cvg::CReal   # Convergence threshold for SCP

    # >> Other <<
    τ_max::Int     # Artificial maximum deferrability
end

"""
    Lander()

Constructor for the quadcopter.

# Returns
- `lander`: the quadcopter definition.
"""
function Lander()::Lander

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
    v_max_V = 0
    v_max_L = 5

    # >> Dynamics <<
    Δt = 0.5
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
    N_targs = [21, 21, 21]
    λ_targs = [1, 2, 3]
    T_targs = 1:n_targs
    ϵ_targs = CVector([0.2, 0.2, 0.2])

    # >> SCP Params <<
    w_buff = 1e3
    w_trust = 1e0
    w_r0 = 1
    w_rf = 1
    sub_iters = 10
    ϵ_cvg = 1e-4

    # >> Other <<
    τ_max = max(N_targs...)
    # τ_max = 5

    # >> Make quadcopter object <<
    lander = Lander(
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
        r0,
        v0,
        n_targs,
        rf_targs,
        vf_targs,
        N_targs,
        λ_targs,
        T_targs,
        ϵ_targs,
        Δt,
        n,
        m,
        A_c,
        B_c,
        p_c,
        w_buff,
        w_trust,
        w_r0,
        w_rf,
        sub_iters,
        ϵ_cvg,
        τ_max,
    )

    return lander
end


# ..:: Solution Structures ::..

"""
`Solution` stores the optimal solution.
"""
mutable struct Solution
    # >> Raw data <<
    t::CVector        # [s] Time vector
    r::CMatrix        # [m] Position trajectory
    v::CMatrix        # [m/s] Velocity trajectory
    T::CMatrix        # [m/s^2] Thrust vector
    Γ::CVector        # [m/s^2] Slack thrust magnitude
    r0_relax::CVector # [m] Initial position relaxation
    rf_relax::CVector # [m] Terminal position relaxation
    cost::CReal       # Optimization's optimal cost

    # >> Processed data <<
    T_nrm::CVector  # [N] Thrust magnitude
    γ::CVector      # [rad] Pointing angle
end

"""
`SolutionSCP` stores the SCP-optimal solution.
"""
mutable struct SolutionSCP
    # >> Raw data <<
    t::CVector        # [s] Time vector
    r::CMatrix        # [m] Position trajectory
    v::CMatrix        # [m/s] Velocity trajectory
    T::CMatrix        # [m/s^2] Thrust vector
    Γ::CVector        # [m/s^2] Slack thrust magnitude
    ν::CMatrix        # [m] Virtual buffer variables for linearized constraints
    μ::CMatrix        # [m] Virtual buffer slack for linearized constraints
    r0_relax::CVector # [m] Initial position relaxation
    rf_relax::CVector # [m] Terminal position relaxation
    cost::CReal       # Optimization's optimal cost

    # >> Processed data <<
    T_nrm::CVector  # [N] Thrust magnitude
    γ::CVector      # [rad] Pointing angle
end

"""
`DDTOSolution` stores the optimal solution from DDTO.
"""
mutable struct DDTOSolution
    targ_sols::Vector{Solution} # Contains the `Solution` to each target
    costs_sol::CVector          # Costs for each target
    cost_dd::CReal              # Cost for deferred decision
    idx_dd::Int                 # Deferred decision branch point index
end

"""
`DDTOSolutionSCP` stores an SCP subproblem solution from DDTO.
"""
mutable struct DDTOSolutionSCP
    targ_sols::Vector{SolutionSCP} # Contains the `SolutionSCP` to each target
    costs_sol::CVector             # Costs for each target
    cost_dd::CReal                 # Cost for deferred decision
    idx_dd::Int                    # Deferred decision branch point index
    η::CVector                     # Trust region variables
end

"""
`BranchSolution` stores the solution of a branch from DDTO
"""
mutable struct BranchSolution
    sol::Solution  # Contains the `Solution` for the branch
    cost_dd::CReal # Cost for deferred decision
    idx_dd::Int    # Deferred decision branch point index
end


# ..:: Solution Initialization Functions ::..

"""
    EmptySolution()

    Constructor for an empty Solution.

# Arguments
- :out sol: empty solution.
"""
function EmptySolution()::Solution

    t = CVector(undef,0)
    r = CMatrix(undef,0,0)
    v = CMatrix(undef,0,0)
    T = CMatrix(undef,0,0)
    Γ = CVector(undef,0)
    r0_relax = CVector(undef,0)
    rf_relax = CVector(undef,0)
    cost = Inf
    T_nrm = CVector(undef,0)
    γ = CVector(undef,0)

    return Solution(t,r,v,T,Γ,r0_relax,rf_relax,cost,T_nrm,γ)
end

"""
    EmptyDDTOSolution()

Constructor for an empty DDTOSolution.

# Arguments
- :in n_targs: Number of targets for the solution
- :out sol: empty solution.
"""
function EmptyDDTOSolution(n_targs)::DDTOSolution

    targ_sols = Vector{Solution}(undef, n_targs)
    costs_sol = CVector(undef, 0)
    cost_dd   = 0
    idx_dd    = 0

    # Initialize each `Solution` with `EmptySolution()`
    for j = 1:n_targs
        targ_sols[j] = EmptySolution()
    end

    return DDTOSolution(targ_sols,costs_sol,cost_dd,idx_dd)
end

"""
    EmptySolutionSCP()

    Constructor for an empty SolutionSCP.

# Arguments
- :out sol: empty solution.
"""
function EmptySolutionSCP()::SolutionSCP

    t = CVector(undef,0)
    r = CMatrix(undef,0,0)
    v = CMatrix(undef,0,0)
    T = CMatrix(undef,0,0)
    Γ = CVector(undef,0)
    ν = CMatrix(undef,0,0)
    μ = CMatrix(undef,0,0)
    r0_relax = CVector(undef,0)
    rf_relax = CVector(undef,0)
    cost = Inf
    T_nrm = CVector(undef,0)
    γ = CVector(undef,0)

    return SolutionSCP(t,r,v,T,Γ,ν,μ,r0_relax,rf_relax,cost,T_nrm,γ)
end

"""
    EmptyDDTOSolutionSCP()

Constructor for an empty DDTOSolutionSCP.

# Arguments
- :in n_targs: Number of targets for the solution
- :out sol: empty solution.
"""
function EmptyDDTOSolutionSCP(n_targs)::DDTOSolutionSCP

    targ_sols = Vector{SolutionSCP}(undef, n_targs)
    costs_sol = CVector(undef, 0)
    cost_dd   = 0
    idx_dd    = 0
    η         = CVector(undef,0)

    # Initialize each `SolutionSCP` with `EmptySolutionSCP()`
    for j = 1:n_targs
        targ_sols[j] = EmptySolutionSCP()
    end
    
    return DDTOSolutionSCP(targ_sols,costs_sol,cost_dd,idx_dd,η)
end

"""
    EmptyBranchSolution()

    Constructor for an empty BranchSolution.

# Arguments
- :out sol: empty solution.
"""
function EmptyBranchSolution()::BranchSolution
    sol = EmptySolution()
    cost_dd = 0
    idx_dd = 1
    return BranchSolution(sol, cost_dd, idx_dd)
end