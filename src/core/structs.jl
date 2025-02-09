# ..:: Solution Structures ::..
"""
`Solution` stores the optimal solution.
`DDTOSolution` stores the optimal DDTO bundled solution.
`AlgorithmParams` stores algorithm parameters.
"""

mutable struct Solution
    t::CVector        # [s] Time vector
    x::CMatrix        # State trajectory
    u::CMatrix        # Control signal
    cost::CReal       # Optimal cost
end

mutable struct DDTOSolution
    targs::Vector{Solution} # Contains the `Solution` to each target
end

mutable struct AlgorithmParams
    # >> Base traj opt parameters <<
    z0::Vector{CReal}              # Initial state (inf = no constraint)
    u0::Vector{CReal}              # Initial input (inf = no constraint)
    nx::Int                        # Number of states
    nu::Int                        # Number of controls

    # >> DDTO target conditions <<
    n_targs::Int                   # Current number of targets
    zf_targs::Matrix{CReal}        # Terminal state of each target (inf = no constraint)
    uf_targs::Matrix{CReal}        # Terminal input of each target (inf = no constraint)
    λ_targs::Vector{Int}           # Order of target rejection
    J_targs::Vector{Int}           # Target indexing set
    ID_targs::Vector{Int}          # ID associated with each target (for tracking purposes)
    τ_targs::Vector{Int}           # Deferrability index allocation (in order specified by λ_targs) -- set automatically in `solve_tree_ddto`
    α_targs::Vector{CReal}         # Relative weight for deferrability of each target
    ϵ_targs::Vector{CReal}         # Optimality tolerances

    # >> SCP Params <<
    ctcs_enabled::Bool             # Determines if Continuous-Time Constraint Satisfaction (CTCS) should be used
    warmstart_method::String       # Determines what warmstart method we should use (types: linear, single, ddto)
    use_suboptimality::Bool        # Determines if we should compute reference solutions and apply a suboptimality constraint
    use_single_cvx::Bool           # Determines if we should use single-target cvx-optimized solution instead of a SCP-optimized solution for single-target SCP
    w_obj_sing::CReal              # Objective penalty weight (Single-Target)
    w_obj_ddto::CReal              # Objective penalty weight (DDTO)
    w_ctrl::CReal                  # Virtual control penalty weight
    w_buff::CReal                  # Virtual buffer penalty weight
    w_trust::CReal                 # Trust region penalty weight
    ϵ_ctrl::CReal                  # Convergence threshold for virtual control penalty
    ϵ_buff::CReal                  # Convergence threshold for virtual buffer penalty
    ϵ_trust::CReal                 # Convergence threshold for trust region penalty
    ϵ_ctcs::CReal                  # Relaxation tolerance for CTCS violation constraint
    scp_iters::Int                 # Number of SCP subproblem iterations

    # >> Time dilation & discretization <<
    N::Int                         # Number of nodes (for all targets)
    Δt_min::CReal                  # [s] Minimum wall time step
    Δt_max::CReal                  # [s] Maximum wall time step\
    ToF_min::CReal                 # [s] Minimum physical time-of-flight for all targets
    ToF_max::CReal                 # [s] Maximum physical time-of-flight for all targets    
    disc::Int                      # Discretization hold order (currently can either choose 0 or 1)
    N_msi::Int                     # Number of multiple shooting integration steps (per node interval)
    N_sim::Int                     # Number of post-processing simulation steps (per node interval)
    differentiator::String         # Type of differentiation scheme (types: sympy, forwarddiff)

    # DDTO-CVX specific
    gss_cvx::Bool                  # Determine if golden section search should be used to find optimal `N_cvx``
    Δt_cvx::CReal                  # [s] Time step

    # >> Affine scaling parameters <<
    Sx::Matrix{CReal}              # Scaling transformation matrix for state "x"
    sx::Vector{CReal}              # Scaling affine vector for state "x"
    Su::Matrix{CReal}              # Scaling transformation matrix for state "u"
    su::Vector{CReal}              # Scaling affine vector for state "u"
end

# ..:: Constructors for structs ::..

function EmptySolution()::Solution

    t = CVector(undef,0)
    x = CMatrix(undef,0,0)
    u = CMatrix(undef,0,0)
    cost = Inf

    return Solution(t,x,u,cost)
end

function EmptyDDTOSolution(n_targs)::DDTOSolution

    targs = Vector{Solution}(undef, n_targs)
    for j = 1:n_targs
        targs[j] = EmptySolution()
    end

    return DDTOSolution(targs)
end

function AlgorithmParams()::AlgorithmParams
    # >> Base traj opt parameters <<
    z0 = CVector(undef,0)
    u0 = CVector(undef,0)
    nx = 0
    nu = 0

    # >> DDTO target conditions <<
    n_targs = 0
    zf_targs = CMatrix(undef,0,0)
    uf_targs = CMatrix(undef,0,0)
    λ_targs = Array{Int}(undef,0)
    J_targs = Array{Int}(undef,0)
    ID_targs = Array{Int}(undef,0)
    τ_targs = Array{Int}(undef,0)
    α_targs = CVector(undef,0)
    ϵ_targs = CVector(undef,0)

    # >> SCP Params <<
    ctcs_enabled = true
    ddto_warmstart = "linear"
    use_suboptimality = true
    use_single_cvx = false
    w_obj_sing = .01
    w_obj_ddto = .01
    w_ctrl = 50.
    w_buff = 50.
    w_trust = 1.
    ϵ_ctrl = 1e-3
    ϵ_buff = 1e-3
    ϵ_trust = 1e-3
    ϵ_ctcs = 1e-4
    scp_iters = 10

    # >> Time dilation & discretization <<
    N = 11
    Δt_min = 0.01
    Δt_max = 2.
    ToF_min = 0.
    ToF_max = 10.
    disc = 1
    gss_cvx = true
    Δt_cvx = (Δt_min + Δt_max)/2
    N_msi = 10.
    N_sim = 40.
    differentiator = "forwarddiff"

    # >> Affine scaling parameters <<
    Sx = zeros(nx,nx)
    Su = zeros(nu,nu)
    sx = zeros(nx)
    su = zeros(nu)

    return AlgorithmParams(
        z0,
        u0,
        nx,
        nu,
        n_targs,
        zf_targs,
        uf_targs,
        λ_targs,
        J_targs,
        ID_targs,
        τ_targs,
        α_targs,
        ϵ_targs,
        ctcs_enabled,
        ddto_warmstart,
        use_suboptimality,
        use_single_cvx,
        w_obj_sing,
        w_obj_ddto,
        w_ctrl,
        w_buff,
        w_trust,
        ϵ_ctrl,
        ϵ_buff,
        ϵ_trust,
        ϵ_ctcs,
        scp_iters,
        N,
        Δt_min,
        Δt_max,
        ToF_min,
        ToF_max,
        disc,
        N_msi,
        N_sim,
        differentiator,
        gss_cvx,
        Δt_cvx,
        Sx,
        sx,
        Su,
        su,
    )
end

# ..:: Template functions ::..
# These functions are problem-specific and defined for a specific `params` object in other folders besides `core`,
# with the `prob.jl` and `dynamics.jl` files.

function prob_constraints(mdl::JuMP.Model, x::JuMP.VariableRef, u::JuMP.VariableRef, params::Any, ref_traj::Solution)
    return 0
end

function prob_cost(mdl::JuMP.Model, x::JuMP.VariableRef, u::JuMP.VariableRef, params::Any)
    return 0
end

function param_update_law!(params::Any)
    return 0
end

function dynamics_linearized(t_ref::CReal, x_ref::CVector, ν_ref::CVector, params::Any)
    return 0
end

function dynamics_nonlinear(t::CReal, x::CVector, ν::CVector, params::Any)
    return 0
end

function dynamics_linear(params::Any)
    return 0
end