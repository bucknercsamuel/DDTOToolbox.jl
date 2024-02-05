#= DDTO for double integrator -- Parameter Structures and Functions.

Author: Samuel Buckner (UW-ACL)
=#

# ..:: Double Integrator Object ::..

"""
`DIntegrator2DoFParams` holds the quadcopter parameters.
"""
mutable struct DIntegrator2DoFParams

    # >> Constraint parameters <<
    u_max::CReal   # [N] Maximum acceleration
    z0::CVector    # [m] Initial state
    nx::Int        # [-] Number of states
    nu::Int        # [-] Number of controls

    # >> DDTO target conditions <<
    n_targs::Int          # Current number of targets
    zf_targs::CMatrix     # [m] Terminal state of each target
    λ_targs::Vector{Int}  # Order of target rejection
    T_targs::Vector{Int}  # Tag for each target (deprecated)
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

# ..:: Default DIntegrator2DoFParams Constructor ::..

function DIntegrator2DoFParams()::DIntegrator2DoFParams

    # >> Constraint parameters <<
    u_max = 20
    nx = 4
    nu = 3

    # Boundary parameters
    e_x_R2 = [1;0]
    e_y_R2 = [0;1]
    n_targs = 4
    z0 = zeros(4)
    rf_targs = hcat(
        +35.55*e_x_R2 + 5.63*e_y_R2,
        +25.45*e_x_R2 + 25.45*e_y_R2,
        +0*e_x_R2     + 36*e_y_R2,
        +35.55*e_x_R2 - 5.63*e_y_R2
    )
    vf_targs = zeros(2,n_targs)
    zf_targs = vcat(rf_targs,vf_targs)
    λ_targs = [3,2,1,4]
    T_targs = [1,2,3,4]
    τ_targs = zeros(n_targs)
    α_targs = [0,0,1,0]
    ϵ_targs = 0.7*ones(n_targs)
    
    # >> SCP Params <<
    w_obj = 1e3
    w_ctrl = 1e4
    w_buff = 0
    w_trust = 1e2
    ϵ_ctrl = 1e-2
    ϵ_buff = 1e-2
    ϵ_trust = 1e-2
    scp_iters = 100

    # >> Discretization & time dilation <<
    # for DDTO-cvx: Δt = (Δt_min+Δt_max)/2
    N = 11
    Δt_min = 1e-6
    Δt_max = 1
    ToF_max = (N-1)*(Δt_min+Δt_max)/2
    disc = 1

    # >> Affine scaling parameters <<
    rmin = z0[1:2]
    rmax = [max(zf_targs[k,:]...) for k∈1:2]
    Δτ = 1/(N-1)
    Sx,sx = scaling_matrices([rmin; -ones(2)],    [rmax; ones(2)])
    Su,su = scaling_matrices([-u_max*ones(2); 0], [u_max*ones(2); Δt_max/Δτ])

    # >> Make params object <<
    params = DIntegrator2DoFParams(
        u_max,
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

# ..:: Post-processed Solution Structure ::..
mutable struct DIntegrator2DoFSolution
    τ::CVector       # [s] Dilated Time Vector
    t::CVector       # [s] Wall-clock Time vector
    r::CMatrix       # [m] Position trajectory
    v::CMatrix       # [m/s] Velocity trajectory
    a::CMatrix       # [m/s^2] Acceleration trajectory
    s::CVector       # Time dilation factor (if using free-final-time)
    cost::CReal      # Optimization's optimal cost
end

mutable struct DIntegrator2DoFDDTOSolution
    targs::Vector{DIntegrator2DoFSolution}  # Contains the `DIntegrator2DoFSolution` for each target
end

# ..:: Constructors for empty `*Solution` structs ::..

function EmptyDIntegrator2DoFSolution()::DIntegrator2DoFSolution

    τ = CVector(undef,0)
    t = CVector(undef,0)
    r = CMatrix(undef,0,0)
    v = CMatrix(undef,0,0)
    a = CMatrix(undef,0,0)
    s = CVector(undef,0)
    cost = Inf

    return DIntegrator2DoFSolution(τ,t,r,v,a,s,cost)
end

function EmptyDIntegrator2DoFDDTOSolution(n_targs)::DIntegrator2DoFDDTOSolution
    targs = Vector{DIntegrator2DoFSolution}(undef, n_targs)
    for j = 1:n_targs
        targs[j] = EmptyDIntegrator2DoFSolution()
    end
    return DIntegrator2DoFDDTOSolution(targs)
end

# ..:: Function to convert raw `Solution` data for each branch to a `DIntegrator2DoFSolution` ::..

function process_solutions(solution::DDTOSolution, params::DIntegrator2DoFParams)::DIntegrator2DoFDDTOSolution
    solution_proc = EmptyDIntegrator2DoFDDTOSolution(params.n_targs)
    for k = 1:params.n_targs
        # Obtain raw data from solution
        cost = solution.targs[k].cost
        x = solution.targs[k].x
        u = solution.targs[k].u

        # Determine if time dilation was used
        if size(u)[1] < params.nu
            time_dilation = false # ddto-cvx
        else
            time_dilation = true # ddto-scp
        end

        # Handle time based on time dilation flag
        if time_dilation
            τ = solution.targs[k].t # recorded solution time was dilated!
            if ~isempty(u)
                t = time_dilation_control_to_wall_clock_time(u[end,:], τ, params.disc)
            else
                t = 0
            end
        else
            t = solution.targs[k].t
            τ = range(0,1,length(solution.targs[k].t)) |> collect
        end

        # Post-processing
        r = x[1:2,:]
        v = x[3:4,:]
        a = u[1:2,:]
        if time_dilation
            s = u[3,:]
        else
            s = zeros(params.N)
        end
        solution_proc.targs[k] = DIntegrator2DoFSolution(τ,t,r,v,a,s,cost)
    end

    return solution_proc
end