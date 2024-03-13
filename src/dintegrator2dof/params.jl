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
    
    # >> Algorithm parameters <<
    a::AlgorithmParams
end

# ..:: Default DIntegrator2DoFParams Constructor ::..

function DIntegrator2DoFParams(;autogen_targs=false, autogen_targ_count=2)::DIntegrator2DoFParams

    # >> Constraint parameters <<
    u_max = 20

    # >> Algorithm parameters <<
    a = AlgorithmParams()
    a.nx = 5 # (pos, vel, accel integral)
    a.nu = 2 # (accel)

    # Boundary parameters
    e_x_R2 = [1;0]
    e_y_R2 = [0;1]
    a.n_targs = 4
    a.z0 = zeros(5)
    rf_targs = hcat(
        +35.55*e_x_R2 + 5.63*e_y_R2,
        +25.45*e_x_R2 + 25.45*e_y_R2,
        +0*e_x_R2     + 36*e_y_R2,
        +35.55*e_x_R2 - 5.63*e_y_R2
    )
    vf_targs = zeros(2,a.n_targs)
    a.zf_targs = vcat(rf_targs,vf_targs,Inf*ones(1,a.n_targs))
    a.λ_targs = [3,2,1,4]
    a.T_targs = [1,2,3,4]
    a.τ_targs = zeros(a.n_targs)
    a.α_targs = [1,1,1,1]
    # a.α_targs = [0,0,1,0]
    a.ϵ_targs = 0.7*ones(a.n_targs)
    a.u0 = Inf*ones(a.nu)
    a.uf_targs = Inf*ones(a.nu,a.n_targs)
    
    # >> SCP Params <<
    a.ctcs_enabled = false
    a.ddto_warmstart = false
    a.use_suboptimality = true
    a.w_obj_sing = 1
    a.w_obj_ddto = 10*a.w_obj_sing
    a.w_ctrl = 50
    a.w_buff = 0
    a.w_trust = 10
    a.ϵ_ctrl = 1e-2
    a.ϵ_buff = 1e-2
    a.ϵ_trust = 1e-2
    a.ϵ_ctcs = 1e-4
    a.scp_iters = 100

    # >> Discretization & time dilation <<
    # for DDTO-cvx: Δt = (Δt_min+Δt_max)/2
    a.N = 11
    a.Δt_min = 1e-6
    a.Δt_max = 1
    a.ToF_max = (a.N-1)*(a.Δt_min+a.Δt_max)/2
    a.ToF_min = a.ToF_max
    a.disc = 1

    # >> Affine scaling parameters <<
    rmin = [min([a.z0[k,:]; a.zf_targs[k,:]]...) for k∈1:2]
    rmax = [max([a.z0[k,:]; a.zf_targs[k,:]]...) for k∈1:2]
    a.Sx,a.sx = scaling_matrices([rmin; -ones(2); 0], [rmax; ones(2); u_max*a.ToF_max])
    a.Su,a.su = scaling_matrices(-u_max*ones(2), u_max*ones(2))

    # Overwrite targets if we are using automatic target generation
    if autogen_targs
        # Generate targets
        a.n_targs = autogen_targ_count
        rf_targs = generate_random_targets(a.n_targs, 20, [30;30])

        # Update target parameters
        a.zf_targs = vcat(rf_targs, zeros(2,a.n_targs), Inf*ones(1,a.n_targs))
        a.λ_targs = collect(1:a.n_targs)
        a.T_targs = collect(1:a.n_targs)
        a.τ_targs = zeros(a.n_targs)
        a.α_targs = ones(a.n_targs)
        a.ϵ_targs = a.ϵ_targs[1]*ones(a.n_targs)
        a.u0 = Inf*ones(a.nu)
        a.uf_targs = Inf*ones(a.nu,a.n_targs)
    end

    # >> Make params object <<
    params = DIntegrator2DoFParams(
        u_max,
        a
    )

    return params
end

function generate_random_targets(N::Int, radius, vertex)
    rf_targs = Matrix(undef, 2, N)
    for j = 1:N
        r_targ = radius * rand(Float64)
        θ_targ = 2 * pi * rand(Float64)
        rf_targs[:,j] = vertex + r_targ*cos(θ_targ)*[1;0] + r_targ*sin(θ_targ)*[0;1]
    end
    return rf_targs
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
    solution_proc = EmptyDIntegrator2DoFDDTOSolution(params.a.n_targs)
    for k = 1:params.a.n_targs
        # Obtain raw data from solution
        cost = solution.targs[k].cost
        x = solution.targs[k].x
        u = solution.targs[k].u

        # Determine if time dilation was used
        if params.a.nu == 2
            time_dilation = false # ddto-cvx
        else
            time_dilation = true # ddto-scp
        end

        # Handle time based on time dilation flag
        if time_dilation
            τ = solution.targs[k].t # recorded solution time was dilated!
            if ~isempty(u)
                t = time_dilation_control_to_wall_clock_time(u[end,:], τ, params.a.disc)
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
            s = zeros(params.a.N)
        end
        solution_proc.targs[k] = DIntegrator2DoFSolution(τ,t,r,v,a,s,cost)
    end

    return solution_proc
end