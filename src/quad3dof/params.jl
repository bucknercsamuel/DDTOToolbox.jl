#= Shared 3-DOF quadcopter -- Parameter Structures and Functions.

Author: Samuel Buckner (UW-ACL)
=#

Quad3DoFParams = Union{Quad3DoFCageParams,Quad3DoFHaloParams}

# ..:: Post-processed Solution Structure ::..
mutable struct Quad3DoFSolution
    τ::CVector       # [s] Dilated Time Vector
    t::CVector       # [s] Wall-clock Time vector
    x::CMatrix       # [-] State trajectory
    u::CMatrix       # [-] Input trajectory
    r::CMatrix       # [m] Position trajectory
    v::CMatrix       # [m/s] Velocity trajectory
    T::CMatrix       # [m/s^2] Thrust vector
    s::CVector       # Time dilation factor (if using free-final-time)
    T_nrm::CVector   # [N] Thrust magnitude
    ∫T::CVector      # [N] Cumulative thrust
    γ::CVector       # [rad] Pointing angle
    cost::CReal      # Optimization's optimal cost
end

mutable struct Quad3DoFDDTOSolution
    targs::Vector{Quad3DoFSolution}  # Contains the `Quad3DoFSolution` for each target
end

# ..:: Constructors for empty `*Solution` structs ::..

function EmptyQuad3DoFSolution()::Quad3DoFSolution

    τ = CVector(undef,0)
    t = CVector(undef,0)
    x = CMatrix(undef,0,0)
    u = CMatrix(undef,0,0)
    r = CMatrix(undef,0,0)
    v = CMatrix(undef,0,0)
    T = CMatrix(undef,0,0)
    s = CVector(undef,0)
    T_nrm = CVector(undef,0)
    ∫T = CVector(undef,0)
    γ = CVector(undef,0)
    cost = Inf

    return Quad3DoFSolution(τ,t,x,u,r,v,T,s,T_nrm,∫T,γ,cost)
end

function EmptyQuad3DoFDDTOSolution(n_targs)::Quad3DoFDDTOSolution
    targs = Vector{Quad3DoFSolution}(undef, n_targs)
    for j = 1:n_targs
        targs[j] = EmptyQuad3DoFSolution()
    end
    return Quad3DoFDDTOSolution(targs)
end

# ..:: Function to convert raw `Solution` data for each branch to a `Quad3DoFSolution` ::..

function process_solutions(solution::DDTOSolution, params::Quad3DoFParams)::Quad3DoFDDTOSolution
    solution_proc = EmptyQuad3DoFDDTOSolution(params.a.n_targs)
    for k = 1:params.a.n_targs
        # Obtain raw data from solution
        cost = solution.targs[k].cost
        τ = solution.targs[k].t
        x = solution.targs[k].x
        u = solution.targs[k].u
        if ~isempty(u)
            t = time_dilation_control_to_wall_clock_time(u[end,:], τ, params.a.disc)
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
        solution_proc.targs[k] = Quad3DoFSolution(τ,t,x,u,r,v,T,s,T_nrm,∫T,γ,cost)
    end

    return solution_proc
end

# ..:: Extra custom functions ::..

function custom_scaling!(params::Quad3DoFParams)
    if typeof(params) == Quad3DoFCageParams && params.cage_bounds_enabled
        rmin = [params.x_arena_lims[1]; params.y_arena_lims[1]; params.z_arena_lims[1]]
        rmax = [params.x_arena_lims[2]; params.y_arena_lims[2]; params.z_arena_lims[2]]
    else
        rmin = [min([params.a.z0[k,:]; params.a.zf_targs[k,:]]...) for k∈1:3]
        rmax = [max([params.a.z0[k,:]; params.a.zf_targs[k,:]]...) for k∈1:3]
    end
    params.a.Sx,params.a.sx = scaling_matrices([rmin; -params.v_max_L*ones(3); 0], [rmax; params.v_max_L*ones(3); params.a.ToF_max*params.ρ_max])
    params.a.Su,params.a.su = scaling_matrices(-params.ρ_max*ones(3), params.ρ_max*ones(3))
end