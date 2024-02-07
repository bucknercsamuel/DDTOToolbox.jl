"""
NOTE: This function contains the core problem constraints shared by `solve_subproblem_decoupled` and `solve_subproblem_ddto`
"""
function core_problem(
        mdl::JuMP.Model, 
        x::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        u::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        params::Quad3DoFCageParams, 
        ref_traj::Solution
    )

    # ..:: Setup ::..
    # Extract state and control elements
    r = x[1:3,:]
    v = x[4:6,:]
    ∫T = x[7,:]
    T = u[1:3,:]
    s = u[4,:]
    N = size(x,2)
    N_ctrl = size(u,2)

    # Virtual buffers
    ν_thrust = @variable(mdl, [1:N_ctrl])
    ν_obs = @variable(mdl, [1:params.n_obstacles,1:N])

    # Extract reference trajectory elements
    x_ref = ref_traj.x
    u_ref = ref_traj.u
    r_ref = x_ref[1:3,:]
    v_ref = x_ref[4:6,:]
    ∫T_ref = x_ref[7,:]
    T_ref = u_ref[1:3,:]
    s_ref = u_ref[4,:]
    
    # ..:: Constraints ::..
    # Constant altitude constraint
    @constraint(mdl, [k=1:N-1], r[3,k+1] == r[3,k])

    # Thrust bounds
    Χ(k) = normalize(T_ref[:,k])
    @constraint(mdl, [k=1:N_ctrl], vcat(params.ρ_max, T[:,k]) in SecondOrderCone())
    @constraint(mdl, [k=1:N_ctrl], params.ρ_min - dot(Χ(k),T[:,k]) <= ν_thrust[k])
    # @constraint(mdl, [k=1:N_ctrl], params.ρ_min - dot(Χ(k),T[:,k]) <= 0)

    # Attitude pointing constraint
    @constraint(mdl, [k=1:N_ctrl], vcat(dot(T[:,k],e_z)/cos(params.γ_p), T[:,k]) in SecondOrderCone())

    # Velocity upper bound
    @constraint(mdl, [k=1:N], vcat(params.v_max_L,v[1:2,k]) in SecondOrderCone())

    # Cage bounds
    # @constraint(mdl, [k=1:N], r[1,k] >= params.x_arena_lims[1])
    # @constraint(mdl, [k=1:N], r[1,k] <= params.x_arena_lims[2])
    # @constraint(mdl, [k=1:N], r[2,k] >= params.y_arena_lims[1])
    # @constraint(mdl, [k=1:N], r[2,k] <= params.y_arena_lims[2])
    # @constraint(mdl, [k=1:N], r[3,k] >= params.z_arena_lims[1])
    # @constraint(mdl, [k=1:N], r[3,k] <= params.z_arena_lims[2])

    # Obstacle constraints
    for o = 1:params.n_obstacles
        H = params.H_obstacles[o]
        for k = 1:N
            Δr = r_ref[:,k] - params.p_obstacles[:,o]
            δr = r[:,k] - r_ref[:,k]
            ξ  = max(norm(H*Δr,2),1e-2)
            ζ  = transpose(H)*H*Δr / ξ
            @constraint(mdl, ξ + dot(ζ,δr) >= params.R_obstacles[o] + ν_obs[o,k])
            # @constraint(mdl, ξ + dot(ζ,δr) >= params.R_obstacles[o])
        end
    end

    # Process variables
    J_running, J_term = objective_function(mdl, x, u, params)
    ν_buff = [vec(ν_obs);vec(ν_thrust)]

    return J_running, J_term, ν_buff
end

function path_constraint_eval(x::Vector,u::Vector,params::Quad3DoFCageParams)

    # ..:: Setup ::..
    # Extract state and control elements
    r = x[1:3,:]
    v = x[4:6,:]
    T = u[1:3,:]
    N = size(x,2)
    N_ctrl = size(u,2)

    # Vectors for equality & inequality constraints
    h = []
    g = []

    # >> Equality constraints <<
    append!(h, r[3] - params.z0[3]) # Zero altitude in plane (TODO: make sure this doesn't create bugs!!)

    # >> Inequality constraints <<
    append!(g, norm(T) - params.ρ_max) # Thrust upper bound
    append!(g, params.ρ_min - norm(T)) # Thrust lower bound
    append!(g, norm(T) - dot(T,e_z)/cos(params.γ_p)) # Attitude pointing
    append!(g, norm(v[1:2]) - params.v_max_L) # Lateral velocity
    # append!(g, +(r[1] - params.x_arena_lims[2])) # Cage bound
    # append!(g, -(r[1] - params.x_arena_lims[1])) # Cage bound
    # append!(g, +(r[2] - params.y_arena_lims[2])) # Cage bound
    # append!(g, -(r[2] - params.y_arena_lims[1])) # Cage bound
    # append!(g, +(r[3] - params.z_arena_lims[2])) # Cage bound
    # append!(g, -(r[3] - params.z_arena_lims[1])) # Cage bound

    # Obstacle constraints
    for o = 1:params.n_obstacles
        H = params.H_obstacles[o]
        p = params.p_obstacles[:,o]
        R = params.R_obstacles[o]
        append!(g, R - norm(H*(r-p)))
    end

    # >> Determine the integrand ξ for CTCS
    n_g = length(g)
    ξ = sum(([max(0,g[j])^2 for j∈1:n_g])) + sum(h.^2)

    return ξ,g,h
end

function objective_function(
        mdl::JuMP.Model, 
        x::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        u::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        params::Quad3DoFCageParams
    )
    if params.ctcs_enabled
        ∫T = x[end-1,:]
    else
        ∫T = x[end,:]
    end

    J_running = 0
    J_term = ∫T[end]
    return J_running, J_term
end
