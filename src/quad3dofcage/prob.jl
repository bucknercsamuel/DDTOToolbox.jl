function prob_cost(
        mdl::JuMP.Model, 
        x::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        u::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        params::Quad3DoFCageParams;
        nonconvex::Bool = true
    )
    J_running = 0
    J_term = 0

    # If we are using a nonconvex model (SCP), use the thrust norm integral state
    if nonconvex
        ∫T = params.a.ctcs_enabled ? x[end-1,:] : x[end,:]
        J_term = ∫T[end]

    # If we are using convex model, just use the sum of thrust norm directly
    else
        N_ctrl = size(u,2)
        μ = @variable(mdl, [1:N_ctrl])
        @constraint(mdl, [k=1:N_ctrl], vcat(μ[k], u[:,k]) in SecondOrderCone())
        J_running = μ / fill(params.ρ_max, N_ctrl)
    end

    return J_running, J_term
end

function prob_constraints(
        mdl::JuMP.Model, 
        x::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        u::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        params::Quad3DoFCageParams, 
        ref_traj::Solution;
        obstacles::Bool = true,
        nonconvex::Bool = true
    )

    # ..:: Setup ::..
    # Extract state and control elements
    r = x[1:3,:]
    v = x[4:6,:]
    ∫T = x[7,:]
    T = u[1:3,:]
    N = size(x,2)
    N_ctrl = size(u,2)

    if nonconvex
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
    end
    
    # ..:: Constraints ::..
    # Constant altitude constraint
    @constraint(mdl, [k=1:N-1], r[3,k+1] == r[3,k])

    # Thrust bounds
    @constraint(mdl, [k=1:N_ctrl], vcat(params.ρ_max, T[:,k]) in SecondOrderCone())
    if nonconvex
        Χ(k) = normalize(T_ref[:,k])
        @constraint(mdl, [k=1:N_ctrl], params.ρ_min - dot(Χ(k),T[:,k]) <= ν_thrust[k])
    end

    # Attitude pointing constraint
    @constraint(mdl, [k=1:N_ctrl], vcat(dot(T[:,k],e_z)/cos(params.γ_p), T[:,k]) in SecondOrderCone())

    # Velocity upper bound
    @constraint(mdl, [k=1:N], vcat(params.v_max_L,v[1:2,k]) in SecondOrderCone())

    # Cage bounds
    if params.cage_bounds_enabled
        @constraint(mdl, [k=1:N], r[1,k] >= params.x_arena_lims[1])
        @constraint(mdl, [k=1:N], r[1,k] <= params.x_arena_lims[2])
        @constraint(mdl, [k=1:N], r[2,k] >= params.y_arena_lims[1])
        @constraint(mdl, [k=1:N], r[2,k] <= params.y_arena_lims[2])
        @constraint(mdl, [k=1:N], r[3,k] >= params.z_arena_lims[1])
        @constraint(mdl, [k=1:N], r[3,k] <= params.z_arena_lims[2])
    end

    # Obstacle constraints
    if obstacles && nonconvex
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
    end

    if nonconvex
        ν_buff = [vec(ν_obs);vec(ν_thrust)]
        return ν_buff
    end
end

function prob_constraints_eval(
        x::Vector,
        u::Vector,
        params::Quad3DoFCageParams; 
        sympy=false, 
        obstacles=true
    )

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
    append!(h, r[3] - params.h_constant)

    # >> Inequality constraints <<
    append!(g, norm(T) - params.ρ_max) # Thrust upper bound
    append!(g, params.ρ_min - norm(T)) # Thrust lower bound
    append!(g, norm(T) - dot(T,e_z)/cos(params.γ_p)) # Attitude pointing
    append!(g, norm(v[1:2]) - params.v_max_L) # Lateral velocity
    if params.cage_bounds_enabled
        append!(g, +(r[1] - params.x_arena_lims[2]))
        append!(g, -(r[1] - params.x_arena_lims[1]))
        append!(g, +(r[2] - params.y_arena_lims[2]))
        append!(g, -(r[2] - params.y_arena_lims[1]))
        append!(g, +(r[3] - params.z_arena_lims[2]))
        append!(g, -(r[3] - params.z_arena_lims[1]))
    end

    # Obstacle constraints
    if obstacles
        for o = 1:params.n_obstacles
            H = params.H_obstacles[o]
            p = params.p_obstacles[:,o]
            R = params.R_obstacles[o]
            append!(g, R - norm(H*(r-p)))
        end
    end

    # Determine the integrand ξ for CTCS
    n_g = length(g)
    if sympy
        zero = symbols("zero", real=true)
        ξ = sum(([(max(zero,g[j]).subs(zero,0))^2 for j∈1:n_g])) + sum(h.^2)
    else
        ξ = sum(([max(0,g[j])^2 for j∈1:n_g])) + sum(h.^2)
    end

    return ξ,g,h
end