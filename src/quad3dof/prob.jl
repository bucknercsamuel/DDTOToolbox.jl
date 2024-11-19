function prob_cost(
        mdl::JuMP.Model, 
        x::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        u::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        params::Quad3DoFParams;
        nonconvex::Bool = true)
        
    J_running = 0
    J_term = 0

    # If we are using a nonconvex model (SCP), use the thrust norm integral state
    if nonconvex
        ∫T = params.a.ctcs_enabled ? x[end-1,:] : x[end,:]
        J_term = ∫T[end] / params.ρ_max

    # If we are using convex model, just use the sum of thrust norm directly
    else
        N_ctrl = size(u,2)
        μ = @variable(mdl, [1:N_ctrl])
        @constraint(mdl, [k=1:N_ctrl], vcat(μ[k], u[:,k]) in SecondOrderCone())
        J_running = params.a.Δt_cvx * μ ./ params.ρ_max
    end

    return J_running, J_term
end

function prob_constraints(
        mdl::JuMP.Model, 
        x::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        u::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        params::Quad3DoFParams, 
        ref_traj::Solution,
        targ_idx::Int;
        obstacles::Bool = true,
        nonconvex::Bool = true)

    glideslope = true
    if targ_idx == 0
        glideslope = false
    end

    # Scenario-specific constraint management
    hold_altitude = false
    if typeof(params) == Quad3DoFCageParams
        hold_altitude = true
    end

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
    if hold_altitude
        @constraint(mdl, [k=1:N-1], r[3,k+1] == r[3,k])
    end

    # Thrust bounds
    @constraint(mdl, [k=1:N_ctrl], vcat(params.ρ_max, T[:,k]) in SecondOrderCone())
    if nonconvex
        Χ(k) = normalize(T_ref[:,k])
        @constraint(mdl, [k=1:N_ctrl], params.ρ_min - dot(Χ(k),T[:,k]) <= ν_thrust[k])
    end

    # Attitude pointing constraint
    @constraint(mdl, [k=1:N_ctrl], vcat(dot(T[:,k],e_z)/cos(params.γ_p), T[:,k]) in SecondOrderCone())

    # Approach cone / glideslope constraint
    if glideslope
        rf = params.a.zf_targs[1:3,targ_idx]
        @constraint(mdl, [k=1:N], vcat(dot(r[:,k] - rf,e_z)/cos(params.γ_gs), r[:,k] - rf) in SecondOrderCone())
    end

    # # Velocity bounds
    @constraint(mdl, [k=1:N], vcat(params.v_max_L,v[1:2,k]) in SecondOrderCone())
    @constraint(mdl, [k=1:N], v[3,k] <= params.v_max_V)
    @constraint(mdl, [k=1:N], v[3,k] >= params.v_min_V)

    # Cage bounds
    if hasproperty(params, :cage_bounds_enabled)
        if params.cage_bounds_enabled
            @constraint(mdl, [k=1:N], r[1,k] >= params.x_arena_lims[1])
            @constraint(mdl, [k=1:N], r[1,k] <= params.x_arena_lims[2])
            @constraint(mdl, [k=1:N], r[2,k] >= params.y_arena_lims[1])
            @constraint(mdl, [k=1:N], r[2,k] <= params.y_arena_lims[2])
            @constraint(mdl, [k=1:N], r[3,k] >= params.z_arena_lims[1])
            @constraint(mdl, [k=1:N], r[3,k] <= params.z_arena_lims[2])
        end
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
        params::Quad3DoFParams,
        targ_idx::Int; 
        sympy=false, 
        obstacles=true,
        rf_gs=nothing)

    glideslope = true
    if targ_idx == 0
        glideslope = false
    end

    # Scenario-specific constraint management
    hold_altitude = false
    if typeof(params) == Quad3DoFCageParams
        hold_altitude = true
    end

    # ..:: Setup ::..
    # Extract state and control elements
    r = x[1:3,:]
    v = x[4:6,:]
    T = u[1:3,:]
    N = size(x,2)
    N_ctrl = size(u,2)

    # Evaluated quantities for equality & inequality constraints
    h = 0
    g = 0

    # >> Equality constraints <<
    if hold_altitude
        h += augment_equality(r[3] - params.h_constant)
    end

    # >> Inequality constraints <<
    g += augment_inequality(norm(T)/params.ρ_max - 1) # Thrust upper bound
    g += augment_inequality(params.ρ_min/params.ρ_max - norm(T)/params.ρ_max) # Thrust lower bound
    g += augment_inequality(norm(T)/params.ρ_max - dot(T,e_z)/(params.ρ_max*cos(params.γ_p))) # Attitude pointing
    if glideslope
        if isnothing(rf_gs)
            rf_gs = params.a.zf_targs[1:3,targ_idx]
        end
        g += augment_inequality(norm(r-rf_gs) - dot(r-rf_gs,e_z)/cos(params.γ_gs)) # Glideslope
    end
    g += augment_inequality(norm(v[1:2]) - params.v_max_L) # Lateral velocity
    g += augment_inequality( v[3] - params.v_max_V) # Max vertical velocity
    g += augment_inequality( params.v_min_V - v[3]) # Min vertical velocity
    if hasproperty(params, :cage_bounds_enabled)
        if params.cage_bounds_enabled
            g += augment_inequality(+(r[1] - params.x_arena_lims[2]))
            g += augment_inequality(-(r[1] - params.x_arena_lims[1]))
            g += augment_inequality(+(r[2] - params.y_arena_lims[2]))
            g += augment_inequality(-(r[2] - params.y_arena_lims[1]))
            g += augment_inequality(+(r[3] - params.z_arena_lims[2]))
            g += augment_inequality(-(r[3] - params.z_arena_lims[1]))
        end
    end

    # Obstacle constraints
    if obstacles
        for o = 1:params.n_obstacles
            H = params.H_obstacles[o]
            p = params.p_obstacles[:,o]
            R = params.R_obstacles[o]
            g += augment_inequality(R - norm(H*(r-p)))
        end
    end

    # Set integrand ξ for CTCS
    ξ = g + h
    return ξ,g,h
end

function relu_huber_slope1(x)
    if x <= 0
        return 0
    elseif x <= 1/2
        return x^2
    else
        return x - 1/4
    end
end

function huber_slope1(x)
    if abs(x) <= 1/2
        return x^2
    elseif x < -1/2
        return -x - 1/4
    else
        return x - 1/4
    end
end

function augment_inequality(g::Union{Float64,ForwardDiff.Dual})
    return relu_huber_slope1(g)
end

function augment_equality(h)
    return huber_slope1(h)
end

# function augment_inequality(g::Sym)
#     zero = symbols("zero", real=true)
#     return max(zero,g).subs(zero,0)^augment_power
# end