#=
3-DOF quadcopter optimal-control cost and constraint transcription (JuMP),
plus CTCS constraint-augmentation helpers used by nonlinear dynamics.
=#

"""
    prob_cost(mdl, x, u, params::Quad3DoFParams; nonconvex=true) -> (J_running, J_term)

Build the 3-DOF thrust fuel cost: integral state terminal cost for SCP
(`nonconvex=true`) or SOC-epigraph thrust sum for convex models.

# Arguments
- `mdl`: JuMP model receiving any auxiliary epigraph variables.
- `x`: state variables or affine expressions along the horizon.
- `u`: control variables or affine expressions along the horizon.
- `params`: 3-DOF scenario parameters (thrust bounds, CTCS layout).
- `nonconvex`: if `true`, use the thrust-norm integral state terminal cost; if
  `false`, use convex SOC epigraph terms (default `true`).

# Returns
- `J_running`: per-knot running cost terms (nonzero only in the convex model).
- `J_term`: terminal fuel cost term (nonzero only in the SCP model).
"""
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

"""
    prob_constraints(mdl, x, u, params::Quad3DoFParams, ref_traj, targ_idx; obstacles=true, nonconvex=true) -> ν_buff

Impose 3-DOF path constraints (thrust, pointing, velocity, cage bounds,
linearized obstacles). Returns virtual-buffer variables for nonconvex SCP;
cage scenarios also enforce constant altitude.

# Arguments
- `mdl`: JuMP model receiving constraints and virtual buffers.
- `x`: state trajectory variables or affine expressions.
- `u`: control trajectory variables or affine expressions.
- `params`: 3-DOF scenario parameters (limits, obstacles, cage flags).
- `ref_traj`: reference trajectory for linearized obstacle constraints.
- `targ_idx`: target index (`0` denotes the shared trunk segment).
- `obstacles`: if `true`, include linearized obstacle buffers (default `true`).
- `nonconvex`: if `true`, return virtual-buffer slacks for SCP (default `true`).

# Returns
- `ν_buff`: concatenated virtual-buffer variables for SCP, or an empty vector in
  the convex transcription.
"""
function prob_constraints(
        mdl::JuMP.Model, 
        x::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        u::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        params::Quad3DoFParams, 
        ref_traj::Solution,
        targ_idx::Int;
        obstacles::Bool = true,
        nonconvex::Bool = true)

    # Never impose glideslope constraint for the trunk trajectory segment (idx = 0)
    # may lead to infeasible solutions
    glideslope = true
    if targ_idx == 0
        glideslope = false
    end

    # Scenario-specific constraint management
    # (cage scenario must fly in a plane -- constant altitude)
    hold_altitude = false
    if isa(params, Quad3DoFCageParams)
        hold_altitude = true
        glideslope = false
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
    if hasproperty(params, :γ_p)
        @constraint(mdl, [k=1:N_ctrl], vcat(dot(T[:,k],e_z), T[:,k]*cos(params.γ_p)) in SecondOrderCone())
    end

    # # Approach cone / glideslope constraint
    # if glideslope
    #     rf = params.a.zf_targs[1:3,targ_idx]
    #     @constraint(mdl, [k=1:N], vcat(dot(r[:,k] - rf,e_z), (r[:,k] - rf)*cos(params.γ_gs)) in SecondOrderCone())
    # end

    # # Velocity bounds
    @constraint(mdl, [k=1:N], vcat(params.v_max_L,v[1:2,k]) in SecondOrderCone())
    @constraint(mdl, [k=1:N], v[3,k] <= params.v_max_V)
    if hasproperty(params, :v_min_V)
        @constraint(mdl, [k=1:N], v[3,k] >= params.v_min_V)
    end

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
    else
        ν_buff = []
    end
    return ν_buff
end

"""
    prob_cost_eval(x, u, params::Quad3DoFParams; nonconvex=true) -> (J_running, J_term)

Evaluate the 3-DOF fuel cost on numeric trajectories (no JuMP model).

# Arguments
- `x`: numeric state trajectory matrix.
- `u`: numeric control trajectory matrix.
- `params`: 3-DOF scenario parameters.
- `nonconvex`: if `true`, use integral-state terminal cost; else sum thrust norms
  (default `true`).

# Returns
- `J_running`: scalar or vector running cost contribution.
- `J_term`: terminal fuel cost contribution.
"""
function prob_cost_eval(
        x::Matrix,
        u::Matrix,
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
        μ = [norm(u[:,k]) for k = 1:N_ctrl]
        J_running = params.a.Δt_cvx * sum(μ) / params.ρ_max
    end

    return J_running, J_term
end

"""
    prob_constraints_eval(x, u, params::Quad3DoFParams, targ_idx; sympy=false, obstacles=true) -> (ξ, g, h)

Evaluate CTCS constraint-violation integrands at a single state/control sample:
inequality penalty `g`, equality penalty `h`, and combined `ξ = g + h`.

# Arguments
- `x`: state sample vector (or column matrix treated as a single knot).
- `u`: control sample vector (or column matrix treated as a single knot).
- `params`: 3-DOF scenario parameters defining limits and obstacles.
- `targ_idx`: target index (`0` = trunk; selects glideslope reference when enabled).
- `sympy`: reserved for symbolic evaluation paths (default `false`).
- `obstacles`: if `true`, include obstacle violation penalties (default `true`).

# Returns
- `ξ`: combined CTCS integrand `g + h`.
- `g`: accumulated inequality violation penalty.
- `h`: accumulated equality violation penalty.
"""
function prob_constraints_eval(
        x::Vector,
        u::Vector,
        params::Quad3DoFParams,
        targ_idx::Int; 
        sympy=false, 
        obstacles=true)

    # Never impose glideslope constraint for the trunk trajectory segment (idx = 0)
    glideslope = true
    if targ_idx == 0
        glideslope = false
    end

    # Scenario-specific constraint management
    hold_altitude = false
    if isa(params, Quad3DoFCageParams)
        hold_altitude = true
        glideslope = false
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
    if hasproperty(params, :γ_p)
        g += augment_inequality(norm(T)*cos(params.γ_p)/params.ρ_max - dot(T,e_z)/(params.ρ_max)) # Attitude pointing
    end
    if glideslope
        rf_gs = params.a.zf_targs[1:3,targ_idx]
        g += augment_inequality(norm(r-rf_gs)*cos(params.γ_gs) - dot(r-rf_gs,e_z)) # Glideslope
    end
    g += augment_inequality(norm(v[1:2]) - params.v_max_L) # Lateral velocity
    g += augment_inequality( v[3] - params.v_max_V) # Max vertical velocity
    if hasproperty(params, :v_min_V)
        g += augment_inequality( params.v_min_V - v[3]) # Min vertical velocity
    end
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
            p = @view params.p_obstacles[:,o]
            R = params.R_obstacles[o]
            g += augment_inequality(R - norm(H*(r-p)))
        end
    end

    # Set integrand ξ for CTCS
    ξ = g + h
    return ξ,g,h
end

"""
    relu_huber_slope1(x) -> Number

Huber-smoothed ReLU used to softly penalize inequality violations in CTCS.

# Arguments
- `x`: raw inequality residual (scalar).

# Returns
- Smoothed nonnegative penalty value applied to `x`.
"""
function relu_huber_slope1(x)
    if x <= 0
        return 0
    elseif x <= 1/2
        return x^2
    else
        return x - 1/4
    end
end

"""
    huber_slope1(x) -> Number

Huber loss with unit slope outside the quadratic region; used for equality
penalties in CTCS.

# Arguments
- `x`: raw equality residual (scalar).

# Returns
- Smoothed penalty value applied to `x`.
"""
function huber_slope1(x)
    if abs(x) <= 1/2
        return x^2
    elseif x < -1/2
        return -x - 1/4
    else
        return x - 1/4
    end
end

"""
    augment_inequality(g) -> Number

Map a raw inequality residual to a CTCS penalty via `relu_huber_slope1`.

# Arguments
- `g`: raw inequality constraint residual.

# Returns
- CTCS inequality penalty value.
"""
function augment_inequality(g::Union{Float64,ForwardDiff.Dual})
    return relu_huber_slope1(g)
end

"""
    augment_equality(h) -> Number

Map a raw equality residual to a CTCS penalty via `huber_slope1`.

# Arguments
- `h`: raw equality constraint residual.

# Returns
- CTCS equality penalty value.
"""
function augment_equality(h)
    return huber_slope1(h)
end

# function augment_inequality(g::Sym)
#     zero = symbols("zero", real=true)
#     return max(zero,g).subs(zero,0)^augment_power
# end