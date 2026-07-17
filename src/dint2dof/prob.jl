#=
2-DOF double-integrator cost and constraint transcription for convex and SCP
optimal-control problems.
=#

"""
    prob_cost(mdl, x, u, params::DIntegrator2DoFParams; nonconvex=true) -> (J_running, J_term)

Build acceleration fuel cost: integral-state terminal cost for SCP, or
squared-acceleration epigraph terms for the convex model.

# Arguments
- `mdl`: JuMP model receiving epigraph variables in the convex case.
- `x`: state trajectory variables or affine expressions.
- `u`: control trajectory variables or affine expressions.
- `params`: 2-DOF parameters (`u_max`, CTCS layout).
- `nonconvex`: if `true`, terminal integral-state cost; else convex SOC terms (default `true`).

# Returns
- `J_running`: per-knot running cost (convex model) or zero.
- `J_term`: terminal fuel cost (SCP model) or zero.
"""
function prob_cost(
        mdl::JuMP.Model, 
        x::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        u::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        params::DIntegrator2DoFParams;
        nonconvex::Bool = true
    )
    J_running = 0
    J_term = 0

    # If we are using a nonconvex model (SCP), use the acceleration norm integral state
    if nonconvex
        ∫a = params.a.ctcs_enabled ? x[end-1,:] : x[end,:]
        J_term = ∫a[end] / params.u_max

    # If we are using convex model, just use the sum of acceleration norm directly
    else
        N_ctrl = size(u,2)
        μ = @variable(mdl, [1:N_ctrl])
        @constraint(mdl, [k=1:N_ctrl], vcat(μ[k], u[:,k]) in SecondOrderCone())
        J_running = μ.^2 ./ params.u_max^2
    end

    return J_running, J_term
end

"""
    prob_constraints(mdl, x, u, params::DIntegrator2DoFParams, ref_traj; obstacles=true, nonconvex=true)

Impose the maximum-acceleration second-order-cone constraint for the 2-DOF
double integrator. Returns an empty virtual-buffer list in the nonconvex case.

# Arguments
- `mdl`: JuMP model receiving acceleration constraints.
- `x`: state trajectory variables or affine expressions (unused beyond sizing).
- `u`: control trajectory variables or affine expressions.
- `params`: 2-DOF parameters (`u_max`).
- `ref_traj`: reference trajectory (unused; API compatibility).
- `obstacles`: reserved for obstacle constraints (unused in this benchmark).
- `nonconvex`: if `true`, return an empty virtual-buffer vector (default `true`).

# Returns
- `ν_buff`: empty vector when `nonconvex=true`; otherwise implicit `nothing`.
"""
function prob_constraints(
        mdl::JuMP.Model, 
        x::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        u::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        params::DIntegrator2DoFParams, 
        ref_traj::Solution;
        obstacles::Bool = true,
        nonconvex::Bool = true
    )

    # ..:: Setup ::..
    # Extract necessary state and control elements
    a = u[1:2,:]
    N_ctrl = size(u,2)

    # ..:: Constraints ::..
    @constraint(mdl, [k=1:N_ctrl], vcat(params.u_max, a[:,k]) in SecondOrderCone()) # maximum acceleration norm

    if nonconvex
        ν_buff = []
        return ν_buff
    end
end