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
        J_running = μ.^2 / fill(params.u_max^2, N_ctrl)
    end

    return J_running, J_term
end

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