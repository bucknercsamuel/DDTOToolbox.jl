function core_problem(
        mdl::JuMP.Model, 
        x::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        u::Union{Matrix{JuMP.VariableRef},Matrix{AffExpr}}, 
        params::DIntegrator2DoFParams, 
        ref_traj::Solution; # not used for this problem, purely convex path constraints/obj
        τ::Int = -1 # provides the option to get the cost only up to τ (defaults to not giving the cost if -1)
    )
    """
    NOTE: This function contains the scenario-specific path constraints and objective
    """

    # ..:: Setup ::..
    # Extract state and control elements
    r = x[1:2,:]
    v = x[3:4,:]
    a = u[1:2,:]
    s = u[end,:]
    N = size(x,2)
    N_ctrl = size(u,2)
    nx = size(x,1)

    # ..:: Constraints ::..
    @constraint(mdl, [k=1:N_ctrl], vcat(params.u_max, a[:,k]) in SecondOrderCone()) # maximum acceleration norm

    # ..:: Objective ::..
    if nx == 4 # using ddto-cvx
        μ = @variable(mdl, [1:N_ctrl])
        @constraint(mdl, [k=1:N_ctrl], vcat(μ[k], a[:,k]) in SecondOrderCone())
        J_running = μ.^2 ./ fill(params.u_max^2, N_ctrl)
        J_term = 0
        if τ > 0
            J_τ = μ[1:τ].^2 ./ fill(params.u_max^2, τ)
        end
    elseif nx == 5 # using ddto-scp with extra state for integral of acceleration
        J_running = 0
        J_term = x[end,end] / params.u_max
        if τ > 0
            J_τ = x[end,τ] / params.u_max
        end
    else
        error("Incorrect number of states")
    end

    
    # Concatenate all virtual buffers
    ν_buff = [] # no nonconvex path constraints to be buffered
    if τ > 0
        return J_running, J_term, ν_buff, J_τ
    elseif τ == 0
        return J_running, J_term, ν_buff, 0
    else
        return J_running, J_term, ν_buff
    end
end