function dynamics_linear_noaugment(params::Quad3DoFParams)
    return double_integrator_dynamics(dim=3, mass=params.mass, gravity=params.g)
end

function dynamics_linear(params::Quad3DoFParams)
    return double_integrator_dynamics(dim=3, mass=params.mass, gravity=params.g, augment=true, augment_dim=1)
end

function dynamics_nonlinear_nondilated(
        t,
        x,
        u,
        params::Quad3DoFParams)

    # Compute 3-DOF dynamics
    A,B,p = dynamics_linear_noaugment(params)
    f_3dof = A*x[1:6] + B*u[1:3] + p

    # Add drag term (if enabled)
    if params.drag_term_enabled
        v = x[4:6]
        v_aug = vcat(zeros(3),v)
        f_3dof += params.C_d*params.S_A*params.ρ*norm(v)*v_aug/2
    end

    # Compute additional states (thrust integral)
    ∫T = norm(u)

    # Stack function together and apply time dilation (chain rule)
    f = [f_3dof;∫T]
    
    return f
end

function dynamics_nonlinear(
    t::CReal,
    x::CVector,
    ν::CVector,
    params::Quad3DoFParams)::CVector

    f = dynamics_nonlinear_nondilated(t,x,ν[1:end-1],params)
    s = ν[end]
    z = s*f

    return z
end

function dynamics_nonlinear_nondilated_ctcs(
        t::CReal,
        x::Vector,
        u::Vector,
        params::Quad3DoFParams,
        targ_idx::Int)::Vector

    # Dynamics and CTCS state
    f_3dof = dynamics_nonlinear_nondilated(t,x,u,params)
    ξ,_,_ = prob_constraints_eval(x,u,params,targ_idx) # CTCS violation

    # if TEMP_FLAG
    #     ξ *= 0
    # end

    # Stack function together and apply time dilation (chain rule)
    f = [f_3dof;ξ]

    return f
end

function dynamics_nonlinear_ctcs(
        t::CReal,
        x::Vector,
        ν::Vector,
        params::Quad3DoFParams,
        targ_idx::Int)::Vector

    u = ν[1:end-1]
    s = ν[end]
    f = dynamics_nonlinear_nondilated_ctcs(t,x,u,params,targ_idx)
    z = s*f
    return z
end

function dynamics_linearized_ctcs(
        t_ref::CReal,
        x_ref::CVector,
        ν_ref::CVector,
        params::Quad3DoFParams,
        targ_idx::Int)::Tuple{CMatrix,CMatrix,CVector,CVector}

    # Parse reference control
    u_ref = ν_ref[1:end-1]
    s_ref = ν_ref[end]

    # Obtain nondilated dynamics jacobians
    if params.a.differentiator == "sympy"
        ∂f_∂x,∂f_∂u = evaluate_jacobians_sympy(t_ref,x_ref,u_ref,params,targ_idx)
    elseif params.a.differentiator == "forwarddiff"
        nx = length(x_ref)
        nu = length(u_ref)
        fun(z) = dynamics_nonlinear_nondilated_ctcs(t_ref,z[1:nx],z[nx+1:end],params,targ_idx)
        ∂f_∂z = ForwardDiff.jacobian(fun,vcat(x_ref,u_ref))
        ∂f_∂x = ∂f_∂z[:,1:nx]
        ∂f_∂u = ∂f_∂z[:,nx+1:end]
    else
        error("Please choose a valid differentiator option")
    end

    # ∂f_∂s: Evaluate nondilated nonlinear dynamics
    ∂f_∂s = dynamics_nonlinear_ctcs(t_ref,x_ref,vcat(u_ref,1),params,targ_idx)

    # Package partials as linearized matrices
    A = s_ref*∂f_∂x
    B = Matrix([s_ref*∂f_∂u ∂f_∂s])
    Σ = []
    z = -(s_ref*∂f_∂x*x_ref + s_ref*∂f_∂u*u_ref)

    return(A,B,Σ,z)
end