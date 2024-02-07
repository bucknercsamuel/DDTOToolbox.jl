function dynamics_nonlinear_nondilated_ctcs(
    t::CReal,
    x::CVector,
    u::CVector,
    params::Quad3DoFCageParams)::CVector

    A,B,p = dynamics_linear(params)
    f_3dof = A*x[1:6] + B*u + p

    # Compute additional states
    Tnorm = norm(u) # thrust integral
    ξ,_,_ = path_constraint_eval(x,u,params) # CTCS violation

    # Stack function together and apply time dilation (chain rule)
    f = [f_3dof;Tnorm;ξ]

    return f
end

function dynamics_nonlinear_ctcs(
    t::CReal,
    x::CVector,
    ν::CVector,
    params::Quad3DoFCageParams)::CVector
    u = ν[1:end-1]
    s = ν[end]
    f = dynamics_nonlinear_nondilated_ctcs(t,x,u,params)
    z = s*f
    return z
end

function dynamics_linearized_ctcs(
    t_ref::CReal,
    x_ref::CVector,
    ν_ref::CVector,
    params::Quad3DoFCageParams)::Tuple{CMatrix,CMatrix,CVector,CVector}

    # # Parse reference control
    # u_ref = ν_ref[1:end-1]
    # s_ref = ν_ref[end]
    # nx = length(x_ref)

    # # Evaluate the numerical jacobian at this reference point
    # y = [x_ref; u_ref]
    # f(y) = dynamics_nonlinear_nondilated_ctcs(0.,y[1:nx],y[nx+1:end],params)
    # ∂f_∂y = jacobian(Forward, f, y);
    # ∂f_∂x = ∂f_∂y[:,1:nx]
    # ∂f_∂u = ∂f_∂y[:,nx+1:end]

    # # ∂f_∂s: Evaluate nondilated nonlinear dynamics
    # ∂f_∂s = dynamics_nonlinear_nondilated_ctcs(t_ref,x_ref,u_ref,params)

    # # Package partials as linearized matrices
    # A = s_ref*∂f_∂x
    # B = Matrix([s_ref*∂f_∂u ∂f_∂s])
    # Σ = []
    # z = -(s_ref*∂f_∂x*x_ref + s_ref*∂f_∂u*u_ref)

    dyn_nl = (t,x,ν) -> dynamics_nonlinear_ctcs(t,x,ν,params)
    A,B,z = numerical_jacobian(t_ref,x_ref,ν_ref,dyn_nl)
    Σ = []
    
    return(A,B,Σ,z)
end