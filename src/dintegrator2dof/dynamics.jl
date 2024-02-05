function dynamics_linear(params::DIntegrator2DoFParams)
    A = CMatrix([
        zeros(2,2) I(2);
        zeros(2,2) zeros(2,2)
    ])
    B = CMatrix([
        zeros(2,2);
        I(2)
    ])
    return A,B
end

function dynamics_nonlinear(
    t::CReal,
    x::CVector,
    ν::CVector,
    params::DIntegrator2DoFParams)::CVector

    # Compute 2-DOF non-dilated dynamics
    A,B = dynamics_linear(params)
    u = ν[1:end-1]
    s = ν[end]
    f = A*x + B*u
    z = s*f # dilate dynamics w/ chain rule
    
    return z
end

function dynamics_linearized(
    t_ref::CReal,
    x_ref::CVector,
    ν_ref::CVector,
    params::DIntegrator2DoFParams)::Tuple{CMatrix,CMatrix,CVector,CVector}

    # Parse reference control
    u_ref = ν_ref[1:end-1]
    s_ref = ν_ref[end]

    # Obtain linear matrices as our differentials
    A,B = dynamics_linear(params)
    ∂f_∂x = A
    ∂f_∂u = B

    # ∂f_∂s: Evaluate nondilated nonlinear dynamics
    ∂f_∂s = dynamics_nonlinear(t_ref,x_ref,vcat(u_ref,1),params)

    # Package partials as linearized matrices
    A = s_ref*∂f_∂x
    B = Matrix([s_ref*∂f_∂u ∂f_∂s])
    Σ = []
    z = -(s_ref*∂f_∂x*x_ref + s_ref*∂f_∂u*u_ref)

    return(A,B,Σ,z)
end