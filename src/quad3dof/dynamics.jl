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