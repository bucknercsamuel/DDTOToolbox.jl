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

    # Compute 3-DOF dynamics (with drag term)
    A,B,p = dynamics_linear_noaugment(params)
    f_3dof = A*x[1:6] + B*u[1:3] + p

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

function dynamics_linearized(
    t_ref::CReal,
    x_ref::CVector,
    ν_ref::CVector,
    params::Quad3DoFParams)::Tuple{CMatrix,CMatrix,CVector,CVector}

    # Parse reference control
    u_ref = ν_ref[1:end-1]
    s_ref = ν_ref[end]

    # Necessary reference quantities and constants for partials
    T1,T2,T3 = u_ref
    m = params.mass

    # Matrices to populate
    nx = length(x_ref)
    nu = length(u_ref)
    ∂f_∂x = zeros(nx,nx)
    ∂f_∂u = zeros(nx,nu)

    # ∂f_∂x and ∂f_∂u: copy from `generate_dynamics_partials` output
    ∂f_∂x[1,4] = 1.00000000000000
    ∂f_∂x[2,5] = 1.00000000000000
    ∂f_∂x[3,6] = 1.00000000000000
    ∂f_∂u[4,1] = 1/m
    ∂f_∂u[5,2] = 1/m
    ∂f_∂u[6,3] = 1/m
    ∂f_∂u[7,1] = T1/sqrt(T1^2 + T2^2 + T3^2)
    ∂f_∂u[7,2] = T2/sqrt(T1^2 + T2^2 + T3^2)
    ∂f_∂u[7,3] = T3/sqrt(T1^2 + T2^2 + T3^2)

    # ∂f_∂s: Evaluate nondilated nonlinear dynamics
    ∂f_∂s = dynamics_nonlinear(t_ref,x_ref,vcat(u_ref,1),params)

    # Package partials as linearized matrices
    A = s_ref*∂f_∂x
    B = Matrix([s_ref*∂f_∂u ∂f_∂s])
    Σ = []
    z = -(s_ref*∂f_∂x*x_ref + s_ref*∂f_∂u*u_ref)

    return(A,B,Σ,z)
end

function generate_dynamics_partials(params::Quad3DoFParams)

    # Symbols for differentiable quantities
    r1,r2,r3 = symbols("r1 r2 r3", real=true)
    v1,v2,v3 = symbols("v1 v2 v3", real=true)
    T1,T2,T3 = symbols("T1 T2 T3", real=true)
    intT     = symbols("intT"; real=true)
    
    # Symbols for constants
    g1,g2,g3 = symbols("g1 g2 g3", real=true)
    g = [g1;g2;g3]
    m = symbols("m", real=true)

    # Symbol canonicalization
    x = [r1;r2;r3;v1;v2;v3;intT]
    u = [T1;T2;T3]
    nx,nu = length(x),length(u) 

    # Evaluate nondilated nonlinear dynamics
    f = dynamics_nonlinear_nondilated(0,x,u,params)

    # Print out all partial elements
    print_sympy_partials(f,x,u)
end