function dynamics_nonlinear(
    t::CReal,
    x::CVector,
    ν::CVector,
    params::Quad3DoFCageParams)::CVector

    # Dynamics matrices
    A = CMatrix([
        zeros(3,3) I(3);
        zeros(3,3) zeros(3,3)
    ])
    B = CMatrix([
        zeros(3,3);
        I(3)/params.mass
    ])
    p = CVector(vcat(zeros(3),params.g))

    # Compute nonlinear dynamics
    u = ν[1:end-1]
    s = ν[end]
    f = A*x + B*u + p
    z = s*f
    
    return z
end


function dynamics_linearized(
    t_ref::CReal,
    x_ref::CVector,
    ν_ref::CVector,
    params::Quad3DoFCageParams)::Tuple{CMatrix,CMatrix,CVector}

    # Parse reference control
    u_ref = ν_ref[1:end-1]
    s_ref = ν_ref[end]

    # Necessary constants for partials
    m = params.mass

    # Matrices to populate
    nx = length(x_ref)
    nu = length(u_ref)
    ∂f_∂x = zeros(nx,nx)
    ∂f_∂u = zeros(nx,nu)

    # ∂f_∂x and ∂f_∂u: copy from `generate_dynamics_partials` output
    ∂f_∂x[1,1] = 0
    ∂f_∂x[1,2] = 0
    ∂f_∂x[1,3] = 0
    ∂f_∂x[1,4] = 1.00000000000000
    ∂f_∂x[1,5] = 0
    ∂f_∂x[1,6] = 0
    ∂f_∂u[1,1] = 0
    ∂f_∂u[1,2] = 0
    ∂f_∂u[1,3] = 0
    ∂f_∂x[2,1] = 0
    ∂f_∂x[2,2] = 0
    ∂f_∂x[2,3] = 0
    ∂f_∂x[2,4] = 0
    ∂f_∂x[2,5] = 1.00000000000000
    ∂f_∂x[2,6] = 0
    ∂f_∂u[2,1] = 0
    ∂f_∂u[2,2] = 0
    ∂f_∂u[2,3] = 0
    ∂f_∂x[3,1] = 0
    ∂f_∂x[3,2] = 0
    ∂f_∂x[3,3] = 0
    ∂f_∂x[3,4] = 0
    ∂f_∂x[3,5] = 0
    ∂f_∂x[3,6] = 1.00000000000000
    ∂f_∂u[3,1] = 0
    ∂f_∂u[3,2] = 0
    ∂f_∂u[3,3] = 0
    ∂f_∂x[4,1] = 0
    ∂f_∂x[4,2] = 0
    ∂f_∂x[4,3] = 0
    ∂f_∂x[4,4] = 0
    ∂f_∂x[4,5] = 0
    ∂f_∂x[4,6] = 0
    ∂f_∂u[4,1] = 1/m
    ∂f_∂u[4,2] = 0
    ∂f_∂u[4,3] = 0
    ∂f_∂x[5,1] = 0
    ∂f_∂x[5,2] = 0
    ∂f_∂x[5,3] = 0
    ∂f_∂x[5,4] = 0
    ∂f_∂x[5,5] = 0
    ∂f_∂x[5,6] = 0
    ∂f_∂u[5,1] = 0
    ∂f_∂u[5,2] = 1/m
    ∂f_∂u[5,3] = 0
    ∂f_∂x[6,1] = 0
    ∂f_∂x[6,2] = 0
    ∂f_∂x[6,3] = 0
    ∂f_∂x[6,4] = 0
    ∂f_∂x[6,5] = 0
    ∂f_∂x[6,6] = 0
    ∂f_∂u[6,1] = 0
    ∂f_∂u[6,2] = 0
    ∂f_∂u[6,3] = 1/m

    # ∂f_∂s: Evaluate nondilated nonlinear dynamics
    ∂f_∂s = dynamics_nonlinear(t_ref,x_ref,vcat(u_ref,1),params)

    # Package partials as linearized matrices
    A = s_ref*∂f_∂x
    B = Matrix([s_ref*∂f_∂u ∂f_∂s])
    w = -(s_ref*∂f_∂x*x_ref + s_ref*∂f_∂u*u_ref)

    return(A,B,w)
end

function generate_dynamics_partials(params::Quad3DoFCageParams)

    # Symbols for differentiable quantities
    r1,r2,r3 = symbols("r1 r2 r3", real=true)
    v1,v2,v3 = symbols("v1 v2 v3", real=true)
    T1,T2,T3 = symbols("T1 T2 T3", real=true)

    # Symbols for constants
    g1,g2,g3 = symbols("g1 g2 g3", real=true)
    g = [g1;g2;g3]
    m = symbols("m", real=true)

    # Symbol canonicalization
    x = [r1;r2;r3;v1;v2;v3]
    u = [T1;T2;T3]
    nx,nu = length(x),length(u) 

    # Dynamics matrices
    A = Matrix([
        zeros(3,3) I(3);
        zeros(3,3) zeros(3,3)
    ])
    B = Matrix([
        zeros(3,3);
        I(3)/m
    ])
    p = Vector(vcat(zeros(3),g))

    # Evaluate nondilated nonlinear dynamics
    f = A*x + B*u + p

    # Print out all partial elements
    for i = 1:nx
        for j = 1:nx
            ∂fi_∂xj = diff(f[i],x[j])
            print("∂f_∂x[$(i),$(j)] = $(string(∂fi_∂xj))\n")
        end
        for j = 1:nu
            ∂fi_∂uj = diff(f[i],u[j])
            print("∂f_∂u[$(i),$(j)] = $(string(∂fi_∂uj))\n")
        end
    end
end