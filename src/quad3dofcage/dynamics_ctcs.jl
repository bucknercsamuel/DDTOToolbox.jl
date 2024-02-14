heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), one(x)) # needed for symbolic maximum differentiation

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

    # Parse reference control
    u_ref = ν_ref[1:end-1]
    s_ref = ν_ref[end]

    # Matrices to populate
    nx = length(x_ref)
    nu = length(u_ref)
    ∂f_∂x = zeros(nx,nx)
    ∂f_∂u = zeros(nx,nu)

    # Parameters
    r = x_ref[1:3]
    v = x_ref[4:6]
    intT = x_ref[7]
    intxi = x_ref[8]
    T = u_ref[1:3]
    mass = params.mass
    g = params.g
    rho_min = params.ρ_min
    rho_max = params.ρ_max
    gamma_p = params.γ_p
    v_max_V = params.v_max_V
    v_max_L = params.v_max_L
    h_constant = params.h_constant
    x_arena_lims = params.x_arena_lims
    y_arena_lims = params.y_arena_lims
    z_arena_lims = params.z_arena_lims
    R_obstacles = zeros(params.n_obstacles)
    p_obstacles = zeros(3,params.n_obstacles)
    H_obstacles = [zeros(3,3) for _=1:params.n_obstacles]
    for o = 1:params.n_obstacles
        R_obstacles[o] = params.R_obstacles[o]
        for j = 1:3
            p_obstacles[j,o] = params.p_obstacles[j,o]
            for k = 1:3
                H_obstacles[o][j,k] = params.H_obstacles[o][j,k]
            end
        end
    end

    # ∂f_∂x and ∂f_∂u: copy from `generate_dynamics_partials` output (without obstacles)
    ∂f_∂x[1,4] = 1.00000000000000
    ∂f_∂x[2,5] = 1.00000000000000
    ∂f_∂x[3,6] = 1.00000000000000
    ∂f_∂x[8,3] = -2*h_constant + 2*r[3]
    ∂f_∂x[8,4] = 2*v[1]*heaviside(-v_max_L + sqrt(v[1]^2 + v[2]^2))*max(0, -v_max_L + sqrt(v[1]^2 + v[2]^2))/sqrt(v[1]^2 + v[2]^2)
    ∂f_∂x[8,5] = 2*v[2]*heaviside(-v_max_L + sqrt(v[1]^2 + v[2]^2))*max(0, -v_max_L + sqrt(v[1]^2 + v[2]^2))/sqrt(v[1]^2 + v[2]^2)
    ∂f_∂u[4,1] = 1/mass
    ∂f_∂u[5,2] = 1/mass
    ∂f_∂u[6,3] = 1/mass
    ∂f_∂u[7,1] = T[1]/sqrt(T[1]^2 + T[2]^2 + T[3]^2)
    ∂f_∂u[7,2] = T[2]/sqrt(T[1]^2 + T[2]^2 + T[3]^2)
    ∂f_∂u[7,3] = T[3]/sqrt(T[1]^2 + T[2]^2 + T[3]^2)
    ∂f_∂u[8,1] = 2*T[1]*heaviside(-rho_max + sqrt(T[1]^2 + T[2]^2 + T[3]^2))*max(0, -rho_max + sqrt(T[1]^2 + T[2]^2 + T[3]^2))/sqrt(T[1]^2 + T[2]^2 + T[3]^2) - 2*T[1]*heaviside(rho_min - sqrt(T[1]^2 + T[2]^2 + T[3]^2))*max(0, rho_min - sqrt(T[1]^2 + T[2]^2 + T[3]^2))/sqrt(T[1]^2 + T[2]^2 + T[3]^2) + 2*T[1]*heaviside(-1.0*T[3]/cos(gamma_p) + sqrt(T[1]^2 + T[2]^2 + T[3]^2))*max(0, -1.0*T[3]/cos(gamma_p) + sqrt(T[1]^2 + T[2]^2 + T[3]^2))/sqrt(T[1]^2 + T[2]^2 + T[3]^2)
    ∂f_∂u[8,2] = 2*T[2]*heaviside(-rho_max + sqrt(T[1]^2 + T[2]^2 + T[3]^2))*max(0, -rho_max + sqrt(T[1]^2 + T[2]^2 + T[3]^2))/sqrt(T[1]^2 + T[2]^2 + T[3]^2) - 2*T[2]*heaviside(rho_min - sqrt(T[1]^2 + T[2]^2 + T[3]^2))*max(0, rho_min - sqrt(T[1]^2 + T[2]^2 + T[3]^2))/sqrt(T[1]^2 + T[2]^2 + T[3]^2) + 2*T[2]*heaviside(-1.0*T[3]/cos(gamma_p) + sqrt(T[1]^2 + T[2]^2 + T[3]^2))*max(0, -1.0*T[3]/cos(gamma_p) + sqrt(T[1]^2 + T[2]^2 + T[3]^2))/sqrt(T[1]^2 + T[2]^2 + T[3]^2)
    ∂f_∂u[8,3] = 2*T[3]*heaviside(-rho_max + sqrt(T[1]^2 + T[2]^2 + T[3]^2))*max(0, -rho_max + sqrt(T[1]^2 + T[2]^2 + T[3]^2))/sqrt(T[1]^2 + T[2]^2 + T[3]^2) - 2*T[3]*heaviside(rho_min - sqrt(T[1]^2 + T[2]^2 + T[3]^2))*max(0, rho_min - sqrt(T[1]^2 + T[2]^2 + T[3]^2))/sqrt(T[1]^2 + T[2]^2 + T[3]^2) + 2*(T[3]/sqrt(T[1]^2 + T[2]^2 + T[3]^2) - 1.0/cos(gamma_p))*heaviside(-1.0*T[3]/cos(gamma_p) + sqrt(T[1]^2 + T[2]^2 + T[3]^2))*max(0, -1.0*T[3]/cos(gamma_p) + sqrt(T[1]^2 + T[2]^2 + T[3]^2))

    # Manually add obstacle Jacobian terms (SymPy cannot produce them efficiently)
    for o = 1:params.n_obstacles
        H = params.H_obstacles[o]
        p = params.p_obstacles[:,o]
        R = params.R_obstacles[o]
        ∂f_∂x[8,1:3] += -max(0,2*(R-norm(H*(r-p))))*(H*(H*(r-p))/norm(H*(r-p)))
    end

    # ∂f_∂s: Evaluate nondilated nonlinear dynamics
    ∂f_∂s = dynamics_nonlinear_ctcs(t_ref,x_ref,vcat(u_ref,1),params)

    # Package partials as linearized matrices
    A = s_ref*∂f_∂x
    B = Matrix([s_ref*∂f_∂u ∂f_∂s])
    Σ = []
    z = -(s_ref*∂f_∂x*x_ref + s_ref*∂f_∂u*u_ref)

    return(A,B,Σ,z)
end

function generate_dynamics_partials_ctcs(params::Quad3DoFCageParams)

    # Symbols for differentiable quantities
    r  = [symbols("r[$(j)]", real=true) for j=1:3]
    v  = [symbols("v[$(j)]", real=true) for j=1:3]
    T  = [symbols("T[$(j)]", real=true) for j=1:3]
    ∫T = symbols("intT"; real=true)
    ∫ξ = symbols("intxi"; real=true)
    
    # Symbols for constants
    mass = symbols("mass", real=true)
    g = [symbols("g[$(j)]", real=true) for j=1:3]
    ρ_min = symbols("rho_min", real=true)
    ρ_max = symbols("rho_max", real=true)
    γ_p = symbols("gamma_p", real=true)
    v_max_V = symbols("v_max_V", real=true)
    v_max_L = symbols("v_max_L", real=true)
    h_constant = symbols("h_constant", real=true)
    x_arena_lims = [symbols("x_arena_lims[$(j)]", real=true) for j=1:2]
    y_arena_lims = [symbols("y_arena_lims[$(j)]", real=true) for j=1:2]
    z_arena_lims = [symbols("z_arena_lims[$(j)]", real=true) for j=1:2]
    R_obstacles = Vector(undef,params.n_obstacles)
    p_obstacles = Matrix(undef,3,params.n_obstacles)
    H_obstacles = Vector{Matrix}(undef,params.n_obstacles)
    for o = 1:params.n_obstacles
        R_obstacles[o] = symbols("R_obstacles[$(o)]", real=true)
        H_obstacles[o] = Matrix(undef,3,3)
        for j = 1:3
            p_obstacles[j,o] = symbols("p_obstacles[$(j){}$(o)]", real=true)[1]
            for k = 1:3
                H_obstacles[o][j,k] = symbols("H_obstacles[$(o)][$(j){}$(k)]", real=true)[1]
            end
        end
    end

    # Construct a custom parameter object for these symbols
    # (fill in non-symbolic parameters w/ numerical data from original params)
    params_sympy = Quad3DoFCageParams{Any,Any}(
        g,
        params.ρ,
        params.n_rotor,
        mass,
        ρ_min,
        ρ_max,
        γ_p,
        v_max_V,
        v_max_L,
        h_constant,
        params.n_obstacles,
        R_obstacles,
        p_obstacles,
        H_obstacles,
        x_arena_lims,
        y_arena_lims,
        z_arena_lims,
        params.a
    )

    # Symbol canonicalization
    x = [r;v;∫T;∫ξ]
    u = T
    nx,nu = length(x),length(u) 

    # Evaluate nondilated nonlinear dynamics
    A,B,p = dynamics_linear(params_sympy)
    f_3DoF = A*x[1:6] + B*u + p

    # Additional state derivatives
    Tnorm = norm(u) # thrust integral
    ξ,_,_ = path_constraint_eval(x,u,params_sympy; sympy=true, obstacles=false) # CTCS violation
    f = [f_3DoF;Tnorm;ξ]

    # Print out all partial elements
    for i = 1:nx
        for j = 1:nx
            ∂fi_∂xj = diff(f[i],x[j])
            if ∂fi_∂xj != 0
                print(process_sympy_string("∂f_∂x[$(i),$(j)] = $(string(∂fi_∂xj))\n"))
            end
        end
    end
    for i = 1:nx
        for j = 1:nu
            ∂fi_∂uj = diff(f[i],u[j])
            if ∂fi_∂uj != 0
                print(process_sympy_string("∂f_∂u[$(i),$(j)] = $(string(∂fi_∂uj))\n"))
            end
        end
    end
end