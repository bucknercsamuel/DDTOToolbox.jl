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
    ξ,_,_ = prob_constraints_eval(x,u,params_sympy; sympy=true, obstacles=false) # CTCS violation
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

function generate_dynamics_partials(params::Quad3DoFCageParams)

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
    A,B,p = dynamics_linear_noaugment(params)
    f_3DoF = A*x[1:6] + B*u + p
    f = [f_3DoF; norm(u)]

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