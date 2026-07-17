#=
Hand-coded / SymPy-generated Jacobians of CTCS-augmented HALO-scenario
dynamics (including drag and glideslope terms), plus symbolic codegen helper.
=#

"""
    evaluate_jacobians_sympy(t_ref, x_ref, u_ref, params::Quad3DoFHaloParams, targ_idx) -> (∂f_∂x, ∂f_∂u)

Evaluate nondilated CTCS dynamics Jacobians for the HALO scenario at a
reference point (drag, glideslope, thrust/velocity penalties, obstacles).

# Arguments
- `t_ref`: reference time `[s]` (unused by the autonomous model).
- `x_ref`: reference state vector at the linearization point.
- `u_ref`: reference thrust control vector.
- `params`: HALO scenario parameters (drag, glideslope, obstacles).
- `targ_idx`: target index (`0` uses a trivial glideslope reference for the trunk).

# Returns
- `∂f_∂x`: state Jacobian of nondilated CTCS dynamics.
- `∂f_∂u`: control Jacobian of nondilated CTCS dynamics.
"""
function evaluate_jacobians_sympy(
        t_ref::CReal,
        x_ref::CVector,
        u_ref::CVector,
        params::Quad3DoFHaloParams,
        targ_idx::Int)
        
    # Matrices to populate
    nx = length(x_ref)
    nu = length(u_ref)
    ∂f_∂x = zeros(nx,nx)
    ∂f_∂u = zeros(nx,nu)
    
    # Parameters
    r = x_ref[1:3]
    v = x_ref[4:6]
    intT = x_ref[7]
    T = u_ref[1:3]
    mass = params.mass
    g = params.g
    rho = params.ρ
    C_d = params.C_d
    S_A = params.S_A
    rho_min = params.ρ_min
    rho_max = params.ρ_max
    gamma_p = params.γ_p
    gamma_gs = params.γ_gs
    v_min_V = params.v_min_V
    v_max_V = params.v_max_V
    v_max_L = params.v_max_L
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
    if targ_idx == 0
        # make it trivial to satisfy glideslope constraint
        # by placing rf_gs directly below r in z-axis
        rf_gs = copy(r)
        rf_gs[3] -= 1
    else
        rf_gs = params.a.zf_targs[1:3,targ_idx]
    end

    if typeof(params) == Quad3DoFCageParams
        h_constant = params.h_constant
        ∂f_dx_83 = -2*h_constant + 2*r[3]
    else
        ∂f_dx_83 = 0
    end

    # ∂f_∂x and ∂f_∂u: copy from `generate_dynamics_partials` output (without obstacles)
    ∂f_∂x[1,4] = 1.00000000000000
    ∂f_∂x[2,5] = 1.00000000000000
    ∂f_∂x[3,6] = 1.00000000000000
    ∂f_∂x[4,4] = C_d*S_A*rho*v[1]^2/(2*sqrt(v[1]^2 + v[2]^2 + v[3]^2)) + C_d*S_A*rho*sqrt(v[1]^2 + v[2]^2 + v[3]^2)/2
    ∂f_∂x[4,5] = C_d*S_A*rho*v[1]*v[2]/(2*sqrt(v[1]^2 + v[2]^2 + v[3]^2))
    ∂f_∂x[4,6] = C_d*S_A*rho*v[1]*v[3]/(2*sqrt(v[1]^2 + v[2]^2 + v[3]^2))
    ∂f_∂x[5,4] = C_d*S_A*rho*v[1]*v[2]/(2*sqrt(v[1]^2 + v[2]^2 + v[3]^2))
    ∂f_∂x[5,5] = C_d*S_A*rho*v[2]^2/(2*sqrt(v[1]^2 + v[2]^2 + v[3]^2)) + C_d*S_A*rho*sqrt(v[1]^2 + v[2]^2 + v[3]^2)/2
    ∂f_∂x[5,6] = C_d*S_A*rho*v[2]*v[3]/(2*sqrt(v[1]^2 + v[2]^2 + v[3]^2))
    ∂f_∂x[6,4] = C_d*S_A*rho*v[1]*v[3]/(2*sqrt(v[1]^2 + v[2]^2 + v[3]^2))
    ∂f_∂x[6,5] = C_d*S_A*rho*v[2]*v[3]/(2*sqrt(v[1]^2 + v[2]^2 + v[3]^2))
    ∂f_∂x[6,6] = C_d*S_A*rho*v[3]^2/(2*sqrt(v[1]^2 + v[2]^2 + v[3]^2)) + C_d*S_A*rho*sqrt(v[1]^2 + v[2]^2 + v[3]^2)/2
    ∂f_∂x[8,1] = 2*(r[1] - rf_gs[1])*heaviside(-(1.0*r[3] - 1.0*rf_gs[3])/cos(gamma_gs) + sqrt((r[1] - rf_gs[1])^2 + (r[2] - rf_gs[2])^2 + (r[3] - rf_gs[3])^2))*max(0, -(1.0*r[3] - 1.0*rf_gs[3])/cos(gamma_gs) + sqrt((r[1] - rf_gs[1])^2 + (r[2] - rf_gs[2])^2 + (r[3] - rf_gs[3])^2))/sqrt((r[1] - rf_gs[1])^2 + (r[2] - rf_gs[2])^2 + (r[3] - rf_gs[3])^2)
    ∂f_∂x[8,2] = 2*(r[2] - rf_gs[2])*heaviside(-(1.0*r[3] - 1.0*rf_gs[3])/cos(gamma_gs) + sqrt((r[1] - rf_gs[1])^2 + (r[2] - rf_gs[2])^2 + (r[3] - rf_gs[3])^2))*max(0, -(1.0*r[3] - 1.0*rf_gs[3])/cos(gamma_gs) + sqrt((r[1] - rf_gs[1])^2 + (r[2] - rf_gs[2])^2 + (r[3] - rf_gs[3])^2))/sqrt((r[1] - rf_gs[1])^2 + (r[2] - rf_gs[2])^2 + (r[3] - rf_gs[3])^2)
    ∂f_∂x[8,3] = 2*((r[3] - rf_gs[3])/sqrt((r[1] - rf_gs[1])^2 + (r[2] - rf_gs[2])^2 + (r[3] - rf_gs[3])^2) - 1.0/cos(gamma_gs))*heaviside(-(1.0*r[3] - 1.0*rf_gs[3])/cos(gamma_gs) + sqrt((r[1] - rf_gs[1])^2 + (r[2] - rf_gs[2])^2 + (r[3] - rf_gs[3])^2))*max(0, -(1.0*r[3] - 1.0*rf_gs[3])/cos(gamma_gs) + sqrt((r[1] - rf_gs[1])^2 + (r[2] - rf_gs[2])^2 + (r[3] - rf_gs[3])^2))
    ∂f_∂x[8,4] = 2*v[1]*heaviside(-v_max_L + sqrt(v[1]^2 + v[2]^2))*max(0, -v_max_L + sqrt(v[1]^2 + v[2]^2))/sqrt(v[1]^2 + v[2]^2)
    ∂f_∂x[8,5] = 2*v[2]*heaviside(-v_max_L + sqrt(v[1]^2 + v[2]^2))*max(0, -v_max_L + sqrt(v[1]^2 + v[2]^2))/sqrt(v[1]^2 + v[2]^2)
    ∂f_∂x[8,6] = -2*heaviside(-v[3] + v_min_V)*max(0, -v[3] + v_min_V) + 2*heaviside(v[3] - v_max_V)*max(0, v[3] - v_max_V)
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
    
    return ∂f_∂x,∂f_∂u
end

"""
    generate_dynamics_partials_ctcs(params::Quad3DoFHaloParams)

Symbolically differentiate HALO CTCS dynamics and print nonzero partials for
codegen into `evaluate_jacobians_sympy`. Requires SymPy.

# Arguments
- `params`: HALO scenario parameters used to build the symbolic model.

# Returns
- none; partial derivatives are printed to stdout.
"""
function generate_dynamics_partials_ctcs(params::Quad3DoFHaloParams)

    # Symbols for differentiable quantities
    r  = [symbols("r[$(j)]", real=true) for j=1:3]
    v  = [symbols("v[$(j)]", real=true) for j=1:3]
    T  = [symbols("T[$(j)]", real=true) for j=1:3]
    ∫T = symbols("intT"; real=true)
    ∫ξ = symbols("intxi"; real=true)
    
    # Symbols for constants
    mass = symbols("mass", real=true)
    g = [symbols("g[$(j)]", real=true) for j=1:3]
    ρ = symbols("rho", real=true)
    C_d = symbols("C_d", real=true)
    S_A = symbols("S_A", real=true)
    ρ_min = symbols("rho_min", real=true)
    ρ_max = symbols("rho_max", real=true)
    γ_p = symbols("gamma_p", real=true)
    γ_gs = symbols("gamma_gs", real=true)
    v_min_V = symbols("v_min_V", real=true)
    v_max_V = symbols("v_max_V", real=true)
    v_max_L = symbols("v_max_L", real=true)
    rf_gs = [symbols("rf_gs[$(j)]", real=true) for j=1:3]
    
    # Construct a custom parameter object for these symbols
    # (fill in non-symbolic parameters w/ numerical data from original params)
    params_sympy = Quad3DoFHaloParams{Any,Any}(
        g,
        ρ,
        C_d,
        S_A,
        params.n_rotor,
        mass,
        ρ_min,
        ρ_max,
        params.drag_term_enabled,
        params.ϵ_subopt,
        γ_gs,
        γ_p,
        v_min_V,
        v_max_V,
        v_max_L,
        params.n_obstacles,
        params.R_obstacles,
        params.p_obstacles,
        params.H_obstacles,
        params.n_targs_min,
        params.n_targs_max,
        params.R_targs,
        params.R_targs_min,
        params.p_targs,
        params.w_des,
        params.a,
        params.w_obj_decay_factor
    )

    # Symbol canonicalization
    x = [r;v;∫T;∫ξ]
    u = T

    # Evaluate nondilated nonlinear dynamics and CTCS state
    f_3DoF = dynamics_nonlinear_nondilated(0,x,u,params_sympy)
    ξ,_,_ = prob_constraints_eval(x,u,params_sympy,1;sympy=true,obstacles=false,rf_gs=rf_gs) # CTCS violation
    f = [f_3DoF;ξ]

    # Print out all partial elements
    print_sympy_partials(f,x,u)
end