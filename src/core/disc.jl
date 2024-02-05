"""
Provides a variety of discretization techniques for both LTI and LTV systems.
"""

function c2d_LTI_affine(A_c::CMatrix, B_c::CMatrix, p_c::CVector, Δt::CReal, disc::Int)::Tuple{CMatrix, CMatrix, CMatrix, CVector}
    if disc == 0
        A,Bm,p = c2d_LTI_affine_zoh(A_c,B_c,p_c,Δt)
        Bp = zeros(size(Bm))
    elseif disc == 1
        A,Bm,Bp,p = c2d_LTI_affine_foh(A_c,B_c,p_c,Δt)
    else
        error("Please select a valid discretization hold order.")
    end
    return A,Bm,Bp,p
end

function c2d_LTI_affine_zoh(A_c::CMatrix, B_c::CMatrix, p_c::CVector, Δt::CReal)::Tuple{CMatrix, CMatrix, CVector}
    # Discretize LTI vehicle dynamics at Δt time step using zeroth-order
    # hold (ZOH) of the affine state-space form:
    #       dx/dt = A*x + B*u + p
    #
    # :in A_c: Continuous-time A matrix
    # :in B_c: Continuous-time B matrix
    # :in p_c: Continuous-time p (affine) vector
    # :in Δt: the discretization time step.
    # :out: a tuple (A,B,p) for the discrete-time update equation
    #         x_{k+1} = A*x_k + B*u_k + p

    n,m = size(B_c)

    _M = exp(CMatrix([
        A_c B_c p_c;
        zeros(m+1,n+m+1)
    ])*Δt)

    A = _M[1:n,1:n]
    B = _M[1:n,n+1:n+m]
    p = _M[1:n,n+m+1]
    return (A,B,p)
end

function c2d_LTI_affine_foh(A_c::CMatrix, B_c::CMatrix, p_c::CVector, Δt::CReal)::Tuple{CMatrix, CMatrix, CMatrix, CVector}
    # Discretize LTI vehicle dynamics at Δt time step using first-order
    # hold (FOH) of the affine state-space form:
    #       dx/dt = A*x + B*u + p
    #
    # :in A_c: Continuous-time A matrix
    # :in B_c: Continuous-time B matrix
    # :in p_c: Continuous-time p (affine) vector
    # :in Δt: the discretization time step.
    # :out: a tuple (A,Bm,Bp,p) for the discrete-time update equation
    #         x_{k+1} = A*x_k + Bm*u_k + Bp*u_{k+1} + p

    n,m = size(B_c)

    h   = Δt
    Φ   = I(n) + h*A_c 
    Γ   = (h*I(n) + (h^2/2)*A_c)*B_c
    Γ_1 = ((h^2/2)*I(n) + (h^3/6)*A_c)*B_c

    A  = Φ
    Bm = Γ - (1/Δt)*Γ_1
    Bp = (1/Δt)*Γ_1
    p  = (Δt*I(n) + (Δt^2/2)*A_c)*p_c

    return (A, Bm, Bp, p)
end

function c2d_nonlinear(
        t_ref::Vector,
        x_ref::Array,
        u_ref::Array,
        dyn_nl::Function,
        dyn_lin::Function,
        disc::Int;
        p_ref::Array=[],
        num_disc_steps::Int=10
    )::Tuple{Array,Array,Array,Array,Array,Vector,Bool}
    # Integrate a continuous-time linear-time-varying (CT-LTV) system of the form:
    #     ̇x(t) = A(t)x(t) + B(t)u(t) + Σ(t)p
    # To obtain the DT-LTV discretization:
    #     x(k+1) ≈ A(k)x(k) + B(k)u(k) + Σ(k)p + z(k)
    #
    # :in t_ref: reference time signal
    # :in x_ref: reference state signal
    # :in u_ref: reference control signal
    # :in f = dyn_nl(t,x,u,p): Nonlinear dynamics function
    # :in A,B,Σ,z = dyn_lin(t,x,u,p): Linearized dynamics function
    # :in disc: Discretization hold order (0 = ZOH, 1 = FOH)
    # :in p_ref: reference parameter signal (optional)
    # :in num_disc_steps: Number of integrator discretization steps (optional)

    nx,N = size(x_ref)
    nu = size(u_ref,1)
    np = size(p_ref,1)

    Ak = zeros(nx,nx,N-1)
    Bmk = zeros(nx,nu,N-1)
    Bpk = zeros(nx,nu,N-1) # Only used for disc == 1 (FOH)
    Σk = zeros(nx,np,N-1)
    zk = zeros(nx,N-1)
    δk = zeros(N-1)

    if disc == 0
        f0 = vec(vcat(reshape(I(nx),:,1), zeros(nx*nu), zeros(nx*np,1), zeros(nx))) # vec operation
    elseif disc == 1
        f0 = vec(vcat(reshape(I(nx),:,1), zeros(nx*nu), zeros(nx*nu), zeros(nx*np,1), zeros(nx))) # vec operation
    else
        error("Please select a valid discretization hold order.")
    end

    prop_fun = (t,z,t_span) -> ode_nonlinear(t,z,optimal_controller(t,t_ref,u_ref,disc),p_ref,nx,dyn_nl,dyn_lin,disc;t_span=t_span)

    h_min = 0.0001
    disc_failed = false
    for k = 1:(N-1)

        # Setup 
        _f = vec(vcat(x_ref[:,k], f0))
        t_span = [t_ref[k],t_ref[k+1]]
        Δt_prop = max((1/num_disc_steps)*(t_span[2]-t_span[1]), h_min)

        # Propagate (RK4) and record defect (δk)
        prop_fun_ = (t,f) -> prop_fun(t,f,t_span)
        ~,F = rk4(prop_fun_,_f,t_span[1],t_span[2],Δt_prop)
        f_ = F[:,end]
        δk[k] = norm(f_[1:nx] - x_ref[:,k+1])

        # Construct output matrices for this timestep (de-vec operation)
        Ak_ = zeros(nx,nx)
        try
            Ak_ = inv(reshape(f_[nx+1:nx+nx^2],nx,nx))
        catch e
            println("! Discretization Error: discretization inverse failed!")
            disc_failed = true
        end
        Ak[:,:,k] = Ak_
        if disc == 0
            Bmk[:,:,k] =                   Ak_*reshape(f_[nx+nx^2+1               : nx+nx^2+nx*nu           ],nx,nu)
            Σk[:,:,k]  = !isempty(p_ref) ? Ak_*reshape(f_[nx+nx^2+nx*nu+1         : nx+nx^2+nx*nu+nx*np     ],nx,1) : zeros(nx,np,1)
            zk[:,k]    =                   Ak_*reshape(f_[nx+nx^2+nx*nu+nx+1      : nx+nx^2+nx*nu+nx+nx     ],nx,1)
        elseif disc == 1
            Bmk[:,:,k] =                   Ak_*reshape(f_[nx+nx^2+1               : nx+nx^2+nx*nu           ],nx,nu)
            Bpk[:,:,k] =                   Ak_*reshape(f_[nx+nx^2+nx*nu+1         : nx+nx^2+2*nx*nu         ],nx,nu)
            Σk[:,:,k]  = !isempty(p_ref) ? Ak_*reshape(f_[nx+nx^2+2*nx*nu+1       : nx+nx^2+2*nx*nu+nx*np   ],nx,1) : zeros(nx,np,1)
            zk[:,k]    =                   Ak_*reshape(f_[nx+nx^2+2*nx*nu+nx*np+1 : nx+nx^2+2*nx*nu+nx*np+nx],nx,1)
        end

        # if max(norm.(δk)...) > 1e-4
        #     println("Warning: Integration defect is non-trivial (maxδk = $(max(δk...)))!")
        # end
    end

    isbad(x) = isnan(x) || isinf(x)
    if any(isbad.(Ak)) || any(isbad.(Bmk)) || any(isbad.(Bpk)) || any(isbad.(Σk)) || any(isbad.(zk))
        println("! Discretization Error: Integrated elements contain Inf or NaN!")
        disc_failed = true
    end

    return Ak,Bmk,Bpk,Σk,zk,δk,disc_failed
end

function ode_nonlinear(
        t::Float64,
        f::Vector,
        u::Vector,
        p::Vector,
        nx::Int,
        dyn_nl::Function,
        dyn_lin::Function,
        disc::Int;
        t_span::Array=[0,0]
    )::Vector
    # Obtain the function evaluation of the vector-concatenated
    # integrand used by `c2d_nonlinear`
    #
    # :in t: evaluated time
    # :in x: evaluated state
    # :in u: evaluated control
    # :in p: evaluated parameter
    # :in f = dyn_nl(t,x,u,p): Nonlinear dynamics function
    # :in A,B,Σ,z = dyn_lin(t,x,u,p): Linearized dynamics function
    # :in disc: Discretization hold order (0 = ZOH, 1 = FOH)
    # :in t_span: time span of solution nodal segment

    x = f[1:nx]
    ΦAinv = reshape(f[nx+1:nx+nx^2],nx,nx)

    A,B,Σ,z = dyn_lin(t,x,u,p)
    if disc == 0
        f0 = dyn_nl(t,x,u,p)
        f1 = -ΦAinv*A
        f2 = ΦAinv*B
        f3 = !isempty(p) ? ΦAinv*Σ : []
        f4 = ΦAinv*z
        feval = [vec(f0);vec(f1);vec(f2);vec(f3);vec(f4)]
    elseif disc == 1
        tkm = t_span[1]
        tkp = t_span[2]
        Δt = tkp - tkm
        f0 = dyn_nl(t,x,u,p)
        f1 = -ΦAinv*A
        f2 = ΦAinv*B*(tkp-t)/Δt
        f3 = ΦAinv*B*(t-tkm)/Δt
        f4 = !isempty(p) ? ΦAinv*Σ : []
        f5 = ΦAinv*z
        feval = [vec(f0);vec(f1);vec(f2);vec(f3);vec(f4);vec(f5)]
    end

    return feval
end