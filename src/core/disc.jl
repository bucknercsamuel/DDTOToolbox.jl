"""
Provides a variety of discretization techniques for both LTI and LTV systems.
"""

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

function c2d_nonlinear_varying_zoh(
    ref_traj::Solution,
    dyn_nl::Function,
    dyn_lin::Function
)::Tuple{Array,Array,Array,CVector}

    t_ref = ref_traj.t    
    x_ref = ref_traj.x
    u_ref = ref_traj.u

    n,N = size(x_ref)
    m = size(u_ref,1)

    Ak = zeros(n,n,N-1)
    Bk = zeros(n,m,N-1)
    wk = zeros(n,N-1)
    δk = zeros(N-1)

    ABw0 = vec(vcat(reshape(I(n),:,1), zeros(n*m), zeros(n)))

    h_min = 0.0001
    for k = 1:(N-1)

        _z = vec(vcat(x_ref[:,k], ABw0))
        t_span = [t_ref[k],t_ref[k+1]]
        Δt_prop = max((1/40)*(t_span[2]-t_span[1]), h_min)

        prop_fun = (t,z) -> nonlinear_varying_zoh_ode(t,z,optimal_controller(t,t_ref,u_ref),n,dyn_nl,dyn_lin)
        ~,Z = rk4(prop_fun,_z,t_span[1],t_span[2],Δt_prop)
        z_ = Z[:,end]

        # display(_z[1:n])
        # display(z_[1:n])
        # display(x_ref[:,k+1])
        δk[k] = norm(z_[1:n] - x_ref[:,k+1])

        Ak_ = reshape(z_[n+1:n+n^2],n,n)
        Ak[:,:,k] = Ak_
        Bk[:,:,k] = Ak_ * reshape(z_[n+n^2+1:n+n^2+n*m],n,m)
        wk[:,k] = Ak_ * reshape(z_[n+n^2+n*m+1:n+n^2+n*m+n],n,1)
    end
    
    display(max(δk...))
    return Ak,Bk,wk,δk
end

function nonlinear_varying_zoh_ode(
    t::CReal,
    z::CVector,
    u::CVector,
    n::Int,
    dyn_nl::Function,
    dyn_lin::Function
)::CVector

    x = z[1:n]
    ΦA = reshape(z[n+1:n+n^2],n,n)
    ΦAinv = inv(ΦA)

    A,B,w = dyn_lin(t,x,u)

    f0 = reshape(dyn_nl(t,x,u),:,1)
    f1 = reshape(A*ΦA,:,1)
    f2 = reshape(ΦAinv*B,:,1)
    f3 = reshape(ΦAinv*w,:,1)

    feval = vec(vcat(f0,f1,f2,f3))
    return feval
end