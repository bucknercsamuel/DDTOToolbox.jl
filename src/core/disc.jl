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

function c2d_nonlinear(
    ref_traj::Solution,
    dyn_nl::Function,
    dyn_lin::Function,
    disc::Int;
    max_steps=10
)::Tuple{Array,Array,Array,Array,CVector}

    t_ref = ref_traj.t    
    x_ref = ref_traj.x
    u_ref = ref_traj.u

    n,N = size(x_ref)
    m = size(u_ref,1)

    Ak = zeros(n,n,N-1)
    Bmk = zeros(n,m,N-1)
    Bpk = zeros(n,m,N-1) # Only used for disc == 1 (FOH)
    wk = zeros(n,N-1)
    δk = zeros(N-1)

    if disc == 0
        z0 = vec(vcat(reshape(I(n),:,1), zeros(n*m), zeros(n))) # vec operation
    elseif disc == 1
        z0 = vec(vcat(reshape(I(n),:,1), zeros(n*m), zeros(n*m), zeros(n))) # vec operation       
    else
        error("Please select a valid discretization hold order.")
    end

    prop_fun = (t,z,t_span) -> ode_nonlinear(t,z,optimal_controller(t,t_ref,u_ref,disc),n,dyn_nl,dyn_lin,disc;t_span=t_span)

    h_min = 0.0001
    for k = 1:(N-1)

        # Setup 
        _z = vec(vcat(x_ref[:,k], z0))
        t_span = [t_ref[k],t_ref[k+1]]
        Δt_prop = max((1/max_steps)*(t_span[2]-t_span[1]), h_min)

        # Propagate and record defect
        prop_fun_ = (t,z) -> prop_fun(t,z,t_span)
        ~,Z = rk4(prop_fun_,_z,t_span[1],t_span[2],Δt_prop)
        z_ = Z[:,end]
        δk[k] = norm(z_[1:n] - x_ref[:,k+1])

        # Construct output matrices for this timestep (de-vec operation)
        Ak_ = reshape(z_[n+1:n+n^2],n,n)
        Ak[:,:,k] = Ak_
        if disc == 0
            Bmk[:,:,k] = Ak_*reshape(z_[n+n^2+1     : n+n^2+n*m  ],n,m)
            wk[:,k]    = Ak_*reshape(z_[n+n^2+n*m+1 : n+n^2+n*m+n],n,1)
        elseif disc == 1
            Bmk[:,:,k] = Ak_*reshape(z_[n+n^2+1       : n+n^2+n*m    ],n,m)
            Bpk[:,:,k] = Ak_*reshape(z_[n+n^2+n*m+1   : n+n^2+2*n*m  ],n,m)
            wk[:,k]    = Ak_*reshape(z_[n+n^2+2*n*m+1 : n+n^2+2*n*m+n],n,1)
        end
    end
    
    # display(max(δk...))
    return Ak,Bmk,Bpk,wk,δk
end

function ode_nonlinear(
    t::CReal,
    z::CVector,
    u::CVector,
    n::Int,
    dyn_nl::Function,
    dyn_lin::Function,
    disc::Int;
    t_span::Array=[0,0]
)::CVector

    x = z[1:n]
    ΦA = reshape(z[n+1:n+n^2],n,n)
    ΦAinv = inv(ΦA)

    A,B,w = dyn_lin(t,x,u)
    if disc == 0
        f0 = reshape(dyn_nl(t,x,u),:,1)
        f1 = reshape(A*ΦA,:,1)
        f2 = reshape(ΦAinv*B,:,1)
        f3 = reshape(ΦAinv*w,:,1)
        feval = vec(vcat(f0,f1,f2,f3))
    elseif disc == 1
        tkm = t_span[1]
        tkp = t_span[2]
        Δt = tkp - tkm
        f0 = reshape(dyn_nl(t,x,u),:,1)
        f1 = reshape(A*ΦA,:,1)
        f2 = reshape(ΦAinv*B*(tkp-t)/Δt,:,1)
        f3 = reshape(ΦAinv*B*(t-tkm)/Δt,:,1)
        f4 = reshape(ΦAinv*w,:,1)
        feval = vec(vcat(f0,f1,f2,f3,f4))
    end

    return feval
end