"""
Provides a variety of discretization techniques for both LTI and LTV systems.
"""

<<<<<<< HEAD
function c2d_nonlinear(
        TS_batch::Vector{Tuple{CReal,CReal}},
        X_batch::Vector{Tuple{CVector,CVector}},
        U_batch::Vector, # do not specify type further for generality
        dyn_nl::Function,
        dyn_lin::Function,
        disc::Int;
        p_batch::Vector{CVector}=Vector{CVector}(undef,0),
        cpu_parallel::Bool=false,
        gpu_parallel::Bool=false
    )::Tuple{Vector{CMatrix},Vector{CMatrix},Vector{CMatrix},Vector{CMatrix},Vector{CMatrix},Vector}
    # Integrate a continuous-time linear-time-varying (CT-LTV) system of the form:
    #     ̇x(t) = A(t)x(t) + B(t)u(t) + Σ(t)p
    # To obtain the DT-LTV discretization:
    #     x(k+1) ≈ A(k)x(k) + B(k)u(k) + Σ(k)p + z(k)
    #
    # Uses exact discretization, variational method (inverse-free)
    #
    # :in TS_batch: batch time-spans for each knot interval
    # :in X_batch: batch boundary states for each knot interval
    # :in U_batch: batch control parameters for each knot interval (contains all parameters activated between each knot point interval)
    # :in f(k) = dyn_nl(k)(t,x,u,p): Nonlinear dynamics function for kth batch index
    # :in (A,B,Σ,z)(k) = dyn_lin(k)(t,x,u,p): Linearized dynamics function for kth batch index
    # :in disc: Discretization hold order (0 = ZOH, 1 = FOH)
    # :in p_batch: reference SCP parameter signal (optional) (not the same thing as the ODEProblem parameters)
    # :in num_disc_steps: Number of integrator discretization steps (optional)
    # :in h_min: Minimum time step for integration (optional)

    # Set up empty p_batch if not provided
    if length(p_batch) == 0
        p_batch = Vector{CVector}(undef,length(X_batch))
        for k = 1:length(X_batch)
            p_batch[k] = CVector(undef,0)
        end
    end

    # Sizing variables
    N = length(X_batch)
    nx = length(X_batch[1][1])
    nu = length(U_batch[1][1])
    np = length(p_batch[1])
    nx2 = nx^2
    nxnu = nx*nu
    nxnp = nx*np

    # Initialize output matrices
    Ak = Vector{CMatrix}(undef,N)
    Bmk = Vector{CMatrix}(undef,N)
    Bpk = Vector{CMatrix}(undef,N) # Only used for disc == 1 (FOH)
    Σk = Vector{CMatrix}(undef,N)
    wk = Vector{CMatrix}(undef,N)
    δk = Vector{CReal}(undef,N)

    # Construct initial condition for system matrices (vectorized)
    if disc == 0
        S0 = vec(vcat(reshape(I(nx),:,1), zeros(nxnu), zeros(nxnp,1), zeros(nx))) # vec operation
    elseif disc == 1
        S0 = vec(vcat(reshape(I(nx),:,1), zeros(nxnu), zeros(nxnu), zeros(nxnp,1), zeros(nx))) # vec operation
    else
        error("Please select a valid discretization hold order.")
    end

    # Define the propagation function in terms of:
    #   ODE state vector `z` = [x; ΦA; ΦB; ΦΣ; Φw] (vectorized)
    #   Batch index `k` (for indexing time and control input parameters)
    #   Current time `t`
    function prop_fun!(dz,z,k,t)
       dz[:] = SVector{length(z)}(ode_nonlinear(t,z,optimal_controller(t,TS_batch[k],U_batch[k],disc),p_batch[k],nx,nu,np,nx2,nxnu,nxnp,dyn_nl(k),dyn_lin(k),disc;t_span=TS_batch[k]))
    end

    # Define the ODE problem for each batch index
    u0 = SVector{N}([vec(vcat(X_batch[k][1], S0)) for k=1:N])
    t_span = SVector{N}(TS_batch)
    p = SVector{N}(collect(1:N))
    prob = ODEProblem{true}(prop_fun!,u0[1],t_span[1],p[1])
    prob_func = (prob,k,_) -> remake(prob,u0=u0[k],tspan=t_span[k],p=p[k])
    batchprob = EnsembleProblem(prob,prob_func=prob_func,safetycopy=false)

    # Solve the ODE problem
    if gpu_parallel
        sol = DiffEqGPU.solve(batchprob, GPUTsit5(), EnsembleGPUKernel(CUDA.CUDABackend()), trajectories=N)
    elseif cpu_parallel
        sol = DifferentialEquations.solve(batchprob, Tsit5(), EnsembleThreads(), trajectories=N)
    else
        sol = DifferentialEquations.solve(batchprob, Tsit5(), EnsembleSerial(), trajectories=N)
    end

    # Extract propagated system matrices for the batch
    for k = 1:N
        # Extract final state vector for this time interval
        zf = sol.u[k][:,end]

        # Construct output matrices (de-vec operation)
        if disc == 0
            Ak[k]  = devec(zf,nx,nx,nx2,nx)
            Bmk[k] = devec(zf,nx,nu,nxnu,nx+nx2)
            Σk[k]  = devec(zf,nx,np,nxnp,nx+nx2+nxnu)
            wk[k]  = devec(zf,nx,1,nx,nx+nx2+nxnu+nxnp)
        elseif disc == 1
            Ak[k]  = devec(zf,nx,nx,nx2,nx)
            Bmk[k] = devec(zf,nx,nu,nxnu,nx+nx2)
            Bpk[k] = devec(zf,nx,nu,nxnu,nx+nx2+nxnu)
            Σk[k]  = devec(zf,nx,np,nxnp,nx+nx2+2*nxnu)
            wk[k]  = devec(zf,nx,1,nx,nx+nx2+2*nxnu+nxnp)
        end

        # Record defect
        δk[k] = norm(zf[1:nx] - X_batch[k][2])
    end
    
    return Ak,Bmk,Bpk,Σk,wk,δk
end

function ode_nonlinear(
        t::Float64,
        z::CVector,
        u::CVector,
        p::CVector,
        nx::Int,
        nu::Int,
        np::Int,
        nx2::Int,
        nxnu::Int,
        nxnp::Int,
        dyn_nl::Function,
        dyn_lin::Function,
        disc::Int;
        t_span::Tuple{CReal,CReal}=(0,0)
    )::CVector
    # Obtain the function evaluation of the vector-concatenated
    # integrand used by `c2d_nonlinear`
    #
    # :in t: evaluated time
    # :in z: evaluated integration state
    # :in u: evaluated control
    # :in p: evaluated parameter
    # :in nx,...,nxnp: sizing variables to reduce evaluations
    # :in f = dyn_nl(t,x,u,p): Nonlinear dynamics function
    # :in A,B,Σ,z = dyn_lin(t,x,u,p): Linearized dynamics function
    # :in disc: Discretization hold order (0 = ZOH, 1 = FOH)
    # :in t_span: time span of solution nodal segment

    x = z[1:nx]
    A,B,Σ,w = dyn_lin(t,x,u,p)
    fx = dyn_nl(t,x,u,p)
    if disc == 0
        ΦA = devec(z,nx,nx,nx2,nx)
        ΦB = devec(z,nx,nu,nxnu,nx+nx2)
        ΦΣ = devec(z,nx,np,nxnp,nx+nx2+nxnu)
        Φw = devec(z,nx,1,nx,nx+nx2+nxnu+nxnp)
        fA = A*ΦA
        fB = A*ΦB + B
        fΣ = isempty(p) ? [] : A*ΦΣ + Σ
        fw = A*Φw + w
        feval = [vec(fx);vec(fA);vec(fB);vec(fΣ);vec(fw)]
    elseif disc == 1
        ΦA  = devec(z,nx,nx,nx2,nx)
        ΦBm = devec(z,nx,nu,nxnu,nx+nx2)
        ΦBp = devec(z,nx,nu,nxnu,nx+nx2+nxnu)
        ΦΣ  = devec(z,nx,np,nxnp,nx+nx2+2*nxnu)
        Φw  = devec(z,nx,1,nx,nx+nx2+2*nxnu+nxnp)
        tkm = t_span[1]
        tkp = t_span[2]
        Δt = tkp - tkm
        fA = A*ΦA
        fBm = A*ΦBm + B*(tkp-t)/Δt
        fBp = A*ΦBp + B*(t-tkm)/Δt
        fΣ = isempty(p) ? [] : A*ΦΣ + Σ
        fw = A*Φw + w
        feval = [vec(fx);vec(fA);vec(fBm);vec(fBp);vec(fΣ);vec(fw)]
    end

    return feval
end

function add_traj_to_c2d_batch!(traj::Solution, TS_batch::Vector, X_batch::Vector, U_batch::Vector; disc::Int=0, remove_zeros_from_traj::Bool=true)::Vector{Int}
    # Extract from trajectory object
    t_ref = traj.t
    x_ref = traj.x
    u_ref = traj.u
    N = length(t_ref)

    # Conditionals
    if remove_zeros_from_traj # for numerical stability
        remove_ref_zeros!(x_ref, u_ref)
    end
    if disc == 0
        u_ref = [u_ref u_ref[:,end]]
    end

    # Add to C2D batch (see c2d_nonlinear)
    idx_start = length(TS_batch)
    append!(TS_batch, Vector([(t_ref[k],t_ref[k+1]) for k = 1:N-1]))
    append!(X_batch, Vector([(x_ref[:,k],x_ref[:,k+1]) for k = 1:N-1]))
    append!(U_batch, Vector([(u_ref[:,k],u_ref[:,k+1]) for k = 1:N-1]))

    # Get idxs corresponding to this traj
    idxs = collect(idx_start+1:idx_start+N-1)
    return idxs
end

=======
>>>>>>> a94024612f595e2e498cb1d9dc6cf7a44bd27ec5
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

<<<<<<< HEAD
=======
function c2d_nonlinear(
        t_ref::Vector,
        x_ref::Array,
        u_ref::Array,
        dyn_nl::Function,
        dyn_lin::Function,
        disc::Int;
        p_ref::Array=[],
        num_disc_steps::Int=10
    )::Tuple{Array,Array,Array,Array,Array,Vector}
    # Integrate a continuous-time linear-time-varying (CT-LTV) system of the form:
    #     ̇x(t) = A(t)x(t) + B(t)u(t) + Σ(t)p
    # To obtain the DT-LTV discretization:
    #     x(k+1) ≈ A(k)x(k) + B(k)u(k) + Σ(k)p + z(k)
    #
    # Uses exact discretization, variational method (inverse-free)
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
    nx2 = nx^2
    nxnu = nx*nu
    nxnp = nx*np

    Ak = zeros(nx,nx,N-1)
    Bmk = zeros(nx,nu,N-1)
    Bpk = zeros(nx,nu,N-1) # Only used for disc == 1 (FOH)
    Σk = zeros(nx,np,N-1)
    wk = zeros(nx,N-1)
    δk = zeros(N-1)

    if disc == 0
        f0 = vec(vcat(reshape(I(nx),:,1), zeros(nx*nu), zeros(nx*np,1), zeros(nx))) # vec operation
    elseif disc == 1
        f0 = vec(vcat(reshape(I(nx),:,1), zeros(nx*nu), zeros(nx*nu), zeros(nx*np,1), zeros(nx))) # vec operation
    else
        error("Please select a valid discretization hold order.")
    end

    prop_fun = (t,z,t_span) -> ode_nonlinear(t,z,optimal_controller(t,t_ref,u_ref,disc),p_ref,nx,nu,np,nx2,nxnu,nxnp,dyn_nl,dyn_lin,disc;t_span=t_span)

    t_ref_true = time_dilation_control_to_wall_clock_time(u_ref[end,:], t_ref, 1)
    h_min = 0.0001
    for k = 1:(N-1)

        # println("================ Iteration "*string(k)*" ================")

        # Setup 
        _z = vec(vcat(x_ref[:,k], f0))
        t_span = [t_ref[k],t_ref[k+1]]
        Δt_prop = max((1/num_disc_steps)*(t_span[2]-t_span[1]), h_min)

        # Propagate (RK4) and record defect (δk)
        prop_fun_ = (t,z) -> prop_fun(t,z,t_span)
        ~,Z = rk4_batch(prop_fun_,_z,t_span[1],t_span[2],Δt_prop)
        z = Z[:,end]
        δk[k] = norm(z[1:nx] - x_ref[:,k+1])

        # Construct output matrices for this timestep (de-vec operation)
        if disc == 0
            Ak[:,:,k]  = devec(z,nx,nx,nx2,nx)
            Bmk[:,:,k] = devec(z,nx,nu,nxnu,nx+nx2)
            Σk[:,:,k]  = devec(z,nx,np,nxnp,nx+nx2+nxnu)
            wk[:,k]    = devec(z,nx,1,nx,nx+nx2+nxnu+nxnp)
        elseif disc == 1
            Ak[:,:,k]  = devec(z,nx,nx,nx2,nx)
            Bmk[:,:,k] = devec(z,nx,nu,nxnu,nx+nx2)
            Bpk[:,:,k] = devec(z,nx,nu,nxnu,nx+nx2+nxnu)
            Σk[:,:,k]  = devec(z,nx,np,nxnp,nx+nx2+2*nxnu)
            wk[:,k]    = devec(z,nx,1,nx,nx+nx2+2*nxnu+nxnp)
        end

        x_prop = z[1:nx]
        isbad = x -> any(isinf.(x)) | any(isnan.(x))
        if norm(x_prop) >= 1e5 || isbad(x_prop)
            println("================ Iteration "*string(k)*" ================")
            println("x_k:")
            display(x_ref[:,k])
            println("x_k+1:")
            display(x_ref[:,k+1])
            println("x_prop:")
            display(x_prop)
            println("Z:")
            display(Z[1:nx,:]')
            # println("Time delta:")
            # display(t_ref_true[k+1] - t_ref_true[k])
            println("Ak:")
            display(Ak[:,:,k])
        end

        # if k == 16
        #     println("Iteration "*string(k))
        #     println("xm")
        #     display(x_ref[:,k])
        #     println("xp")
        #     display(x_ref[:,k+1])
        #     println("um")
        #     display(u_ref[:,k])
        #     println("up")
        #     display(u_ref[:,k+1])
        #     println("x_prop")
        #     display(z[1:nx])
        #     println("z")
        #     display(z)
        #     println("Ak")
        #     display(Ak[:,:,k])
        #     println("Bmk")
        #     display(Bmk[:,:,k])
        #     println("Bpk")
        #     display(Bpk[:,:,k])
        #     println("Sk")
        #     display(Σk[:,:,k])
        #     println("wk")
        #     display(wk[:,k])
        # end
    end
    
    # disc_failed = false
    # isbad = x -> isinf.(x) | isnan.(x)
    # if any(isbad.((Ak,Bmk,Bpk,Σk,wk))...)
    #     disc_failed = true
    # end

    return Ak,Bmk,Bpk,Σk,wk,δk
end

function ode_nonlinear(
        t::Float64,
        z::Vector,
        u::Vector,
        p::Vector,
        nx::Int,
        nu::Int,
        np::Int,
        nx2::Int,
        nxnu::Int,
        nxnp::Int,
        dyn_nl::Function,
        dyn_lin::Function,
        disc::Int;
        t_span::Array=[0,0]
    )::Vector
    # Obtain the function evaluation of the vector-concatenated
    # integrand used by `c2d_nonlinear`
    #
    # :in t: evaluated time
    # :in z: evaluated integration state
    # :in u: evaluated control
    # :in p: evaluated parameter
    # :in nx,...,nxnp: sizing variables to reduce evaluations
    # :in f = dyn_nl(t,x,u,p): Nonlinear dynamics function
    # :in A,B,Σ,z = dyn_lin(t,x,u,p): Linearized dynamics function
    # :in disc: Discretization hold order (0 = ZOH, 1 = FOH)
    # :in t_span: time span of solution nodal segment

    x = z[1:nx]
    A,B,Σ,w = dyn_lin(t,x,u,p)
    f0 = dyn_nl(t,x,u,p)
    if disc == 0
        ΦA = devec(z,nx,nx,nx2,nx)
        ΦB = devec(z,nx,nu,nxnu,nx+nx2)
        ΦΣ = devec(z,nx,np,nxnp,nx+nx2+nxnu)
        Φw = devec(z,nx,1,nx,nx+nx2+nxnu+nxnp)
        fA = A*ΦA
        fB = A*ΦB + B
        fΣ = isempty(p) ? [] : A*ΦΣ + Σ
        fw = A*Φw + w
        feval = [vec(f0);vec(fA);vec(fB);vec(fΣ);vec(fw)]
    elseif disc == 1
        ΦA  = devec(z,nx,nx,nx2,nx)
        ΦBm = devec(z,nx,nu,nxnu,nx+nx2)
        ΦBp = devec(z,nx,nu,nxnu,nx+nx2+nxnu)
        ΦΣ  = devec(z,nx,np,nxnp,nx+nx2+2*nxnu)
        Φw  = devec(z,nx,1,nx,nx+nx2+2*nxnu+nxnp)
        tkm = t_span[1]
        tkp = t_span[2]
        Δt = tkp - tkm
        fA = A*ΦA
        fBm = A*ΦBm + B*(tkp-t)/Δt
        fBp = A*ΦBp + B*(t-tkm)/Δt
        fΣ = isempty(p) ? [] : A*ΦΣ + Σ
        fw = A*Φw + w
        feval = [vec(f0);vec(fA);vec(fBm);vec(fBp);vec(fΣ);vec(fw)]
    end

    return feval
end

>>>>>>> a94024612f595e2e498cb1d9dc6cf7a44bd27ec5
function devec(z::Vector,n::Int,m::Int,nm::Int,off::Int)::Matrix
    return reshape(z[off+1:off+nm],n,m)
end