#=
Continuous-to-discrete transcription utilities for LTI and LTV systems,
including exact variational discretization (ZOH/FOH) used by SCP solvers.
=#

"""
    c2d_nonlinear(TS_batch, X_batch, U_batch, dyn_nl, dyn_lin, disc; p_batch, cpu_parallel, gpu_parallel) -> (Ak, Bmk, Bpk, Σk, wk, δk)

Exact variational discretization of a CT-LTV system
``\\dot{x} = A(t)x + B(t)u + Σ(t)p`` over a batch of knot intervals, yielding
discrete updates ``x_{k+1} ≈ A_k x_k + B_k^- u_k + B_k^+ u_{k+1} + Σ_k p + w_k``.

# Arguments
- `TS_batch`: time spans `(t⁻, t⁺)` for each knot interval
- `X_batch`: boundary states `(x⁻, x⁺)` for each interval
- `U_batch`: control endpoint pairs `(u⁻, u⁺)` activated on each interval
- `dyn_nl`: factory `k -> (t,x,u,p) -> f` returning nonlinear dynamics for batch index `k`
- `dyn_lin`: factory `k -> (t,x,u,p) -> (A,B,Σ,w)` returning linearized dynamics for batch index `k`
- `disc::Int`: hold order (`0` = ZOH, `1` = FOH)
- `p_batch`: optional per-interval SCP parameter vectors (empty entries if unused)
- `cpu_parallel::Bool`: if `true`, integrate the ensemble with `EnsembleThreads`
- `gpu_parallel::Bool`: if `true`, integrate with the GPU ensemble backend

# Returns
- `Ak`: discrete state transition matrices per interval
- `Bmk`: discrete control matrices multiplying `u_k` (and ZOH `B`)
- `Bpk`: discrete control matrices multiplying `u_{k+1}` (FOH only; unused for ZOH)
- `Σk`: discrete parameter maps per interval
- `wk`: discrete affine / particular-solution terms per interval
- `δk`: state defects ``\\|x^+_{\\mathrm{prop}} - x^+_{\\mathrm{ref}}\\|`` per interval
"""
function c2d_nonlinear(
        TS_batch::Vector{Tuple{CReal,CReal}},
        X_batch::Vector{Tuple{CVector,CVector}},
        U_batch::Vector, # do not specify type further for generality
        dyn_nl::Function,
        dyn_lin::Function,
        disc::Int;
        p_batch::Vector{CVector}=Vector{CVector}(undef,0),
        cpu_parallel::Bool=true, 
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
    batchprob = EnsembleProblem(prob,prob_func=prob_func)

    # Solve the ODE problem
    if gpu_parallel
        sol = DiffEqGPU.solve(batchprob, GPUTsit5(), EnsembleGPUKernel(CUDA.CUDABackend()), trajectories=N)
    elseif cpu_parallel
        sol = OrdinaryDiffEq.solve(batchprob, Tsit5(), EnsembleThreads(), trajectories=N)
    else
        sol = OrdinaryDiffEq.solve(batchprob, Tsit5(), EnsembleSerial(), trajectories=N)
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

"""
    ode_nonlinear(t, z, u, p, nx, nu, np, nx2, nxnu, nxnp, dyn_nl, dyn_lin, disc; t_span) -> CVector

Evaluate the vectorized variational integrand used by [`c2d_nonlinear`](@ref).

# Arguments
- `t::Float64`: evaluation time
- `z::CVector`: concatenated integration state ``[x; \\mathrm{vec}(Φ_A); \\ldots]``
- `u::CVector`: control at time `t`
- `p::CVector`: parameter vector at time `t`
- `nx`, `nu`, `np`: state, control, and parameter dimensions
- `nx2`, `nxnu`, `nxnp`: precomputed products `nx^2`, `nx*nu`, `nx*np`
- `dyn_nl`: nonlinear dynamics `(t,x,u,p) -> f`
- `dyn_lin`: linearized dynamics `(t,x,u,p) -> (A,B,Σ,w)`
- `disc::Int`: hold order (`0` = ZOH, `1` = FOH)
- `t_span`: knot interval `(t⁻, t⁺)` (required for FOH blending)

# Returns
- `feval::CVector`: time derivative of the concatenated variational state
"""
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

"""
    add_traj_to_c2d_batch!(traj, TS_batch, X_batch, U_batch; disc=0, remove_zeros_from_traj=true) -> Vector{Int}

Append a reference [`Solution`](@ref) to the continuous-to-discrete batch buffers.

# Arguments
- `traj::Solution`: reference trajectory to append
- `TS_batch`: batch of time spans (mutated)
- `X_batch`: batch of boundary states (mutated)
- `U_batch`: batch of control endpoint pairs (mutated)
- `disc::Int`: hold order; if `0`, duplicates the final control for ZOH endpoints
- `remove_zeros_from_traj::Bool`: if `true`, replace exact zeros in `traj` for numerics

# Returns
- `idxs::Vector{Int}`: batch indices corresponding to intervals from this trajectory

# Notes
Mutates `TS_batch`, `X_batch`, and `U_batch`.
"""
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

"""
    c2d_LTI_affine(A_c, B_c, p_c, Δt, disc) -> (A, Bm, Bp, p)

Discretize continuous-time LTI affine dynamics at step `Δt`.

# Arguments
- `A_c::CMatrix`: continuous-time state matrix
- `B_c::CMatrix`: continuous-time control matrix
- `p_c::CVector`: continuous-time affine term
- `Δt::CReal`: discretization step size
- `disc::Int`: hold order (`0` = ZOH, `1` = FOH)

# Returns
- `A`: discrete state matrix
- `Bm`: discrete control matrix multiplying `u_k`
- `Bp`: discrete control matrix multiplying `u_{k+1}` (zeros for ZOH)
- `p`: discrete affine term
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

"""
    c2d_LTI_affine_zoh(A_c, B_c, p_c, Δt) -> (A, B, p)

ZOH discretization of ``\\dot{x} = A_c x + B_c u + p_c`` via matrix exponential.

# Arguments
- `A_c::CMatrix`: continuous-time state matrix
- `B_c::CMatrix`: continuous-time control matrix
- `p_c::CVector`: continuous-time affine term
- `Δt::CReal`: discretization step size

# Returns
- `A`: discrete state matrix
- `B`: discrete control matrix multiplying `u_k`
- `p`: discrete affine term
"""
function c2d_LTI_affine_zoh(A_c::CMatrix, B_c::CMatrix, p_c::CVector, Δt::CReal)::Tuple{CMatrix, CMatrix, CVector}
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

"""
    c2d_LTI_affine_foh(A_c, B_c, p_c, Δt) -> (A, Bm, Bp, p)

FOH discretization of ``\\dot{x} = A_c x + B_c u + p_c``.

# Arguments
- `A_c::CMatrix`: continuous-time state matrix
- `B_c::CMatrix`: continuous-time control matrix
- `p_c::CVector`: continuous-time affine term
- `Δt::CReal`: discretization step size

# Returns
- `A`: discrete state matrix
- `Bm`: discrete control matrix multiplying `u_k`
- `Bp`: discrete control matrix multiplying `u_{k+1}`
- `p`: discrete affine term
"""
function c2d_LTI_affine_foh(A_c::CMatrix, B_c::CMatrix, p_c::CVector, Δt::CReal)::Tuple{CMatrix, CMatrix, CMatrix, CVector}
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

"""
    devec(z, n, m, nm, off) -> Matrix

Reshape a contiguous block from a vectorized state into an `n×m` matrix.

# Arguments
- `z::Vector`: concatenated vector containing the block
- `n::Int`: number of rows of the reshaped block
- `m::Int`: number of columns of the reshaped block
- `nm::Int`: block length (`n * m`)
- `off::Int`: index offset before the block starts (exclusive)

# Returns
- `n×m` matrix formed from `z[off+1:off+nm]`
"""
function devec(z::Vector,n::Int,m::Int,nm::Int,off::Int)::Matrix
    return reshape(z[off+1:off+nm],n,m)
end