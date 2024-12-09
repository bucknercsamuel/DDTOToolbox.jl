using DDTOSCP

# Sample configuration
function SampleConfig()::Quad3DoFHaloParams{CReal,Int}
    # Load default params
    params = Quad3DoFHaloParams()

    # Configurables
    r0 = [0,0,150] # [m] Initial position (NED frame)
    v0 = [0,0,0]   # [m/s] Initial velocity (NED frame)
    params.a.n_targs = 4

    # Set sample boundary conditions for n_targs_max = 4 targets
    params.a.z0 = [
        r0;
        v0;
        0;
        Inf
    ]
    params.a.u0 = [
        0;
        0;
        0;
        Inf
    ]
    
    # Set terminal condition targets to be on a circle of radius 100m, make sure to include position and velocity in R3
    zf_targs = [[100*cos(2*pi*k/params.a.n_targs); 100*sin(2*pi*k/params.a.n_targs); 0; 0; 0; 0; Inf] for k = 1:params.a.n_targs]
    params.a.zf_targs = reshape(hcat(zf_targs...), params.a.nx, params.a.n_targs)

    return params
end
params = SampleConfig()

# Set number of knot points (changes discretization mesh resolution)
params.a.N = 100
N = params.a.N

# Generates a reference trajectory of specified size to be able to test different integration methods
j_targ = 1
ref_traj = generate_initial_guess_scp(params,j_targ)

# Process the trajectory into an integration batch
t_ref = ref_traj.t
x_ref = ref_traj.x
u_ref = ref_traj.u
if params.a.disc == 0
    u_ref = [u_ref u_ref[:,end]]
end

# Build C2D batch configuration (see c2d_nonlinear)
TS_batch = Vector([(t_ref[k],t_ref[k+1]) for k = 1:N-1])
X_batch = Vector([(x_ref[:,k],x_ref[:,k+1]) for k = 1:N-1])
U_batch = Vector([(u_ref[:,k],u_ref[:,k+1]) for k = 1:N-1])

# Dynamics
dynamics_ctcs = DynamicsLinearizedCTCS(params)
dyn_lin = (t,x,u,p) -> dynamics_ctcs(t,x,u,params,j_targ)
dyn_nl  = (t,x,u,p) -> dynamics_nonlinear_ctcs(t,x,u,params,j_targ)

# Other
disc = params.a.disc
p_batch = []

# Sizing variables
N = length(X_batch)
nx = length(X_batch[1][1])
nu = length(U_batch[1][1])
np = size(p_batch,1)
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
    f0 = vec(vcat(reshape(I(nx),:,1), zeros(nxnu), zeros(nxnp,1), zeros(nx))) # vec operation
elseif disc == 1
    f0 = vec(vcat(reshape(I(nx),:,1), zeros(nxnu), zeros(nxnu), zeros(nxnp,1), zeros(nx))) # vec operation
else
    error("Please select a valid discretization hold order.")
end

# Define the propagation function in terms of:
#   ODE state vector `z` = [x; ΦA; ΦB; ΦΣ; Φw] (vectorized)
#   Batch index `k` (for indexing time and control input parameters)
#   Current time `t`
function prop_fun!(dz,z,k,t)
   dz[:] = ode_nonlinear(t,z,optimal_controller(t,TS_batch[k],U_batch[k],disc),p_batch,nx,nu,np,nx2,nxnu,nxnp,dyn_nl,dyn_lin,disc;t_span=TS_batch[k])
end
# prop_fun(z,k,t) = ode_nonlinear(t,z,optimal_controller(t,TS_batch[k],U_batch[k],disc),p_batch,nx,nu,np,nx2,nxnu,nxnp,dyn_nl,dyn_lin,disc;t_span=TS_batch[k])

# @time begin
#     for k = 1:N
#         # Setup 
#         z0 = vec(vcat(X_batch[k][1], f0))
#         t_span = TS_batch[k]

#         # # Propagate (RK4) and record defect (δk)
#         # prop_fun_ = (t,z) -> prop_fun(z,k,t)
#         # Δt_rk4_step = max((1/num_disc_steps)*(t_span[2]-t_span[1]), h_min)
#         # ~,Z = rk4_batch(prop_fun_,z0,t_span[1],t_span[2],Δt_rk4_step)
#         # zf = Z[:,end]

#         # Propagate with DifferentialEquations.jl
#         dz = zeros(length(z0))
#         prob = DifferentialEquations.OrdinaryDiffEq.ODEProblem{true}(prop_fun!,z0,t_span,k) # parameter is the batch index
#         sol = DifferentialEquations.OrdinaryDiffEq.solve(prob, Tsit5())
#         zf = sol.u[end]

#         # Construct output matrices for this timestep (de-vec operation)
#         if disc == 0
#             Ak[k]  = devec(zf,nx,nx,nx2,nx)
#             Bmk[k] = devec(zf,nx,nu,nxnu,nx+nx2)
#             Σk[k]  = devec(zf,nx,np,nxnp,nx+nx2+nxnu)
#             wk[k]  = devec(zf,nx,1,nx,nx+nx2+nxnu+nxnp)
#         elseif disc == 1
#             Ak[k]  = devec(zf,nx,nx,nx2,nx)
#             Bmk[k] = devec(zf,nx,nu,nxnu,nx+nx2)
#             Bpk[k] = devec(zf,nx,nu,nxnu,nx+nx2+nxnu)
#             Σk[k]  = devec(zf,nx,np,nxnp,nx+nx2+2*nxnu)
#             wk[k]  = devec(zf,nx,1,nx,nx+nx2+2*nxnu+nxnp)
#         end

#         # Record defect
#         δk[k] = norm(zf[1:nx] - X_batch[k][2])
#     end
# end

u0 = DDTOSCP.StaticArrays.SVector{N}([vec(vcat(X_batch[k][1], f0)) for k=1:N])
t_span = DDTOSCP.StaticArrays.SVector{N}(TS_batch)
p = DDTOSCP.StaticArrays.SVector{N}(collect(1:N))
prob = DDTOSCP.DifferentialEquations.ODEProblem{true}(prop_fun!,u0[1],t_span[1],p[1])
prob_func = (prob,k,repeat) -> DDTOSCP.DifferentialEquations.remake(prob,u0=u0[k],tspan=t_span[k],p=p[k])
batchprob = DDTOSCP.DifferentialEquations.EnsembleProblem(prob,prob_func=prob_func,safetycopy=false)
# sol = DDTOSCP.DiffEqGPU.solve(batchprob, GPUTsit5(), EnsembleGPUKernel(CUDA.CUDABackend()), trajectories=N)
sol = DDTOSCP.DifferentialEquations.solve(batchprob, DDTOSCP.DifferentialEquations.Tsit5(), DDTOSCP.DifferentialEquations.EnsembleDistributed(), trajectories=N)