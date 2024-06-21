using DDTOSCP
include("sim_landing.jl")
include("plots.jl")

# Initialize the quadcopter vehicle
quad = Quad3DoFHaloParams()

# Initial conditions
r0 = [0,0,150] # [m] Initial position (NED frame)
v0 = [0,0,0]   # [m/s] Initial velocity (NED frame)

# Dynamics
dynamics = (t,x,T,U,quad) -> dynamics_nonlinear_nondilated(t,x,optimal_controller(t,T,U,quad.a.disc),quad)

# Set randomization seed
Random.seed!(0)

# Simulate
results = simulate_halo_landing(r0,v0,greedy=true,greedy_dt=5)

# Plot results
build_plots(results; interactive=true)