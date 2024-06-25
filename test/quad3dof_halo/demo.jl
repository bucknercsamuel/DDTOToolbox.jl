using DDTOSCP
include("sim_landing.jl")
include("plots.jl")

# Initialize the quadcopter vehicle
quad = Quad3DoFHaloParams()

# Initial conditions
r0 = [0,0,300] # [m] Initial position (NED frame)
v0 = [0,0,0]   # [m/s] Initial velocity (NED frame)

# Dynamics
dynamics = (t,x,T,U,quad) -> dynamics_nonlinear_nondilated(t,x,optimal_controller(t,T,U,quad.a.disc),quad)

# Set randomization seed
Random.seed!(123)

# Simulate
# results = simulate_halo_landing(quad,r0,v0,dynamics,greedy=true,greedy_dt=1)
results = simulate_halo_landing(quad,r0,v0,dynamics)

# Plot results
with_theme(theme3d; fontsize=fontsize) do
    screens = (
        plot_3d_trajs(results)
    )
end
hold_interactive(screens)