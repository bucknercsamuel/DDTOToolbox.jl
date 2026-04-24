using DDTOToolbox
include("plots/plot_defaults.jl")
include("plots/plot_3d_trajs.jl")
include("plots/plot_states.jl")

# Initialize the quadcopter vehicle
quad = Quad3DoFHaloParams()

# Initial conditions
r0 = [0,0,150] # [m] Initial position (NED frame)
v0 = [0,0,0]   # [m/s] Initial velocity (NED frame)


# Dynamics
dynamics = (t,x,T,U,quad) -> dynamics_nonlinear_nondilated(t,x,optimal_controller(t,T,U,quad.a.disc),quad)

# Set randomization seed
Random.seed!(123)

# Simulate
# greedy = false
greedy = true
# dt = Inf
dt = 1
# dt = 0.1
results = simulate_halo_landing(quad,r0,v0,dynamics,greedy=greedy,greedy_dt=dt)

# Plot results
screens = []
with_theme(theme3d; fontsize=fontsize) do
    push!(screens, plot_3d_trajs(results))
    push!(screens, plot_states(results, integrated_sim=false))
end
hold_interactive(screens)
;