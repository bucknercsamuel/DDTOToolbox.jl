using DDTOSCP
include("sim_landing.jl")
include("plots.jl")

# # Initialize the quadcopter vehicle
# quad = Quad3DoFHaloParams()

# # Initial conditions
# r0 = [0,0,150] # [m] Initial position (NED frame)
# v0 = [0,0,0]   # [m/s] Initial velocity (NED frame)

# # Dynamics
# dynamics = (t,x,T,U,quad) -> dynamics_nonlinear_nondilated(t,x,optimal_controller(t,T,U,quad.a.disc),quad)

# # Set randomization seed
# Random.seed!(0)

# # Simulate
# greedy_dts = [1,3,5,10,20,Inf]
# results_all = []
# append!(results_all, [simulate_halo_landing(copy(quad),r0,v0,dynamics)])
# for dt in greedy_dts
#     append!(results_all, [simulate_halo_landing(copy(quad),r0,v0,dynamics,greedy=true,greedy_dt=dt)])
# end

# Plot results
screens = ()
with_theme(theme3d; fontsize=fontsize) do
    screens = (
        plot_greedy_compare(results_all)
    )
end
hold_interactive(screens)