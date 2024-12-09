using DDTOSCP
using DataFrames
include("sim_landing.jl")
include("plots.jl")

# Initialize the quadcopter vehicle
quad = Quad3DoFHaloParams()

# Initial conditions
r0 = [0.,0.,150.] # [m] Initial position (NED frame)
v0 = [0.,0.,0.]   # [m/s] Initial velocity (NED frame)

# Dynamics
dynamics = (t,x,T,U,quad) -> dynamics_nonlinear_nondilated(t,x,optimal_controller(t,T,U,quad.a.disc),quad)

# Set randomization seed
Random.seed!(1243)
# Random.seed!(1234)

# Greedy vs DDTO config
<<<<<<< HEAD
greedy = false; dt = -1.;
# greedy = true; dt = 1.;
=======
# greedy = false; dt = -1.;
greedy = true; dt = 1.;
>>>>>>> a94024612f595e2e498cb1d9dc6cf7a44bd27ec5
# greedy = true; dt = Inf;

# Simulate
results = simulate_halo_landing(quad,r0,v0,dynamics,greedy=greedy,greedy_dt=dt,R_ROI=r0[3]/3, n_target_pool=1000, n_obs=10)

# Display results
delta(x::Vector,k) = x[k+1] - x[k]
delta(x::Matrix,k) = x[:,k+1] - x[:,k]
final_time = results["sim_time"][end]
cum_thrust = sum([norm(results["sim_control"][1:3,k],2)*delta(results["sim_time"],k) for k = 1:size(results["sim_control"],2)-1])
cum_jerk   = sum([norm(delta(results["sim_control"][1:3,:],k),2) for k = 1:size(results["sim_control"],2)-1])
final_radii = max(results["targs_radii"][:,end]...)
safe_site = final_radii >= quad.R_targs_min
df = DataFrame(
    FinalTime=final_time,
    CumThrust=cum_thrust,
    CumJerk=cum_jerk,
    FinalSiteRadius=final_radii,
    FinalSiteSafe=safe_site
)
println("")
display(df)

# Plot results
screens = []
with_theme(theme2d; fontsize=fontsize) do
    # push!(screens, paper_plot_trajallocation(quad, results, interactive=false))
    push!(screens, plot_3d_trajs(results, interactive=false))
    push!(screens, plot_2d_trajs_XY(results, interactive=false))
    # push!(screens, plot_states(results, integrated_sim=false))
end
hold_interactive(screens)
;