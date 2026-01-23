using DDTOToolbox
using DataFrames
include("sim_landing.jl")
include("plots.jl")

# Initialize the quadcopter vehicle
quad = Quad3DoFHaloParams()
quad.R_targs_min = 2.
quad.n_targs_min = 2
quad.n_targs_max = 3
quad.γ_gs = 89*3.14159/180
quad.v_min_V = -10.
quad.v_max_L = 10.
quad.a.ToF_min = 10.
quad.a.ToF_max = 60.
quad.a.Δt_min = .5*quad.a.ToF_min/(quad.a.N-1)
quad.a.Δt_max =  2*quad.a.ToF_max/(quad.a.N-1)
quad.a.Δt_cvx = (quad.a.Δt_min + quad.a.Δt_max)/2.
quad.w_obj_decay_factor = 1.6
quad.a.scp_iters = 50
quad.ϵ_subopt = 0.02

# Parameters
r0 = [0.,0.,150.] # [m] Initial position (NED frame)
v0 = [0.,0.,0.]   # [m/s] Initial velocity (NED frame)
seed = 123            # Random seed
n_obs = 8             # Number of obstacles
n_target_pool = 8     # Number of targets in the target pool
greedy_dts = [1.,Inf] # Greedy update timestep test options
R_ROI = r0[3]/3       # [m] Radius of the region of interest for targets
target_noise_std = 0.2 # Standard deviation of target noise
target_noise_crossweight = 0.05 # Cross-weighting factor for target noise

# Set randomization seed
Random.seed!(seed)

# Greedy vs DDTO config
greedy = false; dt = -1.;
# greedy = true; dt = 1.;
# greedy = true; dt = Inf;

# Simulate
results,error_code = simulate_halo_landing(
    quad,r0,v0,
    greedy=greedy,
    greedy_dt=dt,
    R_ROI=R_ROI, 
    n_target_pool=n_target_pool, 
    n_obs=n_obs,
    target_noise_std=target_noise_std,
    target_noise_crossweight=target_noise_crossweight)

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
interactive = true
target_colors = generate_custom_colors(n_target_pool)
with_theme(theme2d; fontsize=fontsize) do
    # push!(screens, paper_plot_trajallocation(quad, results, interactive=interactive))
    push!(screens, plot_3d_trajs(results, interactive=interactive))
    push!(screens, plot_2d_trajs_XY(results, interactive=interactive))
    # push!(screens, plot_states(results, integrated_sim=true, interactive=interactive))
end
hold_interactive(screens)
;