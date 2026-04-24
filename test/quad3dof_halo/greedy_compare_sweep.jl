using DDTOToolbox
using DataFrames
using LinearAlgebra
include("sim_landing.jl")
include("plots/plot_defaults.jl")
include("plots/paper_plot_greedy_compare.jl")

# Initialize the quadcopter vehicle
quad = Quad3DoFHaloParams()
quad.R_targs_min = 1.
quad.n_targs_min = 2
quad.n_targs_max = 3
quad.γ_gs = 89*3.14159/180
quad.v_min_V = -10.
quad.v_max_L = 10.
quad.a.ToF_min = 5.
quad.a.ToF_max = 40.
quad.a.Δt_min = .5*quad.a.ToF_min/(quad.a.N-1)
quad.a.Δt_max =  2*quad.a.ToF_max/(quad.a.N-1)
quad.a.Δt_cvx = (quad.a.Δt_min + quad.a.Δt_max)/2.
quad.w_obj_decay_factor = 1.6
quad.a.scp_iters = 50
quad.ϵ_subopt = 0.01
quad.a.use_single_cvx = true

# Parameters
r0 = [0,0,150]        # [m] Initial position (NED frame)
v0 = [0,0,0]          # [m/s] Initial velocity (NED frame)
n_obs = 8             # Number of obstacles
n_target_pool = 8     # Number of targets in the target pool
greedy_dts = [1.,Inf] # Greedy update timestep test options
R_ROI = r0[3]/3       # [m] Radius of the region of interest for targets
target_noise_std = .3 # Standard deviation of target noise
target_noise_crossweight = .05 # Cross-weighting factor for target noise

# Randomization Seed
# good options: {146542,23}
seed = 51
idxs = []
for k = 150:200
    seed = k

    # Simulate DDTO
    Random.seed!(seed)
    results_all = []
    results,_ = simulate_halo_landing(
        copy(quad),r0,v0,
        n_target_pool=n_target_pool,
        n_obs=n_obs,
        R_ROI=R_ROI,
        target_noise_std=target_noise_std,
        target_noise_crossweight=target_noise_crossweight,
        h_cut=75.,
        h_term=0.)
    append!(results_all, [results])

    # Simulate greedy variants
    for dt in greedy_dts
        Random.seed!(seed)
        quad.a.use_single_cvx = false
        results,_ = simulate_halo_landing(
            copy(quad),r0,v0,
            greedy=true,
            greedy_dt=dt,
            n_target_pool=n_target_pool,
            n_obs=n_obs,
            R_ROI=R_ROI,
            target_noise_std=target_noise_std,
            target_noise_crossweight=target_noise_crossweight)
        append!(results_all, [results])
    end

    # Display results summary
    isint(x) = x - floor(x) == 0
    condint(x) = isint(x) ? Int(x) : x
    strcvt(x) = isinf(x) ? "∞" : string(condint(x))
    delta(x,k) = x[k+1] - x[k]
    delta(x::Matrix,k) = x[:,k+1] - x[:,k]
    labels_greedy = ["Greedy-"*strcvt(dt) for dt in greedy_dts]
    labels = ["DDTO", labels_greedy...]
    final_time = [results["sim_time"][end] for results in results_all]
    cum_thrust = [sum([norm(results["sim_control"][1:3,k],2)*delta(results["sim_time"],k) for k = 1:size(results["sim_control"],2)-1]) for results in results_all]
    cum_jerk   = [sum([norm(delta(results["sim_control"][1:3,:],k),2) for k = 1:size(results["sim_control"],2)-1]) for results in results_all]
    final_radii = [max(results["targs_radii"][:,end]...) for results in results_all]
    avg_solve_time = [mean(results["guid_update_solve_time"]) for results in results_all]
    safe_site = [radii >= quad.R_targs_min for radii in final_radii]
    df = DataFrames.DataFrame(
        Type=labels,
        FinalTime=final_time,
        CumThrust=cum_thrust,
        CumJerk=cum_jerk,
        FinalSiteRadius=final_radii,
        FinalSiteSafe=safe_site,
        AvgSolveTime=avg_solve_time
    )
    display(df)

    if cum_thrust[1] < maximum(cum_thrust[2:end])
        push!(idxs, seed)
    end
end

# screens = []
# interactive = false
# target_colors = generate_custom_colors(n_target_pool)
# with_theme(theme2d; fontsize=fontsize) do
#     push!(screens, paper_plot_greedy_compare(
#         results_all, 
#         azel=(pi/4,pi/8), 
#         interactive=interactive))
# end
# hold_interactive(screens)
# ;