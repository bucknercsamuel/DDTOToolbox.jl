using DDTOSCP
using DataFrames
using LinearAlgebra
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
Random.seed!(12345)

# Simulate
greedy_dts = [0.1,1,5,Inf]
results_all = []
append!(results_all, [simulate_halo_landing(copy(quad),r0,v0,dynamics)])
# for dt in greedy_dts
#     append!(results_all, [simulate_halo_landing(copy(quad),r0,v0,dynamics,greedy=true,greedy_dt=dt)])
# end

# Display results
isint(x) = x - floor(x) == 0
condint(x) = isint(x) ? Int(x) : x
strcvt(x) = isinf(x) ? "∞" : string(condint(x))
labels_greedy = ["Greedy-"*strcvt(dt) for dt in greedy_dts]
labels = ["DDTO", labels_greedy...]
final_time = [results["sim_time"][end] for results in results_all]
cum_thrust = [sum([norm(results["sim_control"][1:3,k],2) for k in size(results["sim_control"],2)]) for results in results_all]
avg_thrust = [thrust/time for (thrust,time) in zip(cum_thrust,final_time)]
final_radii = [max(results["targs_radii"][:,end]...) for results in results_all]
safe_site = [radii >= quad.R_targs_min for radii in final_radii]
df = DataFrame(
    Type=labels,
    FinalTime=final_time,
    CumThrust=cum_thrust,
    AvgThrust=avg_thrust,
    FinalSiteRadius=final_radii,
    FinalSiteSafe=safe_site
)
display(df)

# Plot results
with_theme(theme3d; fontsize=fontsize) do
    screens = [
        plot_greedy_compare(results_all)
    ]
    hold_interactive(screens)
end