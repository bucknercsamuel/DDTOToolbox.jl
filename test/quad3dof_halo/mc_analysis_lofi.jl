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

# Monte-Carlo parameters
n_mc = 50
mc_seeds = 1:n_mc
greedy_dts = [-1,1,Inf] # -1 is DDTO

# Construct labels
isint(x) = x - floor(x) == 0
condint(x) = isint(x) ? Int(x) : x
strcvt(x) = isinf(x) ? "∞" : string(condint(x))
labels_greedy = ["Gr-"*strcvt(dt) for dt in greedy_dts[findall(dt->dt!=-1,greedy_dts)]]
labels = ["DDTO", labels_greedy...]

# Construct data storage container
data_segment = Dict(
    "results" => [],
    "error_code" => [], # 1: successful run, 2: failed run
    "final_time" => [],
    "cum_thrust" => [],
    "avg_thrust" => [],
    "final_radius_truth" => [],
    "safe_site" => []
)
data = Dict()
label_pointer = Dict()
for (k,label) in enumerate(labels)
    data[label] = deepcopy(data_segment)
    label_pointer[greedy_dts[k]] = label
end

# MC Loop
is_greedy = dt -> dt == -1 ? false : true
delta(x,k) = x[k+1] - x[k]
for seed in mc_seeds
    for dt in greedy_dts
        Random.seed!(seed) # same seed for each dt for fair comparison
        label = label_pointer[dt]
        try
            results = simulate_halo_landing(copy(quad),r0,v0,dynamics,greedy=is_greedy(dt),greedy_dt=dt)
            append!(data[label]["results"], results)
            append!(data[label]["error_code"], 1) # successful run
            append!(data[label]["final_time"], results["sim_time"][end])
            append!(data[label]["cum_thrust"], sum([norm(results["sim_control"][1:3,k],2)*delta(results["sim_time"],k) for k = 1:size(results["sim_control"],2)-1]))
            append!(data[label]["avg_thrust"], data[label]["cum_thrust"][end] / data[label]["final_time"][end])
            append!(data[label]["final_radius_truth"], max(results["targs_radii"][:,end]...))
            append!(data[label]["safe_site"], data[label]["final_radius_truth"][end] >= quad.R_targs_min)
        catch e
            append!(data[label]["error_code"], 2) # failed run (set to any value != 1)
            append!(data[label]["results"], [nothing])
            append!(data[label]["final_time"], [nothing])
            append!(data[label]["cum_thrust"], [nothing])
            append!(data[label]["avg_thrust"], [nothing])
            append!(data[label]["final_radius_truth"], [nothing])
            append!(data[label]["safe_site"], [nothing])

        end
    end
end

# Display results
removenothing(x) = x[findall(x->x!=nothing,x)]
mean = x -> sum(removenothing(x)) / length(removenothing(x))
final_time = [mean(data[label]["final_time"]) for label in labels]
cum_thrust = [mean(data[label]["cum_thrust"]) for label in labels]
avg_thrust = [mean(data[label]["avg_thrust"]) for label in labels]
final_radii = [mean(data[label]["final_radius_truth"]) for label in labels]
safe_site = [mean(data[label]["safe_site"]) for label in labels]
df = DataFrame(
    Type=labels,
    FinalTime=final_time,
    CumThrust=cum_thrust,
    AvgThrust=avg_thrust,
    FinalSiteRadius=final_radii,
    FinalSiteSafe=safe_site
)

# Plot results
with_theme(theme2d; fontsize=fontsize) do
    screens = [
        plot_mc_statistics(data, groupings=[(1,),tuple(collect(2:length(data))...)])
    ]
    hold_interactive(screens)
end