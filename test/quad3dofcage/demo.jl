using DDTOSCP
using LinearAlgebra

include("scenarios.jl")
include("plots.jl")

params = scenario_toy1()

scp_solutions, scp_simulations, 
ddtoscp_solutions, ddtoscp_simulations, 
defer_solutions, defer_simulations = solve_ddtoscp(params)

set_fonts()
PyPlot.close("all")
pygui(true)

# ..:: SCP Solutions ::..
plot_parametric_trajectories(
    params,
    scp_solutions,
    scp_simulations;
    display_obstacles=true, 
    fname="decoupled_scp_solutions")
    
if params.free_final_time
    plot_time_dilation(
        params, 
        scp_solutions, 
        scp_simulations;
        fname="plot_time_dilation")
end

plot_thrust_magnitude(
    params, 
    scp_solutions, 
    scp_simulations;
    fname="plot_thrust_magnitude")
    
plot_3vec(
    params, 
    scp_solutions, 
    scp_simulations,
    "r";
    fname="plot_positions")

plot_3vec(
    params, 
    scp_solutions, 
    scp_simulations,
    "v";
    fname="temp")

# ..:: DDTO-SCP Solutions ::..
plot_parametric_trajectories(
    params, 
    ddtoscp_solutions, 
    ddtoscp_simulations;
    defer_solution=defer_solutions,
    defer_simulation=defer_simulations,
    display_obstacles=true,
    fname="ddtoscp_solutions")

if params.free_final_time
    plot_time_dilation(
        params, 
        ddtoscp_solutions, 
        ddtoscp_simulations;
        defer_solution=defer_solutions,
        defer_simulation=defer_simulations,
        fname="plot_time_dilation")
end

plot_thrust_magnitude(
    params, 
    ddtoscp_solutions, 
    ddtoscp_simulations;
    defer_solution=defer_solutions,
    defer_simulation=defer_simulations,
    fname="plot_thrust_magnitude")

plot_3vec(
    params, 
    ddtoscp_solutions, 
    ddtoscp_simulations,
    "r";
    defer_solution=defer_solutions,
    defer_simulation=defer_simulations,
    fname="plot_positions")

plot_3vec(
    params, 
    ddtoscp_solutions, 
    ddtoscp_simulations,
    "v";
    defer_solution=defer_solutions,
    defer_simulation=defer_simulations,
    fname="plot_positions")

;