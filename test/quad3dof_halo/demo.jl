using DDTOSCP
include("sim_landing.jl")

# Initialize the quadcopter vehicle
quad = Quad3DoFHaloParams()

# Initial conditions
r0 = [0,0,150] # [m] Initial position (NED frame)
v0 = [0,0,0]   # [m/s] Initial velocity (NED frame)

# Dynamics
dynamics = (t,x,sol) -> dynamics_nonlinear_nondilated(t,x,optimal_controller(t,sol.t,sol.u,quad.a.disc),quad)

# Set randomization seed
Random.seed!(0)

# Simulate
results = simulate_halo_landing(r0,v0)



# # ..:: Log resulting sim time/state/control into a `Solution` object ::..
# t_sim = results_sim_time
# r_sim = results_sim_state[1:3,:]
# v_sim = results_sim_state[4:6,:]
# T_sim = results_sim_control[1:3,:]
# Γ_sim = results_sim_control[4,:]
# cost_sim = sum(Γ_sim) * Δt_sim
# T_nrm_sim = CVector([norm(T_sim[:,k],2) for k=1:length(Γ_sim)])
# γ_sim     = CVector([acos(dot(T_sim[:,k],e_z)/norm(T_sim[:,k],2)) for k=1:length(Γ_sim)])
# results_sim_sol = Solution(t_sim, r_sim, v_sim, T_sim, Γ_sim, cost_sim, T_nrm_sim, γ_sim)

# # == Plotting ==
# include("plots_utils.jl")
# include("plots_core.jl")
# @pyimport matplotlib.animation as animation
# @pyimport mpl_toolkits as mpl
# set_fonts()
# set_fonts()
# pygui(false)
# replay_addto_landing(quad, 
#                      results_sim_sol, 
#                      results_guid_update_branches, 
#                      results_guid_update_trajs, 
#                      results_guid_update_time, 
#                      results_targs_status, 
#                      results_targs_radii, 
#                      view_az=-45)
# plot_addto_parametric_3D_trajectory(quad, 
#                                     results_sim_sol, 
#                                     results_guid_update_branches, 
#                                     results_guid_update_trajs, 
#                                     view_az=-45)
# plot_addto_signals_subplot(quad, 
#                            results_sim_sol, 
#                            results_targs_radii, 
#                            results_guid_update_time, 
#                            guid_lock_time)