"""
    DDTOToolbox

Deferred Decision Trajectory Optimization toolbox for multi-target trajectory
planning. Provides DDTO-SCP (`solve`), DDTO-CVX (`solve_cvx`),
lexicographic DDTO (`solve_lex`), 3-DOF quadcopter and 2-DOF
double-integrator problem formulations, and Adaptive-DDTO closed-loop guidance
utilities.
"""
module DDTOToolbox

using LinearAlgebra
using JuMP, Clarabel
using Random, Noise, Statistics, Distributions
using OrdinaryDiffEq, StaticArrays
using ForwardDiff
using Printf
using JLD2

export
    # Core
    solve,
    solve_cvx,
    solve_lex,
    Solution,
    DDTOSolution,
    EmptySolution,
    EmptyDDTOSolution,
    Quad3DoFSolution,
    Quad3DoFDDTOSolution,
    EmptyQuad3DoFSolution,
    EmptyQuad3DoFDDTOSolution,
    skyenet_ddtoscp_interface,
    generate_dynamics_partials,
    generate_dynamics_partials_ctcs,
    AlgorithmParams,
    Quad3DoFCageParams,
    Quad3DoFCageSampleScenario,
    Quad3DoFHaloParams,
    Quad3DoFParams,
    DIntegrator2DoFParams,
    scaling_matrices,
    generate_initial_guess_scp,
    generate_initial_guess_ddtoscp,
    custom_scaling!,
    dynamics_linear,
    dynamics_nonlinear,
    dynamics_nonlinear_nondilated,
    DynamicsLinearizedCTCS,
    dynamics_ctcs,
    dynamics_nonlinear_ctcs,
    optimal_controller,
    optimal_controller_nondilated,
    rk4_step,
    time_dilation_control_to_wall_clock_time,
    wall_clock_time_to_time_dilation_control,
    quat_to_dcm,
    ode_nonlinear,
    # ADDTO
    compute_ddto_guidance!,
    check_unsafe_targets!,
    check_branch_switch!,
    check_cutoff_altitude!,
    activate_guidance_lock!,
    extract_trunk_segment,
    extract_guid_lock_segment,
    remove_ddto_target!,
    switch_decision,
    log_results!,
    sim_build_target_pool,
    sim_refresh_targets!,
    sim_update_targets!,
    generate_obstacles!,
    sim_generate_random_targets,
    rk4_step_pyjulia,
    reallocate_targ_dims!,
    sort_des_score!,
    setup_addto_dicts,
    save_results,
    set_ddto_subopt!,
    step_halo_sim,
    # Basic
    CReal,
    CVector,
    CMatrix,
    e_x,
    e_y,
    e_z,
    RAD_2_DEG,
    DEG_2_RAD,
    M_2_KM,
    KM_2_M,
    N_2_KN,
    KN_2_N

# >> Core Functionalities <<
include("core/globals.jl")
include("core/linalg.jl")
include("core/structs.jl")
include("core/utils.jl")
include("core/sim.jl")
include("core/dynamics.jl")
include("core/disc.jl")
include("core/opt_base.jl")
include("core/opt_sing_cvx.jl")
include("core/opt_ddto_cvx.jl")
include("core/opt_sing_scp.jl")
include("core/opt_ddto_scp.jl")
include("core/opt_ddto_lex.jl")

# >> Quad 3-DOF Scenario Functionalities <<
include("quad3dof/params_cage.jl")
include("quad3dof/params_halo.jl")
include("quad3dof/params.jl")
include("quad3dof/param_update_law.jl")
include("quad3dof/prob.jl")
include("quad3dof/dynamics.jl")
include("quad3dof/initial_guess.jl")
include("quad3dof/skyenet_interface.jl")

# >> Double Integrator 2-DOF Scenario Functionalities <<
include("dint2dof/params.jl")
include("dint2dof/prob.jl")
include("dint2dof/dynamics.jl")
include("dint2dof/initial_guess.jl")

# >> Adaptive-DDTO Functionalities <<
include("core/adapt_ddto/algorithm.jl")
include("core/adapt_ddto/sim.jl")
include("core/adapt_ddto/utils.jl")

end # module
