module DDTOSCP

using LinearAlgebra
using Random, Noise
using JuMP, MosekTools, ECOS, OSQP
using Printf
using SymPy

export
    # Core
    solve,
    solve_cvx,
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
    Quad3DoFCageParams,
    Quad3DoFCageSampleScenario,
    Quad3DoFHaloParams,
    Quad3DoFParams,
    DIntegrator2DoFParams,
    scaling_matrices,
    generate_initial_guess_ddtoscp,
    custom_scaling!,
    dynamics_nonlinear,
    dynamics_nonlinear_nondilated,
    optimal_controller,
    rk4_step,
    # ADDTO
    extract_trunk_segment,
    extract_guid_lock_segment,
    remove_ddto_target!,
    switch_decision,
    sim_acquire_new_targets!,
    sim_update_locked_targets!,
    sim_generate_random_targets,
    rk4_step_pyjulia,
    reallocate_targ_dims!,
    sort_des_score!,
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
include("core/adapt_ddto/utils_addto.jl")
include("core/adapt_ddto/utils_percep.jl")
include("core/adapt_ddto/utils_pyjulia.jl")

# >> Quad 3-DOF Scenario Functionalities <<
include("quad3dof/params_cage.jl")
include("quad3dof/params_halo.jl")
include("quad3dof/params.jl")
include("quad3dof/prob.jl")
include("quad3dof/dynamics.jl")
include("quad3dof/dynamics_ctcs.jl")
include("quad3dof/initial_guess.jl")
include("quad3dof/skyenet_interface.jl")

# >> Double Integrator 2-DOF Scenario Functionalities <<
include("dint2dof/params.jl")
include("dint2dof/prob.jl")
include("dint2dof/dynamics.jl")
include("dint2dof/initial_guess.jl")

end # module
