module DDTOSCP

using LinearAlgebra
using Random
using JuMP, MosekTools, ECOS
using Statistics
using Printf

export 
    solve,
    solve_cvx,
    skyenet_ddtoscp_interface,
    generate_dynamics_partials,
    generate_dynamics_partials_ctcs,
    Quad3DoFCageParams,
    Quad3DoFCageSampleScenario,
    DIntegrator2DoFParams,
    scaling_matrices,
    generate_initial_guess_ddtoscp,
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
include("core/sim.jl")
include("core/disc.jl")
include("core/opt_base.jl")
include("core/opt_sing_cvx.jl")
include("core/opt_ddto_cvx.jl")
include("core/opt_sing_scp.jl")
include("core/opt_ddto_scp.jl")

# >> Quad 3-DOF Cage Scenario Functionalities <<
include("quad3dofcage/params.jl")
include("quad3dofcage/prob.jl")
include("quad3dofcage/dynamics.jl")
include("quad3dofcage/dynamics_ctcs.jl")
include("quad3dofcage/initial_guess.jl")
include("quad3dofcage/skyenet_interface.jl")

# >> Double Intregator 2-DOF Scenario Functionalities <<
include("dintegrator2dof/params.jl")
include("dintegrator2dof/prob.jl")
include("dintegrator2dof/dynamics.jl")
include("dintegrator2dof/initial_guess.jl")

end # module
