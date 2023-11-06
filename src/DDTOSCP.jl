module DDTOSCP

using LinearAlgebra
using Random
using JuMP, MosekTools, ECOS
using Statistics
using Printf
using SymPy
using PackageCompiler

export 
    solve,
    skyenet_ddtoscp_interface,
    generate_dynamics_partials,
    Quad3DoFCageParams,
    Quad3DoFCageSampleScenario,
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
include("core/opt.jl")

# >> Quad 3-DOF Cage Scenario Functionalities <<
include("quad3dofcage/params.jl")
include("quad3dofcage/dynamics.jl")
include("quad3dofcage/initial_guess.jl")
include("quad3dofcage/prob.jl")
include("quad3dofcage/skyenet_interface.jl")

end # module
