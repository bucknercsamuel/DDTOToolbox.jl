module DDTOSCP

using LinearAlgebra
using Random
using JuMP, MosekTools, ECOS
using Statistics
using Printf

export 
    Quad3DoFCageParams,
    solve_ddtoscp,
    solve_ddto,
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


include("core/globals.jl")
include("core/structs.jl")
include("core/sim.jl")
include("core/disc.jl")
include("core/opt.jl")
include("core/utils_ddto.jl")
include("core/utils_ddtoscp.jl")

include("quad3dofcage/params.jl")
# include("quad3dofcage/parsers.jl")
include("quad3dofcage/dynamics.jl")
include("quad3dofcage/initial_guess.jl")
include("quad3dofcage/prob_ddto.jl")
include("quad3dofcage/prob_ddtoscp.jl")

end # module
