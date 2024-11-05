# Data types
const CReal = Float64
const CVector = Vector{CReal}
const CMatrix = Matrix{CReal}

# Standard basis vectors
const e_x = CVector([1,0,0])
const e_y = CVector([0,1,0])
const e_z = CVector([0,0,1])

# Unit conversion
const RAD_2_DEG = 180 / π
const DEG_2_RAD = π / 180
const M_2_KM = 1/1000
const KM_2_M = 1000
const N_2_KN = 1/1000
const KN_2_N = 1000

# Colors for print statements
const BOLD = "\u001b[1m"
const GRAY = "\u001b[38;5;248m"
const CYAN = "\u001b[36m"
const RED = "\u001b[31m"
const GREEN = "\u001b[32m"
const YELLOW = "\u001b[33m"
const RESET = "\u001b[0m"

# Set solver
# Current options: {"Clarabel", "ECOS", "MOSEK", "OSQP"}
SOLVER_CTCS_DISABLED = "Clarabel"
SOLVER_CTCS_ENABLED = "Clarabel"

# Set verbose option for each algorithm
VERB_OPT = true # Choose whether to print internal updates for the optimal solution bracket searches
VERB_DDTO = true # Choose whether to print internal updates for the DDTO solution branches

# Configure copy method to work with structures
Base.copy(t::T) where T = T([deepcopy(getfield(t,k)) for k ∈ fieldnames(T)]...)

# Set randomization seed
Random.seed!(12345)