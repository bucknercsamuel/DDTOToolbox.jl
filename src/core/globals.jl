#=
Package-wide constants, type aliases, solver/verbosity settings, and utility
overrides used throughout DDTOToolbox.
=#

# Custom data types
"""Scalar floating-point type used throughout the toolbox."""
const CReal = Float64

"""Column vector of [`CReal`](@ref) values."""
const CVector = Vector{CReal}

"""Matrix of [`CReal`](@ref) values."""
const CMatrix = Matrix{CReal}

# Standard basis vectors
"""Unit vector along the inertial ``x``-axis."""
const e_x = CVector([1,0,0])

"""Unit vector along the inertial ``y``-axis."""
const e_y = CVector([0,1,0])

"""Unit vector along the inertial ``z``-axis."""
const e_z = CVector([0,0,1])

# Unit conversion
"""Conversion factor from radians to degrees."""
const RAD_2_DEG = 180 / π

"""Conversion factor from degrees to radians."""
const DEG_2_RAD = π / 180

"""Conversion factor from meters to kilometers."""
const M_2_KM = 1/1000

"""Conversion factor from kilometers to meters."""
const KM_2_M = 1000

"""Conversion factor from newtons to kilonewtons."""
const N_2_KN = 1/1000

"""Conversion factor from kilonewtons to newtons."""
const KN_2_N = 1000

# Colors for print statements
const BOLD = "\u001b[1m"
const GRAY = "\u001b[38;5;248m"
const CYAN = "\u001b[36m"
const RED = "\u001b[31m"
const GREEN = "\u001b[32m"
const YELLOW = "\u001b[33m"
const ORANGE = "\u001b[38;5;208m"
const RESET = "\u001b[0m"

# Hack: define empty type for SymPy's "Sym" if not already defined (SymPy not imported)
if !@isdefined(Sym)
    """Placeholder type used when SymPy is unavailable."""
    struct Sym end
end

# Set solver
# Current options: {"Clarabel", "ECOS", "MOSEK", "OSQP"}
"""Default JuMP solver name used by optimization routines."""
SOLVER = "Clarabel"
SOLVER_CTCS_DISABLED = SOLVER
SOLVER_CTCS_ENABLED = SOLVER

# Set verbose option for each algorithm
"""Print internal updates for optimal-solution bracket searches when `true`."""
VERB_OPT = true

"""Print internal updates for DDTO branch solves when `true`."""
VERB_DDTO = true

# Configure copy method to work with structures
"""
    Base.copy(t::T) where T

Deep-copy all fields of a mutable struct into a new instance of the same type.

# Arguments
- `t::T`: struct instance to copy

# Returns
- new instance of type `T` whose fields are deep copies of those in `t`
"""
Base.copy(t::T) where T = T([deepcopy(getfield(t,k)) for k ∈ fieldnames(T)]...)

# Set randomization seed
Random.seed!(12345)
