# ..:: Solution Structures ::..
"""
`Solution` stores the optimal solution.
`DDTOSolution` stores the optimal solution from DDTO-SCP.
"""

mutable struct Solution
    t::CVector        # [s] Time vector
    x::CMatrix        # State trajectory
    u::CMatrix        # Control signal
    cost::CReal       # Optimal cost
end

mutable struct DDTOSolution
    targs::Vector{Solution} # Contains the `Solution` to each target
end


# ..:: Constructors for empty `*Solution` structs ::..

function EmptySolution()::Solution

    t = CVector(undef,0)
    x = CMatrix(undef,0,0)
    u = CMatrix(undef,0,0)
    cost = Inf

    return Solution(t,x,u,cost)
end

function EmptyDDTOSolution(n_targs)::DDTOSolution

    targs = Vector{Solution}(undef, n_targs)
    for j = 1:n_targs
        targs[j] = EmptySolution()
    end

    return DDTOSolution(targs)
end

# ..:: Template functions ::..
# These functions are problem-specific and defined for a specific `params` object in other folders besides `core`,
# with the `prob.jl` and `dynamics.jl` files.

function core_problem(mdl::JuMP.Model, x::JuMP.VariableRef, u::JuMP.VariableRef, params::Nothing, ref_traj::Solution)
    return 0
end

function dynamics_linearized(t_ref::CReal, x_ref::CVector, ν_ref::CVector, params::Nothing)
    return 0
end

function dynamics_nonlinear(t::CReal, x::CVector, ν::CVector, params::Nothing)
    return 0
end

function dynamics_linear(params::Nothing)
    return 0
end