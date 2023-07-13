# ..:: Solution Structures ::..
"""
`Solution` stores the optimal solution.
`DDTOSolution` stores the optimal solution from DDTO.
`BranchSolution` stores the trajectory solution of a branch from DDTO
"""

mutable struct Solution
    t::CVector        # [s] Time vector
    x::CMatrix        # State trajectory
    u::CMatrix        # Control signal
    cost::CReal       # Optimal cost
end

mutable struct DDTOSolution
    targ_sols::Vector{Solution} # Contains the `Solution` to each target
    costs_sol::CVector          # Costs for each target
    cost_dd::CReal              # Cost for deferred decision
    idx_dd::Int                 # Deferred decision branch point index
end

mutable struct BranchSolution
    sol::Solution  # Contains the `Solution` for the branch
    cost_dd::CReal # Cost for deferred decision
    idx_dd::Int    # Deferred decision branch point index
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

    targ_sols = Vector{Solution}(undef, n_targs)
    costs_sol = CVector(undef, n_targs)
    cost_dd   = 0
    idx_dd    = 0

    for j = 1:n_targs
        targ_sols[j] = EmptySolution()
    end

    return DDTOSolution(targ_sols,costs_sol,cost_dd,idx_dd)
end

function EmptyBranchSolution()::BranchSolution
    sol = EmptySolution()
    cost_dd = 0
    idx_dd = 1
    return BranchSolution(sol, cost_dd, idx_dd)
end










# """
# `SolutionSCP` stores the SCP-optimal solution.
# """
# mutable struct SolutionSCP
#     # >> Raw data <<
#     t::CVector        # [s] Time vector
#     x::CMatrix        # State trajectory
#     u::CMatrix        # Control signal
#     ν::CMatrix        # [m] Virtual buffer variables for linearized constraints
#     μ::CMatrix        # [m] Virtual buffer slack for linearized constraints
#     cost::CReal       # Optimization's optimal cost
# end

# """
# `DDTOSolutionSCP` stores an SCP subproblem solution from DDTO.
# """
# mutable struct DDTOSolutionSCP
#     targ_sols::Vector{SolutionSCP} # Contains the `SolutionSCP` to each target
#     costs_sol::CVector             # Costs for each target
#     cost_dd::CReal                 # Cost for deferred decision
#     idx_dd::Int                    # Deferred decision branch point index
#     η::CVector                     # Trust region variables
# end

# """
#     EmptySolutionSCP()

#     Constructor for an empty SolutionSCP.

# # Arguments
# - :out sol: empty solution.
# """
# function EmptySolutionSCP()::SolutionSCP

#     t = CVector(undef,0)
#     x = CMatrix(undef,0,0)
#     u = CMatrix(undef,0,0)
#     ν = CMatrix(undef,0,0)
#     μ = CMatrix(undef,0,0)
#     cost = Inf

#     return SolutionSCP(t,x,u,ν,μ,cost)
# end

# """
#     EmptyDDTOSolutionSCP()

# Constructor for an empty DDTOSolutionSCP.

# # Arguments
# - :in n_targs: Number of targets for the solution
# - :out sol: empty solution.
# """
# function EmptyDDTOSolutionSCP(n_targs)::DDTOSolutionSCP

#     targ_sols = Vector{SolutionSCP}(undef, n_targs)
#     costs_sol = CVector(undef, 0)
#     cost_dd   = 0
#     idx_dd    = 0
#     η         = CVector(undef,0)

#     # Initialize each `SolutionSCP` with `EmptySolutionSCP()`
#     for j = 1:n_targs
#         targ_sols[j] = EmptySolutionSCP()
#     end
    
#     return DDTOSolutionSCP(targ_sols,costs_sol,cost_dd,idx_dd,η)
# end