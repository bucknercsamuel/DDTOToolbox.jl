# ..:: Lexicographic DDTO-SCP Solver Functions ::..
# Based on Elango et al. 2025 "Deferred Decision Trajectory Optimization", algorithm 2

function solve_lex(params; single_iter::Bool=false, ref_trajs::Any=nothing, simulate_solutions::Bool=true, process_the_solutions::Bool=true)
    # ..:: Problem setup ::..
    custom_scaling!(params)
    param_augmentation!(params)

    # ..:: Solve for DDTO branching solutions to ALL targets ::..
    ref_trajs_ddtoscp, scp_solutions, scp_costs, scp_converged, elapsed_solver_time = warmstart_ddtoscp(params, ref_trajs; single_iter=single_iter)
    if params.a.n_targs > 1
        time_ddto = @elapsed begin
            ddtoscp_solutions, ddtoscp_converged = solve_tree_ddtolex(params, scp_costs, ref_trajs_ddtoscp; single_iter=single_iter)
            println("\n Solve time for generating DDTO branch solutions to all targets:")
        end
        elapsed_solver_time += time_ddto
    else
        ddtoscp_solutions = copy(scp_solutions)
        ddtoscp_converged = copy(scp_converged)
    end
    println("Total DDTO solve time: ", elapsed_solver_time)

    # ..:: Simulate each target solution from I.C. to T.C. ::..
    if simulate_solutions
        @time begin
            if params.a.ctcs_enabled
                dynamics = (t,x,sol) -> dynamics_nonlinear_ctcs(t,x,optimal_controller(t,sol.t,sol.u,params.a.disc),params,0)
            else
                dynamics = (t,x,sol) -> dynamics_nonlinear(t,x,optimal_controller(t,sol.t,sol.u,params.a.disc),params)
            end
            scp_simulations = simulate(scp_solutions, dynamics, params.a.disc; max_steps=params.a.N_sim)
            ddtoscp_simulations = simulate(ddtoscp_solutions, dynamics, params.a.disc; max_steps=params.a.N_sim)
            println("\n Solve time for forward simulation:")
        end
    end

    # ..:: Post-process each target solution (problem-specific) ::..
    if process_the_solutions
        @time begin
            scp_solutions = process_solutions(scp_solutions, params)
            ddtoscp_solutions = process_solutions(ddtoscp_solutions, params)
            if simulate_solutions
                scp_simulations = process_solutions(scp_simulations, params)
                ddtoscp_simulations = process_solutions(ddtoscp_simulations, params)
            end
            println("\n Solve time for post-processing:")
        end
    end

    param_deaugmentation!(params)
    converged = scp_converged && ddtoscp_converged ? true : false
    if simulate_solutions
        return (
            scp_solutions, 
            scp_simulations, 
            ddtoscp_solutions, 
            ddtoscp_simulations,
            converged,
            elapsed_solver_time)
    else
        return (
            scp_solutions, 
            ddtoscp_solutions,
            converged,
            elapsed_solver_time)
    end
end

function solve_concatenated_ddtolex_subproblem(params, params_con, ref_traj::Solution, scp_iter::Int, ref_costs::CVector, cost_dd::CReal)
    # Given a set of reference trajectories with trunk indexed first, form the concatenated subproblem for the DDTO-LEX problem

    # Force CTCS for this method
    if !params.a.ctcs_enabled
        error("CTCS must be enabled for this method")
    end

    # Obtain parameters associated with state-space
    nx = params.a.nx
    nu = params.a.nu
    N = params.a.N

    # Create indexing functions enumerated over n+1 segments
    # 0 -> 1:nx, 1 -> nx+1:2nx, ..., n -> (n-1)nx+1:(n-1)nx+nx, n+1 -> (n-1)nx+nx+1:N
    J_nx(j) = j*nx+1:(j+1)*nx
    J_nu(j) = j*nu+1:(j+1)*nu

    # # Create virtual params object for the concatenated problem
    nx_con = nx*(params.a.n_targs+1)
    nu_con = nu*(params.a.n_targs+1)

    # Create new boundary conditions for the concatenated problem (Inf = no constraint)
    z0_con = Inf*ones(nx_con)
    u0_con = Inf*ones(nu_con)
    zf_con = Inf*ones(nx_con)
    uf_con = Inf*ones(nu_con)
    z0_con[J_nx(0)] = params.a.z0 # initial condition for trunk
    u0_con[J_nu(0)] = params.a.u0 # initial condition for trunk
    # Terminal conditions for branches
    for j = 1:params.a.n_targs
        zf_con[J_nx(j)] = params.a.zf_targs[:,j]
        uf_con[J_nu(j)] = params.a.uf_targs[:,j]
    end

    # Define concatenated dynamics functions
    function dyn_nl(t,x,u,p)
        z = vcat([dynamics_nonlinear_ctcs(t,x[J_nx(j)],u[J_nu(j)],params,j) for j = 0:params.a.n_targs]...)
        return z
    end
    function dyn_lin(t,x,u,p)
        dynamics_linearized_ctcs = DynamicsLinearizedCTCS(params)
        A = zeros(nx_con, nx_con)
        B = zeros(nx_con, nu_con)
        Σ = [] # not in use currently
        z = zeros(nx_con)
        for j = 0:params.a.n_targs
            A_,B_,Σ_,z_ = dynamics_linearized_ctcs(t,x[J_nx(j)],u[J_nu(j)],params,j)
            A[J_nx(j),J_nx(j)] = A_
            B[J_nx(j),J_nu(j)] = B_
            z[J_nx(j)] = z_
        end
        return A,B,Σ,z
    end

    # Define cost function as sum of time dilation terms for trunk (max deferral)
    function prob_cost_(mdl,x,u)
        s_trunk_idx = J_nu(0)[end]
        s_trunk = u[s_trunk_idx,:]
        J_term = 0
        J_running = sum(s_trunk)
        return J_running, J_term
    end

    # Define constraints function
    # since we are using CT-SCvx, don't impose any regular path constraints, only necessary extra constraints
    function prob_constraints_(mdl,x,u,ref_traj)
        # Must impose equality in state and control between end of trunk and start of branches
        for j = 1:params.a.n_targs
            @constraint(mdl, x[J_nx(0),end] == x[J_nx(j),1])
            @constraint(mdl, u[J_nu(0),end] == u[J_nu(j),1])
        end

        # Get trunk cost + branch cost and make suboptimality constraint
        J_running_trunk,_ = prob_cost(mdl,x[J_nx(0),:],u[J_nu(0),:],params)
        for j = 1:params.a.n_targs
            J_running_branch,J_term_branch = prob_cost(mdl,x[J_nx(j),:],u[J_nu(j),:],params)
            @constraint(mdl, (cost_dd + sum(J_running_trunk) + sum(J_running_branch) + J_term_branch) <= ((1 + params.a.ϵ_targs[j]) * ref_costs[j]))
        end

        ν_buff = []
        return ν_buff
    end

    # Solve the concatenated problem
    return solve_ctscvx_subproblem(
        params_con, 
        ref_traj, 
        z0_con, 
        zf_con, 
        u0_con, 
        uf_con, 
        dyn_nl,
        dyn_lin,
        prob_cost_,
        prob_constraints_,
        scp_iter;
        CTCS_idxs = [J_nx(j)[end] for j = 0:params.a.n_targs],
        dilation_idxs = [J_nu(j)[end] for j = 0:params.a.n_targs],
        TOF_idxs = [] # do not impose TOF constraints for this problem
    )
end

function unconcatenate_ddtolex_solution(ddtolex_sol::Solution, params)
    # Given a concatenated DDTO solution, unconcatenate it into a vector of solutions
    # for each target
    # Obtain parameters associated with state-space
    nx = params.a.nx
    nu = params.a.nu
    N = params.a.N

    # Create indexing functions enumerated over n+1 segments
    # 0 -> 1:nx, 1 -> nx+1:2nx, ..., n -> (n-1)nx+1:(n-1)nx+nx, n+1 -> (n-1)nx+nx+1:N
    J_nx(j) = j*nx+1:(j+1)*nx
    J_nu(j) = j*nu+1:(j+1)*nu

    # Unconcatenate the solution into n+1 segments (trunk + branches)
    ddto_sol_segmented = EmptyDDTOSolution(params.a.n_targs+1)
    for j = 1:params.a.n_targs+1
        ddto_sol_segmented.targs[j].t = ddtolex_sol.t # same time vector for all segments
        ddto_sol_segmented.targs[j].x = ddtolex_sol.x[J_nx(j-1),:]
        ddto_sol_segmented.targs[j].u = ddtolex_sol.u[J_nu(j-1),:]
    end

    # Obtain the cost of each segment independently
    for j = 1:params.a.n_targs+1
        J_running, J_term = prob_cost_eval(ddto_sol_segmented.targs[j].x, ddto_sol_segmented.targs[j].u, params)
        ddto_sol_segmented.targs[j].cost = J_running + J_term
    end

    # Obtain DDTO-format solutions by concatenating the trunk segment with the 2,...,N knots of the branches
    ddto_sol = EmptyDDTOSolution(params.a.n_targs)
    for j = 1:params.a.n_targs
        ddto_sol.targs[j].t = vcat(ddto_sol_segmented.targs[1].t, ddto_sol_segmented.targs[1].t[end] .+ ddto_sol_segmented.targs[j+1].t[2:end])
        ddto_sol.targs[j].x = hcat(ddto_sol_segmented.targs[1].x, ddto_sol_segmented.targs[j+1].x[:,2:end])
        ddto_sol.targs[j].u = hcat(ddto_sol_segmented.targs[1].u, ddto_sol_segmented.targs[j+1].u[:,2:end])
        ddto_sol.targs[j].cost = ddto_sol_segmented.targs[j+1].cost + ddto_sol_segmented.targs[1].cost
    end

    return ddto_sol_segmented, ddto_sol
end

function solve_tree_ddtolex(params, scp_costs, ref_trajs::DDTOSolution; single_iter=false)::Tuple{DDTOSolution,Bool}

    # Initialization
    cost_dd = 0. # running deferred-decision (DD) trajectory segment cost sum
    ref_costs = scp_costs
    params.a.τ_targs = zeros(params.a.n_targs)

    # Define a vertically concatenated reference trajectory in concatenated form
    ref_traj = ref_trajs.targs[1] # trunk initial guess
    for j = 1:length(ref_trajs.targs)
        ref_traj.x = vcat(ref_traj.x, ref_trajs.targs[j].x)
        ref_traj.u = vcat(ref_traj.u, ref_trajs.targs[j].u)
    end

    # Define concatenated params object (virtual, not kept)
    n_segments = params.a.n_targs+1
    params_con = deepcopy(params)
    params_con.a.nx = params.a.nx*n_segments
    params_con.a.nu = params.a.nu*n_segments
    
    # Construct Sx,sx,Su,su for the concatenated problem by just repeating them over n_segments
    # Noting Sx and Su are square matrices, while sx and su are vectors
    params_con.a.Sx = kron(params.a.Sx, I(n_segments))
    params_con.a.sx = kron(params.a.sx, ones(n_segments))
    params_con.a.Su = kron(params.a.Su, I(n_segments))
    params_con.a.su = kron(params.a.su, ones(n_segments))

    # Perform branching in the order of preference
    # for k = 1:(params.a.n_targs-1)
    λ_targ = params.a.λ_targs[1]
    VERB_DDTO && @printf("\n========= Solving DDTO-LEX Stage Problem for Deferred Target #%i =========\n", λ_targ)
    subproblem_ = (ref_traj, k) -> solve_concatenated_ddtolex_subproblem(params, params_con, ref_traj, k, ref_costs, cost_dd)
    (solution, feas_status, scp_converged) = solve_ctscvx_iteration(params_con, ref_traj, subproblem_; single_iter=single_iter)
    ddto_sol_segmented, ddto_sol = unconcatenate_ddtolex_solution(solution, params)

    # Perform parameter updates
    for j in params.a.J_targs
        params.a.τ_targs[j] += params.a.N
    end
    cost_dd += ddto_sol_segmented.targs[1].cost # update with cost of trunk segment
    # end

    return ddto_sol, scp_converged
end






# function solve_tree_ddto_lex(params, ref_costs::CVector; single_iter=false, ref_trajs=nothing)::Tuple{DDTOSolution,Bool}
#     # Top-level lexicographic DDTO-SCP solver for all branch points
#     # Solves the problem lexicographically where targets are dropped one by one
#     # and maximal time to maneuver is solved for each new target set with the trunk trajectory building up
#     #
#     # :in params: The params object
#     # :in ref_costs: Optimal costs from initial condition
#     # :out ddto_sol: Vectorized container for all DDTO branch solutions

#     # Define container for each DDTO branch solution
#     ddto_sol = EmptyDDTOSolution(params.a.n_targs)

#     # Define running deferred-decision (DD) trajectory segment cost sum
#     cost_dd = 0.

#     # Initialization
#     n_targs_total = copy(params.a.n_targs)
#     params.a.τ_targs = zeros(n_targs_total) # initialization
#     ref_initial_control = zeros(params.a.nu)
#     ddto_branch_sol = isnothing(ref_trajs) ? generate_initial_guess_scp(params) : ref_trajs
#     params_ = copy(params) # Temp object to be mutated through DDTO loop
#     find_J_elem(J_targs,j) = findfirst(τ->τ==j, J_targs)
#     J_targs_old = copy(params.a.J_targs)
#     idx_dd = 1
#     τ_opt = 0

#     # Perform branching in the order of preference
#     for k = 1:(n_targs_total-1)
#         λ_targ = params_.a.λ_targs[1]
#         VERB_DDTO && @printf("\n========= Solving DDTO-SCP Stage Problem for Deferred Target #%i =========\n", λ_targ)

#         # Directly solve for maximum deferral time using SCP with concatenated trajectory
#         prev_sol = copy(ddto_branch_sol)
#         prev_τ = copy(τ_opt)
#         ddto_branch_sol,τ_opt,Δcost_dd,scp_converged = solve_max_deferral_ddto_lex(params_, ref_costs[params_.a.J_targs], cost_dd, ref_initial_control; single_iter=single_iter, ref_trajs=ddto_branch_sol)
#         if τ_opt == 0
#             ddto_branch_sol = EmptyDDTOSolution(params_.a.n_targs)
#             for j ∈ params_.a.J_targs
#                 ddto_branch_sol.targs[find_J_elem(params_.a.J_targs,j)].x = prev_sol.targs[find_J_elem(J_targs_old,j)].x[:,prev_τ+1:end]
#                 ddto_branch_sol.targs[find_J_elem(params_.a.J_targs,j)].u = prev_sol.targs[find_J_elem(J_targs_old,j)].u[:,prev_τ+1:end]
#                 ddto_branch_sol.targs[find_J_elem(params_.a.J_targs,j)].cost = prev_sol.targs[find_J_elem(J_targs_old,j)].cost
#             end
#         end
#         J_targs_old = copy(params_.a.J_targs)

#         count = 1
#         for j ∈ params_.a.J_targs
#             if k == 1
#                 ddto_sol.targs[j].x = ddto_branch_sol.targs[j].x
#                 ddto_sol.targs[j].u = ddto_branch_sol.targs[j].u
#             else
#                 ddto_sol.targs[j].x[:,idx_dd:end] = ddto_branch_sol.targs[count].x
#                 ddto_sol.targs[j].u[:,idx_dd:end] = ddto_branch_sol.targs[count].u
#             end
#             ddto_sol.targs[j].cost = ddto_branch_sol.targs[count].cost
#             count += 1
#         end

#         # Determine target to be removed (first in the current list of λ_targs)
#         deleteat!(params_.a.λ_targs, 1)
#         pop_idx = findfirst(i->i==λ_targ, params_.a.J_targs)

#         # Have to do some slicing magic for matrices
#         matrix_slice = collect(1:params_.a.n_targs)
#         deleteat!(matrix_slice, pop_idx)

#         # Update params_ target and IC properties for next branch iteration
#         idx_dd += τ_opt
#         params_.a.n_targs -= 1
#         deleteat!(params_.a.J_targs, pop_idx)
#         deleteat!(params_.a.ϵ_targs, pop_idx)
#         params_.a.z0 = ddto_branch_sol.targs[1].x[:,τ_opt+1]
#         params_.a.zf_targs = params_.a.zf_targs[:,matrix_slice]

#         # Update original params with the defer node index
#         params.a.τ_targs[k] = idx_dd
#         if k == n_targs_total - 1
#             params.a.τ_targs[k+1] = idx_dd
#         end

#         # Parameter update print statements
#         cost_dd += Δcost_dd
#         ref_initial_control = ddto_branch_sol.targs[pop_idx].u[:,τ_opt+1]
#         if VERB_DDTO && (k < n_targs_total-1)
#             @printf("   Removed target %i for next branch iteration\n", λ_targ)
#         end
#     end

#     # Append time vectors to all solutions
#     for j ∈ params.a.J_targs
#         # Reconstruct time vector from solution
#         s = ddto_sol.targs[j].u[end,:]
#         τ = range(0, stop=1, length=params.a.N) |> CVector
#         ddto_sol.targs[j].t = time_dilation_control_to_wall_clock_time(s, τ, params.a.disc)
#     end

#     # Converged solution data
#     println("\nDDTO solution properties:")
#     for j = 1:params.a.n_targs
#         ϵ_subopt = (ddto_sol.targs[j].cost - ref_costs[j])/ref_costs[j] * 100
#         t_defer = ddto_sol.targs[j].t[params.a.τ_targs[j]]
#         @printf("   Target %i -- %2.1f [s] deferred, % 2.1f [%%] suboptimal.\n", j, t_defer, ϵ_subopt)
#     end 

#     return ddto_sol, true # Assume converged for now
# end