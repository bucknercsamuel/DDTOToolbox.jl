#=
Single-target (decoupled) SCP / CT-SCvx solvers used as reference solutions and
warmstarts for DDTO-SCP.
=#

"""
    solve_tree_decoupled(params; single_iter=false, ref_trajs=nothing) -> (solutions, all_scp_solutions_converged)

Solve an independent SCP problem for every target. If `ref_trajs` is omitted,
builds linear initial guesses via `generate_initial_guess_scp`.

# Arguments
- `params`: problem parameters.
- `single_iter`: if `true`, run only one CT-SCvx subproblem iteration per target.
- `ref_trajs`: optional warmstart DDTO solution; generated when omitted.

# Returns
- `solutions`: decoupled SCP solutions, one solution per target.
- `all_scp_solutions_converged`: `true` only if every target subproblem converged.
"""
function solve_tree_decoupled(params; single_iter=false, ref_trajs=nothing)::Tuple{DDTOSolution,Bool}
    # Define container for each `solve_optimal_target` solution
    solutions = EmptyDDTOSolution(params.a.n_targs)

    # ..:: Define initial guess reference trajectories using linear interpolations ::..
    if isnothing(ref_trajs)
        ref_trajs = EmptyDDTOSolution(params.a.n_targs)
        for j = 1:params.a.n_targs
            ref_trajs.targs[j] = generate_initial_guess_scp(params,j)
        end
    end

    # ..:: SCP Iteration ::..
    VERB_OPT && println("\n=== Decoupled SCP solutions for each target ===")
    all_scp_solutions_converged = true
    for j = 1:params.a.n_targs
        VERB_OPT && @printf("DDTO Warmstart: target %i\n", params.a.J_targs[j])
        ref_traj = ref_trajs.targs[j]
        subproblem_ = (params_, ref_traj, k) -> solve_subproblem_decoupled_target(params_, ref_traj, j, k)
        (solution, feas_status, scvx_converged) = solve_ctscvx_iteration(params, ref_traj, subproblem_; single_iter=single_iter)
        solutions.targs[j] = solution
        all_scp_solutions_converged = all_scp_solutions_converged && scvx_converged
    end

    return solutions, all_scp_solutions_converged
end

"""
    solve_subproblem_decoupled_target(params, ref_traj, j_targ, scp_iter)

Solve one CT-SCvx subproblem for target `j_targ` about `ref_traj`, wiring
scenario dynamics/cost/constraints (with CTCS when enabled).

# Arguments
- `params`: problem parameters.
- `ref_traj`: reference solution linearized about in this SCP iteration.
- `j_targ`: target index selecting terminal conditions and dynamics variant.
- `scp_iter`: current SCP/PTR iteration index (passed to the subproblem solver).

# Returns
- wrapped solution of `solve_ctscvx_subproblem`.
"""
function solve_subproblem_decoupled_target(params, ref_traj::Solution, j_targ::Int, scp_iter::Int)

    # Define target-specific dynamics and constraints based on CTCS enablement
    if params.a.ctcs_enabled
        dynamics_linearized_ctcs = DynamicsLinearizedCTCS(params)
        dyn_nl = (t,x,u,p) -> dynamics_nonlinear_ctcs(t,x,u,params,j_targ)
        dyn_lin = (t,x,u,p) -> dynamics_linearized_ctcs(t,x,u,params,j_targ)
        prob_constraints_ = (mdl,x,u,ref_traj) -> prob_constraints(mdl,x,u,params,ref_traj,0;nonconvex=false) # only impose convex constraints
    else
        dyn_nl = (t,x,u,p) -> dynamics_nonlinear(t,x,u,params)
        dyn_lin = (t,x,u,p) -> dynamics_linearized(t,x,u,params)
        prob_constraints_ = (mdl,x,u,ref_traj) -> prob_constraints(mdl,x,u,params,ref_traj,0)
    end
    prob_cost_ = (mdl,x,u) -> prob_cost(mdl,x,u,params)

    # Wrapper function to solve the subproblem for target j_targ
    return solve_ctscvx_subproblem(
        params, 
        ref_traj, 
        params.a.z0, 
        params.a.zf_targs[:,j_targ], 
        params.a.u0, 
        params.a.uf_targs[:,j_targ],
        dyn_nl,
        dyn_lin,
        prob_cost_,
        prob_constraints_,
        scp_iter)
end