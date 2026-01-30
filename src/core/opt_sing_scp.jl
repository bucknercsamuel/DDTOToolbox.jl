# ..:: Single-Target (Decoupled) Solver Functions ::..

function solve_tree_decoupled(params; single_iter=false, ref_trajs=nothing)::Tuple{DDTOSolution,Bool}
    # Solve the OPC for a given set of params and all targets independently
    #
    # :in params: The params object
    # :out solutions: Vectorized container for all single-target solutions

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
        subproblem_ = (ref_traj, k) -> solve_subproblem_decoupled_target(params, ref_traj, j, k)
        (solution, feas_status, scvx_converged) = solve_ctscvx_iteration(params, ref_traj, subproblem_; single_iter=single_iter)
        solutions.targs[j] = solution
        all_scp_solutions_converged = all_scp_solutions_converged && scvx_converged
    end

    return solutions, all_scp_solutions_converged
end

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