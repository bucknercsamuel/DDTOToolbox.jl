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
        VERB_OPT && @printf("Target: %i\n", params.a.T_targs[j])
        feas_status = undef
        solution = ref_trajs.targs[j]
        scp_converged = false

        for k = 1:params.a.scp_iters 
            # Solve SCP subproblem
            (solution, feas_status, scp_converged) = solve_subproblem_decoupled(params, solution, j, k)

            if single_iter
                break # skip all convergence criterion, only going to run a single (potentially-infeasible) iterate!
            end

            if feas_status != MOI.OPTIMAL && feas_status != MOI.ALMOST_OPTIMAL
                @printf("   > SCP subproblem is infeasible (MOI status: %s), exiting subproblem iteration.\n", feas_status)
                break
            end
            if scp_converged
                @printf("   > Convergence condition has been reached, exiting subproblem iteration.\n")
                break
            end
        end
        if !scp_converged
            all_scp_solutions_converged = false
        end
        solutions.targs[j] = solution
        VERB_OPT && @printf("   > Total cost: %.3f\n\n", solution.cost)
    end

    return solutions, all_scp_solutions_converged
end

function solve_subproblem_decoupled(params, ref_traj::Solution, j_targ::Int, scp_iter::Int)::Tuple{Solution, MOI.TerminationStatusCode, Bool}

    # ..:: Setup ::..
    # Optimizer configuration
    if params.a.ctcs_enabled
        mdl, solver_type = solver_setup(SOLVER_CTCS_ENABLED)
    else
        mdl, solver_type = solver_setup(SOLVER_CTCS_DISABLED)
    end
    trust_region_type = solver_type
    # trust_region_type = "QP"

    # Sizing parameters
    nx = params.a.nx
    nu = params.a.nu
    N  = params.a.N
    if params.a.disc == 0
        N = N-1
    elseif params.a.disc == 1
        N = N
    end

    # Param check(s)
    if params.a.disc != 0 && params.a.disc != 1
        error("Please select a valid discretization hold order.")
    end
    
    # ..:: Reference trajectory ::..
    t_ref = ref_traj.t
    x_ref = ref_traj.x
    u_ref = ref_traj.u
    x_ref, u_ref = remove_ref_zeros(x_ref, u_ref)

    # ..:: Optimization variables ::..
    # Unscaled variables
    x_us = @variable(mdl, [1:nx,1:N])
    u_us = @variable(mdl, [1:nu,1:N])

    # Apply affine scaling
    x = params.a.Sx*x_us .+ repeat(params.a.sx, 1, N)
    u = params.a.Su*u_us .+ repeat(params.a.su, 1, N)

    # SCP-specific
    ν_ctrl = @variable(mdl, [1:nx,1:(N-1)]) # virtual control
    if solver_type != "QP"
        η = @variable(mdl, [1:N]) # trust region
    end
    @variable(mdl, μ_ctrl) # virtual control objective slack
    @variable(mdl, μ_buff) # virtual buffer objective slack
    @variable(mdl, η_s) # trust region objective slack

    # ..:: Make the optimization problem ::..

    # Path constraints (problem-specific)
    J_running,J_term = prob_cost(mdl,x,u,params)
    if !params.a.ctcs_enabled
        ν_buff = prob_constraints(mdl,x,u,params,ref_traj)
    else
        ν_buff = []
    end

    # Dynamics
    X(k) = x[:,k]
    U(k) = u[:,k]
    if params.a.ctcs_enabled
        dyn_lin = (t,x,u,p) -> dynamics_linearized_ctcs(t,x,u,params,j_targ)
        dyn_nl  = (t,x,u,p) -> dynamics_nonlinear_ctcs(t,x,u,params,j_targ)
    else
        dyn_lin = (t,x,u,p) -> dynamics_linearized(t,x,u,params)
        dyn_nl  = (t,x,u,p) -> dynamics_nonlinear(t,x,u,params)
    end
    Ak,Bmk,Bpk,_,wk,_,_ = c2d_nonlinear(t_ref,x_ref,u_ref,dyn_nl,dyn_lin,params.a.disc)
    SxInv = inv(params.a.Sx)
    SuInv = inv(params.a.Su)
    if params.a.disc == 0
        @constraint(mdl, [k=1:N-1], SxInv*X(k+1) .== SxInv*(Ak[:,:,k]*X(k) + Bmk[:,:,k]*U(k) + wk[:,k]) + ν_ctrl[:,k])
    elseif params.a.disc == 1
        @constraint(mdl, [k=1:N-1], SxInv*X(k+1) .== SxInv*(Ak[:,:,k]*X(k) + Bmk[:,:,k]*U(k) + Bpk[:,:,k]*U(k+1) + wk[:,k]) + ν_ctrl[:,k])
    end

    # CTCS violation
    if params.a.ctcs_enabled
        @constraint(mdl, [k=1:N], x[end,k]/params.a.ϵ_ctcs <= 1)
        @constraint(mdl, [k=1:N], x[end,k] >= 0)
    end

    # Time dilation
    s = u[end,:]
    t = time_dilation_control_to_wall_clock_time(s, ref_traj.t, params.a.disc)
    Δt = diff(t)
    @constraint(mdl, params.a.ToF_min/params.a.ToF_max <= t[end]/params.a.ToF_max <= 1)
    @constraint(mdl, [k=1:N-1], params.a.Δt_min/params.a.Δt_max <= Δt[k]/params.a.Δt_max <= 1)
    @constraint(mdl, [k=1:N], s[k] >= 0)

    # State boundary conditions
    z0 = params.a.z0
    zf = params.a.zf_targs[:,j_targ]
    for k = 1:nx
        if ~isinf(z0[k])
            @constraint(mdl, SxInv[k,k]*x[k,1] == SxInv[k,k]*z0[k])
        end
        if ~isinf(zf[k])
            @constraint(mdl, SxInv[k,k]*x[k,N] == SxInv[k,k]*zf[k])
        end
    end

    # Input boundary conditions
    u0 = params.a.u0
    uf = params.a.uf_targs[:,j_targ]
    for k = 1:nu
        if ~isinf(u0[k])
            @constraint(mdl, SuInv[k,k]*u[k,1] == SuInv[k,k]*u0[k])
        end
        if ~isinf(uf[k])
            @constraint(mdl, SuInv[k,k]*u[k,N] == SuInv[k,k]*uf[k])
        end
    end

    # Trust region constraints
    δX(k) = SxInv*(X(k) .- x_ref[:,k])
    δU(k) = SuInv*(U(k) .- u_ref[:,k])
    if trust_region_type == "QP"
        η_s = sum([δX(k)'*δX(k) + δU(k)'*δU(k) for k=1:N])
    else
        @constraint(mdl, [k=1:N], δX(k)'*δX(k) + δU(k)'*δU(k) <= η[k])
        @constraint(mdl, vcat(η_s, η) in SecondOrderCone())
        @constraint(mdl, η_s >= 0)
    end

    # Virtualization constraints
    @constraint(mdl, vcat(μ_ctrl, vec(ν_ctrl)) in MOI.NormOneCone(length(vec(ν_ctrl))+1))
    if length(ν_buff) > 0
        @constraint(mdl, vcat(μ_buff, vec(ν_buff)) in MOI.NormOneCone(length(vec(ν_buff))+1))
        @constraint(mdl, μ_buff >= 0)
    else
        @constraint(mdl, μ_buff == 0)
    end
    @constraint(mdl, μ_ctrl >= 0)

    # ..:: Solve the problem and save the solution ::..
    # Cost function
    obj_scale = 1/sqrt(max(params.a.w_obj_sing, params.a.w_trust, params.a.w_buff, params.a.w_ctrl))
    J_cost = sum(J_running) + J_term
    @objective(mdl, Min, 
        (params.a.w_obj_sing * J_cost 
      + params.a.w_trust * η_s 
      + params.a.w_buff * μ_buff
      + params.a.w_ctrl * μ_ctrl)*obj_scale)

    optimize!(mdl)
    feas_status = JuMP.termination_status(mdl)
    if feas_status == MOI.OPTIMAL || feas_status == MOI.ALMOST_OPTIMAL
        solve_status = "Feasible"
    else
        solve_status = "Not Feasible"
        return (EmptySolution(), feas_status, false)
    end

    # Package the solution
    x = value.(x)
    u = value.(u)
    cost = value.(J_cost)
    sol = Solution(ref_traj.t,x,u,cost)

    # Obtain evaluation penalties
    μ_buff_pen = value.(μ_buff)
    μ_ctrl_pen = value.(μ_ctrl)
    η_pen = value.(η_s)

    # Determine convergence based on SCP penalties
    if (μ_ctrl_pen <= params.a.ϵ_ctrl) && (μ_buff_pen <= params.a.ϵ_buff) && (η_pen <= params.a.ϵ_trust)
        scp_sub_cvged = true
    else
        scp_sub_cvged = false
    end
    @printf("   SCP Iter: %2.i | Status: %s | Cost = % .2e | μ_ctrl_pen = % .2e | μ_buff_pen = % .2e | η_pen = % .2e\n", scp_iter, solve_status, cost, μ_ctrl_pen, μ_buff_pen, η_pen)
    flush(stdout)

    return (sol, feas_status, scp_sub_cvged)
end