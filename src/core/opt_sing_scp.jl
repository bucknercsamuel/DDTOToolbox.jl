# ..:: Single-Target (Decoupled) Solver Functions ::..

function solve_tree_decoupled(params; single_iter=false, ref_trajs=nothing)::DDTOSolution
    # Solve the OPC for a given set of params and all targets independently
    #
    # :in params: The params object
    # :out solutions: Vectorized container for all single-target solutions

    # Define container for each `solve_optimal_target` solution
    solutions = EmptyDDTOSolution(params.n_targs)

    # ..:: Define initial guess reference trajectories using linear interpolations ::..
    if isnothing(ref_trajs)
        ref_trajs = EmptyDDTOSolution(params.n_targs)
        for j = 1:params.n_targs
            ref_trajs.targs[j] = generate_initial_guess_scp(params,j)
        end
    end

    # ..:: SCP Iteration ::..
    VERB_OPT && println("\n=== Decoupled SCP solutions for each target ===")
    for j = 1:params.n_targs
        VERB_OPT && @printf("Target: %i\n", params.T_targs[j])
        feas_status = undef
        solution = ref_trajs.targs[j]
        scp_converged = false

        for k = 1:params.scp_iters 
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
        solutions.targs[j] = solution
        VERB_OPT && @printf("   > Total cost: %.3f\n\n", solution.cost)
    end

    return solutions
end

function solve_subproblem_decoupled(params, ref_traj::Solution, j_targ::Int, scp_iter::Int)::Tuple{Solution, MOI.TerminationStatusCode, Bool}

    # ..:: Setup ::..
    # Optimizer configuration
    if SOLVER == "ECOS"
        mdl = Model(optimizer_with_attributes(ECOS.Optimizer, "verbose" => 0, "max_iters" => 1000))
    elseif SOLVER == "MOSEK"
        mdl = Model(Mosek.Optimizer)
        JuMP.set_optimizer_attribute(mdl, "LOG",  0) # disable debugging
        JuMP.set_optimizer_attribute(mdl, "MAX_NUM_WARNINGS", 0) # disable warnings
    else
        error("SOLVER is invalid, please select either ECOS or MOSEK")
    end

    # Sizing parameters
    nx = params.nx
    nu = params.nu
    N  = params.N
    if params.disc == 0
        N = N-1
    elseif params.disc == 1
        N = N
    end

    # Param check(s)
    if params.disc != 0 && params.disc != 1
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
    x = params.Sx*x_us .+ repeat(params.sx, 1, N)
    u = params.Su*u_us .+ repeat(params.su, 1, N)

    # SCP-specific
    ν_ctrl = @variable(mdl, [1:nx,1:(N-1)]) # virtual control
    η = @variable(mdl, [1:N]) # trust region
    @variable(mdl, μ_ctrl) # virtual control slack
    @variable(mdl, μ_buff) # virtual buffer slack
    @variable(mdl, η_s) # trust region slack

    # ..:: Make the optimization problem ::..

    # Path constraints (problem-specific)
    if !params.ctcs_enabled
        J_running,J_term,ν_buff = core_problem(mdl,x,u,params,ref_traj)
    else
        J_running,J_term = objective_function(mdl,x,u,params)
        ν_buff = []
    end

    # Dynamics
    X(k) = x[:,k]
    U(k) = u[:,k]
    if params.ctcs_enabled
        dyn_lin = (t,x,u,p) -> dynamics_linearized_ctcs(t,x,u,params)
        dyn_nl  = (t,x,u,p) -> dynamics_nonlinear_ctcs(t,x,u,params)
    else
        dyn_lin = (t,x,u,p) -> dynamics_linearized(t,x,u,params)
        dyn_nl  = (t,x,u,p) -> dynamics_nonlinear(t,x,u,params)
    end
    Ak,Bmk,Bpk,_,wk,_,_ = c2d_nonlinear(t_ref,x_ref,u_ref,dyn_nl,dyn_lin,params.disc)
    SxInv = inv(params.Sx)
    SuInv = inv(params.Su)
    if params.disc == 0
        @constraint(mdl, [k=1:N-1], SxInv*X(k+1) .== SxInv*(Ak[:,:,k]*X(k) + Bmk[:,:,k]*U(k) + wk[:,k]) + ν_ctrl[:,k])
    elseif params.disc == 1
        @constraint(mdl, [k=1:N-1], SxInv*X(k+1) .== SxInv*(Ak[:,:,k]*X(k) + Bmk[:,:,k]*U(k) + Bpk[:,:,k]*U(k+1) + wk[:,k]) + ν_ctrl[:,k])
    end

    # CTCS violation
    if params.ctcs_enabled
        @constraint(mdl, [k=1:N], x[end,k]/params.ϵ_ctcs <= 1)
        @constraint(mdl, [k=1:N], x[end,k] >= 0)
    end

    # Time dilation
    s = u[end,:]
    t = time_dilation_control_to_wall_clock_time(s, ref_traj.t, params.disc)
    Δt = diff(t)
    @constraint(mdl, t[end]/params.ToF_max <= 1)
    @constraint(mdl, [k=1:N-1], params.Δt_min/params.Δt_max <= Δt[k]/params.Δt_max <= 1)
    @constraint(mdl, [k=1:N], s[k] >= 0)

    # Boundary conditions
    z0 = params.z0
    zf = params.zf_targs[:,j_targ]
    nbd = params.ctcs_enabled ? nx-1 : nx # no boundary conditions to apply for CTCS state
    for k = 1:nbd # inf = no boundary condition to be applied
        if ~isinf(z0[k])
            @constraint(mdl, SxInv[k,k]*x[k,1] == SxInv[k,k]*z0[k])
        end
        if ~isinf(zf[k])
            @constraint(mdl, SxInv[k,k]*x[k,N] == SxInv[k,k]*zf[k])
        end
    end

    # Trust region constraints
    δX(k) = SxInv*(X(k) .- x_ref[:,k])
    δU(k) = SuInv*(U(k) .- u_ref[:,k])
    @constraint(mdl, [k=1:N], δX(k)'*δX(k) + δU(k)'*δU(k) <= η[k])
    if N < N
        @constraint(mdl, δX(N)'*δX(N) <= η[N])
    end

    # SCP slack constraints
    @constraint(mdl, vcat(μ_ctrl, vec(ν_ctrl)) in MOI.NormOneCone(length(vec(ν_ctrl))+1))
    if length(ν_buff) > 0
        @constraint(mdl, vcat(μ_buff, vec(ν_buff)) in MOI.NormOneCone(length(vec(ν_buff))+1))
        @constraint(mdl, μ_buff >= 0)
    else
        @constraint(mdl, μ_buff == 0)
    end
    @constraint(mdl, vcat(η_s, η) in SecondOrderCone())
    @constraint(mdl, μ_ctrl >= 0)
    @constraint(mdl, η_s >= 0)

    # ..:: Solve the problem and save the solution ::..
    # Cost function
    obj_scale = 1/sqrt(max(params.w_obj_sing, params.w_trust, params.w_buff, params.w_ctrl))
    J_cost = sum(J_running) + J_term
    @objective(mdl, Min, 
        (params.w_obj_sing * J_cost 
      + params.w_trust * η_s 
      + params.w_buff * μ_buff
      + params.w_ctrl * μ_ctrl)*obj_scale)

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
    if (μ_ctrl_pen <= params.ϵ_ctrl) && (μ_buff_pen <= params.ϵ_buff) && (η_pen <= params.ϵ_trust)
        scp_sub_cvged = true
    else
        scp_sub_cvged = false
    end
    @printf("   SCP Iter: %2.i | Status: %s | Cost = % .2e | μ_ctrl_pen = % .2e | μ_buff_pen = % .2e | η_pen = % .2e\n", scp_iter, solve_status, cost, μ_ctrl_pen, μ_buff_pen, η_pen)
    flush(stdout)

    return (sol, feas_status, scp_sub_cvged)
end