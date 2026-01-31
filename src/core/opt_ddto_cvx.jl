# ..:: Top-level Solve Function ::..
# Based on Elango et al. 2022 "Deferring Decision in Multitarget Trajectory Optimization", algorithms 1 and 2. 
# Note: deprecated currently in favor of opt_ddto_lex.jl

function solve_cvx(params; simulate_solutions=true, process_the_solutions=true, solve_ddto=true)
    # ..:: Execute solver sequence ::..
    # Apply custom scaling (if not already done)
    custom_scaling!(params)
    
    @time begin
        @time begin
            # ..:: Determine best fixed timestep using bisection-search-wrapped single-shot solutions ::..
            if params.a.gss_cvx
                Δt_opt_targs = zeros(params.a.n_targs)
                for j = 1:params.a.n_targs
                    function bisection_fun(Δt)
                        params.a.Δt_cvx = Δt
                        sol = solve_target_decoupled_cvx(params, j)[1]
                        return sol.cost
                    end
                    ϵ = 1e-3 # numerical protection
                    Δt_opt_targs[j] = (1+ϵ) * bisection_search_min_feasible(bisection_fun, params.a.Δt_min, params.a.Δt_max, verbose=false)[1]
                end
            else
                Δt_opt_targs = fill(params.a.Δt_cvx, params.a.n_targs)
            end

            # ..:: Solve for independently-optimal solutions to each target ::..
            opt_solutions = solve_tree_decoupled_cvx(params, Δt_cvx=Δt_opt_targs)
            opt_costs = CVector(zeros(params.a.n_targs))
            for k = 1:params.a.n_targs
                opt_costs[k] = opt_solutions.targs[k].cost
            end
            println("\n Solve time for generating optimal solutions to each target:")
        end

        if params.a.n_targs > 1 && solve_ddto
            @time begin
                # Compute the fixed dt using a specific update law:
                params.a.Δt_cvx = max(Δt_opt_targs...) * (1 + max(params.a.ϵ_targs...))

                # ..:: Solve for DDTO branching solutions to ALL targets ::..
                ddto_solutions = solve_tree_ddtocvx(params, opt_costs, opt_solutions)
                println("\n Solve time for generating DDTO branch solutions to all targets:")
            end
            println("\n Solve time for the full DDTO solution stack:")
        else
            ddto_solutions = copy(opt_solutions)
        end
    end

    # ..:: Simulate each target solution from I.C. to T.C.
    if simulate_solutions
        @time begin
            dynamics = (t,x,sol) -> dynamics_linear(params)[1]*x + dynamics_linear(params)[2]*optimal_controller(t,sol.t,sol.u,params.a.disc)
            opt_simulations = simulate(opt_solutions, dynamics, params.a.disc, max_steps=params.a.N_sim)
            if solve_ddto
                ddto_simulations = simulate(ddto_solutions, dynamics, params.a.disc, max_steps=params.a.N_sim)
            end
            println("\n Solve time for RK4 simulation:")
        end
    end

    # ..:: Post-processing (problem-specific) ::..
    if process_the_solutions
        @time begin
            opt_solutions    = process_solutions(opt_solutions, params)
            if solve_ddto
                ddto_solutions   = process_solutions(ddto_solutions, params)
            end
            if simulate_solutions
                opt_simulations  = process_solutions(opt_simulations, params)
                if solve_ddto
                    ddto_simulations = process_solutions(ddto_simulations, params)
                end
            end
            println("\n Solve time for post-processing:")
        end
    end

    if simulate_solutions
        if solve_ddto
            return (
                opt_solutions, 
                opt_simulations, 
                ddto_solutions, 
                ddto_simulations)
        else
            return (
                opt_solutions, 
                opt_simulations)
        end
    else
        if solve_ddto
            return (
                opt_solutions, 
                ddto_solutions)
        else
            return opt_solutions
        end
    end
end

# ..:: DDTO-Cvx Solver Functions ::..

function solve_tree_ddtocvx(params, ref_costs::CVector, ref_trajs::DDTOSolution)::DDTOSolution
    # Top-level DDTO solver for all branch points
    #
    # :in params: The params object
    # :in ref_costs: Optimal costs from initial condition
    # :out ddto_sol: Vectorized container for all DDTO branch solutions

    # Define container for each DDTO branch solution
    ddto_sol = EmptyDDTOSolution(params.a.n_targs)

    # Define running deferred-decision (DD) trajectory segment cost sum
    cost_dd = 0.

    # Initialization
    n_targs_total = copy(params.a.n_targs)
    params.a.τ_targs = zeros(n_targs_total) # initialization
    ref_initial_control = zeros(params.a.nu)
    ddto_branch_sol = ref_trajs # initialize to branch solutions
    params_ = copy(params) # Temp object to be mutated through DDTO loop
    find_J_elem(J_targs,j) = findfirst(τ->τ==j, J_targs)
    J_targs_old = copy(params.a.J_targs)
    idx_dd = 1
    τ_opt = 0

    # Perform branching in the order of preference
    for k = 1:(n_targs_total-1)
        λ_targ = params_.a.λ_targs[1]
        VERB_DDTO && @printf("\n========= Solving DDTO Stage Problem for Deferred Target #%i =========\n", λ_targ)

        # Obtain Bisection-optimal DDTO solution for this branch
        prev_sol = copy(ddto_branch_sol)
        prev_τ = copy(τ_opt)
        ddto_branch_sol,τ_opt,Δcost_dd = solve_bisection_ddtocvx(params_, ref_costs[params_.a.J_targs], cost_dd, ref_initial_control)
        if τ_opt == 0
            ddto_branch_sol = EmptyDDTOSolution(params_.a.n_targs)
            for j ∈ params_.a.J_targs
                ddto_branch_sol.targs[find_J_elem(params_.a.J_targs,j)].x = prev_sol.targs[find_J_elem(J_targs_old,j)].x[:,prev_τ+1:end]
                ddto_branch_sol.targs[find_J_elem(params_.a.J_targs,j)].u = prev_sol.targs[find_J_elem(J_targs_old,j)].u[:,prev_τ+1:end]
                ddto_branch_sol.targs[find_J_elem(params_.a.J_targs,j)].cost = prev_sol.targs[find_J_elem(J_targs_old,j)].cost
            end
        end
        J_targs_old = copy(params_.a.J_targs)

        count = 1
        for j ∈ params_.a.J_targs
            if k == 1
                ddto_sol.targs[j].x = ddto_branch_sol.targs[j].x
                ddto_sol.targs[j].u = ddto_branch_sol.targs[j].u
            else
                ddto_sol.targs[j].x[:,idx_dd:end] = ddto_branch_sol.targs[count].x
                ddto_sol.targs[j].u[:,idx_dd:end] = ddto_branch_sol.targs[count].u
            end
            ddto_sol.targs[j].cost = ddto_branch_sol.targs[count].cost
            count += 1
        end

        # Determine target to be removed (first in the current list of λ_targs)
        deleteat!(params_.a.λ_targs, 1)
        pop_idx = findfirst(i->i==λ_targ, params_.a.J_targs)

        # Have to do some slicing magic for matrices
        matrix_slice = collect(1:params_.a.n_targs)
        deleteat!(matrix_slice, pop_idx)

        # Update params_ target and IC properties for next branch iteration
        idx_dd += τ_opt
        params_.a.n_targs -= 1
        deleteat!(params_.a.J_targs, pop_idx)
        deleteat!(params_.a.ϵ_targs, pop_idx)
        params_.a.N -= τ_opt
        params_.a.z0 = ddto_branch_sol.targs[1].x[:,τ_opt+1]
        params_.a.zf_targs = params_.a.zf_targs[:,matrix_slice]

        # Update original params with the defer node index
        params.a.τ_targs[k] = idx_dd
        if k == n_targs_total - 1
            params.a.τ_targs[k+1] = idx_dd
        end

        # Parameter update print statements
        cost_dd += Δcost_dd
        ref_initial_control = ddto_branch_sol.targs[pop_idx].u[:,τ_opt+1]
        if VERB_DDTO && (k < n_targs_total-1)
            @printf("   Removed target %i for next branch iteration\n", λ_targ)
        end
    end

    # Append time vectors to all solutions
    Δt = params.a.Δt_cvx
    tf = Δt * (params.a.N-1)
    t  = CVector(range(0, stop=tf, length=params.a.N))
    for j ∈ params.a.J_targs
        ddto_sol.targs[j].t = t
    end

    # Converged solution data
    println("\nDDTO solution properties:")
    for j = 1:params.a.n_targs
        ϵ_subopt = (ddto_sol.targs[j].cost - ref_costs[j])/ref_costs[j] * 100
        t_defer = ddto_sol.targs[j].t[params.a.τ_targs[j]]
        @printf("   Target %i -- %2.1f [s] deferred, % 2.1f [%%] suboptimal.\n", j, t_defer, ϵ_subopt)
    end 

    return ddto_sol
end

function solve_bisection_ddtocvx(params, ref_costs::CVector, cost_dd::CReal, ref_initial_control::CVector)::Tuple{DDTOSolution,Int,CReal}
    # Uses bisection search to solve quasiconvex optimization problem 
    # to branch to the next-queued target for rejection.
    #
    # :in params: The params object
    # :in ref_costs: Optimal costs
    # :in cost_dd: Running cost for decision deferral
    # :out ddto_solution: Contains the DDTO solution for this target/branch point

    # Initial search bracket
    τ_min = 0
    τ_max = params.a.N - 2

    # Bisection search to solve quasiconvex (QCvx) optimization problem
    VERB_DDTO && println("=== Bisection Search for QCvx Optimization ===")
    iter = 1
    while (τ_max - τ_min) > 1
        # Update τ
        τ = Int(ceil(0.5*(τ_max + τ_min)))

        # Compute feasible DDTO
        _,status_feas,_ = solve_feasible_ddtocvx(params, τ, ref_costs, cost_dd, ref_initial_control)

        # Update τ_min or τ_max based on solution convergence
        if status_feas == MOI.OPTIMAL || status_feas == MOI.ALMOST_OPTIMAL
            solve_status = "Feasible"
        else
            solve_status = "Not Feasible"
        end
        VERB_DDTO && @printf("Iteration: %i, τ: %i for τ_min: %i, τ_max: %i -- %s\n", iter, τ, τ_min, τ_max, solve_status)
        if solve_status == "Feasible"
            τ_min = τ
        else
            τ_max = τ
        end

        # Update iteration count
        iter += 1
    end

    # Set optimal τ
    τ_opt = τ_min
    # VERB_DDTO && println("Bisection search terminated -- reached convergence condition (τ_max - τ_min) = 1")
    VERB_DDTO &&  @printf("Bisected τ_opt: %i\n", τ_opt)

    # Compute converged DDTO solution
    if τ_opt > 0
        ddto_solution,_,Δcost_dd = solve_feasible_ddtocvx(params, τ_opt, ref_costs, cost_dd, ref_initial_control)
    else
        ddto_solution = EmptyDDTOSolution(params.a.n_targs)
        Δcost_dd = 0
    end

    return ddto_solution,τ_opt,Δcost_dd
end

function solve_feasible_ddtocvx(params, τ::Int, ref_costs::CVector, cost_dd::CReal, ref_initial_control::CVector)::Tuple{DDTOSolution, MOI.TerminationStatusCode, CReal}
    # Solve the baseline feasibility problem for DDTO.
    #
    # :in params: The params object
    # :in τ: Branch point index
    # :in ref_costs: Optimal costs from `solve_optimal_pdg_all_targets`
    # :in cost_dd: Running cost for decision deferral
    # :out ddto_solution: Contains the DDTO solution for this target/branch point
    # :out feas_status: Feasibility problem solution status code (see MOI.TerminationStatusCode documentation)

    # ..:: Setup ::..
    # Optimizer configuration
    mdl,_ = solver_setup(SOLVER_CTCS_DISABLED)

    # Sizing parameters
    n = params.a.n_targs
    N = params.a.N
    nx = params.a.nx
    nu = params.a.nu
    Δt = params.a.Δt_cvx
    tf = Δt*(params.a.N-1)
    t  = CVector(range(0, stop=tf, length=params.a.N))
    if params.a.disc == 0
        N_ctrl = N-1
    elseif params.a.disc == 1
        N_ctrl = N
    end

    # Param check(s)
    if params.a.disc != 0 && params.a.disc != 1
        error("Please select a valid discretization hold order.")
    end

    # Dynamics
    A_cont,B_cont,p_cont = dynamics_linear(params)
    A,Bm,Bp,p = c2d_LTI_affine(A_cont, B_cont, p_cont, Δt, params.a.disc)

    # ..:: Optimization variables ::..
    # Unscaled variables
    x_us = @variable(mdl, [1:nx,1:N,1:n])
    u_us = @variable(mdl, [1:nu,1:N_ctrl,1:n])

    # Apply affine scaling
    x = Array{JuMP.AffExpr}(undef,nx,N,n)
    u = Array{JuMP.AffExpr}(undef,nu,N_ctrl,n)
    for j = 1:n
        x[:,:,j] = params.a.Sx*x_us[:,:,j] .+ repeat(params.a.sx, 1, N)
        u[:,:,j] = params.a.Su*u_us[:,:,j] .+ repeat(params.a.su, 1, N_ctrl)
    end
    SxInv = inv(params.a.Sx)
    SuInv = inv(params.a.Su)

    # Convenience functions
    X(k,j) = x[:,k,j] # State at time index k and target j
    U(k,j) = u[:,k,j] # Input at time index k and target j

    # ..:: Make the optimization problem ::..
    # Segment constraints
    J_cost = Vector(undef,n)
    J_running = 0
    for j = 1:n
        # Problem-specific construction
        J_running,J_term = prob_cost(mdl,x[:,:,j],u[:,:,j],params;nonconvex=false)
        prob_constraints(mdl,x[:,:,j],u[:,:,j],params,EmptySolution(),j;nonconvex=false)
        
        # Dynamics
        if params.a.disc == 0
            @constraint(mdl, [k=1:N-1], SxInv*X(k+1,j) .==  SxInv*(A*X(k,j) + Bm*U(k,j) + p))
        elseif params.a.disc == 1
            @constraint(mdl, [k=1:N-1], SxInv*X(k+1,j) .==  SxInv*(A*X(k,j) + Bm*U(k,j) + Bp*U(k+1,j) + p))
        end

        # Suboptimality constraint
        J_cost[j] = sum(J_running) + J_term
        if τ > 0
            @constraint(mdl, (cost_dd + J_cost[j]) / ((1 + params.a.ϵ_targs[j]) * ref_costs[j]) <= 1)
        end
    end
    Δcost_dd = sum(J_running[1:τ])

    # Control identicality constraints
    not(x) = x == 0 ? 1 : 0
    for j = 2:n
        for k = 1:(τ + not(N-N_ctrl))
            @constraint(mdl, SuInv*U(k,j) .== SuInv*U(k,j-1))
        end
    end

    # Control continuity constraint (if using FOH)
    # only used if we have already proceeded some amount of the trajectory, i.e. cost_dd > 0
    if params.a.disc == 1 && cost_dd > 0
        for j = 1:n
            @constraint(mdl, SuInv*U(1,j) == SuInv*ref_initial_control)
        end
    end

    # State boundary conditions
    for j = 1:n    
        for k = 1:nx
            if ~isinf(params.a.z0[k])
                @constraint(mdl, SxInv[k,k]*x[k,1,j] == SxInv[k,k]*params.a.z0[k])
            end
            if ~isinf(params.a.zf_targs[k,j])
                @constraint(mdl, SxInv[k,k]*x[k,end,j] == SxInv[k,k]*params.a.zf_targs[k,j])
            end
        end
    end

    # Input boundary conditions
    for j = 1:n    
        for k = 1:nu
            if ~isinf(params.a.u0[k])
                @constraint(mdl, SuInv[k,k]*u[k,1,j] == SuInv[k,k]*params.a.u0[k])
            end
            if ~isinf(params.a.uf_targs[k,j])
                @constraint(mdl, SuInv[k,k]*u[k,end,j] == SuInv[k,k]*params.a.uf_targs[k,j])
            end
        end
    end

    # ..:: Solve the problem and save the solution ::..
    @objective(mdl, Min, sum(J_cost))
    optimize!(mdl)
    feas_status = JuMP.termination_status(mdl)
    if feas_status == MOI.OPTIMAL || feas_status == MOI.ALMOST_OPTIMAL
        solve_status = "Feasible"
    else
        solve_status = "Not Feasible"
        return (EmptyDDTOSolution(n), feas_status, 0)
    end

    # Determine deferred-decision cost (cost up to τ)

    # Package the solution
    ddto_solution = EmptyDDTOSolution(n)
    for j = 1:n
        ddto_solution.targs[j].t = t
        ddto_solution.targs[j].x = value.(x[:,:,j])
        ddto_solution.targs[j].u = value.(u[:,:,j])
        ddto_solution.targs[j].cost = value.(J_cost[j])
    end
    Δcost_dd = value.(Δcost_dd)

    return (ddto_solution, feas_status, Δcost_dd)

end