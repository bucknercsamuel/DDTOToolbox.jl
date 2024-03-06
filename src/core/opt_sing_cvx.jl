# ..:: Single-Target (Decoupled) Solver Functions ::..

function solve_tree_decoupled_cvx(params)::DDTOSolution
    # Solve the OPC for a given set of params and all targets independently
    # using `solve_optimal_target`
    #
    # :in params: The params object
    # :out solutions: Vectorized container for all single-target solutions

    # Define container for each `solve_optimal_target` solution
    solutions = EmptyDDTOSolution(params.a.n_targs)

    # Obtain solutions for each target
    VERB_OPT && println("\n=== Decoupled optimal solutions for each target ===")
    for j = 1:params.a.n_targs
        solutions.targs[j],_ = solve_target_decoupled_cvx(params, params.a.N, j)
        VERB_OPT && @printf("Target: %i, Cost: %.3f\n", params.a.T_targs[j], solutions.targs[j].cost)
    end

    return solutions
end

function solve_target_decoupled_cvx(params, N::Int, j_targ::Int)::Tuple{Solution, MOI.TerminationStatusCode}
    # Solve the optimal landing (PDG) problem for a given params and single target
    # ** (Not DDTO formulation, but used for comparison) **
    #
    # :in params: The params object.
    # :in N: Time horizon
    # :in j_targ: Target index
    # :out sol: Container for solution variables

    # ..:: Setup ::..
    # Optimizer configuration
    mdl,_ = solver_setup(SOLVER_CTCS_DISABLED)

    # Sizing parameters
    nx = params.a.nx
    nu = params.a.nu
    Δt = (params.a.Δt_min + params.a.Δt_max)/2
    tf = Δt * (N-1)
    t  = CVector(range(0, stop=tf, length=N))
    if params.a.disc == 0
        N_ctrl = N-1
    elseif params.a.disc == 1
        N_ctrl = N
    end

    # Param check(s)
    if params.a.disc != 0 && params.a.disc != 1
        error("Please select a valid discretization hold order.")
    end

    # ..:: Optimization variables ::..
    # Unscaled variables
    x_us = @variable(mdl, [1:nx,1:N])
    u_us = @variable(mdl, [1:nu,1:N_ctrl])

    # Apply affine scaling
    x = params.a.Sx*x_us .+ repeat(params.a.sx, 1, N)
    u = params.a.Su*u_us .+ repeat(params.a.su, 1, N_ctrl)
    SxInv = inv(params.a.Sx)
    SuInv = inv(params.a.Su)
    
    # ..:: Make the optimization problem ::..
    # Problem-specific construction
    J_running,J_term = prob_cost(mdl,x,u,params;nonconvex=false)
    prob_constraints(mdl,x,u,params,EmptySolution();nonconvex=false)

    # Dynamics
    X(k) = x[:,k] # State at time index k
    U(k) = u[:,k] # Input at time index k
    A_cont,B_cont,p_cont = dynamics_linear(params)
    A,Bm,Bp,p = c2d_LTI_affine(A_cont, B_cont, p_cont, Δt, params.a.disc)
    if params.a.disc == 0
        @constraint(mdl, [k=1:N-1], SxInv*X(k+1) .==  SxInv*(A*X(k) + Bm*U(k) + p))
    elseif params.a.disc == 1
        @constraint(mdl, [k=1:N-1], SxInv*X(k+1) .==  SxInv*(A*X(k) + Bm*U(k) + Bp*U(k+1) + p))
    end

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

    # ..:: Solve the problem and save the solution ::..
    J_cost = sum(J_running) + J_term
    @objective(mdl, Min, J_cost)
    optimize!(mdl)
    feas_status = JuMP.termination_status(mdl)
    if feas_status == MOI.OPTIMAL || feas_status == MOI.ALMOST_OPTIMAL
        solve_status = "Feasible"
    else
        solve_status = "Not Feasible"
        return (EmptySolution(), feas_status)
    end

    # Package the solution
    x = value.(x)
    u = value.(u)
    cost = value.(J_cost)
    sol = Solution(t,x,u,cost)

    return sol,feas_status
end