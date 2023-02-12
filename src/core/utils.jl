#= LCvx for Quadcopter Landing -- Utility Functions.

Author: Samuel Buckner (UW-ACL)
=#

# ..:: General-Purpose Functions ::..

function optimal_controller(t::CReal, x::CVector, sol::Solution)::CVector
    # Output the interpolated optimal control input at time t.
    # (interpolation based on hold assumption)
    #
    # :in t: the time at which to compute the optimal control.
    # :in x: the current state (not currently used, no feedback control)
    # :in sol: the optimized solution to track.
    # :out u: the optimal input

    # ZOH interpolation
    i = findlast(τ->τ<=t,sol.t)
    if typeof(i)==Nothing || i>=size(sol.u,2)
        u = sol.u[:,end]
    else
        u = sol.u[:,i]
    end

    # (Other interpolation methods not yet implemented)

    return u
end

function c2d_zoh(params::Params, Δt::CReal)::Tuple{CMatrix, CMatrix, CVector}
    # Discretize params dynamics at Δt time step using zeroth-order
    # hold (ZOH).
    #
    # :in params: the params object.
    # :in Δt: the discretization time step.
    # :out : a tuple (A,B,p) for the discrete-time update equation
    #         x_{k+1} = A*x_k + B*u_k + p

    A_c,B_c,p_c,n,m = params.A_c,params.B_c,params.p_c,params.n,params.m

    _M = exp(CMatrix([
        A_c B_c p_c;
        zeros(m+1,n+m+1)
    ])*Δt)

    A = _M[1:n,1:n]
    B = _M[1:n,n+1:n+m]
    p = _M[1:n,n+m+1]
    return (A,B,p)
end

function rk4(f::Function, x0::CVector, Δt::CReal, T::CReal)::Tuple{CVector, CMatrix}
    # Integrate a system of ordinary differential equations (ODE)
    # using RK4.
    #
    # :in f: the function defining the ODE, dx/dt=f(t,x).
    # :in x0: the initial condition.
    # :in Δt: the integration time step.
    # :in T: the integration final time.
    # :out : a vector storing the integration times.
    # :out : a matrix storing in its columns the integrated state
    #        trajectory.

    # ..:: Make time grid ::..
    t = CVector(0.0:Δt:T)
    if (T-t[end])>=√eps()
        push!(t,T)
    end
    N = length(t)

    # ..:: Initialize ::..
    X = CMatrix(undef,length(x0),N)
    X[:,1] = x0

    # ..:: Integrate ::..
    for n = 1:N-1
        y = X[:,n]
        h = t[n+1]-t[n]
        t_ = t[n]
        k1 = f(t_,y)
        k2 = f(t_+h/2,y+h*k1/2)
        k3 = f(t_+h/2,y+h*k2/2)
        k4 = f(t_+h,y+h*k3)
        X[:,n+1] = y+h/6*(k1+2*k2+2*k3+k4)
    end

    return (t,X)
end

function simulate_cont(branch_solutions::Vector{BranchSolution}, x0::CVector, dynamics::Function, Δt::CReal=1e-2)::Vector{BranchSolution}
    # Simulate the dynamics of each branch solution using a predefined control input
    # trajectory in continuous time with RK4 integration.

    branch_simulations = Vector{BranchSolution}(undef, length(branch_solutions))

    for k=1:length(branch_solutions)
        sol = branch_solutions[k].sol
        dynamics_ = (t,x) -> dynamics(t,x,sol)
        tf = sol.t[end]
        t,X = rk4(dynamics_, x0, Δt, tf)
        U = CMatrix(hcat([optimal_controller(t[n],X[:,n],sol) for n = 1:length(t)]...))
        sim = Solution(t,X,U,sol.cost)
        branch_simulations[k] = BranchSolution(sim,branch_solutions[k].cost_dd,branch_solutions[k].idx_dd)
    end

    return branch_simulations
end

# ..:: DDTO Functions ::..

function extract_target_trajectories(params::Params, sols_ddto::Array{DDTOSolution})::Vector{BranchSolution}

    # Obtain full solutions to each target
    DDTO_target_solutions = Vector{BranchSolution}(undef, params.n_targs)
    net_deferral_idx = 1
    net_cost_dd = 0
    leading_cost_traj = 0
    T_targs = copy(params.T_targs)
    x_trunk = CMatrix(undef,params.n,0)
    u_trunk = CMatrix(undef,params.m,0)

    for j in 1:params.n_targs

        # Obtain branch to the desired target
        deferral_idx = sols_ddto[j].idx_dd
        net_deferral_idx += deferral_idx
        λ_targ = params.λ_targs[j]
        rej_idx = findfirst(i->i==λ_targ, T_targs)
        deleteat!(T_targs,rej_idx)
        sol_branch = sols_ddto[j].targ_sols[rej_idx]
        x_branch = sol_branch.x
        u_branch = sol_branch.u

        # Compute costs
        if j < params.n_targs
            leading_cost_traj = net_cost_dd
        end
        total_cost = sols_ddto[j].targ_sols[rej_idx].cost + leading_cost_traj
        net_cost_dd += sols_ddto[j].cost_dd

        # Concatenate to create the solution to the given target
        x_target = hcat(x_trunk, x_branch)
        u_target = hcat(u_trunk, u_branch)
        t_target = collect(0:length(x_target[1,:])-1) * params.Δt
        sol_target = Solution(t_target, x_target, u_target, total_cost)

        # Build the "trunk" to the deferral point for the next solution
        x_trunk = hcat(x_trunk, sols_ddto[j].targ_sols[1].x[:,1:deferral_idx])
        u_trunk = hcat(u_trunk, sols_ddto[j].targ_sols[1].u[:,1:deferral_idx])

        # Add solution
        def_idx = findfirst(i->i==λ_targ, params.T_targs)
        DDTO_target_solutions[def_idx] = EmptyBranchSolution()
        DDTO_target_solutions[def_idx].sol = sol_target
        DDTO_target_solutions[def_idx].cost_dd = net_cost_dd
        DDTO_target_solutions[def_idx].idx_dd = net_deferral_idx
    end

    return DDTO_target_solutions
end
