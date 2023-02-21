#= LCvx for Quadcopter Landing -- Utility Functions.

Author: Samuel Buckner (UW-ACL)
=#

# ..:: General-Purpose Functions ::..

function optimal_controller(t::CReal, T::CVector, U::CMatrix)::CVector
    # Output the interpolated optimal control input at time t.
    # (interpolation based on hold assumption)
    #
    # :in t: the current time
    # :in T: the time signal history
    # :in sol: the input signal history
    # :out u: the interpolated input at time "t"

    # ZOH interpolation
    i = findlast(τ->τ<=t,T)
    if typeof(i)==Nothing || i>=size(U,2)
        u = U[:,end]
    else
        u = U[:,i]
    end

    # (Other interpolation methods not yet implemented)

    return u
end

function rk4(f::Function, x0::CVector, t0::CReal, tf::CReal, Δt::CReal)::Tuple{CVector, CMatrix}
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
    t = CVector(t0:Δt:tf)
    if (tf-t[end])>=√eps()
        push!(t,tf)
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

function simulate_cont(branch_solutions::Vector{BranchSolution}, dynamics::Function)::Vector{BranchSolution}
    # Simulate the dynamics of each branch solution using a predefined control input
    # trajectory in continuous time with RK4 integration.

    branch_simulations = Vector{BranchSolution}(undef, length(branch_solutions))

    for k=1:length(branch_solutions)
        sol = branch_solutions[k].sol
        dynamics_ = (t,x) -> dynamics(t,x,sol)
        n = size(sol.x,1)
        m = size(sol.u,1)
        T = CVector(undef,0)
        X = CMatrix(undef,n,0)
        U = CMatrix(undef,m,0)
        h_min = 1e-4
        for k = 1:(length(sol.t)-1)
            if k == 1
                x0 = sol.x[:,1]
            else
                x0 = X[:,end]
            end
            Δt_prop = max((1/40)*(sol.t[k+1] - sol.t[k]), h_min)
            T_,X_ = rk4(dynamics_, x0, sol.t[k], sol.t[k+1], Δt_prop)
            U_ = CMatrix(hcat([optimal_controller(T_[n],sol.t,sol.u) for n = 1:length(T_)]...))
            T = vcat(T,T_)
            X = hcat(X,X_)
            U = hcat(U,U_)
        end
        sim = Solution(T,X,U,sol.cost)
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
