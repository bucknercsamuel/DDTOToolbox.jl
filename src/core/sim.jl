#= LCvx for Quadcopter Landing -- Utility Functions.

Author: Samuel Buckner (UW-ACL)
=#

# ..:: General-Purpose Functions ::..

function optimal_controller(t::CReal, T::CVector, U::CMatrix, disc::Int)::CVector
    # Output the interpolated optimal control input at time t.
    # (interpolation based on hold assumption)
    #
    # :in t: the current time
    # :in T: the time signal history
    # :in sol: the input signal history
    # :out u: the interpolated input at time "t"

    i = findlast(τ->τ<=t,T)
    if typeof(i)==Nothing || i>=size(U,2)
        u = U[:,end]
    else
        # ZOH interpolation
        if disc == 0
            u = U[:,i]
        # FOH interpolation
        elseif disc == 1
            i_ = i
            _i = i_ + 1
            t_ = T[i_]
            _t = T[_i]
            u_ = U[:,i_]
            _u = U[:,_i]
            u  = u_ + (t - t_)/(_t - t_)*(_u - u_)
        else
            error("Please select a valid discretization hold order.")
        end
    end

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

function simulate_cont(sol::Solution, dyn::Function, disc::Int; max_steps::Int=40)::Solution
    # Simulate the dynamics of the solution using a predefined control input
    # trajectory in continuous time with RK4 integration.

    dyn_ = (t,x) -> dyn(t,x,sol)
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
        Δt_prop = max((1/max_steps)*(sol.t[k+1] - sol.t[k]), h_min)
        T_,X_ = rk4(dyn_, x0, sol.t[k], sol.t[k+1], Δt_prop)
        U_ = CMatrix(hcat([optimal_controller(T_[n],sol.t,sol.u,disc) for n = 1:length(T_)]...))
        T = vcat(T,T_)
        X = hcat(X,X_)
        U = hcat(U,U_)
    end
    sim = Solution(T,X,U,sol.cost)

    return sim
end

function simulate_branches(branch_solutions::Vector{BranchSolution}, dynamics::Function, disc::Int)::Vector{BranchSolution}
    # Run `simulate_cont` for each branch of the provided solution set
    branch_simulations = Vector{BranchSolution}(undef, length(branch_solutions))
    for k=1:length(branch_solutions)
        sim = simulate_cont(branch_solutions[k].sol, dynamics, disc)
        branch_simulations[k] = BranchSolution(sim,branch_solutions[k].cost_dd,branch_solutions[k].idx_dd)
    end

    return branch_simulations
end

function extract_target_trajectories(params, sols_ddto::Array{DDTOSolution}; SCP=false)::Tuple{Vector{BranchSolution},Vector{BranchSolution}}

    # Obtain full solutions to each target
    DDTO_target_solutions = Vector{BranchSolution}(undef, params.n_targs)
    net_deferral_idx = 1
    net_cost_dd = 0
    leading_cost_traj = 0
    T_targs = copy(params.T_targs)
    n_ = size(sols_ddto[1].targ_sols[1].x,1)
    m_ = size(sols_ddto[1].targ_sols[1].u,1)
    if !SCP t_trunk = CVector(undef,0) end
    x_trunk = CMatrix(undef,n_,0)
    u_trunk = CMatrix(undef,m_,0)
    t_offset = 0

    for j in 1:params.n_targs

        # Obtain branch to the desired target
        deferral_idx = sols_ddto[j].idx_dd
        net_deferral_idx += deferral_idx
        λ_targ = params.λ_targs[j]
        rej_idx = findfirst(i->i==λ_targ, T_targs)
        deleteat!(T_targs,rej_idx)
        sol_branch = sols_ddto[j].targ_sols[rej_idx]
        t_branch = deepcopy(sol_branch.t)
        x_branch = deepcopy(sol_branch.x)
        u_branch = deepcopy(sol_branch.u)
        if !SCP
            t_offset_  = copy(t_branch[deferral_idx+1])
            t_branch .+= t_offset
            t_offset  += copy(t_offset_)
        end

        # Compute costs
        if j < params.n_targs
            leading_cost_traj = net_cost_dd
        end
        total_cost = sols_ddto[j].targ_sols[rej_idx].cost + leading_cost_traj
        net_cost_dd += sols_ddto[j].cost_dd

        # Concatenate to create the solution to the given target
        if !SCP 
            t_target = vcat(t_trunk, t_branch[1:end]) 
        else
            t_target = params.τ
        end
        x_target = hcat(x_trunk, x_branch[:,1:end])
        u_target = hcat(u_trunk, u_branch[:,1:end])

        # Build the "trunk" to the deferral point for the next solution
        if !SCP t_trunk = vcat(t_trunk, t_branch[1:deferral_idx]) end
        x_trunk = hcat(x_trunk, x_branch[:,1:deferral_idx])
        u_trunk = hcat(u_trunk, u_branch[:,1:deferral_idx])

        # Add solution
        def_idx = findfirst(i->i==λ_targ, params.T_targs)
        sol_target = Solution(t_target, x_target, u_target, total_cost)
        DDTO_target_solutions[def_idx] = EmptyBranchSolution()
        DDTO_target_solutions[def_idx].sol = sol_target
        DDTO_target_solutions[def_idx].cost_dd = net_cost_dd
        DDTO_target_solutions[def_idx].idx_dd = net_deferral_idx
    end

    if SCP t_trunk = params.τ[1:size(x_trunk,2)] end
    DDTO_trunk = Vector{BranchSolution}(undef, 1)
    DDTO_trunk_sol = Solution(t_trunk, x_trunk, u_trunk, net_cost_dd)
    DDTO_trunk[1] = EmptyBranchSolution()
    DDTO_trunk[1].sol = DDTO_trunk_sol
    DDTO_trunk[1].cost_dd = -1
    DDTO_trunk[1].idx_dd = -1

    return (DDTO_target_solutions, DDTO_trunk)
end

function time_dilation_control_to_wall_clock_time(∂t_∂τ::Vector, dτ_grid::CReal, disc::Int)
    # Converts time-dilation control to wall-clock time based on discretization method
    if length(∂t_∂τ) > 1
        Δt = Vector(undef,length(∂t_∂τ)-1)
        for k=1:length(Δt)
            if disc == 0
                Δt[k] = dτ_grid * ∂t_∂τ[k]
            elseif disc == 1
                Δt[k] = (1/2) * dτ_grid * (∂t_∂τ[k] + ∂t_∂τ[k+1])
            end
        end
        t = cumsum([0.;Δt])
    else
        t = [0.]
    end
    return t
end

function wall_clock_time_to_time_dilation_control(t::Vector, dτ_grid::CReal, disc::Int)
    # Converts wall-clock time to time dilation control based on discretization method
    N = length(t)
    if disc == 0
        N_ctrl = N-1
    elseif disc == 1
        N_ctrl = N
    end
    Δt = diff(t)
    ∂t_∂τ = Vector(undef,N_ctrl)
    if disc == 0
        for k=1:N_ctrl
            ∂t_∂τ[k] = Δt[k] / dτ_grid
        end
    elseif disc == 1
        n = length(Δt)
        ∂t_∂τ[1] = sum([Δt[k] / dτ_grid for k=1:n])/length(Δt) # boundary condition chosen for numerical properties, but is technically arbitrary!
        for k=1:N_ctrl-1
            ∂t_∂τ[k+1] = 2 * Δt[k] / dτ_grid - ∂t_∂τ[k]
        end
    end
    return ∂t_∂τ
end