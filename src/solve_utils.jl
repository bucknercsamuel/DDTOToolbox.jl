#= LCvx for Quadcopter Landing -- Utility Functions.

Author: Samuel Buckner (UW-ACL)
=#

# ..:: General-Purpose Functions ::..

function optimal_controller(t::CReal, x::CVector, sol::Solution)::CVector
    # Output the optimal control input (quadcopter thrust vector) at time t.
    #
    # :in t: the time at which to compute the optimal control.
    # :in x: the current state (not currently used, no feedback control)
    # :in sol: the optimized solution to track.
    # :out u: the optimal input for the state-space form of the quadcopter
    #         dynamics.

    # Get current optimal thrust (ZOH interpolation)
    i = findlast(τ->τ<=t,sol.t)
    if typeof(i)==Nothing || i>=size(sol.T,2)
        T = sol.T[:,end]
    else
        T = sol.T[:,i]
    end
   
    # Create the input vector for the state-space dynamics
    U = CVector(vcat(T,norm(T,2)))

    return U
end

function c2d_zoh(lander::Lander, Δt::CReal)::Tuple{CMatrix, CMatrix, CVector}
    # Discretize lander dynamics at Δt time step using zeroth-order
    # hold (ZOH).
    #
    # :in lander: the lander object.
    # :in Δt: the discretization time step.
    # :out : a tuple (A,B,p) for the discrete-time update equation
    #         x_{k+1} = A*x_k + B*u_k + p

    A_c,B_c,p_c,n,m = lander.A_c,lander.B_c,lander.p_c,lander.n,lander.m

    _M = exp(CMatrix([
        A_c B_c p_c;
        zeros(m+1,n+m+1)
    ])*Δt)

    A = _M[1:n,1:n]
    B = _M[1:n,n+1:n+m]
    p = _M[1:n,n+m+1]
    return (A,B,p)
end

function rk4_step(x_cur::CVector, f::Function, t_cur::CReal, Δt::CReal)::CVector
    # Integrate a system of ordinary differential equations (ODE)
    # one time-step forward using RK4 (updates x_cur in place)
    #
    # :in x_cur: the current state
    # :in f: the function defining the ODE, dx/dt=f(t,x).
    # :in t_cur: the current time (in DDTO solution)
    # :in Δt: the integration time step.
    # :out x_new: the new state

    # ..:: Integrate one time-step forward ::..
    y = x_cur
    h = Δt
    t_ = t_cur
    k1 = f(t_,y)
    k2 = f(t_+h/2,y+h*k1/2)
    k3 = f(t_+h/2,y+h*k2/2)
    k4 = f(t_+h,y+h*k3)
    x_cur = y+h/6*(k1+2*k2+2*k3+k4)

    return x_cur
end

# ..:: DDTO Functions ::..

function extract_target_trajectories(lander::Lander, sols_ddto::Array{DDTOSolution})::Vector{BranchSolution}

    # Obtain FULL solutions to each target
    DDTO_target_solutions = Vector{BranchSolution}(undef, lander.n_targs)
    net_deferral_idx = 1
    net_cost_dd = 0
    leading_cost_traj = 0
    T_targs = copy(lander.T_targs)
    r_trunk = CMatrix(undef, 3, 0)
    v_trunk = CMatrix(undef, 3, 0)
    T_trunk = CMatrix(undef, 3, 0)
    Γ_trunk = CVector(undef, 0)

    for j in 1:lander.n_targs

        # Obtain branch to the desired target
        deferral_idx = sols_ddto[j].idx_dd
        net_deferral_idx += deferral_idx
        λ_targ = lander.λ_targs[j]
        rej_idx = findfirst(i->i==λ_targ, T_targs)
        deleteat!(T_targs,rej_idx)
        sol_branch = sols_ddto[j].targ_sols[rej_idx]
        r_branch = sol_branch.r
        v_branch = sol_branch.v
        T_branch = sol_branch.T
        Γ_branch = sol_branch.Γ
        r0_relax_branch = sol_branch.r0_relax
        rf_relax_branch = sol_branch.rf_relax

        # Compute costs
        if j < lander.n_targs
            leading_cost_traj = net_cost_dd
        end
        total_cost = sols_ddto[j].targ_sols[rej_idx].cost + leading_cost_traj
        net_cost_dd += sols_ddto[j].cost_dd

        # Concatenate to create the solution to the given target
        r_target = hcat(r_trunk, r_branch)
        v_target = hcat(v_trunk, v_branch)
        T_target = hcat(T_trunk, T_branch)
        Γ_target = vcat(Γ_trunk, Γ_branch)
        t_target = collect(0:length(r_target[1,:])-1) * lander.Δt
        T_nrm_target = CVector([norm(T_target[:,i],2) for i=1:length(Γ_target)])
        γ_target = CVector([acos(dot(T_target[:,k],e_z)/norm(T_target[:,k],2)) for k=1:length(Γ_target)])
        sol_target = Solution(t_target, r_target, v_target, T_target, Γ_target, r0_relax_branch, rf_relax_branch, total_cost, T_nrm_target, γ_target)

        # Build the "trunk" to the deferral point for the next solution
        r_trunk = hcat(r_trunk, sols_ddto[j].targ_sols[1].r[:,1:deferral_idx])
        v_trunk = hcat(v_trunk, sols_ddto[j].targ_sols[1].v[:,1:deferral_idx])
        T_trunk = hcat(T_trunk, sols_ddto[j].targ_sols[1].T[:,1:deferral_idx])
        append!(Γ_trunk, sols_ddto[j].targ_sols[1].Γ[1:deferral_idx])

        # Add solution
        def_idx = findfirst(i->i==λ_targ, lander.T_targs)
        DDTO_target_solutions[def_idx] = EmptyBranchSolution()
        DDTO_target_solutions[def_idx].sol = sol_target
        DDTO_target_solutions[def_idx].cost_dd = net_cost_dd
        DDTO_target_solutions[def_idx].idx_dd = net_deferral_idx
    end

    return DDTO_target_solutions
end
