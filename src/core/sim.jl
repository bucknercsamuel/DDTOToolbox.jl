#= LCvx for Quadcopter Landing -- Simulation Utility Functions.

Author: Samuel Buckner (UW-ACL)
=#

# ..:: General-Purpose Functions ::..

function optimal_controller(t::CReal, T::CVector, U::CMatrix, disc::Int)::CVector
    # Output the interpolated optimal control input at time t.
    # (interpolation based on hold assumption)
    #
    # :in t: the current time
    # :in T: the time signal history
    # :in U: the input signal history
    # :out u: the interpolated input at time "t"

    i = findlast(τ->τ<=t,T)
    if typeof(i)==Nothing || i>=size(U,2)
        u = U[:,end]
    else
        u = CVector([interpolate(t,T[i:i+1],U[k,i:i+1],disc) for k=1:size(U,1)])
    end
    return u
end

function optimal_controller(t::CReal, t_span::Tuple{CReal,CReal}, U::Tuple{CVector,CVector}, disc::Int)::CVector
    # Output the interpolated optimal control input at time t.
    # (interpolation based on hold assumption)
    # Note: this variant is used for upgraded batch C2D function.
    #
    # :in t: the current time
    # :in t_span: knot point interval's time span
    # :in U: knot point interval's activated control parameters
    # :out u: the interpolated input at time "t"

    nu = length(U[1])
    if t >= t_span[2]
        u = U[end]
    elseif t <= t_span[1]
        u = U[1]
    else
        T_vec = CVector([t_span[1],t_span[2]])
        U_vec(k) = CVector([U[1][k],U[2][k]])
        u = CVector([interpolate(t,T_vec,U_vec(k),disc) for k=1:nu])
    end
    return u
end

function interpolate(x::CReal, X::CVector, Y::CVector, disc::Int)
    if disc == 0 # ZOH interpolation
        y = Y[1]
    elseif disc == 1 # FOH interpolation
        y  = Y[1] + (x - X[1])/(X[2] - X[1])*(Y[2] - Y[1])
    else
        error("Please select a valid discretization hold order.")
    end
    return y
end

function rk4_batch(f::Function, x0::CVector, t0::CReal, tf::CReal, Δt::CReal)::Tuple{CVector, CMatrix}
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
        t_cur = t[n]
        x_cur = X[:,n]
        Δt = t[n+1]-t[n]
        X[:,n+1] = rk4_step(x_cur, f, t_cur, Δt)
    end

    return (t,X)
end

function rk4_step(x_cur::CVector, f::Function, t_cur::CReal, Δt::CReal)::CVector
    """
    Integrate a system of ordinary differential equations (ODE)
    one time-step forward using RK4 (updates x_cur in place).

    Args:
        x_cur (CVector): the current state.
        f (Function): the function defining the ODE, dx/dt=f(t,x).
        t_cur (CReal): the current time (in DDTO solution).
        Δt (CReal): the integration time step.

    Returns:
        x_new (CVector): the new state.
    """

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

function simulate(sol::Solution, dyn::Function, disc::Int; max_steps::Int=40, h_min::Float64=1e-4)::Solution
    # Simulate the dynamics of the solution using a predefined control input
    # trajectory in continuous time with RK4 integration.

    dyn_ = (t,x) -> dyn(t,x,sol)
    n = size(sol.x,1)
    m = size(sol.u,1)
    T = CVector(undef,0)
    X = CMatrix(undef,n,0)
    U = CMatrix(undef,m,0)
    x0 = sol.x[:,1]
    for k = 1:(length(sol.t)-1)
        idx_cat = k == (length(sol.t)-1) ? 0 : 1
        Δt_prop = max((1/max_steps)*(sol.t[k+1] - sol.t[k]), h_min)
        T_,X_ = rk4_batch(dyn_, x0, sol.t[k], sol.t[k+1], Δt_prop)
        U_ = CMatrix(hcat([optimal_controller(T_[n],sol.t,sol.u,disc) for n = 1:length(T_)]...))
        T = vcat(T,T_[1:end-idx_cat])
        X = hcat(X,X_[:,1:end-idx_cat])
        U = hcat(U,U_[:,1:end-idx_cat])
        x0 = X_[:,end]
    end
    sim = Solution(T,X,U,sol.cost)
    return sim
end

function simulate(solution::DDTOSolution, dyn::Function, disc::Int; max_steps::Int=40, h_min::Float64=1e-4)::DDTOSolution
    # Run `simulate` for each branch of the provided solution set
    n = length(solution.targs)
    simulation = EmptyDDTOSolution(n)
    for k=1:n
        simulation.targs[k] = simulate(solution.targs[k], dyn, disc; max_steps=max_steps, h_min=h_min)
    end

    return simulation
end