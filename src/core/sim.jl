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

function simulate(sol::Solution, dyn::Function, disc::Int; max_steps::Int=40)::Solution
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

function simulate(solution::DDTOSolution, dyn::Function, disc::Int)::DDTOSolution
    # Run `simulate` for each branch of the provided solution set
    n = length(solution.targs)
    simulation = EmptyDDTOSolution(n)
    for k=1:n
        simulation.targs[k] = simulate(solution.targs[k], dyn, disc)
    end

    return simulation
end

function time_dilation_control_to_wall_clock_time(∂t_∂τ::Vector, τ::Vector, disc::Int)
    # Converts time-dilation control to wall-clock time based on discretization method
    Δτ = diff(τ)
    if length(∂t_∂τ) > 1
        Δt = Vector(undef,length(∂t_∂τ)-1)
        for k=1:length(Δt)
            if disc == 0
                Δt[k] = Δτ[k] * ∂t_∂τ[k]
            elseif disc == 1
                Δt[k] = (1/2) * Δτ[k] * (∂t_∂τ[k] + ∂t_∂τ[k+1])
            end
        end
        t = cumsum([0.;Δt])
    else
        t = [0.]
    end
    return t
end

function wall_clock_time_to_time_dilation_control(t::Vector, τ::Vector, disc::Int)
    # Converts wall-clock time to time dilation control based on discretization method
    N = length(t)
    if disc == 0
        N_ctrl = N-1
    elseif disc == 1
        N_ctrl = N
    end
    Δt = diff(t)
    Δτ = diff(τ)
    ∂t_∂τ = Vector(undef,N_ctrl)
    if disc == 0
        for k=1:N_ctrl
            ∂t_∂τ[k] = Δt[k] / Δτ[k]
        end
    elseif disc == 1
        n = length(Δt)
        ∂t_∂τ[1] = sum([Δt[k] / Δτ[k] for k=1:n])/length(Δt) # boundary condition chosen for numerical properties, but is technically arbitrary!
        for k=1:N_ctrl-1
            ∂t_∂τ[k+1] = 2 * Δt[k] / Δτ[k] - ∂t_∂τ[k]
        end
    end
    return ∂t_∂τ
end

"""
Compute Jacobians with direct numerical differentiation (simple central differencing) on nonlinear dynamics
"""
function numerical_jacobian(t_ref, x_ref, u_ref, dyn_nl; pert=1e-4)
    # Setup
    if pert < 1e-10
        error("Required perturbation is too small")
    end
    nx = length(x_ref)
    nu = length(u_ref)
    A = zeros(nx,nx)
    B = zeros(nx,nu)

    # Numerical A
    pertI = pert*I(nx)
    for k=1:nx
        fp = dyn_nl(t_ref, x_ref + pertI[:,k], u_ref)
        fm = dyn_nl(t_ref, x_ref - pertI[:,k], u_ref)
        A[:,k] = (fp-fm) / (2*pert)
    end

    # Numerical B
    pertI = pert*I(nu)
    for k=1:nu
        fp = dyn_nl(t_ref, x_ref, u_ref + pertI[:,k])
        fm = dyn_nl(t_ref, x_ref, u_ref - pertI[:,k])
        B[:,k] = (fp-fm) / (2*pert)
    end

    # Numerical z
    z = dyn_nl(t_ref, x_ref, u_ref) - (A*x_ref + B*u_ref)

    return A,B,z
end