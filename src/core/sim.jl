#=
Trajectory simulation utilities: control interpolation under ZOH/FOH, RK4
integration, and forward simulation of single-target and DDTO solutions.
=#

# ..:: General-Purpose Functions ::..

"""
    optimal_controller(t, T, U, disc) -> CVector

Interpolate an open-loop control history at a query time.

# Arguments
- `t::CReal`: query time at which to evaluate the control
- `T::CVector`: time grid associated with the control history
- `U::CMatrix`: control history with one column per time sample
- `disc::Int`: hold order (`0` = zero-order hold, `1` = first-order hold)

# Returns
- `u::CVector`: interpolated control at time `t` (clamped to the final sample if `t` is past the grid)
"""
function optimal_controller(t::CReal, T::CVector, U::CMatrix, disc::Int)::CVector
    i = findlast(τ->τ<=t,T)
    if typeof(i)==Nothing || i>=size(U,2)
        u = U[:,end]
    else
        u = CVector([interpolate(t,T[i:i+1],U[k,i:i+1],disc) for k=1:size(U,1)])
    end
    return u
end

"""
    optimal_controller(t, t_span, U, disc) -> CVector

Interpolate control over a single knot interval. Used by batch continuous-to-discrete transcription.

# Arguments
- `t::CReal`: query time within (or outside) the knot interval
- `t_span::Tuple{CReal,CReal}`: interval endpoints `(t⁻, t⁺)`
- `U::Tuple{CVector,CVector}`: control values `(u⁻, u⁺)` at the interval endpoints
- `disc::Int`: hold order (`0` = ZOH, `1` = FOH)

# Returns
- `u::CVector`: interpolated control at time `t`
"""
function optimal_controller(t::CReal, t_span::Tuple{CReal,CReal}, U::Tuple{CVector,CVector}, disc::Int)::CVector
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

"""
    interpolate(x, X, Y, disc) -> Number

Scalar interpolation of samples under ZOH or FOH.

# Arguments
- `x::CReal`: query abscissa
- `X::CVector`: sample abscissae (length 2 for FOH)
- `Y::CVector`: sample ordinates matching `X`
- `disc::Int`: hold order (`0` = ZOH uses `Y[1]`, `1` = linear FOH)

# Returns
- interpolated scalar value at `x`
"""
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

"""
    rk4_batch(f, x0, t0, tf, Δt) -> (t, X)

Integrate ``\\dot{x} = f(t,x)`` from `t0` to `tf` with classical RK4.

# Arguments
- `f::Function`: ODE right-hand side with signature `f(t, x)`
- `x0::CVector`: initial state
- `t0::CReal`: integration start time
- `tf::CReal`: integration final time
- `Δt::CReal`: nominal integration step (final step may be shorter to hit `tf`)

# Returns
- `t::CVector`: integration time grid
- `X::CMatrix`: state history with one column per time sample
"""
function rk4_batch(f::Function, x0::CVector, t0::CReal, tf::CReal, Δt::CReal)::Tuple{CVector, CMatrix}
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

"""
    rk4_step(x_cur, f, t_cur, Δt) -> CVector

Advance ``\\dot{x} = f(t,x)`` one classical RK4 step.

# Arguments
- `x_cur::CVector`: state at the beginning of the step
- `f::Function`: ODE right-hand side with signature `f(t, x)`
- `t_cur::CReal`: time at the beginning of the step
- `Δt::CReal`: step size

# Returns
- `x_new::CVector`: state after advancing by `Δt`
"""
function rk4_step(x_cur::CVector, f::Function, t_cur::CReal, Δt::CReal)::CVector
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

"""
    simulate(sol::Solution, dyn, disc; max_steps=40, h_min=1e-4) -> Solution

Forward-simulate a single-target open-loop solution with RK4 on each nodal interval.

# Arguments
- `sol::Solution`: reference trajectory whose control is applied open-loop
- `dyn::Function`: dynamics with signature `dyn(t, x, sol)`
- `disc::Int`: hold order used to interpolate `sol.u`
- `max_steps::Int`: maximum RK4 steps per nodal interval
- `h_min::Float64`: minimum integration step size

# Returns
- `sim::Solution`: densely sampled simulated trajectory (same cost as `sol`)
"""
function simulate(sol::Solution, dyn::Function, disc::Int; max_steps::Int=40, h_min::Float64=1e-4)::Solution
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

"""
    simulate(solution::DDTOSolution, dyn, disc; max_steps=40, h_min=1e-4) -> DDTOSolution

Forward-simulate every branch of a DDTO solution bundle.

# Arguments
- `solution::DDTOSolution`: multi-target solution whose branches are simulated independently
- `dyn::Function`: dynamics with signature `dyn(t, x, sol)` for each branch
- `disc::Int`: hold order used to interpolate each branch control
- `max_steps::Int`: maximum RK4 steps per nodal interval
- `h_min::Float64`: minimum integration step size

# Returns
- `simulation::DDTOSolution`: bundle of simulated branch trajectories
"""
function simulate(solution::DDTOSolution, dyn::Function, disc::Int; max_steps::Int=40, h_min::Float64=1e-4)::DDTOSolution
    n = length(solution.targs)
    simulation = EmptyDDTOSolution(n)
    for k=1:n
        simulation.targs[k] = simulate(solution.targs[k], dyn, disc; max_steps=max_steps, h_min=h_min)
    end

    return simulation
end
