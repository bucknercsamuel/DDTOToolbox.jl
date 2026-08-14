#=
General-purpose utilities for DDTOToolbox: JuMP solver setup, free-final-time
conversion between wall-clock time and time-dilation control, numerical /
SymPy Jacobian helpers, and colored console formatting.
=#

"""
    solver_setup(solver::String) -> (mdl, type)

Construct a JuMP model for the named convex solver.

# Arguments
- `solver::String`: solver name; one of `\"Clarabel\"`, `\"ECOS\"`, `\"MOSEK\"`, or `\"OSQP\"`

# Returns
- `mdl`: configured `JuMP.Model` with verbose logging disabled
- `type::String`: cone class of the solver (`\"QP\"` or `\"SOCP\"`), used to choose trust-region forms
"""
function solver_setup(solver::String)
    type = ""
    if solver == "Clarabel"
        mdl = Model(optimizer_with_attributes(Clarabel.Optimizer,
            "verbose" => 0))
            type = "QP"
    elseif solver == "ECOS"
        mdl = Model(optimizer_with_attributes(ECOS.Optimizer, 
            "verbose" => 0, 
            "max_iters" => 1000))
        type = "SOCP"
    elseif solver == "MOSEK"
        mdl = Model(Mosek.Optimizer)
        JuMP.set_optimizer_attribute(mdl, "LOG",  0) # disable debugging
        JuMP.set_optimizer_attribute(mdl, "MAX_NUM_WARNINGS", 0) # disable warnings
        type = "SOCP"
    elseif solver == "OSQP"
        # mdl = Model(OSQP.Optimizer)
        # JuMP.set_optimizer_attribute(mdl, "LOG",  0) # disable debugging
        # JuMP.set_optimizer_attribute(mdl, "MAX_NUM_WARNINGS", 0) # disable warnings
        mdl = Model(optimizer_with_attributes(OSQP.Optimizer,
            "verbose" => 0, 
            "max_iters" => 10000, 
            "eps_abs" => 1e-8,
            "eps_rel" => 1e-8))
        type = "QP"
    else
        error("solver choice is invalid.")
    end
    return mdl, type
end

"""
    time_dilation_control_to_wall_clock_time(∂t_∂τ, τ, disc) -> Vector

Convert a time-dilation control signal into cumulative wall-clock time.

# Arguments
- `∂t_∂τ::Vector`: time-dilation control samples ``∂t/∂τ`` on the dilated grid
- `τ::Vector`: dilated-time grid corresponding to `∂t_∂τ`
- `disc::Int`: hold order (`0` = ZOH, `1` = FOH)

# Returns
- `t::Vector`: cumulative wall-clock time starting at `0`
"""
function time_dilation_control_to_wall_clock_time(∂t_∂τ::Vector, τ::Vector, disc::Int)
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

"""
    wall_clock_time_to_time_dilation_control(t, τ, disc) -> Vector

Recover time-dilation control from wall-clock and dilated time grids.

# Arguments
- `t::Vector`: wall-clock time samples
- `τ::Vector`: dilated-time samples on the same nodes as `t`
- `disc::Int`: hold order (`0` = ZOH, `1` = FOH)

# Returns
- `∂t_∂τ::Vector`: time-dilation control samples compatible with `disc`
"""
function wall_clock_time_to_time_dilation_control(t::Vector, τ::Vector, disc::Int)
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
    numerical_jacobian(t_ref, x_ref, u_ref, dyn_nl; pert=1e-4) -> (A, B, z)

Approximate continuous-time Jacobians of nonlinear dynamics by central differencing.

# Arguments
- `t_ref`: reference time
- `x_ref`: reference state
- `u_ref`: reference control
- `dyn_nl`: nonlinear dynamics with signature `dyn_nl(t, x, u)`
- `pert`: finite-difference perturbation size (must be ≥ `1e-10`)

# Returns
- `A`: approximate state Jacobian ``∂f/∂x`` at the reference
- `B`: approximate control Jacobian ``∂f/∂u`` at the reference
- `z`: affine remainder such that ``f ≈ A x + B u + z`` evaluated at the reference
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

"""
    process_sympy_string(str::String) -> String

Rewrite a SymPy expression string into Julia-friendly syntax.

# Arguments
- `str::String`: raw SymPy string (may contain `{}`, `Max`, `Heaviside`)

# Returns
- Julia-compatible expression string (commas for matrices, `max`, `heaviside`)
"""
function process_sympy_string(str::String)::String
    str = replace(str, "{}" => ",") # apply commas for matrices (simple hack due to symbols() limitations)
    str = replace(str, "Max" => "max") # necessary for julia convention
    str = replace(str, "Heaviside" => "heaviside") # necessary for julia convention
    return str
end

"""
    print_sympy_partials(f, x, u)

Print nonzero symbolic partial derivatives for pasting into Jacobian evaluators.

# Arguments
- `f`: symbolic dynamics residual vector
- `x`: symbolic state vector
- `u`: symbolic control vector

# Returns
- nothing; prints ``∂f_i/∂x_j`` and ``∂f_i/∂u_j`` assignment lines to stdout
"""
function print_sympy_partials(f,x,u)
    nx,nu = length(x),length(u) 
    for i = 1:nx
        for j = 1:nx
            ∂fi_∂xj = diff(f[i],x[j])
            if ∂fi_∂xj != 0
                print(process_sympy_string("∂f_∂x[$(i),$(j)] = $(string(∂fi_∂xj))\n"))
            end
        end
    end
    for i = 1:nx
        for j = 1:nu
            ∂fi_∂uj = diff(f[i],u[j])
            if ∂fi_∂uj != 0
                print(process_sympy_string("∂f_∂u[$(i),$(j)] = $(string(∂fi_∂uj))\n"))
            end
        end
    end
end

"""
    convert_to_colored_string(value::Float64, tolerance::Float64; specifier=\"% .2e\") -> String

Format a numeric residual with ANSI color by magnitude relative to a tolerance.

# Arguments
- `value::Float64`: residual or penalty value to display
- `tolerance::Float64`: success threshold (green if `value ≤ tolerance`)
- `specifier`: printf-style format string for the numeric value

# Returns
- ANSI-colored formatted string (green / orange within ``10×`` tolerance / red)
"""
function convert_to_colored_string(value::Float64, tolerance::Float64; specifier="% .2e")
    if value <= tolerance
        COLOR = GREEN
    elseif value <= 10*tolerance # within one order of magnitude
        COLOR = ORANGE
    else
        COLOR = RED
    end
    return Printf.format(Printf.Format(COLOR * specifier * RESET), value);
end

"""
    convert_to_colored_string(string::String, success_set::Tuple{String}) -> String

Color a status string according to membership in a success set.

# Arguments
- `string::String`: status label to display (e.g. `\"Feasible\"`)
- `success_set::Tuple{String}`: labels treated as successful (colored green)

# Returns
- ANSI-colored status string (green if in `success_set`, otherwise red)
"""
function convert_to_colored_string(string::String, success_set::Tuple{String})
    if in(string,success_set)
        COLOR = GREEN
    else
        COLOR = RED
    end
    return Printf.format(Printf.Format(COLOR * "%s" * RESET), string);
end

"""
    set_fixed_tof!(params, tf)

Lock every target to the same physical time-of-flight `tf` and disable
golden-section timestep search so DDTO-CVX uses `Δt = tf / (N-1)`.
"""
function set_fixed_tof!(params, tf::Real)
    params.a.ToF_min = tf
    params.a.ToF_max = tf
    params.a.gss_cvx = false
    params.a.Δt_cvx = tf / (params.a.N - 1)
    return params
end

"""
    set_free_tof!(params, tf_min, tf_max)

Allow free final time in `[tf_min, tf_max]` for every target and enable
golden-section search over `Δt` for DDTO-CVX. Per-step bounds are aligned
to that window so SCP / lex / qcvx share the same ToF range.
"""
function set_free_tof!(params, tf_min::Real, tf_max::Real)
    tf_max > tf_min || error("set_free_tof! requires tf_max > tf_min")
    params.a.ToF_min = tf_min
    params.a.ToF_max = tf_max
    params.a.gss_cvx = true
    params.a.Δt_min = tf_min / (params.a.N - 1)
    params.a.Δt_max = tf_max / (params.a.N - 1)
    params.a.Δt_cvx = (params.a.Δt_min + params.a.Δt_max) / 2
    return params
end

"""
    set_objective_ub!(params, J_ub)

Disable ε-suboptimality and apply a fixed, reference-independent upper bound
`J_ub` on the base objective (normalized total thrust) for every target.
"""
function set_objective_ub!(params, J_ub::Real)
    params.a.use_suboptimality = false
    params.a.J_ub_targs = fill(CReal(J_ub), params.a.n_targs)
    return params
end

"""
    base_objective_ub(params, j, ref_cost) -> CReal

Right-hand side of the per-target base-objective constraint: `(1+ε)*J_ref`
when `use_suboptimality`, otherwise the fixed `J_ub_targs[j]`.
"""
function base_objective_ub(params, j::Int, ref_cost)
    if params.a.use_suboptimality
        return (1 + params.a.ϵ_targs[j]) * ref_cost
    end
    if j <= length(params.a.J_ub_targs)
        return params.a.J_ub_targs[j]
    end
    return Inf
end

"""
    constrain_base_objective!(mdl, J_expr, j, params, ref_cost)

Impose the per-target base-objective upper bound on `J_expr`.
"""
function constrain_base_objective!(mdl, J_expr, j::Int, params, ref_cost)
    J_ub = base_objective_ub(params, j, ref_cost)
    if isfinite(J_ub)
        @constraint(mdl, J_expr <= J_ub)
    end
    return nothing
end

"""
    deferrability_objective(α, t_defer) -> CReal

Scalar DDTO deferrability objective ``∑ α_j t_j / max(α)`` (larger is better).
"""
function deferrability_objective(α, t_defer)
    isempty(t_defer) && return 0.0
    αmax = maximum(α)
    αmax == 0 && return sum(t_defer)
    return sum(α[j] * t_defer[j] for j in eachindex(t_defer)) / αmax
end
