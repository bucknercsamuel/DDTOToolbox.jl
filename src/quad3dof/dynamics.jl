#=
3-DOF quadcopter continuous-time dynamics: linear double-integrator forms,
nonlinear (optionally drag-augmented) models, CTCS-augmented dynamics, and
linearization via SymPy or ForwardDiff.
=#

"""
    dynamics_linear_noaugment(params::Quad3DoFParams) -> (A, B, p)

Return unaugmented 3-DOF double-integrator affine dynamics (position/velocity).

# Arguments
- `params`: 3-DOF parameters supplying mass and gravity.

# Returns
- `A`: continuous-time state matrix.
- `B`: continuous-time input matrix.
- `p`: affine drift vector.
"""
function dynamics_linear_noaugment(params::Quad3DoFParams)
    return double_integrator_dynamics(dim=3, mass=params.mass, gravity=params.g)
end

"""
    dynamics_linear(params::Quad3DoFParams) -> (A, B, p)

Return 3-DOF double-integrator dynamics with one integral state augmentation
(thrust-norm integral channel).

# Arguments
- `params`: 3-DOF parameters supplying mass and gravity.

# Returns
- `A`: augmented continuous-time state matrix.
- `B`: augmented continuous-time input matrix.
- `p`: affine drift vector.
"""
function dynamics_linear(params::Quad3DoFParams)
    return double_integrator_dynamics(dim=3, mass=params.mass, gravity=params.g, augment=true, augment_dim=1)
end

"""
    dynamics_nonlinear_nondilated(t, x, u, params::Quad3DoFParams) -> Vector

Evaluate nondilated nonlinear 3-DOF dynamics including optional quadratic drag
and thrust-norm integral rate.

# Arguments
- `t`: time `[s]` (unused by the autonomous model).
- `x`: state vector (position, velocity, ∫T).
- `u`: thrust control vector.
- `params`: 3-DOF scenario parameters (drag, mass, gravity).

# Returns
- Nondilated state derivative vector including `d(∫T)/dt = ‖u‖`.
"""
function dynamics_nonlinear_nondilated(
        t,
        x,
        u,
        params::Quad3DoFParams)

    # Compute 3-DOF dynamics
    A,B,p = dynamics_linear_noaugment(params)
    f_3dof = A*x[1:6] + B*u[1:3] + p

    # Add drag term (if enabled)
    if params.drag_term_enabled
        # Heuristic: do not continue propagating drag term if the norm value has gotten unreasonably large (past maximum constrained value) to avoid integration blowup
        max_vel_mag = sqrt(max(abs(params.v_min_V),abs(params.v_max_V))^2 + params.v_max_L^2)
        v = x[4:6]
        if norm(v) <= max_vel_mag
            v_aug = vcat(zeros(3),v)
            f_3dof .+= params.C_d*params.S_A*params.ρ*norm(v)*v_aug/2
        end
    end

    # Compute additional states (thrust integral)
    ∫T = norm(u)

    # Stack function together and apply time dilation (chain rule)
    f = [f_3dof;∫T]
    
    return f
end

"""
    dynamics_nonlinear(t, x, ν, params::Quad3DoFParams) -> CVector

Time-dilated nonlinear dynamics ``\\dot{x} = s f(x,u)`` where `ν = [u; s]`.

# Arguments
- `t`: time `[s]`.
- `x`: state vector.
- `ν`: augmented control `[u; s]` with time-dilation factor `s`.
- `params`: 3-DOF scenario parameters.

# Returns
- Time-dilated state derivative `s * f(x,u)`.
"""
function dynamics_nonlinear(
    t::CReal,
    x::CVector,
    ν::CVector,
    params::Quad3DoFParams)::CVector

    f = dynamics_nonlinear_nondilated(t,x,ν[1:end-1],params)
    s = ν[end]
    z = s*f

    return z
end

"""
    dynamics_nonlinear_nondilated_ctcs(t, x, u, params, targ_idx) -> Vector

Nondilated dynamics augmented with the CTCS constraint-violation state for
target index `targ_idx` (`0` = trunk).

# Arguments
- `t`: time `[s]`.
- `x`: state vector (including CTCS violation state when enabled).
- `u`: thrust control vector.
- `params`: 3-DOF scenario parameters.
- `targ_idx`: target index selecting constraint set for CTCS evaluation.

# Returns
- Nondilated augmented derivative with CTCS integrand rate appended.
"""
function dynamics_nonlinear_nondilated_ctcs(
        t::CReal,
        x::Vector,
        u::Vector,
        params::Quad3DoFParams,
        targ_idx::Int)::Vector

    # Dynamics and CTCS state
    f_3dof = dynamics_nonlinear_nondilated(t,x,u,params)
    ξ,_,_ = prob_constraints_eval(x,u,params,targ_idx) # CTCS violation

    # Stack function together and apply time dilation (chain rule)
    f = [f_3dof;ξ]

    return f
end

"""
    dynamics_nonlinear_ctcs(t, x, ν, params, targ_idx) -> Vector

Time-dilated CTCS-augmented nonlinear dynamics with control `ν = [u; s]`.

# Arguments
- `t`: time `[s]`.
- `x`: CTCS-augmented state vector.
- `ν`: augmented control `[u; s]`.
- `params`: 3-DOF scenario parameters.
- `targ_idx`: target index for constraint evaluation.

# Returns
- Time-dilated CTCS-augmented state derivative.
"""
function dynamics_nonlinear_ctcs(
        t::CReal,
        x::Vector,
        ν::Vector,
        params::Quad3DoFParams,
        targ_idx::Int)::Vector

    u = ν[1:end-1]
    s = ν[end]
    f = dynamics_nonlinear_nondilated_ctcs(t,x,u,params,targ_idx)
    z = s*f
    return z
end

"""
    DynamicsLinearizedCTCS
    DynamicsLinearizedCTCS(params::Quad3DoFParams) -> DynamicsLinearizedCTCS

Callable linearization cache for CTCS-augmented 3-DOF dynamics. Stores a
Jacobian workspace `∂f_∂z` reused by ForwardDiff.

# Arguments
- `params`: scenario parameters defining augmented state dimension `nx` and control
  dimension `nu` used to size the Jacobian buffer.

# Returns
- Initialized cache ready for repeated ForwardDiff or SymPy linearizations.

# Fields
- `∂f_∂z`: preallocated Jacobian buffer for nondilated CTCS dynamics.
"""
@kwdef mutable struct DynamicsLinearizedCTCS
    ∂f_∂z::CMatrix
    DynamicsLinearizedCTCS(params::Quad3DoFParams) = DynamicsLinearizedCTCS(
        ∂f_∂z = Matrix{Float64}(undef, params.a.nx, params.a.nx+params.a.nu-1))
    DynamicsLinearizedCTCS(args...) = new(args...)
end

"""
    (d::DynamicsLinearizedCTCS)(t_ref, x_ref, ν_ref, params, targ_idx) -> (A, B, Σ, z)

Linearize time-dilated CTCS dynamics about `(x_ref, ν_ref)` using SymPy or
ForwardDiff according to `params.a.differentiator`.

# Arguments
- `d`: Jacobian workspace cache storing `∂f_∂z`.
- `t_ref`: reference time `[s]`.
- `x_ref`: reference state vector.
- `ν_ref`: reference augmented control `[u; s]`.
- `params`: 3-DOF scenario parameters.
- `targ_idx`: target index for CTCS constraint evaluation.

# Returns
- `A`: state Jacobian of the dilated dynamics.
- `B`: control Jacobian including the dilation channel.
- `Σ`: placeholder empty list (unused).
- `z`: affine drift term completing the linearization.
"""
function (d::DynamicsLinearizedCTCS)(
        t_ref::CReal,
        x_ref::CVector,
        ν_ref::CVector,
        params::Quad3DoFParams,
        targ_idx::Int)::Tuple{CMatrix,CMatrix,CVector,CVector}

    # Parse reference control
    u_ref = ν_ref[1:end-1]
    s_ref = ν_ref[end]

    # Obtain nondilated dynamics jacobians
    if params.a.differentiator == "sympy"
        ∂f_∂x,∂f_∂u = evaluate_jacobians_sympy(t_ref,x_ref,u_ref,params,targ_idx)
    elseif params.a.differentiator == "forwarddiff"
        nx = length(x_ref)
        nu = length(u_ref)
        fun(z) = dynamics_nonlinear_nondilated_ctcs(t_ref,z[1:nx],z[nx+1:end],params,targ_idx)
        ForwardDiff.jacobian!(d.∂f_∂z, fun,vcat(x_ref,u_ref))
        ∂f_∂x = @view d.∂f_∂z[:,1:nx]
        ∂f_∂u = @view d.∂f_∂z[:,nx+1:end]
    else
        error("Please choose a valid differentiator option")
    end

    # ∂f_∂s: Evaluate nondilated nonlinear dynamics
    ∂f_∂s = dynamics_nonlinear_ctcs(t_ref,x_ref,vcat(u_ref,1),params,targ_idx)

    # Package partials as linearized matrices
    A = s_ref*∂f_∂x
    B = Matrix([s_ref*∂f_∂u ∂f_∂s])
    Σ = []
    z = -(s_ref*∂f_∂x*x_ref + s_ref*∂f_∂u*u_ref)

    return(A,B,Σ,z)
end