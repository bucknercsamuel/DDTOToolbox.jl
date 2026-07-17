#=
2-DOF double-integrator continuous-time dynamics, analytic linearization, and
SymPy partial-print helper for codegen.
=#

"""
    dynamics_linear_noaugment(params::DIntegrator2DoFParams) -> (A, B, p)

Return unaugmented planar double-integrator affine dynamics.

# Arguments
- `params`: 2-DOF parameters (unused beyond type dispatch).

# Returns
- `A`: continuous-time state matrix.
- `B`: continuous-time input matrix.
- `p`: affine drift vector (zero for this model).
"""
function dynamics_linear_noaugment(params::DIntegrator2DoFParams)
    return double_integrator_dynamics(dim=2)
end

"""
    dynamics_linear(params::DIntegrator2DoFParams) -> (A, B, p)

Return planar double-integrator dynamics with one integral-state augmentation.

# Arguments
- `params`: 2-DOF parameters (unused beyond type dispatch).

# Returns
- `A`: augmented continuous-time state matrix.
- `B`: augmented continuous-time input matrix.
- `p`: affine drift vector.
"""
function dynamics_linear(params::DIntegrator2DoFParams)
    return double_integrator_dynamics(dim=2, augment=true, augment_dim=1)
end

"""
    dynamics_nonlinear(t, x, ν, params::DIntegrator2DoFParams) -> CVector

Time-dilated nonlinear 2-DOF dynamics with acceleration-norm integral rate;
`ν = [u; s]`.

# Arguments
- `t`: time `[s]`.
- `x`: state vector (position, velocity, ∫a).
- `ν`: augmented control `[u; s]`.
- `params`: 2-DOF scenario parameters.

# Returns
- Time-dilated state derivative vector.
"""
function dynamics_nonlinear(
    t::CReal,
    x::CVector,
    ν::CVector,
    params::DIntegrator2DoFParams)::CVector

    # Compute 2-DOF non-dilated dynamics
    A,B,p = dynamics_linear_noaugment(params)
    u = ν[1:end-1]
    s = ν[end]
    f_2dof = A*x[1:end-1] + B*u + p
    f = [f_2dof; norm(u)]
    z = s*f # dilate dynamics w/ chain rule
    
    return z
end

"""
    dynamics_linearized(t_ref, x_ref, ν_ref, params::DIntegrator2DoFParams) -> (A, B, Σ, z)

Analytic linearization of time-dilated 2-DOF dynamics about a reference.

# Arguments
- `t_ref`: reference time `[s]`.
- `x_ref`: reference state vector.
- `ν_ref`: reference augmented control `[u; s]`.
- `params`: 2-DOF scenario parameters.

# Returns
- `A`: state Jacobian of the dilated dynamics.
- `B`: control Jacobian including the dilation channel.
- `Σ`: placeholder empty list (unused).
- `z`: affine drift term completing the linearization.
"""
function dynamics_linearized(
    t_ref::CReal,
    x_ref::CVector,
    ν_ref::CVector,
    params::DIntegrator2DoFParams)::Tuple{CMatrix,CMatrix,CVector,CVector}

    # Parse reference control
    u_ref = ν_ref[1:end-1]
    s_ref = ν_ref[end]
    u1,u2 = u_ref

    # Matrices to populate
    nx = length(x_ref)
    nu = length(u_ref)
    ∂f_∂x = zeros(nx,nx)
    ∂f_∂u = zeros(nx,nu)

    # Printed out partials
    ∂f_∂x[1,1] = 0
    ∂f_∂x[1,2] = 0
    ∂f_∂x[1,3] = 1.00000000000000
    ∂f_∂x[1,4] = 0
    ∂f_∂x[1,5] = 0
    ∂f_∂u[1,1] = 0
    ∂f_∂u[1,2] = 0
    ∂f_∂x[2,1] = 0
    ∂f_∂x[2,2] = 0
    ∂f_∂x[2,3] = 0
    ∂f_∂x[2,4] = 1.00000000000000
    ∂f_∂x[2,5] = 0
    ∂f_∂u[2,1] = 0
    ∂f_∂u[2,2] = 0
    ∂f_∂x[3,1] = 0
    ∂f_∂x[3,2] = 0
    ∂f_∂x[3,3] = 0
    ∂f_∂x[3,4] = 0
    ∂f_∂x[3,5] = 0
    ∂f_∂u[3,1] = 1.00000000000000
    ∂f_∂u[3,2] = 0
    ∂f_∂x[4,1] = 0
    ∂f_∂x[4,2] = 0
    ∂f_∂x[4,3] = 0
    ∂f_∂x[4,4] = 0
    ∂f_∂x[4,5] = 0
    ∂f_∂u[4,1] = 0
    ∂f_∂u[4,2] = 1.00000000000000
    ∂f_∂x[5,1] = 0
    ∂f_∂x[5,2] = 0
    ∂f_∂x[5,3] = 0
    ∂f_∂x[5,4] = 0
    ∂f_∂x[5,5] = 0
    ∂f_∂u[5,1] = u1/sqrt(u1^2 + u2^2)
    ∂f_∂u[5,2] = u2/sqrt(u1^2 + u2^2)

    # ∂f_∂s: Evaluate nondilated nonlinear dynamics
    ∂f_∂s = dynamics_nonlinear(t_ref,x_ref,vcat(u_ref,1),params)

    # Package partials as linearized matrices
    A = s_ref*∂f_∂x
    B = Matrix([s_ref*∂f_∂u ∂f_∂s])
    Σ = []
    z = -(s_ref*∂f_∂x*x_ref + s_ref*∂f_∂u*u_ref)

    return(A,B,Σ,z)
end

"""
    generate_dynamics_partials(params::DIntegrator2DoFParams)

Symbolically differentiate nondilated 2-DOF dynamics and print all partials
for codegen. Requires SymPy.

# Arguments
- `params`: 2-DOF parameters used to build the symbolic double-integrator model.

# Returns
- none; partial derivatives are printed to stdout.
"""
function generate_dynamics_partials(params::DIntegrator2DoFParams)

    # Symbols for differentiable quantities
    r1,r2 = symbols("r1 r2", real=true)
    v1,v2 = symbols("v1 v2", real=true)
    u1,u2 = symbols("u1 u2", real=true)
    intu  = symbols("intu"; real=true)

    # Symbol canonicalization
    x = [r1;r2;v1;v2;intu]
    u = [u1;u2]
    nx,nu = length(x),length(u) 

    # Evaluate nondilated nonlinear dynamics
    A,B,p = dynamics_linear_noaugment(params)
    f_2DoF = A*x[1:end-1] + B*u + p
    f = [f_2DoF; norm(u)]

    # Print out all partial elements
    for i = 1:nx
        for j = 1:nx
            ∂fi_∂xj = diff(f[i],x[j])
            print("∂f_∂x[$(i),$(j)] = $(string(∂fi_∂xj))\n")
        end
        for j = 1:nu
            ∂fi_∂uj = diff(f[i],u[j])
            print("∂f_∂u[$(i),$(j)] = $(string(∂fi_∂uj))\n")
        end
    end
end