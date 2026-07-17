#=
Generic continuous-time double-integrator dynamics helpers shared by scenario
modules (3-DOF quadcopter, 2-DOF double integrator, etc.).
=#

"""
    double_integrator_dynamics(; dim=3, mass=1, gravity=nothing, augment=false, augment_dim=1) -> (A, B, p)

Build continuous-time affine double-integrator state-space matrices
``\\dot{x} = A x + B u + p``.

# Keyword Arguments
- `dim::Int`: spatial dimension of the double integrator (typically `2` or `3`)
- `mass`: vehicle mass used to scale the control matrix `B`
- `gravity`: gravity acceleration in ``\\mathbb{R}^{dim}``; defaults to zeros if `nothing`
- `augment::Bool`: if `true`, append `augment_dim` integrator states to `(A, B, p)`
- `augment_dim::Int`: number of appended integral states when `augment` is `true`

# Returns
- `A`: continuous-time state matrix
- `B`: continuous-time control matrix
- `p`: continuous-time affine forcing vector (gravity and zeros for integral channels)
"""
function double_integrator_dynamics(;
        dim = 3,
        mass = 1,
        gravity = nothing,
        augment = false,
        augment_dim = 1,
    )
    if isnothing(gravity)
        gravity = zeros(dim)
    end

    A = Matrix([
        zeros(dim,dim) I(dim);
        zeros(dim,dim) zeros(dim,dim)
    ])
    B = Matrix([
        zeros(dim,dim);
        I(dim)/mass
    ])
    p = Vector(vcat(zeros(dim),gravity))

    if augment
        A = Matrix([
            A zeros(2*dim,augment_dim);
            zeros(augment_dim,2*dim) I(augment_dim)
        ])
        B = Matrix([
            B;
            zeros(augment_dim,dim)
        ])
        p = vcat(p,zeros(augment_dim))
    end

    return A,B,p
end
