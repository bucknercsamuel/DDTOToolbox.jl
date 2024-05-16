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