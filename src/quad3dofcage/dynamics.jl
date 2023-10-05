function dyn_lin(
    t_ref::CReal,
    x_ref::CVector,
    ν_ref::CVector,
    params::Quad3DoFCageParams)::Tuple{CMatrix,CMatrix,CVector}

    # Obtain original affine system matrices
    A,B,p = params.A_c,params.B_c,params.p_c
    n,m = size(B)

    # Obtain reference time dilation factor
    u_ref = ν_ref[1:end-1]
    s_ref = ν_ref[end]

    # Augment A,B to account for affine term
    A_aff = CMatrix([
        A p;
        zeros(1,n+1)
    ])
    B_aff = CMatrix([
        B;
        zeros(1,m)
    ])
    
    # Compute linearization derivatives
    df_dx = s_ref * A_aff
    df_du = s_ref * B_aff
    df_ds = A_aff * x_ref + B_aff * u_ref

    # Compute linearized A,B,w matrices
    A_ = df_dx
    B_ = CMatrix([df_du df_ds])
    w_ = -(s_ref * B_aff * u_ref + s_ref * A_aff * x_ref)

    return(A_,B_,w_)
end

function dyn_nl(
    t::CReal,
    x::CVector,
    ν::CVector,
    params::Quad3DoFCageParams)::CVector

    # Obtain original affine system matrices
    A,B,p = params.A_c,params.B_c,params.p_c
    n,m = size(B)

    # Augment A,B to account for affine term
    A_aff = CMatrix([
        A p;
        zeros(1,n+1)
    ])
    B_aff = CMatrix([
        B;
        zeros(1,m)
    ])

    u = ν[1:end-1]
    s = ν[end]
    z = s*A_aff*x + s*B_aff*u
    
    return z
end