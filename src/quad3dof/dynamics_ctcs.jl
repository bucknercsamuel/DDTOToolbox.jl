function dynamics_nonlinear_nondilated_ctcs(
    t::CReal,
    x::CVector,
    u::CVector,
    params::Quad3DoFParams,
    targ_idx::Int)::CVector

    # Dynamics and CTCS state
    f_3dof = dynamics_nonlinear_nondilated(t,x,u,params)
    ξ,_,_ = prob_constraints_eval(x,u,params,targ_idx) # CTCS violation

    # Stack function together and apply time dilation (chain rule)
    f = [f_3dof;ξ]

    return f
end

function dynamics_nonlinear_ctcs(
    t::CReal,
    x::CVector,
    ν::CVector,
    params::Quad3DoFParams,
    targ_idx::Int)::CVector
    u = ν[1:end-1]
    s = ν[end]
    f = dynamics_nonlinear_nondilated_ctcs(t,x,u,params,targ_idx)
    z = s*f
    return z
end