#= Adaptive-DDTO -- Utility functions for PyJulia interface (only used for AirSim code).
Author: Samuel Buckner (UW-ACL)
=#

function rk4_step_pyjulia(params, guid_traj::Solution, x_cur::CVector, t_cur::CReal, Δt::CReal)::CVector
    """
    Necessary wrapper for rk4_step where the dynamics function is created here
    (since there is no way to pass an anonymous function from Python to Julia AFAIK)

    Args:
        params (any): The params object.
        guid_traj (Solution): The guidance trajectory.
        x_cur (CVector): Current state.
        t_cur (CReal): Current time.
        Δt (CReal): Time step.

    Returns:
        x_cur (CVector): Updated state.
    """

    # Get current dynamics with controller tracking `guid_traj`
    ct_dyn = (t,x) -> params.A_c*x + params.B_c*optimal_controller(t,x,guid_traj) + params.p_c # Continuous-time dynamics

    # Call `rk4_step`
    x_cur = rk4_step(x_cur, ct_dyn, t_cur, Δt)

    return x_cur
end

function reallocate_targ_dims!(params)
    """
    Reallocate the target dimensions to the maximum number of targets.
    * NOTE: This function will modify the params object.

    Args:
        params (any): the params object.
    """

    params.n_targs = params.n_targs_max
    params.R_targs = CVector(undef, params.n_targs)
    params.rf_targs = CMatrix(undef, 3, params.n_targs)
    params.vf_targs = zeros(3, params.n_targs)
    params.N_targs = Vector{Int}(undef, params.n_targs)
    params.λ_targs = Vector{Int}(undef, params.n_targs)
    params.T_targs = 1:params.n_targs
    params.ϵ_targs = CVector(undef, params.n_targs)
    for (key,~) in params.p_targs
        params.p_targs[key] = CVector(undef, params.n_targs)
    end
end

function sort_des_score!(params)
    """
    Sort the targets by descending desirability score.
    * NOTE: This function will modify the params object.

    Args:
        params (any): the params object.
    """

    des_score = zeros(params.n_targs)
    for j = 1 : params.n_targs
        des_score[j] = 
            params.p_targs["pcd"][j] * params.w_des[1] + 
            params.p_targs["prox_veh"][j] * params.w_des[2] + 
            params.p_targs["prox_clust"][j] * params.w_des[3] + 
            params.p_targs["µ_99"][j] * params.w_des[4] + 
            params.R_targs[j] * params.w_des[5]
    end
    params.λ_targs = sortperm(des_score)

    # Sort the last two targets to always be increasing since there is no rejection preference between the two
    params.λ_targs[end-1:end] = sort(params.λ_targs[end-1:end])
end