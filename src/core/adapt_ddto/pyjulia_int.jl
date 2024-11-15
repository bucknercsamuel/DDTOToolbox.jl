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

    params.a.n_targs = params.n_targs_max
    params.a.zf_targs = vcat(zeros(params.a.nx-1, params.a.n_targs), Inf * ones(1, params.a.n_targs))
    params.a.uf_targs = Inf * ones(params.a.nu,params.a.n_targs)
    params.a.λ_targs = Vector{Int}(undef, params.a.n_targs)
    params.a.ID_targs = Vector{Int}(undef, params.a.n_targs)
    params.a.J_targs = 1:params.a.n_targs
    params.a.τ_targs = zeros(params.a.n_targs)
    params.a.α_targs = ones(params.a.n_targs)
    params.a.ϵ_targs = fill(params.ϵ_subopt, params.a.n_targs)
    params.a.w_obj_ddto = params.a.w_obj_sing / params.a.n_targs
    params.R_targs = CVector(undef, params.a.n_targs)
    for (key,~) in params.p_targs
        params.p_targs[key] = CVector(undef, params.a.n_targs)
    end
end

function sort_des_score!(params)
    """
    Sort the targets by descending desirability score.
    * NOTE: This function will modify the params object.

    Args:
        params (any): the params object.
    """

    des_score = zeros(params.a.n_targs)
    for j = 1 : params.a.n_targs
        des_score[j] = 
            params.p_targs["pcd"][j] * params.w_des[1] + 
            params.p_targs["prox_veh"][j] * params.w_des[2] + 
            params.p_targs["prox_clust"][j] * params.w_des[3] + 
            params.p_targs["µ_99"][j] * params.w_des[4] + 
            params.R_targs[j] * params.w_des[5]
    end
    params.a.λ_targs = sortperm(des_score)
end

function remove_infeasible_targets!(params; pre_compute::Bool=false)
    """
    Remove targets that intersect with obstacles laterally
    * NOTE: This function will modify the params object.

    Args:
        params (any): the params object.
    """
    zf_targs = params.a.zf_targs
    p_obstacles = params.p_obstacles
    R_targs = params.R_targs
    R_obstacles = params.R_obstacles
    J_targs = params.a.J_targs
    for j = 1 : params.a.n_targs
        for k = 1 : params.n_obstacles
            if norm(zf_targs[1:2,j] - p_obstacles[1:2,k]) < (R_targs[j] + R_obstacles[k])
                remove_ddto_target!(params, J_targs[j])
                break
            end
        end
    end
    if pre_compute
        params.a.J_targs = Vector(1:params.a.n_targs)
        sort_des_score!(params)
    end
end

function configure_greedy!(params)
    params.n_targs_min = 1
    params.n_targs_max = 1
end

function save_results(path, results)
    save(path, Dict("data"=>results))
end