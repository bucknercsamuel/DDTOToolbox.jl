function param_update_law!(params::Quad3DoFParams)
    params.a.w_obj_sing /= params.w_obj_decay_factor
    params.a.w_obj_ddto /= params.w_obj_decay_factor
end