#=
Per-SCP-iteration parameter update law for 3-DOF quadcopter problems
(objective weight decay).
=#

"""
    param_update_law!(params::Quad3DoFParams)

Decay single-target and DDTO objective weights by `w_obj_decay_factor` after
each SCP / PTR iteration.

# Arguments
- `params`: 3-DOF parameters whose objective weights are decayed.

# Returns
- none

# Notes
Mutates `params.a.w_obj_sing` and `params.a.w_obj_ddto`.
"""
function param_update_law!(params::Quad3DoFParams)
    params.a.w_obj_sing /= params.w_obj_decay_factor
    params.a.w_obj_ddto /= params.w_obj_decay_factor
end
