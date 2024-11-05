"""
In this script, I investigate a problem scenario that is known to produce a non-convergent SCvx solution (single-target).
This is investigated for understanding fundamental algorithm improvements.
"""

using DDTOSCP
include("plots.jl")

# Obtain default params that were used for this scenario
function DefaultParams()::Quad3DoFHaloParams{CReal,Int}
    # >> Environmental parameters <<
    g = -9.807*e_z
    ρ = 1.225

    # Rotor Parameters (not stored)
    # (See line 21 of https://github.com/microsoft/AirSim/blob/master/AirLib/include/vehicles/multirotor/RotorParams.hpp)
    C_T = 0.109919
    RPM_max = 6396.667
    d_prop = 0.2286

    # >> Vehicle parameters <<
    C_D = 1.3/4.0
    S_A = .18*.11 # overhead rectangular area assuming vehicle's velocity is mostly aligned with body -Z, not including arms
    n_rotor = 4
    mass = 1
    ρ_min = 5                                             # [N] AirSim throttle lower bound default
    ρ_max = n_rotor * C_T * ρ * (RPM_max/60)^2 * d_prop^4 # [N] Max physical thrust of single engine
    drag_term_enabled = true

    # >> Constraint parameters <<
    ϵ_subopt = 0
    γ_gs = 89 * DEG_2_RAD
    γ_p = 89 * DEG_2_RAD
    v_max_V = 1e-3
    v_min_V = -5
    v_max_L = 5

    # Obstacle and boundary parameters 
    # (defaults to empty, scenario-specific)
    n_obstacles = 0
    R_obstacles = CVector(undef,0)
    p_obstacles = CMatrix(undef,0,0)
    H_obstacles = Vector{CMatrix}(undef,0)

    # >> Algorithm parameters <<
    a = AlgorithmParams()
    a.nx = 7 # (position, velocity, thrust 2-norm)
    a.nu = 3 # (thrust)
    a.z0 = Inf * ones(a.nx) # empty initial state (to be populated with current state)
    a.u0 = Inf * ones(a.nu) # empty initial control (to be populated with current control)
    w_obj_decay_factor = 1.2

    # SCP parameters
    a.ctcs_enabled = true
    a.warmstart_method = "single" # types: (linear, single, ddto)
    a.w_obj_sing = .01
    a.w_ctrl = 50.
    a.w_trust = 10.
    a.w_buff = a.w_ctrl
    a.ϵ_ctrl = 5e-3
    a.ϵ_buff = 5e-3
    a.ϵ_trust = 5e-3
    a.scp_iters = 50

    # Time dilation & discretization
    a.N = 20
    a.N_msi = 10
    a.ToF_min = 10.
    a.ToF_max = 120.
    a.Δt_min = .5*a.ToF_min/(a.N-1)
    a.Δt_max =  2*a.ToF_max/(a.N-1)
    a.gss_cvx = true
    a.Δt_cvx = (a.Δt_min + a.Δt_max)/2.
    a.differentiator = "forwarddiff"

    # >> HALO-specific parameters <<
    n_targs_min = 2
    n_targs_max = 7
    R_targs_min = 1.
    R_targs = CVector(undef, a.n_targs)
    p_targs = Dict(
        "pcd" => CVector(undef, a.n_targs),         # Point cloud density
        "prox_veh" => CVector(undef, a.n_targs),    # Proximity of landing site to vehicle
        "prox_clust" => CVector(undef, a.n_targs),  # Proximity to other landing sites ("cluster proximity")
        "µ_99" => CVector(undef, a.n_targs),        # 99th percentile uncertainty
    )
    w_des = [0,0,1,0,0] # Weights for: [pcd, prox_veh, prox_clust, µ_99, R_targs]

    params = Quad3DoFHaloParams{CReal,Int}(
        g,
        ρ,
        C_D,
        S_A,
        n_rotor,
        mass,
        ρ_min,
        ρ_max,
        drag_term_enabled,
        ϵ_subopt,
        γ_gs,
        γ_p,
        v_min_V,
        v_max_V,
        v_max_L,
        n_obstacles,
        R_obstacles,
        p_obstacles,
        H_obstacles,
        n_targs_min,
        n_targs_max,
        R_targs,
        R_targs_min,
        p_targs,
        w_des,
        a,
        w_obj_decay_factor,
    )

    return params
end

# Set scenario-specific data that is known to fail (will fail on 4th iteration of first target)
function FailedParams()::Quad3DoFHaloParams{CReal,Int}
    # Load default params
    params = DefaultParams()

    # Configure for single-target problem
    params.n_targs_min = 1
    params.n_targs_max = 1

    # Resize all parameters that depend on number of max targets
    reallocate_targ_dims!(params)

    # Set boundary conditions for single target
    params.a.z0 = [
        -2.041531619726133;
        -14.72081896693884;
        108.31963744244631;
        -0.5931192362368635;
        -3.3885055525279126;
        -4.847433846477435;
        0.0;
        Inf
    ]
    params.a.u0 = [
        0.6833988245278033;
        -1.4231318020669774;
        10.412767234739698;
        Inf
    ]
    params.a.zf_targs = reshape([
        35.740962162615084;
        135.9817457589394;
        1.0;
        0.0;
        0.0;
        0.0;
        Inf
    ], params.a.nx,1)

    return params
end

# Test
params = FailedParams()
solve(params)
;