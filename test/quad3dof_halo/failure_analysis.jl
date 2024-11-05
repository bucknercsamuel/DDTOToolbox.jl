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
    ϵ_subopt = 0.01
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
    params.n_targs_min = 2
    params.n_targs_max = 7

    # Resize all parameters that depend on number of max targets
    reallocate_targ_dims!(params)

    # Set boundary conditions for single target
    params.a.z0 = [
        -2.23346529721832;
        -0.6415195976948731;
        143.83925087892652;
        0.24608853589483431;
        -0.09865393164977841;
        -4.947810039372024;
        0.0;
        Inf
    ]
    params.a.u0 = [
        3.478020317230397;
        0.44693325387338045;
        9.498009359831745;
        Inf;
        Inf
    ]
    params.a.zf_targs = [
        121.314    16.8385  -63.15    -22.4603   19.2965  -140.803   -118.133;
        31.8854  -97.8167   84.9552  -97.8337  -11.7335     1.3898   -61.3404;
         1.0       1.0       1.0       1.0       1.0        1.0        1.0;
         0.0       0.0       0.0       0.0       0.0        0.0        0.0;
         0.0       0.0       0.0       0.0       0.0        0.0        0.0;
         0.0       0.0       0.0       0.0       0.0        0.0        0.0;
        Inf       Inf       Inf       Inf       Inf        Inf        Inf
    ]

    return params
end

# Test
params = FailedParams()
solve(params)
;