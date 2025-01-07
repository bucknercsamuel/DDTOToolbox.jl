using DDTOSCP
using Random
using Printf
using Debugger
using JLD2

function simulate_halo_landing(
        quad,               # Quad object
        r0,                 # [m] Initial position (NED frame)
        v0;                 # [m/s] Initial velocity (NED frame)
        Δt_sim    = 0.01,   # [s] Simulation integration time-step
        Δt_print  = 1.,     # [s] Simulation printing update time-step
        R_ROI     = 150.,   # [m] Radius of the region of interest for targets
        h_cut     = 50.,    # [m] Altitude condition to commit to best target
        h_term    = 1.,     # [m] Altitude condition to terminate descent phase
        h_eps     = 1.,     # [m] Acceptable altitude error in termination condition
        greedy    = false,  # Select if we should use greedy method instead of DDTO
        greedy_dt = 5,      # Greedy update timestep
        n_target_pool = 10, # Number of targets in the global pool
        n_obs = 0,          # Number of obstacles (default to 0)
    )

    # Modifications if using greedy single-target method
    if greedy
        quad.n_targs_min = 1
        quad.n_targs_max = 1
    end

    # Override the quadcopter's initial conditions
    quad.a.z0 = vcat(r0,v0,0)

    # Initialize thrust control to vertical at hover condition
    init_thrust = -quad.mass*quad.g

    # Build the target pool
    target_pool = sim_build_target_pool(n_target_pool, R_ROI, min_radius=quad.R_targs_min)

    # Build the obstacle pool
    # obs_rad_position = R_ROI/2
    # obs_rad = 5.
    # generate_obstacles!(quad, n_obs, obs_rad_position, obs_rad)
    generate_obstacles!(quad, n_obs, (3,12), (-R_ROI,R_ROI), (-R_ROI,R_ROI), 0)

    # Simulation status
    sim_cur_iter    = 0
    sim_cur_time    = 0.0
    sim_cur_state   = quad.a.z0
    sim_cur_control = init_thrust
    sim_num_ddto    = 0
    error_code      = 0

    # Other variables
    guid,flags,results = setup_addto_dicts(quad)
    save_param_checkpts = false
    time_last_print = 0.0

    # Initial print statements
    println("=== Beginning Simulation ===")
    @printf("Time: %.2f s, Alt: %.2f m, Number of targets: %i\n", sim_cur_time, sim_cur_state[3], quad.n_targs_max)

    # ..:: Main Sim Loop ::..
    max_iter = 1e6
    while !flags["descent_complete"]

        # Greedy update if not using Adaptive-DDTO
        if greedy && guid["cur_time"] > greedy_dt && !flags["guid_lock_activated"]
            flags["update_ddto"] = true
            remove_ddto_target!(quad,1)
        end

        # Execute Adaptive-DDTO algorithm pipeline
        if flags["update_ddto"] && !flags["guid_lock_activated"] && !flags["guid_lock_staged"]
            sim_refresh_targets!(quad, target_pool)
            if save_param_checkpts
                save("quad3dof_halo/tmp/params.jld","params",quad)
            end
            compute_ddto_guidance!(quad, guid, flags, sim_cur_state, sim_cur_control, sim_cur_time)
            sim_num_ddto += 1
        end
        if !flags["guid_lock_activated"]
            check_unsafe_targets!(quad, guid, flags, sim_cur_time) # checks if there are too many unsafe targets in the current target pool
            check_branch_switch!(quad, guid, flags, sim_cur_state, sim_cur_time) # determines if we should switch to deferred target or remain on deferred segment
            if !flags["guid_lock_staged"]
                check_cutoff_altitude!(sim_cur_time, sim_cur_state[3], h_cut, flags) # stages guidance lock if conditions are met
            end
        end
        sim_update_targets!(quad, target_pool)
        if flags["guid_lock_staged"]
            activate_guidance_lock!(quad, guid, flags, sim_cur_time)
        end
        
        # Log results
        sim_cur_control = optimal_controller(guid["cur_time"], guid["cur_traj_sim"].t, guid["cur_traj_sim"].u, quad.a.disc)
        log_results!(quad, results, guid, flags, sim_cur_state, sim_cur_control, sim_cur_time, target_pool=target_pool)

        # Print sim status update(s)
        if (sim_cur_time - time_last_print) >= Δt_print
            defer_time_remaining = guid["defer_time"] - guid["cur_time"]
            if !flags["guid_lock_activated"]
                @printf("Time: %.2f s, Alt: %.2f m, Number of targets: %i, Next deferred target: %i (%.2f s remaining)\n", sim_cur_time, sim_cur_state[3], quad.a.n_targs, guid["defer_targ"], defer_time_remaining)
            else
                @printf("Time: %.2f s, Alt: %.2f m, Guidance locked to target %i!\n", sim_cur_time, sim_cur_state[3], guid["defer_targ"])
            end
            # @printf("   ﹂ Debug: Applying control [%.2f, %.2f, %.2f]\n", sim_cur_control[1], sim_cur_control[2], sim_cur_control[3])
            time_last_print = sim_cur_time
        end

        # Terminate sim if we reach the phase completion condition
        if sim_cur_state[3] <= (h_term + h_eps)
            flags["descent_complete"] = true
            error_code = 1
            @printf("   ﹂ UPDATE [%.2f s]: Terminal altitude condition reached -- landing successful!\n", sim_cur_time)
        end
        if sim_cur_time >= quad.a.ToF_max
            display("Simulation ran for too long, exiting...")
            error_code = 2
            flags["descent_complete"] = true
        end

        # Integrate the currently-tracked guidance trajectory for one time-step
        sim_cur_time, sim_cur_state = step_halo_sim(quad, sim_cur_time, sim_cur_state, guid, Δt_sim=Δt_sim)
    end

    return results,error_code
end