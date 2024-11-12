using DDTOSCP
using Random
using Printf
using Debugger
using JLD2

function simulate_halo_landing(
        quad,              # Quad object
        r0,                # [m] Initial position (NED frame)
        v0,                # [m/s] Initial velocity (NED frame)
        dynamics;          # Dynamics function
        Δt_sim    = 0.01,  # [s] Simulation integration time-step
        Δt_print  = 1.,    # [s] Simulation printing update time-step
        R_ROI     = 150.,   # [m] Radius of the region of interest for targets
        h_cut     = 50.,   # [m] Altitude condition to commit to best target
        h_term    = 1.,    # [m] Altitude condition to terminate descent phase
        h_eps     = 1.,    # [m] Acceptable altitude error in termination condition
        greedy    = false, # Select if we should use greedy method instead of DDTO
        greedy_dt = 5,      # Greedy update timestep
        n_target_pool = 10 # Number of targets in the global pool
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
    target_pool = sim_build_target_pool(n_target_pool, R_ROI, min_radius=quad.R_targs_min, max_radius=5*quad.R_targs_min)

    # Simulation status
    sim_cur_iter    = 0
    sim_cur_time    = 0.0
    sim_cur_state   = quad.a.z0
    sim_cur_control = init_thrust
    sim_num_ddto    = 0

    # Other variables
    guid,flags,results = setup_addto_dicts(quad)
    save_param_checkpts = false
    time_last_print = 0.0
    t_fine = nothing
    u_fine = nothing
    τ_fine = nothing

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
            τ_fine = CVector(range(start=guid["cur_traj"].τ[1],stop=guid["cur_traj"].τ[end],length=1001))
            u_fine = CMatrix(hcat([optimal_controller(τ_fine[n],guid["cur_traj"].τ,guid["cur_traj"].u,quad.a.disc) for n = 1:length(τ_fine)]...))
            t_fine = CVector(time_dilation_control_to_wall_clock_time(u_fine[end,:], τ_fine, quad.a.disc))
        end
        if !flags["guid_lock_activated"]
            check_unsafe_targets!(quad, guid, flags, sim_cur_time)
            check_branch_switch!(quad, guid, flags, sim_cur_state, sim_cur_time)
            check_cutoff_altitude!(sim_cur_state, sim_cur_time, h_cut, flags)
        end
        sim_update_targets!(quad, target_pool)
        if flags["guid_lock_staged"]
            activate_guidance_lock!(quad, guid, flags, sim_cur_time)
            τ_fine = CVector(range(start=guid["cur_traj"].τ[1],stop=guid["cur_traj"].τ[end],length=1001))
            u_fine = CMatrix(hcat([optimal_controller(τ_fine[n],guid["cur_traj"].τ,guid["cur_traj"].u,quad.a.disc) for n = 1:length(τ_fine)]...))
            t_fine = CVector(time_dilation_control_to_wall_clock_time(u_fine[end,:], τ_fine, quad.a.disc))
        end

        # Log results
        sim_cur_control = optimal_controller(guid["cur_time"], t_fine, u_fine, quad.a.disc)
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
            @printf("   ﹂ UPDATE [%.2f s]: Terminal altitude condition reached -- landing successful!\n", sim_cur_time)
        end
        if sim_cur_time >= quad.a.ToF_max
            display("Simulation ran for too long, exiting...")
            flags["descent_complete"] = true
        end

        # Integrate the currently-tracked guidance trajectory for one time-step
        # and go to next iteration
        cur_ct_dyn = (t,x) -> dynamics(t,x,t_fine,u_fine,quad)
        sim_cur_state = rk4_step(sim_cur_state, cur_ct_dyn, guid["cur_time"], Δt_sim)
        sim_cur_time     += Δt_sim
        guid["cur_time"] += Δt_sim
    end

    return results
end