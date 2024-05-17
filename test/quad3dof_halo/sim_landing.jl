using DDTOSCP
using Random
using Printf

function simulate_halo_landing(
        r0,              # [m] Initial position (NED frame)
        v0;              # [m/s] Initial velocity (NED frame)
        Δt_sim   = 0.01, # [s] Simulation integration time-step
        Δt_print = 1.,   # [s] Simulation printing update time-step
        R_ROI    = 50.,  # [m] Radius of the region of interest for targets
        h_cut    = 50.,  # [m] Altitude condition to commit to best target
        h_term   = 1.,   # [m] Altitude condition to terminate descent phase
        h_eps    = 1.    # [m] Acceptable altitude error in termination condition
    )

    # Override the quadcopter's initial conditions
    quad.a.z0 = vcat(r0,v0,0)

    # Obtain N_max landing targets from perception stack
    sim_acquire_new_targets!(quad, R_ROI)

    # Simulation status
    sim_cur_iter    = 0
    sim_cur_time    = 0.0
    sim_cur_state   = quad.a.z0
    sim_cur_control = zeros(quad.a.nu+1)
    sim_num_ddto    = 0

    # Guidance
    guid_cur_opt    = EmptyQuad3DoFDDTOSolution(quad.a.n_targs) # Most recently-computed optimal solution set
    guid_cur_ddto   = EmptyQuad3DoFDDTOSolution(quad.a.n_targs) # Most recently-computed DDTO solution set
    guid_cur_traj   = EmptyQuad3DoFSolution() # Current guidance solution to track
    guid_cur_time   = 0.0 # Current time in guidance solution
    guid_prev_ddto  = EmptyQuad3DoFDDTOSolution(quad.a.n_targs)
    guid_defer_targ = -1 # Next deferred target in consideration (tag number)
    guid_defer_time = 1.e6 # Time until branch point to next deferred target
    guid_lock_time  = 1.e6 # Time at which guidance lock was activated

    # Flags
    flag_update_ddto         = true
    flag_ddto_converged      = false
    flag_log_ddto_results    = false # If set to true, log DDTO results
    flag_guid_lock_activated = false # If set to true, Adaptive-DDTO will be disabled and guidance will fix to the best target at the current time
    flag_descent_complete    = false # If set to true, signals the end of the simulation/descent phase
    flag_guid_lock_staged    = false # If set to true, stage a guidance lock

    # Other variables
    time_last_print = 0.0
    λ_targs_org = quad.a.λ_targs

    # Helper functions
    τ_lu(j) = quad.a.τ_targs[findfirst(i->i==j, quad.a.λ_targs)] # obtain the deferrability index in the trunk of the j-th target

    # Initialize results storage containers
    results_guid_update_ddto_bundles = Array{Quad3DoFDDTOSolution}(undef,0)
    results_guid_update_trajs        = Array{Quad3DoFSolution}(undef, 0)
    results_guid_update_time         = CVector(undef, 0)
    results_sim_time                 = CVector(undef, 0)
    results_sim_state                = CMatrix(undef, quad.a.nx, 0)
    results_sim_control              = CMatrix(undef, quad.a.nu+1, 0)
    results_targs_radii              = CMatrix(undef, quad.n_targs_max, 0)
    results_targs_status             = CMatrix(undef, quad.n_targs_max, 0)
    results_targs_positions          = Array{CMatrix}(undef, 0)

    # Initial print statements
    println("=== Beginning Simulation ===")
    @printf("Time: %.2f s, Alt: %.2f m, Number of targets: %i\n", sim_cur_time, sim_cur_state[3], quad.a.n_targs)

    # ..:: Main Sim Loop ::..
    max_iter = 1e6
    while !flag_descent_complete

        # Compute new DDTO guidance tree (if staged to do so)
        if flag_update_ddto && !flag_guid_lock_activated

            # Obtain (N_max - N_current) new targets from perception stack
            sim_acquire_new_targets!(quad, R_ROI)

            # Set guidance initial conditions as current sim state
            quad.a.z0 = sim_cur_state

            # Guidance solving
            flag_ddto_converged = false
            # try
                guid_cur_opt,guid_cur_ddto,flag_ddto_converged = solve(quad, simulate_solutions=false) # Compute DDTO solution
            # catch e
            #     @printf("---> DDTO ERROR [%.2f s]: %s\n", sim_cur_time, e)
            # end
            if !flag_ddto_converged
                @printf("---> UPDATE [%.2f s]: Guidance lock staged [DDTO computation unsuccessful -- contingency activated!]\n", sim_cur_time)
                guid_cur_ddto = guid_prev_ddto
                flag_guid_lock_staged = true
            end
            guid_prev_ddto = copy(guid_cur_ddto)
            sim_num_ddto += 1

            if !flag_guid_lock_staged
                guid_cur_traj = extract_trunk_segment(quad, guid_cur_ddto) # Track the trunk of DDTO by default
                @printf("---> UPDATE [%.2f s]: DDTO solution successfully recomputed [tracking trunk segment]\n", sim_cur_time)
    
                # Parameter updates
                guid_cur_time = 0.0 # Reset guidance time to zero
                guid_defer_targ = quad.a.λ_targs[1]
                guid_defer_time = guid_cur_ddto.targs[guid_defer_targ].t[τ_lu(guid_defer_targ)]
                λ_targs_org = quad.a.λ_targs

                # If trunk segment has zero length (no deferring could take place),
                # lock guidance to the best target at the current point in time (last index of last DDTO branch solution)
                # as a contingency measure
                if length(guid_cur_traj.t) == 0
                    @printf("---> UPDATE [%.2f s]: Guidance lock staged [DDTO deferral was not possible -- contingency activated!]\n", sim_cur_time)
                    flag_guid_lock_staged = true
                end
            end

            # Flag updates
            flag_update_ddto = false
            flag_log_ddto_results = true
        end
        
        # Check for unsafe targets (radii check)
        if !flag_guid_lock_activated
            cur_targs = copy(quad.a.T_targs)
            for targ in cur_targs
                targ_idx = findfirst(i->i==targ, quad.a.T_targs)

                # Remove target if unsafe
                if quad.R_targs[targ_idx] <= quad.R_targs_min
                    @printf("---> UPDATE [%.2f s]: Removing target %i [bounding radius below the minimum threshold]\n", sim_cur_time, targ)
                    remove_ddto_target!(quad, targ)

                    # If this target was queued for deferral, move to next target for deferral
                    if targ == guid_defer_targ
                        guid_defer_targ = quad.a.λ_targs[1] # Add the next target in the queue to consideration for deferral
                        guid_defer_time = guid_cur_ddto.targs[guid_defer_targ].t[τ_lu(guid_defer_targ)]
                    end
                end

                # Reached minimum target threshold
                if (quad.a.n_targs <= 2) || (quad.a.n_targs <= quad.n_targs_min)
                    @printf("---> UPDATE [%.2f s]: DDTO recomputation staged [target set count below the minimum threshold]\n", sim_cur_time)
                    flag_update_ddto = true
                    break
                end
            end
        end

        # Check for branch switching
        if (guid_cur_time >= guid_defer_time) && !flag_guid_lock_activated
            while guid_cur_time >= guid_defer_time # Ready to determine switch

                # Determine if we should switch or not
                switch_branch = switch_decision(quad, guid_defer_targ)

                # Engage switch by staging DDTO update
                if switch_branch
                    @printf("---> UPDATE [%.2f s]: DDTO recomputation staged [chose to defer to target %i]\n", sim_cur_time, guid_defer_targ)

                    # Remove all targets except for switch target (`guid_defer_targ`)
                    other_targs = copy(quad.a.T_targs)
                    deleteat!(other_targs, findfirst(i->i==guid_defer_targ, other_targs))
                    for targ in other_targs
                        remove_ddto_target!(quad, targ)
                    end

                    flag_update_ddto = true
                    break

                # Remove the current target for deferral and go to the next one
                else
                    @printf("---> UPDATE [%.2f s]: Removing target %i [chose to stay on trunk segment]\n", sim_cur_time, guid_defer_targ)
                    remove_ddto_target!(quad, guid_defer_targ) # Remove the target that was in consideration for deferral
                    guid_defer_targ = quad.a.λ_targs[1] # Add the next target in the queue to consideration for deferral
                    guid_defer_time = guid_cur_ddto.targs[guid_defer_targ].t[τ_lu(guid_defer_targ)]
                end

                # Reached minimum target threshold
                if (quad.a.n_targs <= 2) || (quad.a.n_targs <= quad.n_targs_min)
                    @printf("---> UPDATE [%.2f s]: DDTO recomputation staged [target set count below the minimum threshold]\n", sim_cur_time)
                    flag_update_ddto = true
                    break
                end
            end
        end

        # Update DDTO-locked target parameters
        sim_update_locked_targets!(quad)

        # Check if we have reached the cutoff altitude
        if sim_cur_state[3] <= h_cut && !flag_guid_lock_activated
            @printf("---> UPDATE [%.2f s]: Guidance lock staged [Cutoff altitude reached!]\n", sim_cur_time)
            flag_guid_lock_staged = true
        end

        # Lock guidance to best current target if necessary
        if flag_guid_lock_staged
            
            # Determine the current "best" target in terms of radius and obtain the corresponding trajectory
            targ_best_idx = argmax(quad.R_targs)
            targ_best = quad.a.T_targs[targ_best_idx]
            guid_cur_traj = extract_guid_lock_segment(quad, guid_cur_ddto, targ_best, λ_targs_org)
            
            # Parameter updates
            guid_defer_targ = targ_best
            guid_defer_time = 1e6
            flag_guid_lock_activated = true
            flag_guid_lock_staged = false
            flag_log_ddto_results = true
            guid_lock_time = sim_cur_time
            @printf("---> UPDATE [%.2f s]: Guidance locked to target %i\n", sim_cur_time, guid_defer_targ)

            # Remove all targets except for locked target (`guid_defer_targ`)
            other_targs = copy(quad.a.T_targs)
            deleteat!(other_targs, findfirst(i->i==guid_defer_targ, other_targs))
            for targ in other_targs
                remove_ddto_target!(quad, targ)
            end
        end

        # Integrate the currently-tracked guidance trajectory for one time-step
        cur_ct_dyn = (t,x) -> dynamics(t,x,guid_cur_traj)
        sim_cur_state = rk4_step(sim_cur_state, cur_ct_dyn, guid_cur_time, Δt_sim)
        sim_cur_time  += Δt_sim
        guid_cur_time += Δt_sim

        # Obtain current control for logging
        if flag_ddto_converged
            sim_cur_control = optimal_controller(guid_cur_time, guid_cur_traj.t, guid_cur_traj.u, quad.a.disc)
        else
            sim_cur_control = zeros(quad.a.nu+1)
        end

        # Print sim status update(s)
        if (sim_cur_time - time_last_print) >= Δt_print
            defer_time_remaining = guid_defer_time - guid_cur_time
            if !flag_guid_lock_activated
                @printf("Time: %.2f s, Alt: %.2f m, Number of targets: %i, Next deferred target: %i (%.2f s remaining)\n", sim_cur_time, sim_cur_state[3], quad.a.n_targs, guid_defer_targ, defer_time_remaining)
            else
                @printf("Time: %.2f s, Alt: %.2f m, Guidance locked to target %i!\n", sim_cur_time, sim_cur_state[3], guid_defer_targ)
            end
            # @printf("---> Debug: Applying control [%.2f, %.2f, %.2f]\n", sim_cur_control[1], sim_cur_control[2], sim_cur_control[3])
            time_last_print = sim_cur_time
        end

        # Log continuous sim results
        results_sim_state   = hcat(results_sim_state, sim_cur_state)
        results_sim_control = hcat(results_sim_control, sim_cur_control)
        append!(results_sim_time, sim_cur_time)

        # Log current target radii (if a target index is unallocated, insert -Inf)
        sim_cur_radii = fill(-Inf, quad.n_targs_max)
        sim_cur_radii[quad.a.T_targs] = quad.R_targs
        results_targs_radii = hcat(results_targs_radii, sim_cur_radii)
        
        # Log current target positions (if a target index is unallocated, insert -Inf)
        sim_cur_targ_pos = -Inf * ones(3, quad.n_targs_max)
        sim_cur_targ_pos[:,quad.a.T_targs] = quad.a.zf_targs[1:3,:]
        append!(results_targs_positions, [sim_cur_targ_pos])

        # Log target status (1 = valid, 0 = lost)
        targs_status = zeros(quad.n_targs_max)
        for k=1:quad.n_targs_max
            if k in quad.a.T_targs
                targs_status[k] = 1
            end
        end
        results_targs_status = hcat(results_targs_status, targs_status)
        
        # Log conditional sim results (DDTO)
        if flag_log_ddto_results
            append!(results_guid_update_ddto_bundles, [guid_cur_ddto])
            append!(results_guid_update_trajs, [guid_cur_traj])
            append!(results_guid_update_time, sim_cur_time)
            flag_log_ddto_results = false
        end

        # Terminate sim if we reach the phase completion condition
        if sim_cur_state[3] <= (h_term + h_eps)
            flag_descent_complete = true
            @printf("---> UPDATE [%.2f s]: Terminal altitude condition reached -- landing successful!\n", sim_cur_time)
        end
        if sim_cur_time >= 100
            display("Simulation ran for too long, exiting...")
            flag_descent_complete = true
        end
    end
end