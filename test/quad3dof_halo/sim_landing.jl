using DDTOSCP
using Random
using Printf
using Debugger

function simulate_halo_landing(
        r0,              # [m] Initial position (NED frame)
        v0;              # [m/s] Initial velocity (NED frame)
        Δt_sim   = 0.01, # [s] Simulation integration time-step
        Δt_print = 1.,   # [s] Simulation printing update time-step
        R_ROI    = 50.,  # [m] Radius of the region of interest for targets
        h_cut    = 50.,  # [m] Altitude condition to commit to best target
        h_term   = 1.,   # [m] Altitude condition to terminate descent phase
        h_eps    = 1.;   # [m] Acceptable altitude error in termination condition
        greedy   = true, # Select if we should use greedy method instead of DDTO
        greedy_dt = 5    # Greedy update timestep
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
    guid = Dict()
    guid["cur_opt"]      = EmptyQuad3DoFDDTOSolution(quad.a.n_targs) # Most recently-computed optimal solution set
    guid["cur_ddto"]     = EmptyQuad3DoFDDTOSolution(quad.a.n_targs) # Most recently-computed DDTO solution set
    guid["cur_ddto_sim"] = EmptyQuad3DoFDDTOSolution(quad.a.n_targs) # Most recently-computed DDTO simulation set
    guid["cur_traj"]     = EmptyQuad3DoFSolution() # Current guidance solution to track
    guid["cur_time"]     = 0.0 # Current time in guidance solution
    guid["defer_targ"]   = -1 # Next deferred target in consideration (tag number)
    guid["defer_time"]   = 1.e6 # Time until branch point to next deferred target
    guid["lock_time"]    = 1.e6 # Time at which guidance lock was activated
    guid["λ_targs_org"]  = quad.a.λ_targs # Stores initial preference ordering
    guid["comp_params"]  = Quad3DoFHaloParams()

    # Flags
    flags = Dict()
    flags["update_ddto"]         = true
    flags["ddto_converged"]      = false
    flags["log_ddto_results"]    = false # If set to true, log DDTO results
    flags["guid_lock_activated"] = false # If set to true, Adaptive-DDTO will be disabled and guidance will fix to the best target at the current time
    flags["descent_complete"]    = false # If set to true, signals the end of the simulation/descent phase
    flags["guid_lock_staged"]    = false # If set to true, stage a guidance lock

    # Results (to be logged)
    results = Dict()
    results["guid_update_ddto_params"]       = Array{Quad3DoFHaloParams}(undef,0)
    results["guid_update_ddto_bundles"]      = Array{Quad3DoFDDTOSolution}(undef,0)
    results["guid_update_ddto_bundles_sims"] = Array{Quad3DoFDDTOSolution}(undef,0)
    results["guid_update_trajs"]             = Array{Quad3DoFSolution}(undef, 0)
    results["guid_update_time"]              = CVector(undef, 0)
    results["sim_time"]                      = CVector(undef, 0)
    results["sim_state"]                     = CMatrix(undef, quad.a.nx, 0)
    results["sim_control"]                   = CMatrix(undef, quad.a.nu+1, 0)
    results["targs_radii"]                   = CMatrix(undef, quad.n_targs_max, 0)
    results["targs_status"]                  = CMatrix(undef, quad.n_targs_max, 0)
    results["targs_positions"]               = Array{CMatrix}(undef, 0)

    # Other variables
    time_last_print = 0.0
    t_fine = nothing
    u_fine = nothing
    τ_fine = nothing

    # Initial print statements
    println("=== Beginning Simulation ===")
    @printf("Time: %.2f s, Alt: %.2f m, Number of targets: %i\n", sim_cur_time, sim_cur_state[3], quad.a.n_targs)

    # ..:: Main Sim Loop ::..
    max_iter = 1e6
    while !flags["descent_complete"]

        # Execute Adaptive-DDTO algorithm pipeline
        if flags["update_ddto"] && !flags["guid_lock_activated"]
            sim_acquire_new_targets!(quad, R_ROI)
            compute_ddto_guidance!(quad, guid, flags, sim_cur_state, sim_cur_time)
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
        sim_update_locked_targets!(quad)
        if flags["guid_lock_staged"]
            activate_guidance_lock!(quad, guid, flags, sim_cur_time)
            τ_fine = CVector(range(start=guid["cur_traj"].τ[1],stop=guid["cur_traj"].τ[end],length=1001))
            u_fine = CMatrix(hcat([optimal_controller(τ_fine[n],guid["cur_traj"].τ,guid["cur_traj"].u,quad.a.disc) for n = 1:length(τ_fine)]...))
            t_fine = CVector(time_dilation_control_to_wall_clock_time(u_fine[end,:], τ_fine, quad.a.disc))
        end

        # Integrate the currently-tracked guidance trajectory for one time-step
        cur_ct_dyn = (t,x) -> dynamics(t,x,t_fine,u_fine,quad)
        sim_cur_state = rk4_step(sim_cur_state, cur_ct_dyn, guid["cur_time"], Δt_sim)
        sim_cur_time     += Δt_sim
        guid["cur_time"] += Δt_sim

        # Log results
        log_results!(quad, results, guid, flags, sim_cur_state, sim_cur_time)

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
        if sim_cur_time >= .2*quad.a.ToF_max
            display("Simulation ran for too long, exiting...")
            flags["descent_complete"] = true
        end
    end

    return results
end