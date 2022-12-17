function standalone_interface(scenario::String="default", method::String="SCP")
    ## Define the lander vehicle and scenario parameters
    quad = Lander()

    if scenario=="default"
        ~ # Do nothing
    elseif scenario=="toy1"
        scenario_toy1!(quad)
    elseif scenario=="onr_demo"
        scenario_onr_demo!(quad)
    end

    ## Call DDTO
    execute_ddto_solution(quad, method)

    ## Plot solutions
    # Setup
    set_fonts()
    PyPlot.close("all")
    pygui(false)

    # Plot functionality
    plot_parametric_optimal_trajectories(quad, sols_optimal)
    plot_parametric_ddto_trajectories(quad, DDTO_target_solutions)
    plot_states(quad, DDTO_target_solutions)
end

function execute_ddto_solution(quad::Lander, method::String)
    if method == "Baseline"
        @time begin
            @time begin
                # ..:: Solve for independently-optimal solutions to each target ::..
                (sols_optimal) = solve_optimal_pdg_all_targets(quad)
                costs_optimal = CVector(zeros(quad.n_targs))
                for k = 1:quad.n_targs
                costs_optimal[k] = sols_optimal[k].cost
                end
                println("\n Solve time for generating optimal solutions to each target:")
            end
    
            @time begin
                # ..:: Solve for DDTO branching solutions to ALL targets ::..
                sols_ddto = solve_ddto_pdg(quad, costs_optimal)
                println("\n Solve time for generating DDTO branch solutions to all targets:")
            end
            println("\n Solve time for the full DDTO solution stack:")
        end
        DDTO_target_solutions = extract_target_trajectories(quad, sols_ddto)
    
    elseif method == "SCP"
        @time begin
            @time begin
                # ..:: Solve for independently-optimal solutions to each target ::..
                (sols_optimal) = solve_optimal_pdg_all_targets(quad)
                costs_optimal = CVector(zeros(quad.n_targs))
                for k = 1:quad.n_targs
                costs_optimal[k] = sols_optimal[k].cost
                end
                println("\n Solve time for generating optimal solutions to each target:")
            end
    
            @time begin
                # ..:: Solve for DDTO branching solutions to ALL targets ::..
                sols_ddto = solve_ddto_scp(quad, costs_optimal, sols_optimal)
                println("\n Solve time for generating DDTO branch solutions to all targets:")
            end
            println("\n Solve time for the full DDTO solution stack:")
        end
        DDTO_target_solutions = extract_target_trajectories(quad, sols_ddto)
        
    else
        println("Selection invalid.")
    end
end