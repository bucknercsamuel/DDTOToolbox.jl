function skyenet_interface(
        num_targs::Int;
        r0::CVector;
        v0::CVector;
        rf::CVector;
        vf::CVector;
        K::Int;
        tf::CReal;
        a_min::CReal;
        a_max::CReal;
        v_max::CReal;
        theta_max::CReal;
        ri_relax::CVector;
        rf_relax::CVector;
        subopt_tol::CReal;
        w_buff::CReal;
        w_trust::CReal;
        scp_iters::CReal;
        tau_max::CReal;
        eps_cvg::CReal;
        n::Int;
        c_x::CVector;
        c_y::CVector;
        R::CVector;
        M0::CVector;
        M1::CVector;
    )::Tuple{
        CVector,
        CVector,
        CVector,
        CVector,
        CVector,
        CVector,
        CVector,
        CReal,
        CReal
    }
    
    ## Define the lander vehicle and scenario parameters
    quad = Lander()
    height = -1

    # >> Vehicle parameters <<
    quad.ρ_min = a_min
    quad.ρ_max = a_max

    # >> Constraint parameters <<
    quad.γ_p = theta_max
    quad.v_max_L = v_max

    # >> Dynamics <<
    quad.Δt = tf / (K-1)

    # >> Obstacle parameters <<
    quad.n_obstacles = n
    quad.R_obstacles = R
    quad.p_obstacles = hcat((c_x,c_y,height*ones(n)))
    quad.H_obstacles = vcat((M0,M1))

    # >> Boundary conditions <<
    quad.r0 = r0
    quad.v0 = v0
    quad.rf_targs = []
    quad.vf_targs = []
    for i = 1:3
        for j = 1:num_targs
            quad.rf_targs.append(rf[j+3*(i-1)])
            quad.vf_targs.append(vf[j+3*(i-1)])
        end
    end

    # >> Target conditions <<
    quad.n_targs = num_targs
    quad.N_targs = repeat(K, num_targs)
    quad.λ_targs = collect(1:num_targs)
    quad.T_targs = collect(1:num_targs)
    quad.ϵ_targs = repeat(subopt_tol, num_targs)

    # >> SCP Params <<
    quad.w_buff = w_buff
    quad.w_trust = w_trust
    quad.w_r0 = ri_relax
    quad.w_rf = rf_relax
    quad.sub_iters = scp_iters
    quad.ϵ_cvg = eps_cvg
    
    # >> Other <<
    quad.τ_max = tau_max
    method = "SCP"

    ## Call DDTO
    ~, DDTO_target_solutions = execute_ddto_solution(quad, method)

    ## Package solutions into a usable format
    t = DDTO_target_solutions[0].sol.t
    s = zeros(K)
    r = []
    v = []
    a = []
    r0_relax = DDTO_target_solutions[0].sol.r0_relax
    rf_relax = []
    tau = tf
    dtau = Δt
    for i = 1:3
        for j = 1:K
            for k = 1:num_targs
                r.append(DDTO_target_solutions[k].sol.r[i,j])
                v.append(DDTO_target_solutions[k].sol.v[i,j])
                a.append(DDTO_target_solutions[k].sol.T[i,j] / quad.mass)
            end
        end
    end
    for i = 1:num_targs
        rf_relax.push(DDTO_target_solutions[i].sol.rf_relax)
    end

    return (t,s,r,v,a,r0_relax,rf_relax,tau,dtau)
end

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
    sols_optimal, DDTO_target_solutions = execute_ddto_solution(quad, method)

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

function execute_ddto_solution(quad::Lander, method::String)::Tuple{Vector{Solution},Vector{BranchSolution}}
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
    
    return DDTO_target_solutions
end