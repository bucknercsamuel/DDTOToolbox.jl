#= DDTO SCP Formulation.

Author: Samuel Buckner (UW-ACL)
=#

function solve_ddto_scp(lander::Lander, costs_optimal::CVector, sols_optimal::Vector{Solution})::Vector{DDTOSolution}
    # Top-level DDTO solver for all branch points
    #
    # :in lander: The lander object
    # :in costs_optimal: Optimal costs from initial condition
    # :out ddto_branch_sols: Vectorized container for all DDTO branch solutions

    # Define container for each DDTO branch solution
    ddto_branch_sols = Vector{DDTOSolution}(undef, lander.n_targs)
    for k = 1:(lander.n_targs)
        ddto_branch_sols[k] = EmptyDDTOSolution(lander.n_targs-k+1)
    end

    # Define running deferred-decision (DD) trajectory segment cost sum
    cost_dd_sum = 0.

    # Define initial guess/reference traj as optimal solution
    ref_trajs = sols_optimal

    # Perform branching in the order of preference
    n_targs_total = copy(lander.n_targs)
    lander_ = copy(lander) # Temp object to be mutated through DDTO loop
    for k = 1:(n_targs_total-1)

        if VERB_DDTO
            specifiers = repeat("%.3f, ", lander_.n_targs)
            specifiers = specifiers[1:end-2] # Remove string and comma at end
            format_string = "   Chosen suboptimality tolerances: {"*specifiers*"}\n"

            @printf("\n========= Solving DDTO for Branch #%i =========\n", k)
            @eval @printf($format_string, $lander_.ϵ_targs...)
        end

        # Obtain Bisection-optimal DDTO solution for this branch
        ddto_branch_sols[k] = solve_bisection_qcvx_ddto_scp(lander_, costs_optimal, cost_dd_sum, ref_trajs)

        # Determine target to be removed (first in the current list of λ_targs)
        λ_targ = lander_.λ_targs[1]
        deleteat!(lander_.λ_targs, 1)
        pop_idx = findfirst(i->i==λ_targ, lander_.T_targs)

        # Have to do some slicing magic for matrices
        matrix_slice = collect(1:lander_.n_targs)
        deleteat!(matrix_slice, pop_idx)

        # Update lander_ target and IC properties for next branch iteration
        lander_.n_targs -= 1
        deleteat!(lander_.T_targs, pop_idx)
        deleteat!(lander_.N_targs, pop_idx)
        deleteat!(lander_.ϵ_targs, pop_idx)
        lander_.N_targs .-= ddto_branch_sols[k].idx_dd
        lander_.r0 = ddto_branch_sols[k].targ_sols[1].r[:,ddto_branch_sols[k].idx_dd+1]
        lander_.v0 = ddto_branch_sols[k].targ_sols[1].v[:,ddto_branch_sols[k].idx_dd+1]
        lander_.rf_targs = lander_.rf_targs[:,matrix_slice]
        lander_.vf_targs = lander_.vf_targs[:,matrix_slice]

        # Update ref traj using solution from bisection search
        ref_trajs = Vector{Solution}(undef, lander_.n_targs)
        for j = 1:lander_.n_targs
            ref_trajs[j]   = EmptySolution()
            ref_trajs[j].r = ddto_branch_sols[k].targ_sols[j].r[:,(ddto_branch_sols[k].idx_dd+1):end]
            ref_trajs[j].v = ddto_branch_sols[k].targ_sols[j].v[:,(ddto_branch_sols[k].idx_dd+1):end]
            ref_trajs[j].T = ddto_branch_sols[k].targ_sols[j].T[:,(ddto_branch_sols[k].idx_dd+1):end]
            ref_trajs[j].Γ = ddto_branch_sols[k].targ_sols[j].Γ[(ddto_branch_sols[k].idx_dd+1):end]
        end

        # Update deferred-decision (DD) cost for next branch iteration
        cost_dd_sum += ddto_branch_sols[k].cost_dd

        # Parameter update print statements
        if VERB_DDTO && (k < n_targs_total-1)
            @printf("   Removed target %i for next branch iteration\n", λ_targ)
        end
    end

    # Add a final element to the branch solutions for the final target
    if lander.λ_targs[end-1] > lander.λ_targs[end]
        final_idx = 1
    else
        final_idx = 2
    end
    ddto_branch_sols[end].targ_sols = Vector{Solution}(undef, lander.n_targs)
    ddto_branch_sols[end].costs_sol = [ddto_branch_sols[end-1].costs_sol[final_idx]]
    ddto_branch_sols[end].idx_dd    = 0
    ddto_branch_sols[end].cost_dd   = 0

    # Remove deferred states/controls from previous solution final target
    for j = 1:lander.n_targs
        idx_dd = ddto_branch_sols[end-1].idx_dd
        ddto_branch_sols[end].targ_sols[j]       = EmptySolution()
        ddto_branch_sols[end].targ_sols[j].t     = ddto_branch_sols[end-1].targ_sols[final_idx].t[idx_dd+1:end]
        ddto_branch_sols[end].targ_sols[j].r     = ddto_branch_sols[end-1].targ_sols[final_idx].r[:,idx_dd+1:end]
        ddto_branch_sols[end].targ_sols[j].v     = ddto_branch_sols[end-1].targ_sols[final_idx].v[:,idx_dd+1:end]
        ddto_branch_sols[end].targ_sols[j].T     = ddto_branch_sols[end-1].targ_sols[final_idx].T[:,idx_dd+1:end]
        ddto_branch_sols[end].targ_sols[j].Γ     = ddto_branch_sols[end-1].targ_sols[final_idx].Γ[idx_dd+1:end]
        ddto_branch_sols[end].targ_sols[j].T_nrm = ddto_branch_sols[end-1].targ_sols[final_idx].T_nrm[idx_dd+1:end]
        ddto_branch_sols[end].targ_sols[j].γ     = ddto_branch_sols[end-1].targ_sols[final_idx].γ[idx_dd+1:end]
    end

    return ddto_branch_sols
end

function solve_bisection_qcvx_ddto_scp(lander::Lander, costs_optimal::CVector, cost_dd::CReal, ref_trajs::Vector{Solution})::DDTOSolution
    # Uses bisection search to solve quasiconvex optimization problem 
    # to branch to the next-queued target for rejection.
    #
    # :in lander: The lander object
    # :in costs_optimal: Optimal costs from `solve_optimal_pdg_all_targets`
    # :in cost_dd: Running cost for decision deferral
    # :out ddto_solution: Contains the DDTO solution for this target/branch point

    # Initial search bracket
    τ_min = 0
    τ_max = min(min(lander.N_targs...) - 2, lander.τ_max)

    # Bisection search to solve quasiconvex (QCvx) optimization problem
    VERB_DDTO && println("=== Bisection Search for QCvx Optimization ===")
    iter = 1
    while (τ_max - τ_min) > 1
        # Update τ
        τ = Int(ceil(0.5*(τ_max + τ_min)))

        # Compute SCP solution
        (~, ref_traj_update, status_feas) = solve_scp_subproblem(lander, τ, costs_optimal, cost_dd, ref_trajs)

        # Update τ_min or τ_max based on solution convergence
        if status_feas == MOI.OPTIMAL
            τ_min = τ
            solve_status = "Feasible"
            ref_trajs = ref_traj_update # Update ref traj
        else
            τ_max = τ
            solve_status = "Not Feasible"
        end
        VERB_DDTO && @printf("Bisection Iteration #%i -- τ_min: %i, τ_max: %i, status: %s\n", iter, τ_min, τ_max, solve_status)

        # Update iteration count
        iter += 1
    end

    # Set optimal τ
    τ_opt = τ_min
    VERB_DDTO && println("Bisection search terminated -- reached convergence condition (τ_max - τ_min) = 1")

    # Compute converged DDTO solution SCP iteration
    (ddto_solution_scp, ref_traj_update, status_feas) = solve_scp_subproblem(lander, τ_opt, costs_optimal, cost_dd, ref_trajs)
    ddto_solution_scp.idx_dd = τ_opt

    # Port sols back to `DDTOSolution`
    ddto_solution = EmptyDDTOSolution(lander.n_targs)
    ddto_solution.costs_sol = ddto_solution_scp.costs_sol
    ddto_solution.cost_dd   = ddto_solution_scp.cost_dd
    ddto_solution.idx_dd    = ddto_solution_scp.idx_dd
    for j = 1:lander.n_targs
        ddto_solution.targ_sols[j].t     = ddto_solution_scp.targ_sols[j].t
        ddto_solution.targ_sols[j].r     = ddto_solution_scp.targ_sols[j].r
        ddto_solution.targ_sols[j].v     = ddto_solution_scp.targ_sols[j].v
        ddto_solution.targ_sols[j].T     = ddto_solution_scp.targ_sols[j].T
        ddto_solution.targ_sols[j].Γ     = ddto_solution_scp.targ_sols[j].Γ
        ddto_solution.targ_sols[j].cost  = ddto_solution_scp.targ_sols[j].cost
        ddto_solution.targ_sols[j].T_nrm = ddto_solution_scp.targ_sols[j].T_nrm
        ddto_solution.targ_sols[j].γ     = ddto_solution_scp.targ_sols[j].γ
    end

    # Determine solution convergence
    if status_feas == MOI.OPTIMAL
        ref_trajs = ref_traj_update
        @printf("Bisection search successful -- τ_opt: %i\n", τ_opt)
    else
        error("Bisection search unsuccessful. Problem is unsolved.")
    end
    VERB_DDTO && println("Updated costs to each remaining target from initial condition:")
    for j = 1:lander.n_targs
        VERB_DDTO && @printf("   Target: %i, Cost: %.3f\n", lander.T_targs[j], ddto_solution.costs_sol[j] + cost_dd)
    end

    # Remove excess state/control nodes from solution
    for j = 1:length(ddto_solution.targ_sols)
        N_targ = lander.N_targs[j]
        N_targ_ctrl = N_targ - 1
        ddto_solution.targ_sols[j].t     = ddto_solution.targ_sols[j].t[1:N_targ]
        ddto_solution.targ_sols[j].r     = ddto_solution.targ_sols[j].r[:,1:N_targ]
        ddto_solution.targ_sols[j].v     = ddto_solution.targ_sols[j].v[:,1:N_targ]
        ddto_solution.targ_sols[j].T     = ddto_solution.targ_sols[j].T[:,1:N_targ_ctrl]
        ddto_solution.targ_sols[j].Γ     = ddto_solution.targ_sols[j].Γ[1:N_targ_ctrl]
        ddto_solution.targ_sols[j].T_nrm = ddto_solution.targ_sols[j].T_nrm[1:N_targ_ctrl]
        ddto_solution.targ_sols[j].γ     = ddto_solution.targ_sols[j].γ[1:N_targ_ctrl]
    end

    return ddto_solution

end

function solve_scp_subproblem(lander::Lander, τ::Int, costs_optimal::CVector, cost_dd::CReal, ref_trajs::Vector{Solution})::Tuple{DDTOSolutionSCP, Vector{Solution}, MOI.TerminationStatusCode}

    # SCP subproblem iteration
    status_feas_sub = undef
    sols_feas_sub = undef
    for k = 1:lander.sub_iters

        # Solve SCP subproblem
        (sols_feas_sub, status_feas_sub) = solve_feasible_ddto_scp(lander, τ, costs_optimal, cost_dd, ref_trajs)
        
        if status_feas_sub == MOI.OPTIMAL
            solve_status_sub = "Feasible"
        else
            solve_status_sub = "Not Feasible"
            @printf("   > SCP subproblem is infeasible, exiting subproblem iteration.\n")
            break
        end

        # Virtual buffer extraction
        μ_L1 = []
        ν_L1 = []
        for o = 1:lander.n_obstacles
            μ_obs = []
            ν_obs = []
            for j = 1:lander.n_targs
                ν = sols_feas_sub.targ_sols[j].ν
                μ = sols_feas_sub.targ_sols[j].μ
                append!(ν_obs, ν[o,:])
                append!(μ_obs, μ[o,:])
            end
            append!(ν_L1, norm(ν_obs,1))
            append!(μ_L1, norm(μ_obs,1))
        end
        ν_L1_max = max(ν_L1...)
        μ_L1_max = max(μ_L1...)

        # @printf("   SCP Iter: %2.i | Status: %s | μ_L1 = [%.1e,%.1e,%.1e,%.1e,%.1e]\n", k, solve_status_sub, μ_L1...)
        @printf("   SCP Iter: %2.i | Status: %s | μ_L1,max = %.1e\n", k, solve_status_sub, μ_L1_max)
        if μ_L1_max <= lander.ϵ_cvg
            @printf("   > Convergence condition has been reached, exiting subproblem iteration.\n")
        end
    end

    # Package a reference trajectory solution format (Vector{Solution})
    ref_traj_update = Vector{Solution}(undef, lander.n_targs)
    for j = 1:lander.n_targs
        ref_traj_update[j]       = EmptySolution()
        ref_traj_update[j].t     = sols_feas_sub.targ_sols[j].t
        ref_traj_update[j].r     = sols_feas_sub.targ_sols[j].r
        ref_traj_update[j].v     = sols_feas_sub.targ_sols[j].v
        ref_traj_update[j].T     = sols_feas_sub.targ_sols[j].T
        ref_traj_update[j].Γ     = sols_feas_sub.targ_sols[j].Γ
        ref_traj_update[j].cost  = sols_feas_sub.targ_sols[j].cost
        ref_traj_update[j].T_nrm = sols_feas_sub.targ_sols[j].T_nrm
        ref_traj_update[j].γ     = sols_feas_sub.targ_sols[j].γ
    end

    return (sols_feas_sub, ref_traj_update, status_feas_sub)
end

function solve_feasible_ddto_scp(lander::Lander, τ::Int, costs_optimal::CVector, cost_dd::CReal, reference_targ_trajs::Vector{Solution})::Tuple{DDTOSolutionSCP, MOI.TerminationStatusCode}
    # Solve the baseline feasibility problem for DDTO.
    #
    # :in lander: The lander object
    # :in τ: Branch point index
    # :in costs_optimal: Optimal costs from `solve_optimal_pdg_all_targets`
    # :in cost_dd: Running cost for decision deferral
    # :out ddto_solution: Contains the DDTO solution for this target/branch point
    # :out feas_status: Feasibility problem solution status code (see MOI.TerminationStatusCode documentation)

    # ..:: Discrete time interval ::..

    N  = max(lander.N_targs...)
    n  = lander.n_targs
    Δt = lander.Δt
    tf = Δt * (N-1)
    N_ctrl = N-1 # Number of nodes to apply control constraints for (N-1 for ZOH)
    A,B,p = c2d_zoh(lander,Δt)

    # ..:: Make the optimization problem ::..

    # >> Optimizer setup <<
    # mdl = Model(optimizer_with_attributes(ECOS.Optimizer, "verbose" => 0))
    mdl = Model(Mosek.Optimizer)
    JuMP.set_optimizer_attribute(mdl, "LOG", 0) # disable debugging

    # >> Base optimization variables <<
    @variable(mdl, r[1:3,1:N,1:n])
    @variable(mdl, v[1:3,1:N,1:n])
    @variable(mdl, T[1:3,1:N,1:n])
    @variable(mdl, Γ[1:N,1:n])

    # >> SCP variables <<
    # Boundary condition relaxations
    @variable(mdl, r0_relax[1:3])
    @variable(mdl, rf_relax[1:3,1:n])

    # Virtual buffers
    @variable(mdl, ν[1:lander.n_obstacles,1:N,1:n])
    @variable(mdl, μ[1:lander.n_obstacles,1:N,1:n])

    # Trust region variables
    @variable(mdl, η_x[1:N])
    @variable(mdl, η_u[1:N_ctrl])

    # >> Expression holders <<
    subopt = Array{QuadExpr}(undef, N_ctrl, n)

    # >> Convenience functions <<
    X = (k,j) -> [r[:,k,j]; v[:,k,j]] # State at time index k and target j
    U = (k,j) -> [T[:,k,j]; Γ[k,j]]   # Input at time index k and target j

    # ..:: Constraints ::..

    # >> Iterate through targets <<
    for j = 1:n

        # Target N
        N_targ = lander.N_targs[j]
        N_targ_ctrl = N_targ - 1

        # Slice indexing to n without current target j
        J = collect(1:n)
        deleteat!(J, j)

        # >> Dynamics <<
        @constraint(mdl, [k=1:N_targ-1], X(k+1,j) .== A*X(k,j) + B*U(k,j) + p)

        # >> Constant altitude constraint <<
        @constraint(mdl, [k=1:N_targ-1], r[3,k+1,j] == r[3,k,j])

        # >> Thrust bounds <<
        @constraint(mdl, [k=1:N_targ_ctrl], Γ[k,j] >= lander.ρ_min)
        @constraint(mdl, [k=1:N_targ_ctrl], Γ[k,j] <= lander.ρ_max)
        @constraint(mdl, [k=1:N_targ_ctrl], vcat(Γ[k,j], T[:,k,j]) in MOI.SecondOrderCone(4))

        # >> Attitude pointing constraint <<
        @constraint(mdl, [k=1:N_targ_ctrl], dot(T[:,k,j],e_z) >= Γ[k,j]*cos(lander.γ_p))

        # >> Velocity upper bound <<
        # @constraint(mdl, [k=1:N_targ], vcat(lander.v_max_V,v[3,k,j])   in MOI.SecondOrderCone(2))
        @constraint(mdl, [k=1:N_targ], vcat(lander.v_max_L,v[1:2,k,j]) in MOI.SecondOrderCone(3))

        # >> Cage bounds <<
        @constraint(mdl, [k=1:N_targ], r[1,k,j] >= lander.x_arena_lims[1])
        @constraint(mdl, [k=1:N_targ], r[1,k,j] <= lander.x_arena_lims[2])
        @constraint(mdl, [k=1:N_targ], r[2,k,j] >= lander.y_arena_lims[1])
        @constraint(mdl, [k=1:N_targ], r[2,k,j] <= lander.y_arena_lims[2])
        @constraint(mdl, [k=1:N_targ], r[3,k,j] >= lander.z_arena_lims[1])
        @constraint(mdl, [k=1:N_targ], r[3,k,j] <= lander.z_arena_lims[2])

        # >> Identicality <<
        for k = 1:N_targ_ctrl
            if τ > 0
                if k <= τ
                    for l = 1:n-1
                        @constraint(mdl, U(k,j) .== U(k,J[l]))
                    end
                end
            end
        end

        # >> Suboptimality <<
        for k = 1:N_ctrl
            if k <= N_targ_ctrl
                subopt[k,j] = @expression(mdl, Δt*Γ[k,j])
            else
                subopt[k,j] = @expression(mdl, 0.0)
            end
        end

        # >> Zero out state/control nodes from N_targ+1 to N <<
        @constraint(mdl, [k=N_targ+1:N],           X(k,j) .== zeros(lander.n,1))
        @constraint(mdl, [k=N_targ_ctrl+1:N_ctrl], U(k,j) .== zeros(lander.m,1))

        # >> Boundary conditions << 
        @constraint(mdl, r[:,1,j]      .== lander.r0 + r0_relax)
        @constraint(mdl, v[:,1,j]      .== lander.v0)
        @constraint(mdl, r[:,N_targ,j] .== lander.rf_targs[:,j] + rf_relax[:,j])
        @constraint(mdl, v[:,N_targ,j] .== lander.vf_targs[:,j])

        # >> Suboptimality <<
        @constraint(mdl, sum(subopt[:,j]) + cost_dd .<= (1 + lander.ϵ_targs[j]) * costs_optimal[j])

        # >> SCP constraints <<
        # Extract reference trajectory for target j
        r_ref = reference_targ_trajs[j].r
        v_ref = reference_targ_trajs[j].v
        T_ref = reference_targ_trajs[j].T
        Γ_ref = reference_targ_trajs[j].Γ
        x_ref = vcat(r_ref, v_ref)
        u_ref = vcat(T_ref, reshape(Γ_ref, 1, length(Γ_ref)))

        # Linearization constraints
        for o = 1:lander.n_obstacles
            H = lander.H_obstacles[o]
            for k = 1:N
                Δr = r_ref[:,k] - lander.p_obstacles[:,o]
                δr = r[:,k,j] - r_ref[:,k]
                ξ  = norm(H*Δr,2)
                ζ  = transpose(H)*H*Δr / ξ
                @constraint(mdl, ξ + dot(ζ,δr) >= lander.R_obstacles[o] + ν[o,k,j])
                # @constraint(mdl, ν[o,k,j] >= 0)
                @constraint(mdl, vcat(μ[o,k,j], ν[o,k,j]) in MOI.NormOneCone(2))
            end
        end

        # Trust region constraints
        @constraint(mdl, [k=1:N_targ],      vcat(η_x[k], X(k,j) - x_ref[:,k]) in MOI.SecondOrderCone(lander.n+1))
        @constraint(mdl, [k=1:N_targ_ctrl], vcat(η_u[k], U(k,j) - u_ref[:,k]) in MOI.SecondOrderCone(lander.m+1))
    end

    # >> Cost function <<
    @objective(mdl, Min, 
            sum(subopt) + 
            lander.w_buff * sum(μ.^2) + 
            lander.w_trust * (sum(η_x.^2) + sum(η_u.^2)) +
            lander.w_r0 * sum(r0_relax.^2) + 
            lander.w_rf * sum(rf_relax.^2))


    optimize!(mdl)
    feas_status = JuMP.termination_status(mdl)

    r = value.(r)
    v = value.(v)
    T = value.(T)
    Γ = value.(Γ)
    ν = value.(ν)
    μ = value.(μ)
    η_x = value.(η_x)
    η_u = value.(η_u)
    η = [η_x;η_u]

    # ..:: Determine optimal cost and deferred-decision (DD) cost ::..

    costs_sol = CVector(zeros(n))
    cost_dd  = 0
    for j = 1:n
        N_targ = lander.N_targs[j]
        for k = 1:N_targ-1
            costs_sol[j] += Δt*Γ[k,j]
            if k==τ && j==1
                cost_dd = costs_sol[j]
            end
        end
    end

    # ..:: Package the DDTO Solution ::..

    ddto_solution = EmptyDDTOSolutionSCP(n)
    for j = 1:n
        # Raw data
        ddto_solution.targ_sols[j].t = CVector(range(0, stop=tf, length=lander.N_targs[j]))
        ddto_solution.targ_sols[j].r = r[:,:,j]
        ddto_solution.targ_sols[j].v = v[:,:,j]
        ddto_solution.targ_sols[j].T = T[:,:,j]
        ddto_solution.targ_sols[j].Γ = Γ[:,j]
        ddto_solution.targ_sols[j].ν = ν[:,:,j]
        ddto_solution.targ_sols[j].μ = μ[:,:,j]
        ddto_solution.targ_sols[j].r0_relax = r0_relax[:]
        ddto_solution.targ_sols[j].rf_relax = rf_relax[:,j]
        ddto_solution.targ_sols[j].cost = costs_sol[j]

        # Processed data
        ddto_solution.targ_sols[j].T_nrm = CVector([norm(T[:,k,j],2) for k=1:N_ctrl])
        ddto_solution.targ_sols[j].γ     = CVector([acos(dot(T[:,k,j],e_z)/norm(T[:,k,j],2)) for k=1:N_ctrl])
    end
    ddto_solution.costs_sol = costs_sol
    ddto_solution.cost_dd   = cost_dd
    ddto_solution.η         = η

    return (ddto_solution, feas_status)

end