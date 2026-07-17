#=
Shared optimization primitives for DDTOToolbox: affine variable scaling,
line-search helpers, and continuous-time successive convexification (CT-SCvx)
subproblem iteration / transcription.
=#

# ..:: Numerical Scaling ::..

"""
    scaling_matrices(xmin, xmax) -> (S, s)

Build diagonal affine scaling ``x = S x_s + s`` that maps a roughly unit box
onto bounds `[xmin, xmax]`.

# Arguments
- `xmin`: elementwise lower bounds on the physical variable
- `xmax`: elementwise upper bounds on the physical variable

# Returns
- `S`: diagonal (or scalar) scaling matrix
- `s`: affine offset ``(xmin + xmax)/2``
"""
function scaling_matrices(xmin, xmax)
    make_diagonal(x) = Diagonal(x)
    make_diagonal(x::Number) = x
    s = (xmin + xmax) / 2
    S = make_diagonal(max.(1.0, abs.((xmax - xmin) / 2)))
    return S,s
end

"""
    unscale(xs, xmin, xmax)

Map scaled variables back to physical coordinates.

# Arguments
- `xs`: scaled variable array (vector or leading-dimension state array)
- `xmin`: physical lower bounds used to build the scaling
- `xmax`: physical upper bounds used to build the scaling

# Returns
- `x`: physical-coordinate array with the same shape as `xs`
"""
function unscale(xs, xmin, xmax)
    S,s = scaling_matrices(xmin, xmax)
    dims = size(xs)
    if length(dims) == 2
        x = S*xs .+ s
    else
        xs_reshape = reshape(xs, dims[1], prod(dims[2:end]))
        x_reshape = S*xs_reshape .+ s
        x = reshape(x_reshape, dims...)
    end
    return x
end

"""
    remove_ref_zeros!(x_ref, u_ref; ϵ_small=1e-6)

Replace exact zeros in reference trajectories for numerical stability.

# Arguments
- `x_ref`: reference state trajectory (mutated in place)
- `u_ref`: reference control trajectory (mutated in place)
- `ϵ_small`: replacement value written over exact zeros

# Returns
- nothing

# Notes
Mutates `x_ref` and `u_ref`.
"""
function remove_ref_zeros!(x_ref, u_ref; ϵ_small=1e-6)
    x_ref[x_ref .== 0] .= ϵ_small
    u_ref[u_ref .== 0] .= ϵ_small
end

# ..:: Line Search Optimization ::..

"""
    bisection_search_min_feasible(fun, τ_min::Int, τ_max::Int; ϵ_tol=1, verbose=true) -> Int

Integer bisection that finds the *minimum* feasible `τ` for which `fun(τ)` is
finite. Distinct from DDTO branch bisection, which maximizes deferral.

# Arguments
- `fun::Function`: maps integer `τ` to a cost (`Inf` means infeasible)
- `τ_min::Int`: lower bound of the search bracket
- `τ_max::Int`: upper bound of the search bracket
- `ϵ_tol::Int`: terminate when `τ_max - τ_min ≤ ϵ_tol`
- `verbose::Bool`: print iteration status when `true`

# Returns
- `τ_opt::Int`: smallest feasible integer `τ` found in the bracket
"""
function bisection_search_min_feasible(fun::Function, τ_min::Int, τ_max::Int; ϵ_tol::Int=1, verbose::Bool=true)::Int
    iter = 1
    while (τ_max - τ_min) > ϵ_tol
        # Update τ
        τ = Int(ceil(0.5*(τ_max + τ_min)))

        # Compute feasible DDTO
        cost = fun(τ)

        # Update τ_max or τ_min based on solution convergence
        if ~isinf(cost)
            τ_max = τ
            solve_status = "Feasible"
        else
            τ_min = τ
            solve_status = "Not Feasible"
        end
        verbose && @printf("Iteration: %i, τ_min: %i, τ_max: %i -- %s\n", iter, τ_min, τ_max, solve_status)

        # Update iteration count
        iter += 1
    end

    # Set optimal τ
    τ_opt = τ_max

    return τ_opt
end

"""
    bisection_search_min_feasible(fun, τ_min::Float64, τ_max::Float64; ϵ_tol=1e-3, verbose=true) -> Float64

Floating-point variant of minimum-feasible bisection search over `τ`.

# Arguments
- `fun::Function`: maps `τ` to a cost (`Inf` means infeasible)
- `τ_min::Float64`: lower bound of the search bracket
- `τ_max::Float64`: upper bound of the search bracket
- `ϵ_tol::Float64`: terminate when `τ_max - τ_min ≤ ϵ_tol`
- `verbose::Bool`: print iteration status when `true`

# Returns
- `τ_opt::Float64`: smallest feasible `τ` found in the bracket
"""
function bisection_search_min_feasible(fun::Function, τ_min::Float64, τ_max::Float64; ϵ_tol::Float64=1e-3, verbose::Bool=true)::Float64
    iter = 1
    while (τ_max - τ_min) > ϵ_tol
        # Update τ
        τ = 0.5*(τ_max + τ_min)

        # Compute feasible DDTO
        cost = fun(τ)

        # Update τ_max or τ_min based on solution convergence
        if ~isinf(cost)
            τ_max = τ
            solve_status = "Feasible"
        else
            τ_min = τ
            solve_status = "Not Feasible"
        end
        verbose && @printf("Iteration: %i, τ_min: %.2f, τ_max: %.2f -- %s\n", iter, τ_min, τ_max, solve_status)

        # Update iteration count
        iter += 1
    end

    # Set optimal τ
    τ_opt = τ_max

    return τ_opt
end

"""
    golden_section(f, a, b; tol=1e-3, get_first_feasible=false, verbose=true) -> (x_sol, f_sol, x_last_feas)

Golden-section search minimizing a unimodal scalar function on `[a, b]`
(Kochenderfer & Wheeler, *Algorithms for Optimization*, 2019).

# Arguments
- `f::Function`: scalar objective `x -> f(x)` to minimize (`Inf` treated as infeasible)
- `a::Float64`: search-domain lower bound
- `b::Float64`: search-domain upper bound
- `tol::Float64`: termination tolerance on the bracket width
- `get_first_feasible::Bool`: if `true`, stop at the first finite evaluation
- `verbose::Bool`: print bracket updates when `true`

# Returns
- `x_sol::Float64`: approximate minimizer
- `f_sol::Float64`: objective value `f(x_sol)`
- `x_last_feas::Float64`: last feasible (`finite`) abscissa seen, or `Inf`
"""
function golden_section(f::Function, a::Float64, b::Float64; tol::Float64=1e-3, get_first_feasible::Bool=false, verbose::Bool=true)::Tuple{Float64, Float64, Float64}
    ϕ = (1+√5)/2
    n = ceil(log((b-a)/tol)/log(ϕ)+1)
    ρ = ϕ-1
    d = ρ*b+(1-ρ)*a
    yd = f(d)
    x_sol_last_feas = Inf
    for ~ = 1:n-1
        c = ρ*a+(1-ρ)*b
        yc = f(c)
        if yc < yd
            b,d,yd = d,c,yc
        else
            a,b = b,c
        end
        bracket = sort([a,b,c,d])
        verbose && @printf("Golden Bracket: [%.3f,%.3f,%.3f,%.3f] -- Loss: %.3f\n", bracket..., yc)
        if get_first_feasible && !isinf(yc)
            verbose && println("Feasible solution found, breaking Golden Section Search.")
            break
        end
        if !isinf(yc)
            x_sol_last_feas = b
        end
        flush(stdout)
    end
    x_sol = b
    sol = (x_sol,f(x_sol),x_sol_last_feas)
    return sol
end

# ..:: Continuous-Time Successive Convexification (CT-SCvx) ::..

"""
    solve_ctscvx_iteration(params, ref_traj, subproblem_; single_iter=false) -> (solution, feas_status, scvx_converged)

Run the outer CT-SCvx / PTR loop: repeatedly solve `subproblem_`, apply
[`param_update_law!`](@ref), and stop on feasibility failure, penalty
convergence, or iteration cap.

# Arguments
- `params`: problem parameter object (copied internally for weight updates)
- `ref_traj::Solution`: initial reference trajectory for the first subproblem
- `subproblem_::Function`: `(params, ref_traj, iter) -> (Solution, TerminationStatusCode, Bool)`
- `single_iter::Bool`: if `true`, run only one subproblem iterate

# Returns
- `solution::Solution`: latest subproblem solution
- `feas_status`: MOI termination status of the last solved subproblem
- `scvx_converged::Bool`: `true` if PTR penalty convergence was declared
"""
function solve_ctscvx_iteration(
        params,
        ref_traj::Solution,
        subproblem_::Function;
        single_iter::Bool=false
    )::Tuple{Solution, MOI.TerminationStatusCode, Bool}
    
    feas_status = undef
    solution = ref_traj
    scvx_converged = false
    params_ = copy(params)
    for k = 1:params.a.scp_iters
        # Solve SCvx subproblem
        (solution, feas_status, scvx_converged) = subproblem_(params_, solution, k)

        # Update problem parameters
        param_update_law!(params_)

        if single_iter
            break # skip all convergence criterion, only going to run a single (potentially-infeasible) iterate!
        end

        if feas_status != MOI.OPTIMAL && feas_status != MOI.ALMOST_OPTIMAL
            @printf("   > SCvx subproblem is infeasible (MOI status: %s), exiting subproblem iteration.\n", feas_status)
            break
        end
        if scvx_converged
            @printf("   > Convergence condition has been reached, exiting subproblem iteration.\n")
            break
        end
    end
    VERB_OPT && @printf("   > Total cost: %.3f\n\n", solution.cost)

    return (solution, feas_status, scvx_converged)
end

"""
    solve_ctscvx_subproblem(params, ref_traj, z0, zf, u0, uf, dyn_nl, dyn_lin, prob_cost_, prob_constraints_, scp_iter; CTCS_idxs, dilation_idxs, TOF_idxs) -> (sol, feas_status, scp_sub_cvged)

Formulate and solve one continuous-time SCvx subproblem about `ref_traj`.

# Arguments
- `params`: problem parameters
- `ref_traj::Solution`: reference trajectory for linearization
- `z0::CVector`: initial-state boundary values (`Inf` = unconstrained entry)
- `zf::CVector`: terminal-state boundary values (`Inf` = unconstrained entry)
- `u0::CVector`: initial-control boundary values (`Inf` = unconstrained entry)
- `uf::CVector`: terminal-control boundary values (`Inf` = unconstrained entry)
- `dyn_nl`: nonlinear dynamics `(t,x,u,p) -> f`
- `dyn_lin`: linearized dynamics `(t,x,u,p) -> (A,B,Σ,w)`
- `prob_cost_`: `(mdl,x,u) -> (J_running, J_term)` cost builder
- `prob_constraints_`: `(mdl,x,u,ref_traj) -> ν_buff` constraint builder
- `scp_iter::Int`: current outer SCP iteration index (for logging)
- `CTCS_idxs`: state indices of CTCS violation channels (default: last state)
- `dilation_idxs`: control indices of time-dilation channels (default: last control)
- `TOF_idxs`: dilation indices that also carry time-of-flight bounds

# Returns
- `sol::Solution`: optimized nodal trajectory (empty on infeasibility)
- `feas_status`: MOI termination status
- `scp_sub_cvged::Bool`: `true` if virtual-control/buffer/trust penalties meet tolerances
"""
function solve_ctscvx_subproblem(
        params, 
        ref_traj::Solution, 
        z0::CVector, 
        zf::CVector, 
        u0::CVector, 
        uf::CVector,
        dyn_nl,
        dyn_lin,
        prob_cost_,
        prob_constraints_, 
        scp_iter::Int;
        CTCS_idxs = nothing,
        dilation_idxs = nothing,
        TOF_idxs = nothing
    )::Tuple{Solution, MOI.TerminationStatusCode, Bool}

    # ..:: Setup ::..
    # Optimizer configuration
    if params.a.ctcs_enabled
        mdl, solver_type = solver_setup(SOLVER_CTCS_ENABLED)
    else
        mdl, solver_type = solver_setup(SOLVER_CTCS_DISABLED)
    end
    trust_region_type = solver_type

    # Sizing parameters
    nx = params.a.nx
    nu = params.a.nu
    N  = params.a.N
    if params.a.disc == 0
        N_ctrl = N-1
    elseif params.a.disc == 1
        N_ctrl = N
    end

    # Param check(s)
    if params.a.disc != 0 && params.a.disc != 1
        error("Please select a valid discretization hold order.")
    end
    
    # ..:: Optimization variables ::..
    # Unscaled variables
    x_us = @variable(mdl, [1:nx,1:N])
    u_us = @variable(mdl, [1:nu,1:N_ctrl])

    # Apply affine scaling
    x = params.a.Sx*x_us .+ repeat(params.a.sx,1,N)
    u = params.a.Su*u_us .+ repeat(params.a.su,1,N_ctrl)

    # SCP-specific
    ν_ctrl = @variable(mdl, [1:nx,1:(N-1)]) # virtual control
    η = nothing
    if solver_type != "QP"
        η = @variable(mdl, [1:N]) # trust region penalty terms
    end
    @variable(mdl, μ_ctrl) # virtual control objective slack
    @variable(mdl, μ_buff) # virtual buffer objective slack
    @variable(mdl, η_s) # trust region objective slack

    # ..:: Transcription ::..
    TS_batch = Vector{Tuple{CReal,CReal}}(undef,0)
    X_batch = Vector{Tuple{CVector,CVector}}(undef,0)
    U_batch = Vector{Tuple{CVector,CVector}}(undef,0)

    # Build C2D batch
    idxs_traj = add_traj_to_c2d_batch!(ref_traj, TS_batch, X_batch, U_batch; disc=params.a.disc)

    # Perform linearization and discretization
    result = @timed c2d_nonlinear(TS_batch,X_batch,U_batch,k->dyn_nl,k->dyn_lin,params.a.disc)
    Ak,Bmk,Bpk,_,wk,_ = result[1]
    time_trans = result[2]

    # ..:: Make the optimization problem ::..
    # Objective & constraints
    J_running,J_term = prob_cost_(mdl,x,u)
    ν_buff = prob_constraints_(mdl,x,u,ref_traj)

    # Dynamics
    X(k) = x[:,k]
    U(k) = u[:,k]
    SxInv = inv(params.a.Sx)
    SuInv = inv(params.a.Su)
    if params.a.disc == 0
        @constraint(mdl, [k∈idxs_traj], SxInv*X(k+1) .== SxInv*(Ak[k]*X(k) + Bmk[k]*U(k) + wk[k]) + ν_ctrl[:,k])
    elseif params.a.disc == 1
        @constraint(mdl, [k∈idxs_traj], SxInv*X(k+1) .== SxInv*(Ak[k]*X(k) + Bmk[k]*U(k) + Bpk[k]*U(k+1) + wk[k]) + ν_ctrl[:,k])
    end

    # CTCS violation
    if isnothing(CTCS_idxs)
        CTCS_idxs = [nx]
    end
    if params.a.ctcs_enabled
        for ctcs_idx in CTCS_idxs
            @constraint(mdl, [k=1:N], x[ctcs_idx,k]/params.a.ϵ_ctcs <= 1)
            @constraint(mdl, [k=1:N], x[ctcs_idx,k] >= 0)
        end
    end

    # Time dilation
    if isnothing(dilation_idxs)
        dilation_idxs = [nu]
    end
    if isnothing(TOF_idxs)
        TOF_idxs = dilation_idxs
    end
    for dilation_idx in dilation_idxs
        s = u[dilation_idx,:]
        t = time_dilation_control_to_wall_clock_time(s, ref_traj.t, params.a.disc)
        Δt = diff(t)
        @constraint(mdl, [k=1:N-1], params.a.Δt_min/params.a.Δt_max <= Δt[k]/params.a.Δt_max <= 1)
        @constraint(mdl, [k=1:N], s[k] >= 0)
        if dilation_idx in TOF_idxs
            @constraint(mdl, params.a.ToF_min/params.a.ToF_max <= t[end]/params.a.ToF_max <= 1)
        end
    end

    # State boundary conditions
    for k = 1:nx
        if ~isinf(z0[k])
            @constraint(mdl, SxInv[k,k]*x[k,1] == SxInv[k,k]*z0[k])
        end
        if ~isinf(zf[k])
            @constraint(mdl, SxInv[k,k]*x[k,N] == SxInv[k,k]*zf[k])
        end
    end

    # Input boundary conditions
    for k = 1:nu
        if ~isinf(u0[k])
            @constraint(mdl, SuInv[k,k]*u[k,1] == SuInv[k,k]*u0[k])
        end
        if ~isinf(uf[k])
            @constraint(mdl, SuInv[k,k]*u[k,N] == SuInv[k,k]*uf[k])
        end
    end

    # Trust region constraints
    δX(k) = SxInv*(X(k) .- ref_traj.x[:,k]) 
    δU(k) = k < N_ctrl ? SuInv*(U(k) .- ref_traj.u[:,k]) : zeros(nu)
    if trust_region_type == "QP"
        η_s = sum([δX(k)'*δX(k) + δU(k)'*δU(k) for k=1:N])
    else
        @constraint(mdl, [k=1:N], δX(k)'*δX(k) + δU(k)'*δU(k) <= η[k])
        @constraint(mdl, vcat(η_s, η) in SecondOrderCone())
        @constraint(mdl, η_s >= 0)
    end

    # Virtualization constraints
    @constraint(mdl, vcat(μ_ctrl, vec(ν_ctrl)) in MOI.NormOneCone(length(vec(ν_ctrl))+1))
    if length(ν_buff) > 0
        @constraint(mdl, vcat(μ_buff, vec(ν_buff)) in MOI.NormOneCone(length(vec(ν_buff))+1))
        @constraint(mdl, μ_buff >= 0)
    else
        @constraint(mdl, μ_buff == 0)
    end
    @constraint(mdl, μ_ctrl >= 0)

    # ..:: Solve the problem and save the solution ::..
    # Cost function
    obj_scale = 1/sqrt(max(params.a.w_obj_sing, params.a.w_trust, params.a.w_buff, params.a.w_ctrl))
    J_cost = sum(J_running) + J_term
    @objective(mdl, Min,
        (params.a.w_obj_sing * J_cost
      + params.a.w_trust * η_s 
      + params.a.w_buff * μ_buff
      + params.a.w_ctrl * μ_ctrl)*obj_scale)

    result = @timed optimize!(mdl)
    time_solve = result[2]
    feas_status = JuMP.termination_status(mdl)
    if feas_status == MOI.OPTIMAL || feas_status == MOI.ALMOST_OPTIMAL
        solve_status = "Feasible"
    else
        solve_status = "Not Feasible"
        return (EmptySolution(), feas_status, false)
    end

    # Package the solution
    x = value.(x)
    u = value.(u)
    cost = value.(J_cost)
    sol = Solution(ref_traj.t,x,u,cost)

    # ..:: Determine if PTR subproblem has converged ::..
    μ_buff_pen = value.(μ_buff)
    μ_ctrl_pen = value.(μ_ctrl)
    η_pen = value.(η_s)
    if (μ_ctrl_pen <= params.a.ϵ_ctrl) && (μ_buff_pen <= params.a.ϵ_buff) && (η_pen <= params.a.ϵ_trust)
        scp_sub_cvged = true
    else
        scp_sub_cvged = false
    end

    # Print update
    if scp_iter == 1
        VERB_OPT && @printf("   |------------------------------------- SCP Subproblem ------------------------------------|\n")
        VERB_OPT && @printf("   | Iter |  Status  | Trs [ms] | Slv [ms] |   Cost    | μ_ctrl_pen | μ_buff_pen |   η_pen   |\n")
        VERB_OPT && @printf("   |-----------------------------------------------------------------------------------------|\n")
    end
    VERB_OPT && @printf("   |  %2.i  | %s |   % 4.f   |   % 4.f   | % .2e | %s  | %s  | %s |\n", 
        scp_iter, 
        convert_to_colored_string(solve_status,("Feasible",)),
        time_trans*1e3,
        time_solve*1e3,
        cost,
        convert_to_colored_string(μ_ctrl_pen,params.a.ϵ_ctrl),
        convert_to_colored_string(μ_buff_pen,params.a.ϵ_buff),
        convert_to_colored_string(η_pen,params.a.ϵ_trust))
    if scp_iter == params.a.scp_iters || scp_sub_cvged
        VERB_OPT && @printf("   |-----------------------------------------------------------------------------------------|\n")
    end
    flush(stdout)
    
    return (sol, feas_status, scp_sub_cvged)
end

# ..:: Other ::..

"""
    heaviside(x::AbstractFloat) -> Number

Unit step function used by CTCS / symbolic Jacobian expressions.

# Arguments
- `x::AbstractFloat`: scalar argument

# Returns
- `0` if `x < 0`, otherwise `1` (same type as `x`)
"""
heaviside(x::AbstractFloat) = ifelse(x < 0, zero(x), one(x))