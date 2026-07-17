#=
Initial-guess generators for 2-DOF double-integrator SCP and DDTO-SCP.
=#

"""
    generate_initial_guess_scp(params::DIntegrator2DoFParams, j::Int) -> Solution

Linear state interpolation from IC to target `j` with near-zero acceleration
and uniform time-of-flight dilation control.

# Arguments
- `params`: 2-DOF parameters (IC, terminal state for target `j`, horizon).
- `j`: target index selecting `zf_targs[:,j]`.

# Returns
- `Solution` warmstart with dilated time, interpolated state, and augmented control.
"""
function generate_initial_guess_scp(params::DIntegrator2DoFParams, j::Int)::Solution
    N = params.a.N

    # Uniform mapping from 0 to some known maximum time-of-flight
    τ_ig = range(0, stop=1, length=N) |> CVector
    t_ig = range(0, stop=params.a.ToF_max, length=N) |> CVector

    # Direct linear interpolation from initial to terminal condition
    x_ig = CMatrix(zeros(params.a.nx,N))
    for k = 1:params.a.nx-1
        x_ig[k,:] = CVector(range(params.a.z0[k], stop=params.a.zf_targs[k,j], length=N))
    end
    
    # Use zero acceleration in each axis
    ν_ig = 1e-2*ones(2,N)

    # Augmented control
    s_ig = wall_clock_time_to_time_dilation_control(t_ig, τ_ig, params.a.disc)
    u_ig = vcat(ν_ig, reshape(s_ig,1,length(s_ig)))

    return Solution(τ_ig, x_ig, u_ig, 0)
end

"""
    generate_initial_guess_ddtoscp(params::DIntegrator2DoFParams) -> DDTOSolution

Tree-structured DDTO initial guess using geometric-mean trunk segments and
linear branches to each terminal condition.

# Arguments
- `params`: 2-DOF parameters with deferral ordering `λ_targs`.

# Returns
- Tree-structured `DDTOSolution` warmstart for DDTO-SCP.
"""
function generate_initial_guess_ddtoscp(params::DIntegrator2DoFParams)::DDTOSolution
    N = params.a.N

    # Set node deferrability allocation (may have already been set, overrides)
    set_deferrability_node_allocation!(params)
    
    # Compute tree interpolation recursively 
    # (removing targets by deferrability order one-by-one)
    solution = EmptyDDTOSolution(params.a.n_targs)
    x_ig_trunk = CMatrix(undef, params.a.nx, 0)
    J_rem = Vector(1:params.a.n_targs) # store all remaining targets as we remove them
    x_end_prev = params.a.z0
    iter = 1
    for j ∈ params.a.λ_targs
        τ = params.a.τ_targs[iter]
        
        Δτ = iter == 1 ? params.a.τ_targs[iter] : params.a.τ_targs[iter] - params.a.τ_targs[iter-1]
        idx_cat = iter == 1 ? 1 : 2
        if Δτ > 0
            # Compute geometric mean state between remaining targets and initial condition
            x_mean = (x_end_prev + sum(params.a.zf_targs[:,J_rem], dims=2))/(length(J_rem)+1)

            # Append to the trunk solution using current geometric mean
            x_ig_trunk_new = zeros(params.a.nx,Δτ) |> CMatrix
            for k = 1:params.a.nx-1
                x_ig_trunk_new[k,:] = range(x_end_prev[k], stop=x_mean[k,1], length=Δτ+idx_cat-1)[idx_cat:end] |> CVector
            end
            x_ig_trunk = hcat(x_ig_trunk, x_ig_trunk_new)
        else
            x_mean = x_end_prev
        end

        # Construct the branch segment and concatenate it onto the trunk to form the solution to target j
        x_ig_branch = zeros(params.a.nx,N-τ) |> CMatrix
        for k = 1:params.a.nx-1
            range_ = range(x_mean[k,1], stop=params.a.zf_targs[k,j], length=N-τ+1) |> CVector
            x_ig_branch[k,:] = range_[2:end]
        end
        x_ig = hcat(x_ig_trunk, x_ig_branch)
    
        # Uniform mapping from 0 to some known maximum time-of-flight
        τ_ig = range(0, stop=1, length=N) |> CVector
        t_ig = range(0, stop=params.a.ToF_max, length=N) |> CVector

        # Zero acceleration in each axis
        ν_ig = 1e-2*ones(2,N)

        # Augmented control
        s_ig = wall_clock_time_to_time_dilation_control(t_ig, τ_ig, params.a.disc)
        u_ig = vcat(ν_ig, reshape(s_ig,1,length(s_ig)))

        # Record solution to j-th target
        solution.targs[j] = Solution(τ_ig, x_ig, u_ig, 0)

        # Iteration updates
        iter += 1
        pop_idx = findfirst(i->i==j, J_rem)
        deleteat!(J_rem, pop_idx) # remove j-th target from pool of remaining targets
        x_end_prev = x_ig_trunk[:,end]
    end

    return solution
end