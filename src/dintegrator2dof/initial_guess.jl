function generate_initial_guess_scp(params::DIntegrator2DoFParams, j::Int)::Solution
    N = params.N

    # Uniform mapping from 0 to some known maximum time-of-flight
    τ_ig = range(0, stop=1, length=N) |> CVector
    t_ig = range(0, stop=params.ToF_max, length=N) |> CVector

    # Direct linear interpolation from initial to terminal condition
    x_ig = CMatrix(zeros(params.nx,N))
    for k = 1:params.nx-1
        x_ig[k,:] = CVector(range(params.z0[k], stop=params.zf_targs[k,j], length=N))
    end
    
    # Use zero acceleration in each axis
    ν_ig = zeros(2,N)

    # Augmented control
    s_ig = [0, (diff(t_ig) ./ diff(τ_ig))...]
    u_ig = vcat(ν_ig, reshape(s_ig,1,length(s_ig)))

    return Solution(τ_ig, x_ig, u_ig, 0)
end

function generate_initial_guess_ddtoscp(params::DIntegrator2DoFParams)::DDTOSolution
    N = params.N

    # Set node deferrability allocation (may have already been set, overrides)
    set_deferrability_node_allocation!(params)
    
    # Compute tree interpolation recursively 
    # (removing targets by deferrability order one-by-one)
    solution = EmptyDDTOSolution(params.n_targs)
    x_ig_trunk = CMatrix(undef, params.nx, 0)
    J_rem = Vector(1:params.n_targs) # store all remaining targets as we remove them
    x_end_prev = params.z0
    iter = 1
    for j ∈ params.λ_targs
        τ = params.τ_targs[iter]
        
        Δτ = iter == 1 ? params.τ_targs[iter] : params.τ_targs[iter] - params.τ_targs[iter-1]
        if Δτ > 0
            # Compute geometric mean state between remaining targets and initial condition
            x_mean = (x_end_prev + sum(params.zf_targs[:,J_rem], dims=2))/(length(J_rem)+1)

            # Append to the trunk solution using current geometric mean
            x_ig_trunk_new = zeros(params.nx,Δτ) |> CMatrix
            for k = 1:params.nx
                x_ig_trunk_new[k,:] = range(x_end_prev[k], stop=x_mean[k,1], length=Δτ+1)[2:end] |> CVector
            end
            x_ig_trunk = hcat(x_ig_trunk, x_ig_trunk_new)
        else
            x_mean = x_end_prev
        end

        # Construct the branch segment and concatenate it onto the trunk to form the solution to target j
        x_ig_branch = zeros(params.nx,N-τ) |> CMatrix
        for k = 1:params.nx-1
            range_ = range(x_mean[k,1], stop=params.zf_targs[k,j], length=N-τ+1) |> CVector
            x_ig_branch[k,:] = range_[2:end]
        end
        x_ig = hcat(x_ig_trunk, x_ig_branch)
    
        # Uniform mapping from 0 to some known maximum time-of-flight
        τ_ig = range(0, stop=1, length=N) |> CVector
        t_ig = range(0, stop=params.ToF_max, length=N) |> CVector

        # Zero acceleration in each axis
        ν_ig = zeros(2,N)

        # Augmented control
        s_ig = [0, (diff(t_ig) ./ diff(τ_ig))...]
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