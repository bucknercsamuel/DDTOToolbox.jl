function generate_initial_guess_scp(params::Quad3DoFCageParams, j::Int)::Solution
    N = params.a.N

    # Uniform mapping from 0 to some known maximum time-of-flight
    τ_ig = range(0, stop=1, length=N) |> CVector
    t_ig = range(0, stop=params.a.ToF_max, length=N) |> CVector

    # Direct linear interpolation from initial to terminal condition
    x_ig = CMatrix(zeros(params.a.nx,N))
    for k = 1:params.a.nx-1
        x_ig[k,:] = CVector(range(params.a.z0[k], stop=params.a.zf_targs[k,j], length=N))
    end
    
    # Use average thrust between min and max and note that most force should be along +Z
    ρ_avg = norm(params.g)*params.mass
    ν_ig = CMatrix(zeros(params.a.nu-1,N))
    ν_ig[1,:] = CVector(range(0, stop=0, length=N))
    ν_ig[2,:] = CVector(range(0, stop=0, length=N))
    ν_ig[3,:] = CVector(range(ρ_avg, stop=ρ_avg, length=N))

    # Augmented control
    s_ig = [0, (diff(t_ig) ./ diff(τ_ig))...]
    u_ig = vcat(ν_ig, reshape(s_ig,1,length(s_ig)))

    # Compute cumulative thrust norm for 7th state
    Δt = diff(t_ig)
    x_ig[7,:] = CVector(cumsum([params.a.z0[7],[norm(Δt[k]*ν_ig[:,k]) for k=1:N-1]...]))

    return Solution(τ_ig, x_ig, u_ig, 0)
end

function generate_initial_guess_ddtoscp(params::Quad3DoFCageParams)::DDTOSolution
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
        
        # Compute geometric mean state between remaining targets and initial condition (first 6 states only)
        x_mean = (x_end_prev[1:6] + sum(params.a.zf_targs[1:6,J_rem], dims=2))/(length(J_rem)+1)

        # Append to the trunk solution using current geometric mean
        Δτ = iter == 1 ? params.a.τ_targs[iter] : params.a.τ_targs[iter] - params.a.τ_targs[iter-1]
        x_ig_trunk_new = zeros(params.a.nx,Δτ) |> CMatrix
        for k = 1:params.a.nx-1
            x_ig_trunk_new[k,:] = range(x_end_prev[k], stop=x_mean[k,1], length=Δτ+1)[2:end] |> CVector
        end
        x_ig_trunk = hcat(x_ig_trunk, x_ig_trunk_new)

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

        # Use average thrust between min and max and note that most force should be along +Z
        ρ_avg = (params.ρ_max + params.ρ_min)/2
        ν_ig = zeros(params.a.nu-1,N) |> CMatrix
        ν_ig[1,:] = range(0, stop=0, length=N) |> CVector
        ν_ig[2,:] = range(0, stop=0, length=N) |> CVector
        ν_ig[3,:] = range(ρ_avg, stop=ρ_avg, length=N) |> CVector

        # Augmented control
        s_ig = [0, (diff(t_ig) ./ diff(τ_ig))...]
        u_ig = vcat(ν_ig, reshape(s_ig,1,length(s_ig)))

        # Compute cumulative thrust norm
        x_ig[7,:] = CVector(cumsum([0,[norm(ν_ig[:,k]) for k=1:N-1]...]))

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