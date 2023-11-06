function generate_initial_guess_scp(params::Quad3DoFCageParams, j::Int)::Solution
    N = params.N
    if params.disc == 0
        N_ctrl = N-1
    elseif params.disc == 1
        N_ctrl = N
    end
    
    # Uniform mapping from 0 to some known maximum time-of-flight
    t_ig = CVector(range(0, stop=params.ToF_max, length=N))

    # Direct linear interpolation from initial to terminal condition
    x_ig = CMatrix(zeros(params.nx,N))
    for k = 1:params.nx-1
        x_ig[k,:] = CVector(range(params.z0[k], stop=params.zf_targs[k,j], length=N))
    end
    
    # Use average between min and max and note that most force should be along +Z
    ρ_avg = norm(params.g)*params.mass
    ν_ig = CMatrix(zeros(params.nu-1,N_ctrl))
    ν_ig[1,:] = CVector(range(0, stop=0, length=N_ctrl))
    ν_ig[2,:] = CVector(range(0, stop=0, length=N_ctrl))
    ν_ig[3,:] = CVector(range(ρ_avg, stop=ρ_avg, length=N_ctrl))

    # Use augmented state & control for time dilation
    s_ig = (t_ig[2] - t_ig[1])/params.Δτ
    u_ig = vcat(ν_ig, fill(s_ig,1,N_ctrl))

    # Compute cumulative thrust norm for 7th state
    Δt = diff(t_ig)
    x_ig[7,:] = CVector(cumsum([params.z0[7],[norm(Δt[k]*ν_ig[:,k]) for k=1:N_ctrl-1]...]))

    return Solution(t_ig, x_ig, u_ig, 0)
end

function generate_initial_guess_ddtoscp(τ, params::Quad3DoFCageParams)::Vector{Solution}

    N = params.N
    if params.disc == 0
        N_ctrl = N-1
    elseif params.disc == 1
        N_ctrl = N
    end

    # Compute geometric mean state between all targets and initial condition (first 6 states only)
    zmean = (params.z0[1:6] + sum(params.zf_targs[1:6,:], dims=2))/(params.n_targs+1)
    
    # Compute trunk state interpolation
    x_ig_trunk = CMatrix(zeros(params.nx,τ))
    for k = 1:params.nx-1
        x_ig_trunk[k,:] = range(params.z0[k], stop=zmean[k,1], length=τ) |> CVector
    end
    
    # Compute branch state interpolations for all targets
    solutions = Vector{Solution}(undef, params.n_targs)
    for j = 1:params.n_targs
        x_ig_branch = CMatrix(zeros(params.nx,N-τ))
        for k = 1:params.nx-1
            range_ = range(zmean[k,1], stop=params.zf_targs[k,j], length=N-τ+1) |> CVector
            x_ig_branch[k,:] = range_[2:end]
        end
        x_ig = hcat(x_ig_trunk, x_ig_branch)
    
        # Uniform mapping from 0 to some known maximum time-of-flight
        t_ig = range(0, stop=params.ToF_max, length=N) |> CVector

        # Use average between min and max and note that most force should be along +Z
        ρ_avg = (params.ρ_max + params.ρ_min)/2
        ν_ig = CMatrix(zeros(params.nu-1,N_ctrl))
        ν_ig[1,:] = range(0, stop=0, length=N_ctrl) |> CVector
        ν_ig[2,:] = range(0, stop=0, length=N_ctrl) |> CVector
        ν_ig[3,:] = range(ρ_avg, stop=ρ_avg, length=N_ctrl) |> CVector

        # Convert to dilation-augmented control
        s_ig = (t_ig[2] - t_ig[1])/params.Δτ
        u_ig = vcat(ν_ig, fill(s_ig,1,N_ctrl))

        # Compute cumulative thrust norm
        x_ig[7,:] = CVector(cumsum([0,[norm(ν_ig[:,k]) for k=1:N_ctrl-1]...]))

        solutions[j] = Solution(t_ig, x_ig, u_ig, 0)
    end

    return solutions
end