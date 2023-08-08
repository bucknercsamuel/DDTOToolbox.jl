function generate_initial_guess_scp(params::Quad3DoFCageParams, j::Int)::Solution
    N = params.N_fft
    if params.disc == 0
        N_ctrl = N-1
    elseif params.disc == 1
        N_ctrl = N
    end
    
    # Uniform mapping from 0 to some known maximum time-of-flight
    t_ig = CVector(range(0, stop=params.ToF_max, length=N))

    # Direct linear interpolation from initial to terminal condition
    x_ig = CMatrix(zeros(params.n,N))
    for k = 1:params.n
        x_ig[k,:] = CVector(range(params.z0[k], stop=params.zf_targs[k,j], length=N))
    end
    
    # Use average between min and max and note that most force should be along +Z
    ρ_avg = (params.ρ_max + params.ρ_min)/2
    u_ig = CMatrix(zeros(params.m,N_ctrl))
    u_ig[1,:] = CVector(range(0, stop=0, length=N_ctrl))
    u_ig[2,:] = CVector(range(0, stop=0, length=N_ctrl))
    u_ig[3,:] = CVector(range(ρ_avg, stop=ρ_avg, length=N_ctrl))

    # Use augmented state & control for time dilation
    s_ig = (t_ig[2] - t_ig[1])/(params.τ[2] - params.τ[1])
    x_ig = vcat(x_ig, ones(1,N))
    u_ig = vcat(u_ig, fill(s_ig,1,N_ctrl))

    return Solution(t_ig, x_ig, u_ig, 0)
end

function generate_initial_guess_ddtoscp(params::Quad3DoFCageParams)::Vector{Solution}

    N = 2*params.N_fft
    if params.disc == 0
        N_ctrl = N-1
    elseif params.disc == 1
        N_ctrl = N
    end

    # Compute geometric mean state between all targets and initial condition
    zmean = (params.z0 + sum(params.zf_targs, dims=2))/(params.n_targs+1)

    # Compute trunk state interpolation
    x_ig_trunk = CMatrix(zeros(params.n,params.N_fft))
    for k = 1:params.n
        x_ig_trunk[k,:] = range(params.z0[k], stop=zmean[k,1], length=params.N_fft) |> CVector
    end
    
    # Compute branch state interpolations for all targets
    solutions = Vector{Solution}(undef, params.n_targs)
    for j = 1:params.n_targs
        x_ig_branch = CMatrix(zeros(params.n,params.N_fft))
        for k = 1:params.n
            x_ig_branch[k,:] = range(zmean[k,1], stop=params.zf_targs[k,j], length=params.N_fft) |> CVector
        end
        x_ig = hcat(x_ig_trunk, x_ig_branch)
    
        # Uniform mapping from 0 to some known maximum time-of-flight
        t_ig = CVector(range(0, stop=params.ToF_max, length=N))

        # Use average between min and max and note that most force should be along +Z
        ρ_avg = (params.ρ_max + params.ρ_min)/2
        u_ig = CMatrix(zeros(params.m,N_ctrl))
        u_ig[1,:] = CVector(range(0, stop=0, length=N_ctrl))
        u_ig[2,:] = CVector(range(0, stop=0, length=N_ctrl))
        u_ig[3,:] = CVector(range(ρ_avg, stop=ρ_avg, length=N_ctrl))

        # Use augmented state & control for time dilation
        s_ig = (t_ig[2] - t_ig[1])/(params.τ[2] - params.τ[1])
        x_ig = vcat(x_ig, ones(1,N))
        u_ig = vcat(u_ig, fill(s_ig,1,N_ctrl))

        solutions[j] = Solution(t_ig, x_ig, u_ig, 0)
    end

    return solutions
end