function generate_initial_guess(params::Params, j::Int)::Solution
    if params.free_final_time
        N = params.N_fft
    else
        N = params.N_targs[j]
    end
    if params.disc == 0
        N_ctrl = N-1
    elseif params.disc == 1
        N_ctrl = N
    end
    
    # Uniform mapping from 0 to some known maximum time-of-flight
    if params.free_final_time
        t_ig = CVector(range(0, stop=params.ToF_max, length=N))
    else
        t_ig = CVector(range(0, stop=N*params.Δt, length=N))
    end

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
    u_ig[4,:] = CVector(range(ρ_avg, stop=ρ_avg, length=N_ctrl))

    return Solution(t_ig, x_ig, u_ig, 0)
end