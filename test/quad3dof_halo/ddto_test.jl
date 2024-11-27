using DDTOSCP
include("plots.jl")

function SampleConfig()::Quad3DoFHaloParams{CReal,Int}
    # Load default params
    params = Quad3DoFHaloParams()

    # Configurables
    r0 = [0,0,150] # [m] Initial position (NED frame)
    v0 = [0,0,0]   # [m/s] Initial velocity (NED frame)
    params.a.n_targs = 4

    # Set sample boundary conditions for n_targs_max = 4 targets
    params.a.z0 = [
        r0;
        v0;
        0;
        Inf
    ]
    params.a.u0 = [
        0;
        0;
        0;
        Inf
    ]
    
    # Set terminal condition targets to be on a circle of radius 100m, make sure to include position and velocity in R3
    zf_targs = [[100*cos(2*pi*k/params.a.n_targs); 100*sin(2*pi*k/params.a.n_targs); 0; 0; 0; 0; Inf] for k = 1:params.a.n_targs]
    params.a.zf_targs = reshape(hcat(zf_targs...), params.a.nx, params.a.n_targs)

    return params
end

params = SampleConfig()
solve(params)
;