using DDTOSCP

# Precompiles the DDTO solve function (other functions don't need to be precompiled, mainly just this one)
begin
    function SampleConfig()::Quad3DoFHaloParams{CReal,Int}
        # Load default params
        params = Quad3DoFHaloParams()

        # Setup
        r0 = [0,0,150] # [m] Initial position (NED frame)
        v0 = [0,0,0]   # [m/s] Initial velocity (NED frame)
        params.a.n_targs = 4
        params.n_targs_max = params.a.n_targs
        params.n_targs_min = 2
        
        # Configure param dimensions for this number of targets
        reallocate_targ_dims!(params)

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
        params.a.uf_targs = Inf*ones(params.a.nu, params.a.n_targs)

        # Sort by desirability score
        sort_des_score!(params)

        return params
    end

    params = SampleConfig()
    solve(params)
end