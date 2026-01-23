"""
In this script, I investigate a problem scenario that is known to produce a non-convergent SCvx solution (single-target).
This is investigated for understanding fundamental algorithm improvements.
"""

using DDTOToolbox
include("plots.jl")

function FailedParams()::Quad3DoFHaloParams{CReal,Int}
    # Load default params
    params = Quad3DoFHaloParams()

    # # Configure for single-target problem
    # params.n_targs_min = 1
    # params.n_targs_max = 1
    params.n_targs_min = 2
    params.n_targs_max = 4

    # Resize all parameters that depend on number of max targets
    reallocate_targ_dims!(params)

    # Set boundary conditions for single target
    params.a.z0 = [
        -9.890073931444503
        33.23677487980699
        74.59426231468457
        -3.2592702750600018
         3.0760868384043905
        -2.6733755344933194
         0.0
        Inf
    ]
    params.a.u0 = [
        -3.2481875711891455
        -0.9208103062020808
         7.083114953349579
        Inf
    ]
    params.a.zf_targs = reshape([
        -91.2985   5.12139   9.29315  -13.1941;
        90.8214  17.4355   -5.5555    -8.09507;
         1.0      1.0       1.0        1.0;
         0.0      0.0       0.0        0.0;
         0.0      0.0       0.0        0.0;
         0.0      0.0       0.0        0.0;
        Inf      Inf       Inf        Inf
    ], params.a.nx,params.n_targs_max)

    return params
end

# Test
params = FailedParams()
solve(params)
;