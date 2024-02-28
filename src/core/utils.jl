#= LCvx for Quadcopter Landing -- Extra Utility Functions.

Author: Samuel Buckner (UW-ACL)
=#

function solver_setup(solver::String)
    type = ""
    if solver == "ECOS"
        mdl = Model(optimizer_with_attributes(ECOS.Optimizer, 
            "verbose" => 0, 
            "max_iters" => 1000))
        type = "SOCP"
    elseif solver == "MOSEK"
        mdl = Model(Mosek.Optimizer)
        JuMP.set_optimizer_attribute(mdl, "LOG",  0) # disable debugging
        JuMP.set_optimizer_attribute(mdl, "MAX_NUM_WARNINGS", 0) # disable warnings
        type = "SOCP"
    elseif solver == "OSQP"
        # mdl = Model(OSQP.Optimizer)
        # JuMP.set_optimizer_attribute(mdl, "LOG",  0) # disable debugging
        # JuMP.set_optimizer_attribute(mdl, "MAX_NUM_WARNINGS", 0) # disable warnings
        mdl = Model(optimizer_with_attributes(OSQP.Optimizer,
            "verbose" => 0, 
            "max_iters" => 1000, 
            "eps_abs" => 1e-8,
            "eps_rel" => 1e-8))
        type = "QP"
    else
        error("solver choice is invalid.")
    end
    return mdl, type
end

function time_dilation_control_to_wall_clock_time(∂t_∂τ::Vector, τ::Vector, disc::Int)
    # Converts time-dilation control to wall-clock time based on discretization method
    Δτ = diff(τ)
    if length(∂t_∂τ) > 1
        Δt = Vector(undef,length(∂t_∂τ)-1)
        for k=1:length(Δt)
            if disc == 0
                Δt[k] = Δτ[k] * ∂t_∂τ[k]
            elseif disc == 1
                Δt[k] = (1/2) * Δτ[k] * (∂t_∂τ[k] + ∂t_∂τ[k+1])
            end
        end
        t = cumsum([0.;Δt])
    else
        t = [0.]
    end
    return t
end

function wall_clock_time_to_time_dilation_control(t::Vector, τ::Vector, disc::Int)
    # Converts wall-clock time to time dilation control based on discretization method
    N = length(t)
    if disc == 0
        N_ctrl = N-1
    elseif disc == 1
        N_ctrl = N
    end
    Δt = diff(t)
    Δτ = diff(τ)
    ∂t_∂τ = Vector(undef,N_ctrl)
    if disc == 0
        for k=1:N_ctrl
            ∂t_∂τ[k] = Δt[k] / Δτ[k]
        end
    elseif disc == 1
        n = length(Δt)
        ∂t_∂τ[1] = sum([Δt[k] / Δτ[k] for k=1:n])/length(Δt) # boundary condition chosen for numerical properties, but is technically arbitrary!
        for k=1:N_ctrl-1
            ∂t_∂τ[k+1] = 2 * Δt[k] / Δτ[k] - ∂t_∂τ[k]
        end
    end
    return ∂t_∂τ
end

"""
Compute Jacobians with direct numerical differentiation (simple central differencing) on nonlinear dynamics
"""
function numerical_jacobian(t_ref, x_ref, u_ref, dyn_nl; pert=1e-4)
    # Setup
    if pert < 1e-10
        error("Required perturbation is too small")
    end
    nx = length(x_ref)
    nu = length(u_ref)
    A = zeros(nx,nx)
    B = zeros(nx,nu)

    # Numerical A
    pertI = pert*I(nx)
    for k=1:nx
        fp = dyn_nl(t_ref, x_ref + pertI[:,k], u_ref)
        fm = dyn_nl(t_ref, x_ref - pertI[:,k], u_ref)
        A[:,k] = (fp-fm) / (2*pert)
    end

    # Numerical B
    pertI = pert*I(nu)
    for k=1:nu
        fp = dyn_nl(t_ref, x_ref, u_ref + pertI[:,k])
        fm = dyn_nl(t_ref, x_ref, u_ref - pertI[:,k])
        B[:,k] = (fp-fm) / (2*pert)
    end

    # Numerical z
    z = dyn_nl(t_ref, x_ref, u_ref) - (A*x_ref + B*u_ref)

    return A,B,z
end

function process_sympy_string(str::String)::String
    str = replace(str, "{}" => ",") # apply commas for matrices (simple hack due to symbols() limitations)
    str = replace(str, "Max" => "max") # necessary for julia convention
    str = replace(str, "Heaviside" => "heaviside") # necessary for julia convention
    return str
end