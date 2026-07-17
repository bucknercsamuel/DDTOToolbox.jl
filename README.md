# DDTOToolbox.jl

<video src="https://github.com/user-attachments/assets/601bbc4a-6a70-4a71-8bcc-a6a4f08b8e98" autoplay loop muted playsinline width="100%"></video>

A general-purpose toolbox to construct and solve deferred-decision trajectory optimization (DDTO) problems with high-performance Julia code. These implementations have been used to generate obstacle avoidance maneuvers for quadrotor systems, perform hazard-aware landing tests (with perception-in-the-loop) for quadrotors in the AirSim environment, and perform hardware demonstrations with dynamic ground-based obstacles.

**Features:**  
✔️ A variety of different single-target and deferred-decision trajectory optimizers provided (see below).  
✔️ Customized specification of problem dynamics, constraints, objective and initial guess which can interface with any trajectory optimizer.  
✔️ [JuMP] optimization parsing and batched (parallelized) autodifferentiation using [ForwardDiff] or symbolic differentiation using [SymPy].  
✔️ Closed-loop simulation framework with the Adaptive DDTO algorithm.  
✔️ Support for sysimage and C-interfaceable library generation with [PackageCompiler].  
✔️ ROS2 integration with [AirSim] and [HALSS] perception module (reach out to corresponding authors for more details).  

**Trajectory optimizers supported:**  
✔️ CT-SCvx ([Elango 2025]) → 📁[src/core/opt_sing_scp.jl]  
✔️ Quasiconvex DDTO ([Elango 2022, alg. 2]) → 📁[src/core/opt_ddto_cvx.jl]  
✔️ Lexicographical DDTO ([Elango 2025, alg. 1]) → 📁[src/core/opt_ddto_lex.jl]  
✔️ Adaptive DDTO ([Hayner 2023, alg. 1]) → 📁[src/core/adapt_ddto/algorithm.jl]  
✔️ Graph DDTO (In preparation) → 📁[src/core/opt_ddto_scp.jl]  

## Install



## Run

### Constructing a new problem

To generate a new problem scenario in this toolbox, you must first create a customized parameter object associated with that problem. Currently, we have three different parameter objects associated with three different problem types:
- `DIntegrator2DoFParams`: basic double integrator state-space model in a 2D plane.
- `Quad3DoFCageParams`: quadrotor-specific model in 3D space with cage/arena limits and obstacles.
- `Quad3DoFHaloParams`: quadrotor-specific model in 3D space for autonomous landing scenarios (see the [HALO] framework for more information).

Once a new parameter object is created, you must use it to overload all prototype function templates. Moreover, you must design the following functions for a custom parameter structure `CustomParams`:

```
   prob_constraints(mdl, x, u, params::CustomParams, ref_traj) # JuMP-defined path constraints
   prob_cost(mdl, x, u, params::CustomParams) # JuMP-defined objective function
   dynamics_nonlinear(t, x, ν, params::CustomParams) # nonlinear dynamics function (if nonconvex trajectory optimizer)
   dynamics_linearized(t_ref, x_ref, ν_ref, params::CustomParams) # Linearized dynamics system matrices (if nonconvex trajectory optimizer)
   dynamics_linear(params::CustomParams) # Linear dynamics system matrices (if convex trajectory optimizer)
```
More information on each of these function prototypes is available in 📁[src/core/structs.jl], along with examples for the existing problem types in 📁[src/dint2dof.jl] and 📁[src/quad3dof.jl]. Problems are handled using the mechanisms of time dilation (augmented control) and continuous-time constraint satisfaction (augmented state) as outlined in the [CT-SCvx] seminal paper.

### Running the problem



## Authors

The primary author and maintainer for this repository is [Samuel Buckner] with the [University of Washington Autonomous Controls Laboratory].

[SymPy]: https://www.sympy.org/en/index.html
[AirSim]: https://microsoft.github.io/AirSim/
[JuMP]: https://jump.dev/
[ForwardDiff]: https://github.com/JuliaDiff/ForwardDiff.jl
[PackageCompiler]: https://github.com/Julialang/PackageCompiler.jl
[HALSS]: https://github.com/haynec/HALSS/tree/bafec9ad35f408c2ac76a8239611e45076680eb7
[Elango 2025]: https://www.sciencedirect.com/science/article/abs/pii/S0005109825003589
[Elango 2022, alg. 2]: https://arc.aiaa.org/doi/full/10.2514/6.2022-1583
[Elango 2025, alg. 1]: https://arxiv.org/abs/2502.06623
[Hayner 2023, alg. 1]: https://ieeexplore.ieee.org/abstract/document/10160655
[HALO]: https://haynec.github.io/papers/halo/
[CT-SCvx]: https://www.sciencedirect.com/science/article/abs/pii/S0005109825003589
[src/core/opt_sing_scp.jl]: src/core/opt_sing_scp.jl
[src/core/opt_ddto_cvx.jl]: src/core/opt_ddto_cvx.jl
[src/core/opt_ddto_lex.jl]: src/core/opt_ddto_lex.jl
[src/core/adapt_ddto/algorithm.jl]: src/core/adapt_ddto/algorithm.jl
[src/core/opt_ddto_scp.jl]: src/core/opt_ddto_scp.jl
[src/core/structs.jl]: src/core/structs.jl
[src/dint2dof]: src/dint2dof
[src/quad3dof]: src/quad3dof
[Samuel Buckner]: https://bucknercsamuel.github.io/
[University of Washington Autonomous Controls Laboratory]: https://uwacl.com/