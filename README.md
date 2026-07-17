# DDTOToolbox.jl

<video src="https://github.com/user-attachments/assets/601bbc4a-6a70-4a71-8bcc-a6a4f08b8e98" autoplay loop muted playsinline width="100%"></video>

A general-purpose toolbox to construct and solve deferred-decision trajectory optimization (DDTO) problems with high-performance Julia code. These implementations have been used to generate obstacle avoidance maneuvers for quadrotor systems, perform hazard-aware landing tests (with perception-in-the-loop) for quadrotors in the AirSim environment, and perform hardware demonstrations with dynamic ground-based obstacles.

Features:
✔️ A variety of different single-target and deferred-decision trajectory optimizers provided (see below).
✔️ Customized specification of problem dynamics, constraints, objective and initial guess which can interface with any trajectory optimizer.
✔️ [JuMP] optimization parsing and batched (parallelized) autodifferentiation using [ForwardDiff].
✔️ Closed-loop simulation with the Adaptive DDTO algorithm.
✔️ Support for sysimage and C-interfaceable library generation with [PackageCompiler].
✔️ ROS2 integration with AirSim (reach out to corresponding authors for more detail).

Trajectory optimizers supported:
✔️ CT-SCvx ([Elango 2025]) → 📁[src/core/opt_sing_scp.jl]
✔️ Quasiconvex DDTO ([Elango 2022, alg. 2]) → 📁[src/core/opt_ddto_cvx.jl]
✔️ Lexicographical DDTO ([Elango 2025, alg. 1]) → 📁[src/core/opt_ddto_lex.jl]
✔️ Adaptive DDTO ([Hayner 2023, alg. 1]) → 📁[src/core/adapt_ddto/algorithm.jl]
✔️ Graph DDTO (In preparation) → 📁[src/core/opt_ddto_scp.jl]



[JuMP]: https://jump.dev/
[ForwardDiff]: https://github.com/JuliaDiff/ForwardDiff.jl
[PackageCompiler]: https://github.com/Julialang/PackageCompiler.jl
[Elango 2025]: https://www.sciencedirect.com/science/article/abs/pii/S0005109825003589
[Elango 2022, alg. 2]: https://arc.aiaa.org/doi/full/10.2514/6.2022-1583
[Elango 2025, alg. 1]: https://arxiv.org/abs/2502.06623
[Hayner 2023, alg. 1]: https://ieeexplore.ieee.org/abstract/document/10160655
[src/core/opt_sing_scp.jl]: src/core/opt_sing_scp.jl
[src/core/opt_ddto_cvx.jl]: src/core/opt_ddto_cvx.jl
[src/core/opt_ddto_lex.jl]: src/core/opt_ddto_lex.jl
[src/core/adapt_ddto/algorithm.jl]: src/core/adapt_ddto/algorithm.jl
[src/core/opt_ddto_scp.jl]: src/core/opt_ddto_scp.jl