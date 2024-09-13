# minvol_warren2d
This folder contains work to perform a series of optimizations for a warren truss structure:

$$
\begin{align}
    \min_{x} \quad & V = \rho A^TL\\
  \text{s.t.} \quad & \lvert d_{yi}\rvert \leq d_{\max} \quad \forall i \in 1,\dots, n_n \nonumber\\
  & \lvert\sigma_i \rvert \leq \sigma_{\max} \quad \forall i \in 1,\dots,n_e \nonumber
\end{align}
$$

With constant parameters and initial values defined in `init_problem.jl`.

## Optimization strategy study
File names indicate the type of optimization and the algorithm used: `type_alg.jl`

## Types
- `gb`: gradient based. Our method using AD and custom adjoints.
- `na`: no adjoint. Using out-of-the-box AD.
- `gf`: gradient free.
- `fd`: finite difference (+ a gradient-based optimization problem)

## Algorithms
All optimization is performed using `Nonconvex.jl` and subpackages.

From `NonconvexNLopt.jl`
- `mma`: Method of moving asymptotes (MMA)
- `cobyla`: Constrained Optimization BY Linear Approximation (COBYLA)

From `NonconvexMetaheuristics.jl`
- `ga`: Genetic algorithm (population size of 100 used)

### Stopping criteria
- Maximum number of function evaluations: 25000
- Objective value relative tolerance of 1E-6
- Maximum time of 120s