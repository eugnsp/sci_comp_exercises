# Finite-elements methods

Implement
- the Jacobi,
- the Gauss-Seidel, and
- the SOR

iterative methods for determining the solution of a 2D elliptic boundary value problem

<code>-&Delta;u(x, y) = cos(2x) sin(2y)</code>

in the domain <code>&Omega; = [0, 3] &times; [0, 3]</code> with Dirichlet boundary condition

<code>f(x, y) = 1/8 cos(2x) sin(2y) + 1/10 y</code> at <code>&part;&Omega;</code>.

## Equation solution

Standard FEM with <code>P<sub>4</sub></code> Lagrange elements on coarse mesh:
![Equation solution](/figs/seq/lin/poisson_2d_fem/std.png)

Mixed FEM with <code>P<sub>1</sub></code> Lagrange + <code>DP<sub>0</sub></code> discontinuous Lagrange elements on fine mesh:
![Equation solution](/figs/seq/lin/poisson_2d_fem/mixed.png)
