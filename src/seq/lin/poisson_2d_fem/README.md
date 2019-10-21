# Finite-elements methods

Implement
- the Jacobi,
- the Gauss-Seidel, and
- the SOR

iterative methods for determining the solution of a 2D elliptic boundary value problem

<code>-&Delta;u(x, y) = cos(2x) sin(2y)</code>

in the domain <code>&Omega; = [0, 3] &times; [0, 3]</code> with the Dirichlet boundary conditions `u = 0` at <code>&part;&Omega; &cap; {(x, 0)}</code>, `u = .25` at <code>&part;&Omega; &cap; {(x, 3)}</code>, and zero Neumann boundary condition at <code>&part;&Omega; &cap; ({(0, y)} &cup; {(3, y)})</code>.

## Equation solution

Standard FEM with <code>P<sub>4</sub></code> Lagrange elements on a coarse mesh:
![Equation solution](/figs/seq/lin/poisson_2d_fem/std.png)

Mixed FEM with <code>P<sub>1</sub></code> Lagrange + <code>DP<sub>0</sub></code> discontinuous Lagrange elements on a fine mesh:
![Equation solution](/figs/seq/lin/poisson_2d_fem/mixed.png)
