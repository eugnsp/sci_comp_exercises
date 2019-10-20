# Relaxation methods

Implement
- the Jacobi,
- the Gauss-Seidel, and
- the SOR

iterative methods for determining the solution of a 2D elliptic boundary value problem

<code>-&Delta;u(x, y) = cos(2x) sin(2y)</code>

in the domain <code>&Omega; = [0, 3] &times; [0, 3]</code> with Dirichlet boundary condition

<code>f(x, y) = 1/8 cos(2x) sin(2y) + 1/10 y</code> at <code>&part;&Omega;</code>.

## Equation solution

![Equation solution](/figs/seq/lin/laplace_2d_fdm/relaxation/solution.png)

## Convergence

![Convergence](/figs/seq/lin/laplace_2d_fdm/relaxation/convergence.png)

## Animations

* [The Jacobi method](https://www.youtube.com/watch?v=cEU6AxI-Z_U)
* [The Gauss-Seidel method](https://www.youtube.com/watch?v=YyGgVBPh24Y)
* [The SOR method](https://www.youtube.com/watch?v=Chp4Sif4PN8)
