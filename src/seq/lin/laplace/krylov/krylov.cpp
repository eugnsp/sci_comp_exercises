/*********************************************************************
Iterative methods
-----------------

Implement the Jacobi, the Gauss-Seidel and the SOR iterative methods
for determining the solution of a 2D elliptic boundary value problem.

This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#include "../grid.hpp"
#include "io.hpp"
#include "laplace_cg_solver.hpp"
#include "matrix.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

template<typename T>
T dist_inf(const Matrix<T>& x, const Matrix<T>& y)
{
	assert(x.rows() == y.rows());
	assert(x.cols() == y.cols());

	T dist = 0;
	for (std::size_t col = 0; col < x.cols(); ++col)
		for (std::size_t row = 0; row < x.rows(); ++row)
			dist = std::max(dist, std::abs(x(row, col) - y(row, col)));
	return dist;
}

int main()
{
	const auto rhs_fn = [](double x, double y) { return std::cos(2 * x) * sin(2 * y); };
	const auto sol_fn = [](double x, double y) { return .125 * std::cos(2 * x) * std::sin(2 * y) + .1 * y; };

	const auto n_its = 200u;
	Grid<double> x_grid{0., 3, 50};
	Grid<double> y_grid{0., 3, 50};

	const auto x_labels = [&x_grid](auto i) { return x_grid[i]; };
	const auto y_labels = [&y_grid](auto i) { return y_grid[i]; };

	const auto true_sol = at_grid_pts(x_grid, y_grid, sol_fn);
	write_gnuplot("solution.dat", true_sol, x_labels, y_labels);

	std::cout << "System size: " << x_grid.n << " x " << y_grid.n << std::endl;

	// Conjugate gradient

	Laplace_cg_solver<double> cg_solver(x_grid, y_grid, rhs_fn, sol_fn);
	const auto cg_res = cg_solver.run(n_its);
	const auto& cg_sol = cg_solver.solution();

	std::cout << "Conjugate gradient solver: du = " << dist_inf(cg_sol, true_sol) << std::endl;

	write_vec("convergence.txt", cg_res);

	// Animation of conjugate gradient

	Laplace_cg_solver<double> cg_anim_solver(x_grid, y_grid, rhs_fn, sol_fn);
	cg_anim_solver.run(100, [&](auto it)
	{;
		write_gnuplot("cg_" + std::to_string(it) + ".dat", cg_anim_solver.solution(), x_labels, y_labels);
	});

	return 0;
}
