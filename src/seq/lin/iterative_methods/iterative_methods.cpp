/*********************************************************************
Iterative methods
-----------------

Implement the Jacobi, the Gauss-Seidel and the SOR iterative methods
for determining the solution of a 2D elliptic boundary value problem.

This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#include "grid.hpp"
#include "io.hpp"
#include "laplace_jacobi_solver.hpp"
#include "laplace_gauss_seidel_solver.hpp"
#include "laplace_sor_solver.hpp"
#include "matrix.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <utility>
#include <vector>

template<typename T>
T dist_inf(const Matrix<T>& x, const Matrix<T>& y)
{
	assert(x.rows() == y.rows());
	assert(x.cols() == y.cols());

	T diff = 0;
	for (std::size_t col = 0; col < x.cols(); ++col)
		for (std::size_t row = 0; row < x.rows(); ++row)
			diff = std::max(diff, std::abs(x(row, col) - y(row, col)));
	return diff;
}

int main()
{
	const auto rhs_fn = [](double x, double y) { return std::cos(2 * x) * sin(2 * y); };
	const auto sol_fn = [](double x, double y) { return (1. / 8.) * std::cos(2 * x) * std::sin(2 * y); };

	const auto n_its = 5'000u;
	Grid<double> x_grid{0., 3, 50};
	Grid<double> y_grid{0., 3, 50};

	const auto x_labels = [&x_grid](auto i) { return x_grid[i]; };
	const auto y_labels = [&y_grid](auto i) { return y_grid[i]; };
	const auto true_sol = at_grid_pts(x_grid, y_grid, sol_fn);
	write_gnuplot("solution.dat", true_sol, x_labels, y_labels);

	std::cout << "System size: " << x_grid.n << " x " << y_grid.n << std::endl;

	// Jacobi

	Laplace_jacobi_solver<double> jacobi_solver(x_grid, y_grid, rhs_fn, sol_fn);
	const auto jacobi_du = jacobi_solver.run(n_its);
	const auto& jacobi_sol = jacobi_solver.solution();

	std::cout << "Jacobi solver: du = " << dist_inf(jacobi_sol, true_sol) << std::endl;
	write_gnuplot("jacobi.dat", jacobi_sol, x_labels, y_labels);

	// Gauss-Seidel

	Laplace_gauss_seidel_solver<double> gauss_seidel_solver(x_grid, y_grid, rhs_fn, sol_fn);
	const auto gauss_seidel_du = gauss_seidel_solver.run(n_its);
	const auto& gauss_seidel_sol = gauss_seidel_solver.solution();

	std::cout << "Gauss-Seidel solver: du = " << dist_inf(gauss_seidel_sol, true_sol) << std::endl;
	write_gnuplot("gauss_seidel.dat", gauss_seidel_sol, x_labels, y_labels);

	// SOR

	Laplace_sor_solver<double> sor_solver(x_grid, y_grid, rhs_fn, sol_fn);
	const auto sor_du = sor_solver.run(n_its);
	const auto& sor_sol = sor_solver.solution();

	std::cout << "SOR solver: du = " << dist_inf(sor_sol, true_sol) << std::endl;
	write_gnuplot("sor.dat", sor_sol, x_labels, y_labels);

	write_vec("convergence.txt", jacobi_du, gauss_seidel_du, sor_du);

	return 0;
}
