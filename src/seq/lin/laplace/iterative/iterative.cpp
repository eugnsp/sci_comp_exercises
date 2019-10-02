/*********************************************************************
Iterative methods
-----------------

Implement the Jacobi, the Gauss-Seidel and the SOR iterative methods
for determining the solution of a 2D elliptic boundary value problem.

This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#include "../grid.hpp"
#include "io.hpp"
#include "laplace_even_odd_sor_solver.hpp"
#include "laplace_gauss_seidel_solver.hpp"
#include "laplace_jacobi_solver.hpp"
#include "laplace_sor_solver.hpp"
#include "matrix.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

int main()
{
	const auto rhs_fn = [](double x, double y) { return std::cos(2 * x) * sin(2 * y); };
	const auto sol_fn = [](double x, double y) { return .125 * std::cos(2 * x) * std::sin(2 * y) + .1 * y; };

	const auto n_its = 3'000u;
	Grid<double> x_grid{0., 3, 50};
	Grid<double> y_grid{0., 3, 50};

	const auto x_labels = [&x_grid](auto i) { return x_grid[i]; };
	const auto y_labels = [&y_grid](auto i) { return y_grid[i]; };

	const auto true_sol = at_grid_pts(x_grid, y_grid, sol_fn);
	write_gnuplot("solution.dat", true_sol, x_labels, y_labels);

	std::cout << "System size: " << x_grid.n << " x " << y_grid.n << std::endl;

	// Jacobi

	Laplace_jacobi_solver<double> jacobi_solver(x_grid, y_grid, rhs_fn, sol_fn);
	const auto jacobi_res = jacobi_solver.run(n_its, [&](auto it)
	{
		if (it < 1'000)
			write_gnuplot("jacobi_" + std::to_string(it) + ".dat", jacobi_solver.solution(), x_labels, y_labels);
	});
	const auto& jacobi_sol = jacobi_solver.solution();

	std::cout << "Jacobi solver: du = " << norm_sup(jacobi_sol - true_sol) << std::endl;

	// Gauss-Seidel

	Laplace_gauss_seidel_solver<double> gauss_seidel_solver(x_grid, y_grid, rhs_fn, sol_fn);
	const auto gauss_seidel_res = gauss_seidel_solver.run(n_its, [&](auto it)
	{
		if (it < 500)
			write_gnuplot(
				"gauss_seidel_" + std::to_string(it) + ".dat", gauss_seidel_solver.solution(), x_labels, y_labels);
	});
	const auto& gauss_seidel_sol = gauss_seidel_solver.solution();

	std::cout << "Gauss-Seidel solver: du = " << norm_sup(gauss_seidel_sol - true_sol) << std::endl;

	// SOR

	Laplace_sor_solver<double> sor_solver(x_grid, y_grid, rhs_fn, sol_fn);
	const auto sor_res = sor_solver.run(n_its, [&](auto it)
	{
		if (it < 100)
			write_gnuplot("sor_" + std::to_string(it) + ".dat", sor_solver.solution(), x_labels, y_labels);
	});
	const auto& sor_sol = sor_solver.solution();

	std::cout << "SOR solver: du = " << norm_sup(sor_sol - true_sol) << std::endl;

	// Even-odd SOR

	Laplace_even_odd_sor_solver<double> even_odd_sor_solver(x_grid, y_grid, rhs_fn, sol_fn);
	const auto even_odd_sor_res = even_odd_sor_solver.run(n_its);
	const auto& even_odd_sor_sol = even_odd_sor_solver.solution();

	std::cout << "Even-odd SOR solver: du = " << norm_sup(even_odd_sor_sol - true_sol) << std::endl;

	write_vec("convergence.txt", jacobi_res, gauss_seidel_res, sor_res, even_odd_sor_res);
	return 0;
}
