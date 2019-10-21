/*********************************************************************
Iterative methods
-----------------

Implement the Jacobi, the Gauss-Seidel and the SOR iterative methods
for determining the solution of a 2D elliptic boundary value problem.

This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#include "../grid.hpp"
#include "../relaxation/laplace_even_odd_sor_solver.hpp"
#include "io.hpp"
#include "laplace_cg_solver.hpp"

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
	const auto sol_fn = [](double x, double y) {
		return .125 * std::cos(2 * x) * std::sin(2 * y) + .1 * y;
	};

	const auto n_its = 350u;
	Grid<double> x_grid{0., 3, 50};
	Grid<double> y_grid{0., 3, 50};

	const auto true_sol = fn_sample(x_grid, y_grid, sol_fn);
	write_gnuplot("solution.dat", true_sol, x_grid, y_grid);

	std::cout << "System size: " << x_grid.n << " x " << y_grid.n << std::endl;

	// Even-odd SOR (as a base line)

	Laplace_even_odd_sor_solver<double> even_odd_sor_solver(x_grid, y_grid, rhs_fn, sol_fn);
	const auto even_odd_sor_res = even_odd_sor_solver.run(n_its);
	const auto& even_odd_sor_sol = even_odd_sor_solver.solution();

	std::cout << "Even-odd SOR solver: du = " << norm_sup(even_odd_sor_sol - true_sol) << std::endl;

	// Conjugate gradient

	Laplace_cg_solver<double> cg_solver(x_grid, y_grid, rhs_fn, sol_fn);
	const auto cg_res = cg_solver.run(n_its, [&, it = 0]() mutable {
		if (it < 100)
			write_gnuplot(
				"cg_" + std::to_string(it++) + ".dat", cg_solver.solution(), x_grid, y_grid);
	});
	const auto& cg_sol = cg_solver.solution();

	std::cout << "Conjugate gradient solver: du = " << norm_sup(cg_sol - true_sol) << std::endl;

	write_vec("convergence.txt", even_odd_sor_res, cg_res);
	return 0;
}
