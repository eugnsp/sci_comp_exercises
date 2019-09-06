/*********************************************************************
Jacobi method
-------------
Problems and solutions in scientific computing by W.-H. Steeb et al.
Chapter 9, problem 23

Implement the iterative Jacobi method for determining the solution
of a diagonally dominant system of linear equations.

This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#include "matrix.hpp"
#include "io.hpp"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

template<typename T>
T dist1(const std::vector<T>& x, const std::vector<T>& y)
{
	assert(x.size() == y.size());

	T diff = 0;
	for (std::size_t i = 0; i < x.size(); ++i)
		diff += std::abs(x[i] - y[i]);
	return diff;
}

// Computes y += mat * x
template<typename T>
void mul_add(const Matrix<T>& mat, const std::vector<T>& x, std::vector<T>& y)
{
	assert(mat.rows() == y.size());
	assert(mat.cols() == x.size());

	for (std::size_t j = 0; j < mat.cols(); ++j)
		for (std::size_t i = 0; i < mat.rows(); ++i)
			y[i] += mat(i, j) * x[j];
}

// Computes y -= mat' * x, where mat' is mat with zero diagonal entries
template<typename T>
void nondiag_mul_sub(const Matrix<T>& mat, const std::vector<T>& x, std::vector<T>& y)
{
	assert(mat.rows() == y.size());
	assert(mat.cols() == x.size());

	for (std::size_t j = 0; j < mat.cols(); ++j)
		for (std::size_t i = 0; i < mat.rows(); ++i)
			if (i != j)
				y[i] -= mat(i, j) * x[j];
}

// Solves the linear system mat * sol = rhs using the Jacobi method
template<typename T>
unsigned int jacobi(const Matrix<T>& mat, const std::vector<T>& rhs, std::vector<T>& sol, const double delta = 1e-8)
{
	// To solve A x = b, iterate
	// 		x_{k+1} = B x_k + c,
	// where
	// 		B = -D^{-1} (A - D), c = D^{-1} b
	// with the preconditioner D = diag(A)

	assert(mat.rows() == mat.cols());
	assert(mat.rows() == rhs.size());
	assert(mat.rows() == sol.size());
	const auto n = mat.rows();

	std::vector<T> new_sol;

	constexpr unsigned int max_iters = 1'000u;
	unsigned int iter = 1;
	while (true)
	{
		new_sol = rhs;
		nondiag_mul_sub(mat, sol, new_sol);
		for (std::size_t i = 0; i < n; ++i)
			new_sol[i] /= mat(i, i);

		if (dist1(sol, new_sol) < delta)
			break;

		sol = new_sol;

		if (++iter >= max_iters)
			throw std::runtime_error("Too many iterations");
	}

	return iter;
}

int main()
{
	Matrix<double> mat(4, 4, {7, 1, 1, 2, 1, 6, -2, 3, 1, -2, 7, 4, 2, 3, 4, 6});
	std::vector<double> rhs{1, 2, 3, 4};
	std::vector<double> sol{0, 0, 0, 0};

	std::cout << "System matrix:\n"
			  << mat << '\n'
			  << "Right-hand side:\n"
			  << std::fixed << std::setprecision(4) << rhs << '\n'
			  << "Initial guess:\n"
			  << sol << std::endl;

	const auto iter = jacobi(mat, rhs, sol);

	std::cout << "Solution after " << iter << " iterations:\n" << sol << std::endl;

	std::vector<double> rhs_check{0, 0, 0, 0};
	mul_add(mat, sol, rhs_check);

	std::cout << "Computed right-hand side:\n" << rhs_check << std::endl;

	for (std::size_t i = 0; i < 4; ++i)
		rhs_check[i] -= rhs[i];

	std::cout << "Errors:\n" << std::scientific << rhs_check << std::endl;

	return 0;
}
