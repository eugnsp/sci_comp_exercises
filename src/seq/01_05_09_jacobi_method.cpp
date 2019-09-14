/*********************************************************************
Jacobi method
-------------
Problems and solutions in scientific computing by W.-H. Steeb et al.
Chapter 5, problem 9

Implement the Jacobi eigenvalue algorithm to find eigenvalues and
eigenvectors of a real symmetric matrix.

This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#include "matrix.hpp"
#include "io.hpp"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

template<typename T>
T off(const Matrix<T>& m)
{
	assert(m.rows() == m.cols());
	const auto n = m.rows();

	T res = 0;
	for (std::size_t row = 0; row < n - 1; ++row)
		for (std::size_t col = row + 1; col < n; ++col)
			res += std::abs(m(row, col));

	return res;
}

template<typename T>
void rotate(Matrix<T>& m, const std::size_t row1, const std::size_t col1, const std::size_t row2,
	const std::size_t col2, const T cos, const T sin)
{
	const auto m1 = m(row1, col1);
	const auto m2 = m(row2, col2);
	m(row1, col1) = m1 * cos - m2 * sin;
	m(row2, col2) = m1 * sin + m2 * cos;
}

// Finds eigenvalues and eigenvectors of a real symmetric matrix using the Jacobi method;
// the algorithm is a simplified version of that described in "Linear algebra" (J.H.Wilkinson)
// and "Numerical recipes" (W.H.Press)
template<typename T>
unsigned int jacobi_eigenpairs(
	Matrix<T>& mat, Matrix<double>& vecs, std::vector<double>& vals, const double delta = 1e-8)
{
	assert(mat.rows() == mat.cols());
	const auto n = mat.rows();

	vecs.resize(n, n);
	vals.resize(n);

	vecs.fill(0);
	for (std::size_t i = 0; i < n; i++)
	{
		vals[i] = mat(i, i);
		vecs(i, i) = 1;
	}

	constexpr unsigned int max_iters = 50u;
	unsigned int iter = 1;
	while (true)
	{
		if (off(mat) < delta * n * n)
			break;

		for (std::size_t row = 0; row < n - 1; row++)
			for (std::size_t col = row + 1; col < n; col++)
			{
				if (std::abs(mat(row, col)) < delta)
					continue;

				const auto theta = (vals[col] - vals[row]) / (2 * mat(row, col));
				const auto t = std::copysign(1 / (std::abs(theta) + std::hypot(1, theta)), theta);
				const auto cos = 1 / std::hypot(1, t);
				const auto sin = t * cos;
				const auto h = t * mat(row, col);

				vals[row] -= h;
				vals[col] += h;

				mat(row, col) = 0;
				for (std::size_t i = 0; i < row; i++)
					rotate(mat, i, row, i, col, cos, sin);
				for (std::size_t i = row + 1; i < col; i++)
					rotate(mat, row, i, i, col, cos, sin);
				for (std::size_t i = col + 1; i < n; i++)
					rotate(mat, row, i, col, i, cos, sin);

				for (std::size_t j = 0; j < n; j++)
					rotate(vecs, j, row, j, col, cos, sin);
			}

		if (++iter >= max_iters)
			throw std::runtime_error("Too many iterations");
	}

	return iter;
}

void test(Matrix<double> mat)
{
	std::cout << "Matrix:\n" << mat << std::endl;

	Matrix<double> vecs;
	std::vector<double> vals;
	const auto iter = jacobi_eigenpairs(mat, vecs, vals);

	std::cout << "After " << iter << " ierations:\n\n";

	std::cout << "Eigenvalues:\n" << std::fixed << std::setprecision(7);

	for (auto v : vals)
		std::cout << v << ' ';
	std::cout << '\n' << std::endl;

	std::cout << "Eigenvectors:\n" << vecs << std::endl;
	std::cout << "-------------\n" << std::endl;
}

int main()
{
	test(hilbert_matrix<double>(4));
	test(frank_matrix<double>(4));

	return 0;
}
