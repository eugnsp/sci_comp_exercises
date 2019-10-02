/*********************************************************************
Jacobi method
-------------

Implement the Jacobi eigenvalue algorithm to find eigenvalues and
eigenvectors of a real symmetric matrix.

This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#include "io.hpp"
#include <es_la/dense.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

template<typename T>
T off(const es_la::Matrix_x<T>& m)
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
void rotate(es_la::Matrix_x<T>& m, const std::size_t row1, const std::size_t col1, const std::size_t row2,
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
	es_la::Matrix_x<T> mat, es_la::Matrix_x<T>& vecs, es_la::Vector_x<T>& vals, const T delta)
{
	assert(mat.rows() == mat.cols());
	const auto n = mat.rows();

	vecs.resize(n, n);
	vecs = 0;
	vecs.diag_view() = 1;

	vals.resize(n);
	vals = mat.diag_view();

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

template<typename T>
bool have_same_elements(es_la::Vector_x<double>& vec1, es_la::Vector_x<double>& vec2, const T delta)
{
	if (vec1.size() != vec2.size())
		return false;

	std::sort(vec1.data(), vec1.data() + vec1.size());
	std::sort(vec2.data(), vec2.data() + vec2.size());

	return std::equal(vec1.data(), vec1.data() + vec1.size(), vec2.data(), vec2.data() + vec2.size(),
		[delta](auto val1, auto val2) { return std::abs(val1 - val2) < delta; });
}

bool test(es_la::Matrix_x<double> mat)
{
	std::cout << "Matrix size: " << mat.rows() << " x " << mat.cols() << std::endl;

	es_la::Matrix_x<double> vecs;
	es_la::Vector_x<double> vals;

	const auto delta = 1e-10;
	const auto iters = jacobi_eigenpairs(mat, vecs, vals, delta);
	std::cout << "The Jacobi diagonalization took " << iters << " iterations." << std::endl;

	es_la::Matrix_x<double> vals_diag(vals.size(), vals.size(), 0);
	vals_diag.diag_view() = vals;
	std::cout << "||M * V - V * D|| = " << norm_sup(mat * vecs - vecs * vals_diag) << '\n';

	es_la::Vector_xd true_vals;
	es_la::eigenvalues(mat, true_vals);

	const auto f = have_same_elements(vals, true_vals, delta);
	std::cout << "Eigenvalues are " << (f ? "" : "in") << "correct.\n" << std::endl;
	return f;
}

int main()
{
	if (!test(es_la::hilbert_matrix<double>(5)))
		return -1;
	if (!test(es_la::hilbert_matrix<double>(10)))
		return -1;

	if (!test(es_la::frank_matrix<double>(5)))
		return -1;
	if (!test(es_la::frank_matrix<double>(10)))
		return -1;
	if (!test(es_la::frank_matrix<double>(20)))
		return -1;

	return 0;
}
