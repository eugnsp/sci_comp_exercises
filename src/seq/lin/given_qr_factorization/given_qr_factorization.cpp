// This file is covered by the LICENSE file in the root of this project.

#include "io.hpp"

#include <esl/dense.hpp>
#include <esl/io.hpp>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <utility>

template<typename T>
bool is_eq(const esl::Matrix_x<T>& mat_a, const esl::Matrix_x<T>& mat_b, T tol = 1e-6)
{
	if (mat_a.rows() != mat_b.rows() || mat_a.cols() != mat_b.cols())
		return false;

	for (std::size_t col = 0; col < mat_a.cols(); ++col)
		for (std::size_t row = 0; row < mat_a.rows(); ++row)
			if (std::abs(mat_a(row, col) - mat_b(row, col)) > tol)
				return false;

	return true;
}

template<typename T>
int sign(T val)
{
	return (0 < val) - (val < 0);
}

template<typename T>
std::pair<T, T> givens(const T a, const T b)
{
	assert(b != 0);

	T cos, sin;
	if (std::abs(a) < std::abs(b))
	{
		const auto tau = a / b;
		sin = 1 / std::sqrt(1 + tau * tau);
		cos = sin * tau;
	}
	else
	{
		const auto tau = b / a;
		cos = 1 / std::sqrt(1 + tau * tau);
		sin = cos * tau;
	}
	return {cos, sin};
}

template<typename T>
T rho(const T cos, const T sin)
{
	if (cos == 0)
		return 1;
	if (std::abs(sin) < std::abs(cos))
		return sign(cos) * sin / 2;
	else
		return sign(sin) * 2 / cos;
}

template<typename T>
std::pair<T, T> inv_rho(const T rho)
{
	if (rho == 1)
		return {0, 1};
	if (std::abs(rho) < 1)
	{
		const auto sin = 2 * rho;
		return {std::sqrt(1 - sin * sin), sin};
	}
	else
	{
		const auto cos = 2 / rho;
		return {cos, std::sqrt(1 - cos * cos)};
	}
}

template<typename T>
void rotate(esl::Matrix_x<T>& m, const std::size_t row1, const std::size_t col1,
	const std::size_t row2, const std::size_t col2, const T cos, const T sin)
{
	const auto m1 = m(row1, col1);
	const auto m2 = m(row2, col2);
	m(row1, col1) = m1 * cos + m2 * sin;
	m(row2, col2) = m2 * cos - m1 * sin;
}

// Computes the QR factorization of the given matrix `a`, the matrix R
// is stored in the upper triangular part of `a`, the matrix Q can be
// computed from the lower triangular part
template<typename T>
void givens_qr_factorize(esl::Matrix_x<T>& qr)
{
	for (std::size_t col = 0; col < qr.cols(); ++col)
		for (std::size_t row = qr.rows() - 1; row > col; --row)
		{
			if (qr(row, col) == 0)
				continue;

			const auto [cos, sin] = givens(qr(row - 1, col), qr(row, col));
			for (std::size_t i = col; i < qr.cols(); ++i)
				rotate(qr, row - 1, i, row, i, cos, sin);

			qr(row, col) = rho(cos, sin);
		}
}

template<typename T>
esl::Matrix_x<T> givens_qr_get_q(const esl::Matrix_x<T>& qr)
{
	esl::Matrix_x<T> q = esl::id_matrix<T>(qr.rows());
	for (std::size_t col = 0; col < qr.cols(); ++col)
		for (std::size_t row = qr.rows() - 1; row > col; --row)
		{
			const auto [cos, sin] = inv_rho(qr(row, col));
			for (std::size_t i = 0; i < qr.rows(); ++i)
				rotate(q, i, row - 1, i, row, cos, sin);
		}

	return q;
}

template<typename T>
void test(esl::Matrix_x<T> a)
{
	std::cout << std::fixed;
	std::cout << "Matrix:\n" << printer(a, 5) << std::endl;

	auto a0 = a;
	givens_qr_factorize(a);
	auto q = givens_qr_get_q(a);

	// After computing Q, zero lower triangular part of A
	for (std::size_t col = 0; col < a.cols(); ++col)
		for (std::size_t row = col + 1; row < a.rows(); ++row)
			a(row, col) = 0;

	const auto qr = (q * a).eval();
	std::cout << "Q:\n" << printer(q, 5) << std::endl;
	std::cout << "R:\n" << printer(a, 5) << std::endl;
	std::cout << "Q R:\n" << printer(qr, 5) << std::endl;

	if (is_eq(a0, qr))
		std::cout << "QR is correct.\n";
	else
		std::cout << "QR is incorrect!\n";

	std::cout << "-------------\n" << std::endl;
}

int main()
{
	test(esl::hilbert_matrix<double>(5));
	test(esl::hilbert_matrix<double>(4));

	test(esl::random_real_matrix(4));
	test(esl::random_real_matrix(5));
	test(esl::random_real_matrix(6));

	return 0;
}
