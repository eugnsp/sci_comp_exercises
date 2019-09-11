/*********************************************************************
MPI matrix-matrix multiplication
--------------------------------
Parallel algorithms by H.Casanova et al.
Exercise 4.4

Implement the matrix multiplication algorithm from Sec. 4.2 using MPI.

This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#include "io.hpp"
#include "matrix.hpp"
#include "mpi.hpp"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iostream>
#include <optional>
#include <vector>

#define USE_NON_BLOCKING_COMM
inline constexpr std::size_t size = 2'000;
using M = unsigned long long;

struct Extent
{
	std::size_t first;
	std::size_t last;

	std::size_t size() const
	{
		assert(first <= last);
		return last - first;
	}
};

// Return the row extent for a matrix allocated in the given processor `mpi_rank`
Extent alloc_extent(unsigned int mpi_size, unsigned int mpi_rank)
{
	const auto block_size = size / mpi_size;
	return {mpi_rank * block_size, (mpi_rank + 1 < mpi_size) ? (mpi_rank + 1) * block_size : size};
}

Extent alloc_extent()
{
	return alloc_extent(mpi_comm_size(), mpi_comm_rank());
}

template<typename T>
void init_a_matrix(Matrix<T>& mat)
{
	const auto extent = alloc_extent();
	for (std::size_t col = 0; col < mat.cols(); ++col)
		for (std::size_t row = 0; row < mat.rows(); ++row)
			mat(row, col) = static_cast<T>(row + col + extent.first);
}

template<typename T>
void init_b_matrix(Matrix<T>& mat)
{
	const auto extent = alloc_extent();
	for (std::size_t col = 0; col < mat.cols(); ++col)
		for (std::size_t row = 0; row < mat.rows(); ++row)
			mat(row, col) = static_cast<T>(1 + row + col + extent.first);
}

template<typename T>
bool check_c_matrix(const Matrix<T>& mat)
{
	const auto extent = alloc_extent();
	const auto s_sp1 = size * (size + 1) / 2;
	const auto s_sm1 = (size - 1) * size / 2;
	const auto sm1_s_sp1 = (size - 1) * size * (size + 1) / 3;

	if (mat.rows() != size || mat.cols() != extent.size())
		return false;

	for (std::size_t col = 0; col < mat.cols(); ++col)
		for (std::size_t row = 0; row < mat.rows(); ++row)
		{
			const auto s_col = col + extent.first;
			const auto val = static_cast<T>(sm1_s_sp1 + s_sm1 * s_col + s_sp1 * row + size * row * s_col);
			if (mat(row, col) != val)
				return false;
		}

	return true;
}

template<typename M>
void matrix_mul(Matrix<M>& a, const Matrix<M>& b, Matrix<M>& c)
{
	const auto mpi_size = mpi_comm_size();
	const auto mpi_rank = mpi_comm_rank();
	const auto prev_mpi_rank = (mpi_rank + mpi_size - 1) % mpi_size;
	const auto next_mpi_rank = (mpi_rank + 1) % mpi_size;

	// Buffer for receiving a submatrix of A from the previous processor
	Matrix<M> prev_a_;

	for (unsigned int r = 0; r < mpi_size; ++r)
	{
		const auto curr_p = (mpi_rank + mpi_size - r) % mpi_size;
		const auto prev_p = (mpi_rank + mpi_size - r - 1) % mpi_size;

		prev_a_.resize(size, alloc_extent(mpi_size, prev_p).size());

#ifdef USE_NON_BLOCKING_COMM
		auto send_status = mpi_isend(a.data(), a.size(), next_mpi_rank);
		auto recv_status = mpi_irecv(prev_a_.data(), prev_a_.size(), prev_mpi_rank);
#else
		mpi_send(a.data(), a.size(), next_mpi_rank);
		mpi_recv(prev_a.data(), prev_a.size(), prev_mpi_rank);
#endif

		const auto b_row_shift = alloc_extent(mpi_size, curr_p).first;
		for (std::size_t col = 0; col < c.cols(); ++col)
			for (std::size_t i = 0; i < a.cols(); ++i)
				for (std::size_t row = 0; row < size; ++row)
					c(row, col) += a(row, i) * b(i + b_row_shift, col);

#ifdef USE_NON_BLOCKING_COMM
		mpi_wait(recv_status);
		mpi_wait(send_status);
#endif

		a.swap(prev_a_);
	}
}

template<typename T>
std::optional<bool> check_result(const Matrix<T>& c)
{
	const int check = check_c_matrix(c);
	if (mpi_comm_rank() == 0)
	{
		std::vector<int> checks(mpi_comm_size());
		checks.front() = check;
		mpi_gather_recv(checks.data(), 1);
		return std::all_of(checks.begin(), checks.end(), [](auto c) { return !!c; });
	}
	else
	{
		mpi_gather_send(&check, 1);
		return {};
	}
}

int main(int argc, char* argv[])
{
	Mpi mpi(argc, argv);

	const auto mpi_size = mpi_comm_size();
	const auto mpi_rank = mpi_comm_rank();
	const auto extent = alloc_extent(mpi_size, mpi_rank);

	// Submatrix of the matrix A
	Matrix<M> a(size, extent.size());

	// Submatrix of the matrix B
	Matrix<M> b(size, extent.size());

	// Submatrix of the matrix C
	Matrix<M> c(size, extent.size(), 0);

	init_a_matrix(a);
	init_b_matrix(b);

	if (mpi_rank == 0)
	{
		std::cout << "Task sizes distribution: ";
		for (unsigned int r = 0; r < mpi_size; ++r)
			std::cout << alloc_extent(mpi_size, r).size() << ' ';
		std::cout << '\n' << std::endl;
	}

	matrix_mul(a, b, c);

	if (auto check = check_result(c); check)
	{
		if (*check)
			std::cout << "Matrix C = AB is correct." << std::endl;
		else
			std::cout << "Matrix C = AB is incorrect!" << std::endl;

		if (!*check)
			return 1;
	}

	return 0;
}
