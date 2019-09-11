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

template<typename T>
void init_a_matrix(Matrix<T>& mat, Extent extent)
{
	for (std::size_t col = 0; col < mat.cols(); ++col)
		for (std::size_t row = 0; row < mat.rows(); ++row)
			mat(row, col) = static_cast<T>(row + col + extent.first);
}

template<typename T>
void init_b_matrix(Matrix<T>& mat, Extent extent)
{
	for (std::size_t col = 0; col < mat.cols(); ++col)
		for (std::size_t row = 0; row < mat.rows(); ++row)
			mat(row, col) = static_cast<T>(1 + row + col + extent.first);
}

template<typename T>
bool check_c_matrix(const Matrix<T>& mat, Extent extent)
{
	const auto size = mat.rows();
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

Extent block_extent(unsigned int mpi_size, unsigned int mpi_rank)
{
	const auto block_size = size / mpi_size;
	return {mpi_rank * block_size, (mpi_rank + 1 < mpi_size) ? (mpi_rank + 1) * block_size : size};
}

int main(int argc, char* argv[])
{
	Mpi mpi(argc, argv);

	const auto mpi_size = mpi_comm_size();
	const auto mpi_rank = mpi_comm_rank();
	const auto prev_mpi_rank = (mpi_rank + mpi_size - 1) % mpi_size;
	const auto next_mpi_rank = (mpi_rank + 1) % mpi_size;

	const auto extent = block_extent(mpi_size, mpi_rank);
	if (mpi_rank == 0)
	{
		std::cout << "Task sizes distribution: ";
		for (unsigned int r = 0; r < mpi_size; ++r)
			std::cout << block_extent(mpi_size, r).size() << ' ';
		std::cout << '\n' << std::endl;
	}

	Matrix<M> a(size, extent.size());
	Matrix<M> b(size, extent.size());
	Matrix<M> c(size, extent.size(), 0);
	Matrix<M> prev_a;

	init_a_matrix(a, extent);
	init_b_matrix(b, extent);

	for (unsigned int r = 0; r < mpi_size; ++r)
	{
		const auto curr_p = (mpi_rank + mpi_size - r) % mpi_size;
		const auto prev_p = (mpi_rank + mpi_size - r - 1) % mpi_size;

		prev_a.resize(size, block_extent(mpi_size, prev_p).size());

#ifdef USE_NON_BLOCKING_COMM
		auto send_status = mpi_isend(a.data(), a.size(), next_mpi_rank);
		auto recv_status = mpi_irecv(prev_a.data(), prev_a.size(), prev_mpi_rank);
#else
		mpi_send(a.data(), a.size(), next_mpi_rank);
		mpi_recv(prev_a.data(), prev_a.size(), prev_mpi_rank);
#endif

		const auto b_row_shift = block_extent(mpi_size, curr_p).first;
		for (std::size_t col = 0; col < extent.size(); ++col)
			for (std::size_t i = 0; i < a.cols(); ++i)
				for (std::size_t row = 0; row < size; ++row)
					c(row, col) += a(row, i) * b(i + b_row_shift, col);

#ifdef USE_NON_BLOCKING_COMM
		mpi_wait(recv_status);
		mpi_wait(send_status);
#endif

		a.swap(prev_a);
	}

	// Check result
	if (mpi_rank == 0)
	{
		std::vector<int> checks(mpi_size);
		checks.front() = check_c_matrix(c, extent);
		mpi_gather_recv(checks.data(), 1);

		if (std::all_of(checks.begin(), checks.end(), [](auto c) { return !!c; }))
		{
			std::cout << "Matrix C = AB is correct." << std::endl;
			return 0;
		}
		else
		{
			std::cout << "Matrix C = AB is incorrect!" << std::endl;
			return -1;
		}
	}
	else
	{
		const int check = check_c_matrix(c, extent);
		mpi_gather_send(&check, 1);
		return 0;
	}
}
