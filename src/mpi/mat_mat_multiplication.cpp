/*********************************************************************
Matrix-matrix multiplication
----------------------------

Implement the matrix multiplication algorithm using MPI.

This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#include "io.hpp"
#include "mpi.hpp"
#include <esl/dense.hpp>

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
void init_matrix(esl::Matrix_x<T>& mat, std::size_t offset)
{
	for (std::size_t col = 0; col < mat.cols(); ++col)
		for (std::size_t row = 0; row < mat.rows(); ++row)
			mat(row, col) = static_cast<T>(row + col + offset);
}

template<typename T>
bool is_ab_product(const esl::Matrix_x<T>& mat)
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
			const auto val =
				static_cast<T>(sm1_s_sp1 + s_sm1 * s_col + s_sp1 * row + size * row * s_col);
			if (mat(row, col) != val)
				return false;
		}

	return true;
}

template<typename M>
void mpi_mul_add(esl::Matrix_x<M>& a, const esl::Matrix_x<M>& b, esl::Matrix_x<M>& c)
{
	const auto mpi_size = mpi_comm_size();
	const auto mpi_rank = mpi_comm_rank();
	const auto prev_mpi_rank = (mpi_rank + mpi_size - 1) % mpi_size;
	const auto next_mpi_rank = (mpi_rank + 1) % mpi_size;

	// Buffer for receiving a submatrix of A from the previous processor
	esl::Matrix_x<M> prev_a_;

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
std::optional<bool> check(const esl::Matrix_x<T>& c)
{
	const int check = is_ab_product(c);
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
	esl::Matrix_x<M> a(size, extent.size());

	// Submatrix of the matrix B
	esl::Matrix_x<M> b(size, extent.size());

	// Submatrix of the matrix C
	esl::Matrix_x<M> c(size, extent.size(), 0);

	// Initialize A and B with some fixed values, so that
	// the correctness of the product C = AB can be easily checked
	init_matrix(a, extent.first);
	init_matrix(b, extent.first + 1);

	if (mpi_rank == 0)
	{
		std::cout << "Task sizes distribution: ";
		for (unsigned int r = 0; r < mpi_size; ++r)
			std::cout << alloc_extent(mpi_size, r).size() << ' ';
		std::cout << '\n' << std::endl;
	}

	mpi_mul_add(a, b, c);

	if (const auto ch = check(c); ch)
	{
		if (*ch)
		{
			std::cout << "Matrix C = AB is correct." << std::endl;
			return 0;
		}
		else
		{
			std::cout << "Matrix C = AB is incorrect!" << std::endl;
			return 1;
		}
	}

	return 0;
}
