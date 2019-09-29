/*********************************************************************
Matrix class and auxiliary routines
-----------------------------------

This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#pragma once
#include <cassert>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <random>
#include <utility>
#include <vector>

template<typename T>
class Matrix
{
public:
	using Type = T;

	using Container = std::vector<Type>;
	using Reference = typename Container::reference;
	using Const_reference = typename Container::const_reference;

public:
	Matrix() = default;

	Matrix(std::size_t rows, std::size_t cols) : data_(rows * cols), rows_(rows), cols_(cols)
	{}

	Matrix(std::size_t rows, std::size_t cols, const Type& value) : Matrix(rows, cols)
	{
		fill(value);
	}

	Matrix(std::size_t rows, std::size_t cols, std::initializer_list<Type> init) : data_(init), rows_(rows), cols_(cols)
	{}

	Matrix(const Matrix&) = default;
	Matrix(Matrix&&) = default;

	Matrix& operator=(const Matrix&) = default;
	Matrix& operator=(Matrix&&) = default;

	Reference operator()(std::size_t row, std::size_t col)
	{
		assert(row < rows_ && col < cols_);
		return data_[row + col * rows_];
	}

	Const_reference operator()(std::size_t row, std::size_t col) const
	{
		assert(row < rows_ && col < cols_);
		return data_[row + col * rows_];
	}

	Type* data()
	{
		return data_.data();
	}

	const Type* data() const
	{
		return data_.data();
	}

	std::size_t rows() const
	{
		return rows_;
	}

	std::size_t cols() const
	{
		return cols_;
	}

	std::size_t size() const
	{
		return rows_ * cols_;
	}

	void resize(std::size_t rows, std::size_t cols)
	{
		rows_ = rows;
		cols_ = cols;
		data_.resize(rows_ * cols_);
	}

	void fill(const Type& value)
	{
		data_.assign(data_.size(), value);
	}

	void swap(Matrix& other)
	{
		std::swap(data_, other.data_);
		std::swap(rows_, other.rows_);
		std::swap(cols_, other.cols_);
	}

private:
	Container data_;
	std::size_t rows_ = 0;
	std::size_t cols_ = 0;
};

template<typename T>
void transpose(Matrix<T>& mat)
{
	assert(mat.rows() == mat.cols());
	const auto n = mat.rows();

	for (std::size_t row = 1; row < n; ++row)
		for (std::size_t col = 0; col < row; ++col)
			std::swap(mat(row, col), mat(col, row));
}

// Computes C += A * B
template<typename T>
void mul_add(const Matrix<T>& mat_a, const Matrix<T>& mat_b, Matrix<T>& mat_c)
{
	assert(mat_a.rows() == mat_c.rows());
	assert(mat_a.cols() == mat_b.rows());
	assert(mat_b.cols() == mat_c.cols());

	for (std::size_t col = 0; col < mat_c.cols(); ++col)
		for (std::size_t i = 0; i < mat_a.cols(); ++i)
			for (std::size_t row = 0; row < mat_c.rows(); ++row)
				mat_c(row, col) += mat_a(row, i) * mat_b(i, col);
}

template<typename T>
bool is_eq(const Matrix<T>& mat_a, const Matrix<T>& mat_b, T tol = 1e-6)
{
	if (mat_a.rows() != mat_b.rows() || mat_a.cols() != mat_b.cols())
		return false;

	for (std::size_t col = 0; col < mat_a.cols(); ++col)
		for (std::size_t row = 0; row < mat_a.rows(); ++row)
			if (std::abs(mat_a(row, col) - mat_b(row, col)) > tol)
				return false;

	return true;
}

// Returns the identity matrix
template<typename T>
Matrix<T> id_matrix(const std::size_t rows, const std::size_t cols)
{
	Matrix<T> mat(rows, cols);
	for (std::size_t col = 0; col < cols; ++col)
		for (std::size_t row = 0; row < rows; ++row)
			mat(row, col) = (row == col) ? T{1} : T{0};

	return mat;
}

// Returns the identity matrix
template<typename T>
Matrix<T> id_matrix(const std::size_t size)
{
	return id_matrix<T>(size, size);
}

// Returns the Hilbert matrix
template<typename T>
Matrix<T> hilbert_matrix(const std::size_t rows, const std::size_t cols)
{
	Matrix<T> mat(rows, cols);
	for (std::size_t col = 0; col < cols; ++col)
		for (std::size_t row = 0; row < rows; ++row)
			mat(row, col) = T{1} / (1 + row + col);

	return mat;
}

// Returns the Hilbert matrix
template<typename T>
Matrix<T> hilbert_matrix(const std::size_t size)
{
	return hilbert_matrix<T>(size, size);
}

// Returns the Frank matrix
template<typename T>
Matrix<T> frank_matrix(const std::size_t rows, const std::size_t cols)
{
	Matrix<T> mat(rows, cols);
	for (std::size_t col = 0; col < cols; ++col)
		for (std::size_t row = 0; row < rows; ++row)
			mat(row, col) = 1 + std::min(row, col);

	return mat;
}

// Returns the Frank matrix
template<typename T>
Matrix<T> frank_matrix(const std::size_t size)
{
	return frank_matrix<T>(size, size);
}

// Returns a random matrix
template<typename T>
Matrix<T> random_matrix(const std::size_t rows, const std::size_t cols)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<T> dist;

	Matrix<T> mat(rows, cols);
	for (std::size_t col = 0; col < cols; ++col)
		for (std::size_t row = 0; row < rows; ++row)
			mat(row, col) = dist(gen);

	return mat;
}

// Returns a random matrix
template<typename T>
Matrix<T> random_matrix(const std::size_t size)
{
	return random_matrix<T>(size, size);
}
