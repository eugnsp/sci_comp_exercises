// This file is covered by the LICENSE file in the root of this project.

#pragma once
#include <cassert>
#include <cstddef>
#include <initializer_list>
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

	const Container& data() const
	{
		return data_;
	}

	std::size_t rows() const
	{
		return rows_;
	}

	std::size_t cols() const
	{
		return cols_;
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

private:
	Container data_;
	std::size_t rows_ = 0;
	std::size_t cols_ = 0;
};
