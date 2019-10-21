/*********************************************************************
1D grid class
-------------

This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#pragma once
#include <esl/dense.hpp>

#include <cstddef>

template<typename T>
struct Grid
{
	const T min;
	const T max;
	const std::size_t n;

	T dx() const
	{
		return (max - min) / (n - 1);
	}

	T operator[](std::size_t i) const
	{
		return min + (max - min) * i / (n - 1);
	}

	T operator()(std::size_t i) const
	{
		return (*this)[i];
	}
};

// Samples a function over a 2D grid formed by the direct product
// of grids `x` and `y`, returns a 2D matrix of `fn(x[row], y[col])`
template<typename T, class Fn>
esl::Matrix_x<T> fn_sample(const Grid<T>& x, const Grid<T>& y, Fn fn)
{
	esl::Matrix_x<T> values(x.n, y.n);
	for (std::size_t iy = 0; iy < y.n; ++iy)
		for (std::size_t ix = 0; ix < x.n; ++ix)
			values(ix, iy) = fn(x[ix], y[iy]);

	return values;
}
