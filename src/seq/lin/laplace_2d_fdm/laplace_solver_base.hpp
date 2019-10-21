/*********************************************************************
This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#pragma once
#include "grid.hpp"
#include <esl/dense.hpp>

#include <cstddef>
#include <vector>

template<typename T, class Solver>
class Laplace_solver_base
{
public:
	template<class Rhs_fn, class Bnd_fn>
	Laplace_solver_base(const Grid<T>& x, const Grid<T>& y, Rhs_fn rhs_fn, Bnd_fn bnd_fn) :
		nx_(x.n - 2), ny_(y.n - 2), alpha_x_(1 / (x.dx() * x.dx())),
		alpha_y_(1 / (y.dx() * y.dx())), alpha_(2 * (alpha_x_ + alpha_y_))
	{
		const Grid<T> x_internal{x[1], x[nx_], nx_};
		const Grid<T> y_internal{y[1], y[ny_], ny_};
		rhs_ = fn_sample(x_internal, y_internal, rhs_fn);

		sol_.resize(x.n, y.n);
		sol_ = 0;
		for (std::size_t iy = 0; iy < y.n; ++iy)
		{
			sol_(0, iy) = bnd_fn(x.min, y[iy]);
			sol_(x.n - 1, iy) = bnd_fn(x.max, y[iy]);
		}

		for (std::size_t ix = 0; ix < x.n; ++ix)
		{
			sol_(ix, 0) = bnd_fn(x[ix], y.min);
			sol_(ix, y.n - 1) = bnd_fn(x[ix], y.max);
		}
	}

	std::vector<T> run(unsigned int n_its)
	{
		return run(n_its, [] {});
	}

	template<class Fn>
	std::vector<T> run(unsigned int n_its, Fn&& fn)
	{
		std::vector<T> ress;
		ress.reserve(n_its);

		static_cast<Solver*>(this)->do_run(n_its, ress, fn);
		return ress;
	}

	const esl::Matrix_x<T>& solution() const
	{
		return sol_;
	}

protected:
	template<class Matrix>
	T mul_nondiag_a(const Matrix& mat, std::size_t ix, std::size_t iy) const
	{
		const auto sx = mat(ix - 1, iy) + mat(ix + 1, iy);
		const auto sy = mat(ix, iy - 1) + mat(ix, iy + 1);
		return -(alpha_x_ * sx + alpha_y_ * sy);
	}

	template<class Matrix>
	T mul_a(const Matrix& mat, std::size_t ix, std::size_t iy) const
	{
		return alpha_ * mat(ix, iy) + mul_nondiag_a(mat, ix, iy);
	}

protected:
	esl::Matrix_x<T> rhs_;
	esl::Matrix_x<T> sol_;

	const std::size_t nx_; // Number of free dofs along x
	const std::size_t ny_; // Number of free dofs along y

	const T alpha_x_; // = 1 / (dx)^2
	const T alpha_y_; // = 1 / (dy^2)
	const T alpha_;	  // = 2 * (alpha_x + alpha_y)
};
