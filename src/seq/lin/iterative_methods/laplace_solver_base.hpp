#pragma once
#include "grid.hpp"
#include "matrix.hpp"
#include <cstddef>
#include <vector>

template<typename T>
class Laplace_solver_base
{
public:
	template<class Rhs_fn, class Bnd_fn>
	Laplace_solver_base(const Grid<T>& x, const Grid<T>& y, Rhs_fn rhs_fn, Bnd_fn bnd_fn) :
		nx_(x.n - 2), ny_(y.n - 2), inv_ddx_(1 / (x.dx() * x.dx())), inv_ddy_(1 / (y.dx() * y.dx()))
	{
		const Grid<T> x_internal{x[1], x[nx_], nx_};
		const Grid<T> y_internal{y[1], y[ny_], ny_};
		rhs_ = at_grid_pts(x_internal, y_internal, rhs_fn);

		sol_.resize(x.n, y.n);
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
		std::vector<T> du;
		du.reserve(n_its);

		do_run(n_its, du);
		return du;
	}

	const Matrix<T>& solution() const
	{
		return sol_;
	}

protected:
	virtual void do_run(unsigned int n_its, std::vector<T>& du) = 0;

protected:
	Matrix<T> rhs_;
	Matrix<T> sol_;

	const std::size_t nx_;	// Number of free dofs along x
	const std::size_t ny_;	// Number of free dofs along y

	const T inv_ddx_;			// = 1 / (dx)^2
	const T inv_ddy_;			// = 1 / (dy^2)
};
