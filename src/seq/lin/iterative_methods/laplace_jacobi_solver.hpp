#pragma once
#include "grid.hpp"
#include "laplace_solver_base.hpp"
#include "matrix.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

template<typename T>
class Laplace_jacobi_solver : public Laplace_solver_base<T>
{
private:
	using Base = Laplace_solver_base<T>;

public:
	using Base::Base;

private:
	void do_run(unsigned int n_its, std::vector<T>& du) override
	{
		Matrix<T> temp(nx_, ny_);

		const auto alpha = .5 / (inv_ddx_ + inv_ddy_);

		while (n_its-- > 0)
		{
			for (std::size_t iy = 1; iy <= ny_; ++iy)
				for (std::size_t ix = 1; ix <= nx_; ++ix)
				{
					const auto sx = sol_(ix - 1, iy) + sol_(ix + 1, iy);
					const auto sy = sol_(ix, iy - 1) + sol_(ix, iy + 1);
					temp(ix - 1, iy - 1) = alpha * (rhs_(ix - 1, iy - 1) + inv_ddx_ * sx + inv_ddy_ * sy);
				}

			T dist = 0;
			for (std::size_t iy = 0; iy < ny_; ++iy)
				for (std::size_t ix = 0; ix < nx_; ++ix)
				{
					dist = std::max(dist, std::abs(temp(ix, iy) - sol_(ix + 1, iy + 1)));
					sol_(ix + 1, iy + 1) = temp(ix, iy);
				}

			du.push_back(std::log10(dist));
		}
	}

private:
	using Base::rhs_;
	using Base::sol_;

	using Base::nx_;
	using Base::ny_;

	using Base::inv_ddx_;
	using Base::inv_ddy_;
};
