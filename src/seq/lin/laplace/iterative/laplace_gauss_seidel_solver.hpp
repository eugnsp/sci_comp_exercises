/*********************************************************************
This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#pragma once
#include "../laplace_solver_base.hpp"
#include "matrix.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

template<typename T>
class Laplace_gauss_seidel_solver : public Laplace_solver_base<T, Laplace_gauss_seidel_solver<T>>
{
private:
	using Base = Laplace_solver_base<T, Laplace_gauss_seidel_solver<T>>;

public:
	using Base::Base;

	template<class Fn>
	void do_run(unsigned int n_its, std::vector<T>& ress, Fn&& fn)
	{
		const auto alpha = 2 * (inv_ddx_ + inv_ddy_);
		const auto inv_alpha = 1 / alpha;

		for (auto it = 0u; it < n_its; ++it)
		{
			T res = 0;
			for (std::size_t iy = 1; iy <= ny_; ++iy)
				for (std::size_t ix = 1; ix <= nx_; ++ix)
				{
					const auto sx = sol_(ix - 1, iy) + sol_(ix + 1, iy);
					const auto sy = sol_(ix, iy - 1) + sol_(ix, iy + 1);
					const auto r = rhs_(ix - 1, iy - 1) - alpha * sol_(ix, iy) + inv_ddx_ * sx + inv_ddy_ * sy;
					res = std::max(res, std::abs(r));
				}

			ress.push_back(std::log10(res));

			for (std::size_t iy = 1; iy <= ny_; ++iy)
				for (std::size_t ix = 1; ix <= nx_; ++ix)
				{
					const auto sx = sol_(ix - 1, iy) + sol_(ix + 1, iy);
					const auto sy = sol_(ix, iy - 1) + sol_(ix, iy + 1);
					sol_(ix, iy) = inv_alpha * (rhs_(ix - 1, iy - 1) + inv_ddx_ * sx + inv_ddy_ * sy);
				}

			fn(it);
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
