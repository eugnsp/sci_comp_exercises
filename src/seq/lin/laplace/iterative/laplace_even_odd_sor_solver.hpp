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
class Laplace_even_odd_sor_solver : public Laplace_solver_base<T, Laplace_even_odd_sor_solver<T>>
{
private:
	using Base = Laplace_solver_base<T, Laplace_even_odd_sor_solver<T>>;

public:
	using Base::Base;

	template<class Fn>
	void do_run(unsigned int n_its, std::vector<T>& ress, Fn&& fn)
	{
		Matrix<T> temp;

		const auto rho_jacobi_sq_over_4 = rho_jacobi_sq() / 4;
		T omega = 1;

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

			for (std::size_t ix_s : {1, 2})
			{
				auto ix_f = ix_s;
				for (std::size_t iy = 1; iy <= ny_; ++iy, ix_f = 3 - ix_f)
					for (auto ix = ix_f; ix <= nx_; ix += 2)
					{
						const auto sx = sol_(ix - 1, iy) + sol_(ix + 1, iy);
						const auto sy = sol_(ix, iy - 1) + sol_(ix, iy + 1);
						const auto new_sol = inv_alpha * (rhs_(ix - 1, iy - 1) + inv_ddx_ * sx + inv_ddy_ * sy);
						sol_(ix, iy) += omega * (new_sol - sol_(ix, iy));
					}

				if (it == 0 && ix_s == 1)
					omega = 1 / (1 - 2 * rho_jacobi_sq_over_4);
				else
					omega = 1 / (1 - rho_jacobi_sq_over_4 * omega);
			}

			fn(it);
		}
	}

private:
	T rho_jacobi_sq() const
	{
		const auto alpha = inv_ddy_ / inv_ddx_;
		const auto rho = (std::cos(M_PI / nx_) + alpha * std::cos(M_PI / ny_)) / (1 + alpha);
		return rho * rho;
	}

private:
	using Base::rhs_;
	using Base::sol_;

	using Base::nx_;
	using Base::ny_;

	using Base::inv_ddx_;
	using Base::inv_ddy_;
};
