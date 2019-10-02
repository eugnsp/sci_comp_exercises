/*********************************************************************
This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#pragma once
#include "../laplace_solver_base.hpp"
#include <es_la/dense.hpp>

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
		auto sol = sol_.view(1, nx_, 1, ny_);
		es_la::Fn_matrix a_sol(nx_, ny_, [this, &sol](auto ix, auto iy) { return mul_a(sol, ix, iy); });

		for (auto it = 0u; it < n_its; ++it)
		{
			ress.push_back(std::log10(norm_sup(rhs_ - a_sol)));

			for (std::size_t iy = 0; iy < ny_; ++iy)
				for (std::size_t ix = 0; ix < nx_; ++ix)
					sol(ix, iy) = (rhs_(ix, iy) - mul_nondiag_a(sol, ix, iy)) / alpha_;

			fn(it);
		}
	}

private:
	using Base::mul_a;
	using Base::mul_nondiag_a;

	using Base::rhs_;
	using Base::sol_;

	using Base::nx_;
	using Base::ny_;

	using Base::alpha_x_;
	using Base::alpha_y_;
	using Base::alpha_;
};
