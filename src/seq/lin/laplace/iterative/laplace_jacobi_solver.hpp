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
class Laplace_jacobi_solver : public Laplace_solver_base<T, Laplace_jacobi_solver<T>>
{
private:
	using Base = Laplace_solver_base<T, Laplace_jacobi_solver<T>>;

public:
	using Base::Base;

	template<class Fn>
	void do_run(unsigned int n_its, std::vector<T>& ress, Fn&& fn)
	{
		auto sol = sol_.view(1, nx_, 1, ny_);
		es_la::Fn_matrix nda_sol(nx_, ny_, [this, &sol](auto ix, auto iy) { return mul_nondiag_a(sol, ix, iy); });

		es_la::Matrix_x<T> temp(nx_, ny_);
		for (auto it = 0u; it < n_its; ++it)
		{
			temp = (rhs_ - nda_sol) / alpha_;
			ress.push_back(std::log10(alpha_ * norm_sup(temp - sol)));
			sol = temp;
			fn(it);
		}
	}

private:
	using Base::mul_nondiag_a;

	using Base::rhs_;
	using Base::sol_;

	using Base::nx_;
	using Base::ny_;

	using Base::alpha_;
};
