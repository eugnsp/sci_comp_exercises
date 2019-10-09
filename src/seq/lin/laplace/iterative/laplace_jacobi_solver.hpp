/*********************************************************************
This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#pragma once
#include "../laplace_solver_base.hpp"
#include <esl/dense.hpp>

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
		auto sol = sol_.view(1, 1, nx_, ny_);
		const auto mul_nda_sol = [this, &sol](auto ix, auto iy) { return mul_nondiag_a(sol, ix, iy); };

		esl::Matrix_x<T> temp(nx_, ny_);
		while (n_its-- > 0)
		{
			temp = (rhs_ - esl::Fn_matrix(nx_, ny_, mul_nda_sol)) / alpha_;
			const auto res = alpha_ * norm_sup(temp - sol);
			ress.push_back(std::log10(res));
			sol = temp;
			fn();
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
