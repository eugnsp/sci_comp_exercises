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
class Laplace_sor_solver : public Laplace_solver_base<T, Laplace_sor_solver<T>>
{
private:
	using Base = Laplace_solver_base<T, Laplace_sor_solver<T>>;

public:
	using Base::Base;

	template<class Fn>
	void do_run(unsigned int n_its, std::vector<T>& ress, Fn&& fn)
	{
		auto sol = sol_.view(1, 1, nx_, ny_);
		const auto a_sol = [this, &sol](auto ix, auto iy) { return mul_a(sol, ix, iy); };

		const auto omega = 2 / (1 + std::sqrt(1 - rho_jacobi_sq()));

		while (n_its-- > 0)
		{
			const auto res = norm_sup(rhs_ - esl::Fn_matrix(nx_, ny_, a_sol));
			ress.push_back(std::log10(res));

			for (std::size_t iy = 0; iy < ny_; ++iy)
				for (std::size_t ix = 0; ix < nx_; ++ix)
				{
					const auto new_sol = (rhs_(ix, iy) - mul_nondiag_a(sol, ix, iy)) / alpha_;
					sol(ix, iy) += omega * (new_sol - sol(ix, iy));
				}

			fn();
		}
	}

private:
	T rho_jacobi_sq() const
	{
		const auto alpha = alpha_y_ / alpha_x_;
		const auto rho = (std::cos(M_PI / nx_) + alpha * std::cos(M_PI / ny_)) / (1 + alpha);
		return rho * rho;
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
