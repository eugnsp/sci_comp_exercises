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
class Laplace_even_odd_sor_solver : public Laplace_solver_base<T, Laplace_even_odd_sor_solver<T>>
{
private:
	using Base = Laplace_solver_base<T, Laplace_even_odd_sor_solver<T>>;

public:
	using Base::Base;

	template<class Fn>
	void do_run(unsigned int n_its, std::vector<T>& ress, Fn&& fn)
	{
		auto sol = sol_.view(1, 1, nx_, ny_);
		const auto a_sol = [this, &sol](auto ix, auto iy) { return mul_a(sol, ix, iy); };

		const auto rho_jacobi_sq_over_4 = rho_jacobi_sq() / 4;
		T omega = 1;

		for (auto it = 0u; it < n_its; ++it)
		{
			const auto res = norm_sup(rhs_ - esl::Fn_matrix(nx_, ny_, a_sol));
			ress.push_back(std::log10(res));

			for (auto ix_s : {0u, 1u})
			{
				std::size_t ix_f = ix_s;
				for (std::size_t iy = 0; iy < ny_; ++iy, ix_f = 1 - ix_f)
					for (std::size_t ix = ix_f; ix < nx_; ix += 2)
					{
						const auto new_sol = (rhs_(ix, iy) - mul_nondiag_a(sol, ix, iy)) / alpha_;
						sol(ix, iy) += omega * (new_sol - sol(ix, iy));
					}

				if (it == 0 && ix_s == 1)
					omega = 1 / (1 - 2 * rho_jacobi_sq_over_4);
				else
					omega = 1 / (1 - omega * rho_jacobi_sq_over_4);
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
