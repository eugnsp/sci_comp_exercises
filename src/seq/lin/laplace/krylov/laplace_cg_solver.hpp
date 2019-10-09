#pragma once
#include "../laplace_solver_base.hpp"
#include <esl/dense.hpp>

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

template<typename T>
class Laplace_cg_solver : public Laplace_solver_base<T, Laplace_cg_solver<T>>
{
private:
	using Base = Laplace_solver_base<T, Laplace_cg_solver<T>>;

public:
	using Base::Base;

	template<class Fn>
	void do_run(unsigned int n_its, std::vector<T>& ress, Fn&& fn)
	{
		auto sol = sol_.view(1, 1, nx_, ny_);
		const auto mul_a_sol = [this, &sol](auto ix, auto iy) { return mul_a(sol, ix, iy); };

		esl::Matrix_x<T> pp(nx_ + 2, ny_ + 2, 0);
		auto p = pp.view(1, 1, nx_, ny_);
		const auto mul_a_p = [this, &p](auto ix, auto iy) { return mul_a(p, ix, iy); };

		esl::Matrix_x<T> q(nx_, ny_), r(nx_, ny_);
		p = r = rhs_ - esl::Fn_matrix(nx_, ny_, mul_a_sol);

		auto rho = dot(r, r);
		while (n_its-- > 0)
		{
			// Note: due to round-off errors, res != norm_sup(r)
			const auto res = norm_sup(rhs_ - esl::Fn_matrix(nx_, ny_, mul_a_sol));
			ress.push_back(std::log10(res));

			q = esl::Fn_matrix(nx_, ny_, mul_a_p);

			const auto alpha = rho / dot(p, q);
			sol += alpha * p;
			r -= alpha * q;

			const auto old_rho = std::exchange(rho, dot(r, r));
			const auto beta = rho / old_rho;
			p = beta * p + r;

			fn();
		}
	}

private:
	using Base::mul_a;

	using Base::rhs_;
	using Base::sol_;

	using Base::nx_;
	using Base::ny_;
};
