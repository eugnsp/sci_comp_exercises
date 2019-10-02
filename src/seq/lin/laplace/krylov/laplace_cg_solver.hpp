#pragma once
#include "../laplace_solver_base.hpp"
#include <es_la/dense.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
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
		auto sol = sol_.view(1, nx_, 1, ny_);
		es_la::Fn_matrix mul_a_sol(nx_, ny_, [this, &sol](auto ix, auto iy) { return mul_a(sol, ix, iy); });

		es_la::Matrix_x<T> pp(nx_ + 2, ny_ + 2, 0);
		auto p = pp.view(1, nx_, 1, ny_);
		es_la::Fn_matrix mul_a_p(nx_, ny_, [this, &p](auto ix, auto iy) { return mul_a(p, ix, iy); });

		es_la::Matrix_x<T> ap(nx_, ny_), r(nx_, ny_);

		p = r = rhs_ - mul_a_sol;

		auto rr_old = dot(r, r);
		for (auto it = 0u; it < n_its; ++it)
		{
			ress.push_back(std::log10(norm_sup(rhs_ - mul_a_sol)));

			ap = mul_a_p;

			const auto alpha = rr_old / dot(p, ap);
			sol += alpha * p;
			r -= alpha * ap;

			const auto rr = dot(r, r);
			const auto beta = rr / rr_old;
			p = beta * p + r;
			rr_old = rr;

			fn(it);
		}
	}

private:
	using Base::mul_a;

	using Base::rhs_;
	using Base::sol_;

	using Base::nx_;
	using Base::ny_;

	using Base::alpha_x_;
	using Base::alpha_y_;
	using Base::alpha_;
};
