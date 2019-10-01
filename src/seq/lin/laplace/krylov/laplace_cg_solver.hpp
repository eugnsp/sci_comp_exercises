#pragma once
#include "../laplace_solver_base.hpp"
#include "matrix.hpp"

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
		Matrix<T> r(nx_, ny_);
		Matrix<T> ap(nx_, ny_);

		const auto alpha = 2 * (inv_ddx_ + inv_ddy_);
		const auto inv_alpha = 1 / alpha;

		for (std::size_t iy = 1; iy <= ny_; ++iy)
			for (std::size_t ix = 1; ix <= nx_; ++ix)
			{
				const auto sx = sol_(ix - 1, iy) + sol_(ix + 1, iy);
				const auto sy = sol_(ix, iy - 1) + sol_(ix, iy + 1);
				r(ix - 1, iy - 1) = rhs_(ix - 1, iy - 1) - alpha * sol_(ix, iy) + inv_ddx_ * sx + inv_ddy_ * sy;
			}

		auto p = r;

		for (auto it = 0u; it < n_its; ++it)
		{
			const auto r_dot_r = dot(r, r);

			for (std::size_t iy = 0; iy < ny_; ++iy)
				for (std::size_t ix = 0; ix < nx_; ++ix)
				{
					const auto sx = (ix == 0 ? 0 : p(ix - 1, iy)) + (ix == nx_ - 1 ? 0 : p(ix + 1, iy));
					const auto sy = (iy == 0 ? 0 : p(ix, iy - 1)) + (iy == ny_ - 1 ? 0 : p(ix, iy + 1));
					ap(ix, iy) = alpha * p(ix, iy) - inv_ddx_ * sx - inv_ddy_ * sy;
				}

			T res = 0;
			auto cg_alpha = r_dot_r / dot(p, ap);
			for (std::size_t iy = 1; iy <= ny_; ++iy)
				for (std::size_t ix = 1; ix <= nx_; ++ix)
				{
					res = std::max(res, cg_alpha * std::abs(p(ix - 1, iy - 1)));
					sol_(ix, iy) += cg_alpha * p(ix - 1, iy - 1);
					r(ix - 1, iy - 1) -= cg_alpha * ap(ix - 1, iy - 1);
				}

			ress.push_back(std::log10(res));

			auto cg_beta = dot(r, r) / r_dot_r;
			for (std::size_t iy = 0; iy < ny_; ++iy)
				for (std::size_t ix = 0; ix < nx_; ++ix)
					p(ix, iy) = r(ix, iy) + cg_beta * p(ix, iy);

			fn(it);
		}
	}

private:
	T dot(const Matrix<T>& x, const Matrix<T>& y)
	{
		T dot = 0;
		for (std::size_t iy = 0; iy < ny_; ++iy)
			for (std::size_t ix = 0; ix < nx_; ++ix)
				dot += x(ix, iy) * y(ix, iy);
		return dot;
	}

private:
	using Base::rhs_;
	using Base::sol_;

	using Base::nx_;
	using Base::ny_;

	using Base::inv_ddx_;
	using Base::inv_ddy_;
};
