/*********************************************************************
Exponent Taylor series
----------------------

Write a function to calculate exp(x) using its Taylor series expansion.

This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#include "io.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <utility>

#define DISPLAY_AS_HEX

template<class Fn, typename T>
auto sum(Fn fn, T sum_init, unsigned long long n_init = 1, const unsigned long long n_step = 1)
{
	while (true)
	{
		const auto new_sum = sum_init + fn(n_init);
		if (new_sum == std::exchange(sum_init, new_sum))
			return sum_init;
		n_init += n_step;
	}
}

template<typename T>
T exp(T x)
{
	T term = 1;
	return sum([&](auto n) { return term *= x / n; }, T{1});
}

template<typename T>
T exp_even_part(T x)
{
	T term = 1;
	return sum([&](auto n) { return term *= x * x / ((n - 1) * n); }, term, 2, 2);
}

template<typename T>
T exp_odd_part(T x)
{
	T term = x;
	return sum([&](auto n) { return term *= x * x / ((n - 1) * n); }, term, 3, 2);
}

void test_exp_precision(int x)
{
	const auto xf = static_cast<float>(x);
	const auto xd = static_cast<double>(x);

#ifndef DISPLAY_AS_HEX
	std::cout << std::setprecision(std::numeric_limits<float>::digits10);
#endif

	std::cout << "x = " << x << '\n'
			  << "float:\n"
			  << "  std::exp(x)        = " << std::exp(xf) << '\n'
			  << "  sum(x)             = " << exp(xf) << '\n'
			  << "  sum_even + sum_odd = " << exp_even_part(xf) + exp_odd_part(xf) << '\n'
			  << "  sum_even           = " << exp_even_part(xf) << '\n'
			  << "  sum_odd            = " << exp_odd_part(xf) << '\n';
	if (x < 0)
		std::cout << "  1 / sum(-x)        = " << 1.f / exp(-xf) << '\n';

#ifndef DISPLAY_AS_HEX
	std::cout << std::setprecision(std::numeric_limits<double>::digits10);
#endif

	std::cout << "double:\n"
			  << "  std::exp(x)        = " << std::exp(xd) << '\n'
			  << "  sum(x)             = " << exp(xd) << '\n'
			  << "  sum_even + sum_odd = " << exp_even_part(xd) + exp_odd_part(xd) << '\n'
			  << "  sum_even           = " << exp_even_part(xd) << '\n'
			  << "  sum_odd            = " << exp_odd_part(xd) << '\n';
	if (x < 0)
		std::cout << "  1 / sum(-x)        = " << 1. / exp(-xd) << '\n';

	std::cout << std::endl;
}

int main()
{
#ifdef DISPLAY_AS_HEX
	std::cout << std::hexfloat;
#else
	std::cout << std::fixed;
#endif

	test_exp_precision(10);
	test_exp_precision(5);
	test_exp_precision(-5);
	test_exp_precision(-10);

	return 0;
}
