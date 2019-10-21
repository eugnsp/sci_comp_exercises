// This file is covered by the LICENSE file in the root of this project.

#include <iomanip>
#include <iostream>
#include <limits>
#include <type_traits>
#include <utility>

template<class T, class Fn>
auto plain_sum(T n0, const T n1, Fn fn)
{
	std::invoke_result_t<Fn, T> sum{};
	while (n0 <= n1)
		sum += fn(n0++);

	return sum;
}

template<class T, class Fn>
auto kahan_sum(T n0, const T n1, Fn fn)
{
	std::invoke_result_t<Fn, T> sum{};
	std::invoke_result_t<Fn, T> carry{};
	while (n0 <= n1)
	{
		const auto y = fn(n0++) - carry;
		const auto t = sum + y;
		carry = (t - sum) - y;
		sum = t;
	}

	return sum;
}

template<typename T>
auto fn(const T n)
{
	return 1 / (n * n);
}

int main()
{
	const int n = 100'000;

	std::cout << std::setprecision(std::numeric_limits<float>::digits10)
			  << plain_sum(1, n, [](auto n) { return fn<float>(n); }) << '\n'
			  << kahan_sum(1, n, [](auto n) { return fn<float>(n); }) << std::endl;

	std::cout << std::setprecision(std::numeric_limits<double>::digits10)
			  << plain_sum(1, n, [](auto n) { return fn<double>(n); }) << '\n'
			  << kahan_sum(1, n, [](auto n) { return fn<double>(n); }) << std::endl;

	std::cout << std::setprecision(std::numeric_limits<long double>::digits10)
			  << plain_sum(1, n, [](auto n) { return fn<long double>(n); }) << std::endl;

	return 0;
}
