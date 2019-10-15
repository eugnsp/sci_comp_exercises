/*********************************************************************
IEEE 754 ULPs
-------------

Write a function that computes the number of full floating-point
intervals between any two positive machine numbers.

This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#include "ieee754.hpp"
#include <cassert>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>

template<class T>
void sort2(T& x, T& y)
{
	if (y < x)
		std::swap(x, y);
}

template<class T>
std::size_t ulps_bf(T x, T y)
{
	assert(x >= 0);
	assert(y >= 0);

	sort2(x, y);

	std::size_t n = 0;
	while (x < y)
	{
		x = std::nextafter(x, y);
		++n;
	}

	return n;
}

template<class T>
std::size_t ulps(T x, T y)
{
	assert(x >= 0);
	assert(y >= 0);

	sort2(x, y);

	constexpr auto n_mantissas = mantissa_mask<T> + 1;
	return n_mantissas * (exponent(y) - exponent(x)) + mantissa(y) - mantissa(x);
}

template<typename T>
void test(const T x, const T y)
{
	const auto nu = ulps(x, y);
	std::cout << '[' << x << ", " << y << "]: " << nu
			  << " (" << (100. * nu / n_numbers<T>()) << "%)"
			  << std::endl;

	if (ulps(x, y) != ulps_bf(x, y))
		throw std::logic_error("Bad number of ULPs");
}

template<typename T>
void info()
{
	std::cout << "Number of ULPs in " << (std::is_same_v<T, float> ? "float" : "double") << ":\n"
			  << "2^" << CHAR_BIT * sizeof(T) << "       = " << static_cast<Bits<T>>(-1) << '\n'
			  << "Zeros      = 2" << '\n'
			  << "Normals    = " << n_normals<T>() << '\n'
			  << "Subnormals = " << n_subnormals<T>() << '\n'
			  << "Total      = " << n_numbers<T>() << '\n'
			  << std::endl;
}

int main()
{
	try
	{
		info<float>();
		test(-0.f, 0.f);
		test(0.f, 1.f);
		test(1.f, 2.f);
		test(0.f, 0.1f);
		test(0.1f, 0.25f);
	}
	catch (std::exception& e)
	{
		std::cout << "Exception: " << e.what() << std::endl;
		return 1;
	}

	std::cout << "OK." << std::endl;
	return 0;
}
