/*********************************************************************
Multiresolution analysis
------------------------
An introduction to scientific computing by I.Danaila et al.
Chapter 6

Implement 1D forward and backward wavelet transforms for Haar,
Schauder and Daubechies wavelets.

This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#define _USE_MATH_DEFINES

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

struct Haar_wavelet
{
	static constexpr double forward[2][2] = {{M_SQRT1_2, M_SQRT1_2}, {M_SQRT1_2, -M_SQRT1_2}};

	static constexpr std::ptrdiff_t backward_shift = 0;
	static constexpr double backward_c[2][1] = {{M_SQRT1_2}, {M_SQRT1_2}};
	static constexpr double backward_d[2][1] = {{M_SQRT1_2}, {-M_SQRT1_2}};
};

struct Schauder_wavelet
{
	static constexpr double forward[2][3] = {{M_SQRT2, 0, 0}, {-M_SQRT1_2, M_SQRT2, -M_SQRT1_2}};

	static constexpr std::ptrdiff_t backward_shift = 0;
	static constexpr double backward_c[2][2] = {{M_SQRT1_2, 0}, {M_SQRT1_2 / 2, M_SQRT1_2 / 2}};
	static constexpr double backward_d[2][2] = {{0, 0}, {M_SQRT1_2, 0}};
};

struct Daubechies_wavelet
{
private:
	static constexpr auto sqrt3 = 1.732050807568877;
	static constexpr auto c0 = (1 + sqrt3) / (4 * M_SQRT2);
	static constexpr auto c1 = (3 + sqrt3) / (4 * M_SQRT2);
	static constexpr auto c2 = (3 - sqrt3) / (4 * M_SQRT2);
	static constexpr auto c3 = (1 - sqrt3) / (4 * M_SQRT2);

public:
	static constexpr double forward[2][4] = {{c0, c1, c2, c3}, {c3, -c2, c1, -c0}};

	static constexpr std::ptrdiff_t backward_shift = -1;
	static constexpr double backward_c[2][2] = {{c2, c0}, {c3, c1}};
	static constexpr double backward_d[2][2] = {{c1, c3}, {-c0, -c2}};
};

template<typename T>
bool is_pow2(T x)
{
	return (x & (x - 1)) == 0;
}

template<class Vec1, class Vec2>
bool are_equal(const Vec1& v1, const Vec2& v2, double tol = 1e-6)
{
	if (v1.size() != v2.size())
		return false;

	for (std::size_t i = 0; i < v1.size(); ++i)
		if (std::abs(v1[i] - v2[i]) > tol)
			return false;

	return true;
}

template<class Vec, typename T>
void zero_small_coeffs(Vec& vec, T max_value)
{
	for (auto& v : vec)
		if (std::abs(v) < max_value)
			v = 0;
}

template<class Fn, typename T>
auto fn_sample(Fn&& fn, std::size_t n, T x_min, T x_max)
{
	std::vector<T> xs;
	std::vector<std::invoke_result_t<Fn, T>> values;
	xs.reserve(n);
	values.reserve(n);

	for (std::size_t i = 0; i < n; ++i)
	{
		const auto x = x_min + (x_max - x_min) * i / n;
		xs.push_back(x);
		values.push_back(fn(x));
	}

	return std::make_pair(xs, values);
}

template<class Wavelet, class In, class Out, std::size_t... is>
void wavelet_forward_impl(const In& in, Out& out, std::index_sequence<is...>)
{
	assert(is_pow2(in.size()));
	auto n = in.size();

	std::vector<typename Out::value_type> tmp(n);

	out = in;
	while (n > 1)
	{
		const auto n2 = n / 2;
		for (std::size_t k = 0; k < n2; ++k)
		{
			tmp[k] = ((Wavelet::forward[0][is] * out[(2 * k + is) % n]) + ...);
			tmp[n2 + k] = ((Wavelet::forward[1][is] * out[(2 * k + is) % n]) + ...);
		}
		std::copy_n(tmp.begin(), n, out.begin());
		n = n2;
	}
}

template<class Wavelet, class In, class Out, std::size_t... is>
void wavelet_backward_impl(const In& in, Out& out, std::index_sequence<is...>)
{
	assert(is_pow2(in.size()));
	const auto n = in.size();

	std::vector<typename Out::value_type> tmp(n);

	out = in;
	std::size_t j = 1;
	while (j < n)
	{
		const auto j2 = 2 * j;
		std::copy_n(out.begin(), j2, tmp.begin());
		for (std::size_t k = 0; k < j; ++k)
		{
			out[2 * k] = ((Wavelet::backward_c[0][is] * tmp[(k + is + Wavelet::backward_shift) % j]) + ...) +
						 ((Wavelet::backward_d[0][is] * tmp[j + (k + is + Wavelet::backward_shift) % j]) + ...);
			out[2 * k + 1] = ((Wavelet::backward_c[1][is] * tmp[(k + is + Wavelet::backward_shift) % j]) + ...) +
							 ((Wavelet::backward_d[1][is] * tmp[j + (k + is + Wavelet::backward_shift) % j]) + ...);
		}
		j = j2;
	}
}

template<class Wavelet, class In, class Out>
void wavelet_forward(const In& in, Out& out)
{
	constexpr auto n = std::extent_v<decltype(Wavelet::forward), 1>;
	wavelet_forward_impl<Wavelet>(in, out, std::make_index_sequence<n>{});
}

template<class Wavelet, class In, class Out>
void wavelet_backward(const In& in, Out& out)
{
	constexpr auto n = std::extent_v<decltype(Wavelet::backward_c), 1>;
	wavelet_backward_impl<Wavelet>(in, out, std::make_index_sequence<n>{});
}

template<class Wavelet>
void wavelet_test(std::string file_name, double eps1, double eps2)
{
	const unsigned int j = 10;
	const std::size_t n = 1u << j;

	const double x_min = 0;
	const double x_max = 1;

	const auto f = [](double x) { return std::exp(-x) * std::sin(4 * M_PI * x); };
	const auto c = fn_sample(f, n, x_min, x_max);

	std::vector<double> d;
	wavelet_forward<Wavelet>(c.second, d);

	std::vector<double> c0;
	wavelet_backward<Wavelet>(d, c0);
	assert(are_equal(c.second, c0));

	std::vector<double> c1;
	zero_small_coeffs(d, eps1);
	wavelet_backward<Wavelet>(d, c1);

	std::vector<double> c2;
	zero_small_coeffs(d, eps2);
	wavelet_backward<Wavelet>(d, c2);

	std::ofstream file(file_name);
	for (std::size_t i = 0; i < n; ++i)
		file << c.first[i] << ' ' << c.second[i] << ' ' << c1[i] << ' ' << c2[i] << std::endl;
}

int main()
{
	wavelet_test<Haar_wavelet>("02_06_multiresolution_analysis_haar.txt", .025, .25);
	wavelet_test<Schauder_wavelet>("02_06_multiresolution_analysis_schauder.txt", .1, .5);
	wavelet_test<Daubechies_wavelet>("02_06_multiresolution_analysis_daubechies.txt", .1, .5);

	return 0;
}
