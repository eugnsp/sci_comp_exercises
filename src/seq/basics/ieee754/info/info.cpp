// This file is covered by the LICENSE file in the root of this project.

#include "../ieee754.hpp"
#include <climits>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

template<typename T>
class Binary_prn
{
public:
	Binary_prn(T value, std::size_t n = CHAR_BIT * sizeof(T)) : value_(value), n_(n)
	{}

	void operator()(std::ostream& os) const
	{
		for (std::size_t i = n_; i > 0; --i)
		{
			const auto bit = T{1} << (i - 1);
			os << ((value_ & bit) ? 1 : 0);
		}
	}

private:
	const T value_;
	const std::size_t n_;
};

template<typename T>
std::ostream& operator<<(std::ostream& os, Binary_prn<T> bin)
{
	bin(os);
	return os;
}

////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
int classify(T value)
{
	const auto m = mantissa(value);
	const auto e = exponent(value);
	constexpr auto max_e = exponent_mask<T> >> mantissa_size<T>;

	if (e == 0)
		return m == 0 ? FP_ZERO : FP_SUBNORMAL;
	if (e == max_e)
		return m == 0 ? FP_INFINITE : FP_NAN;
	return FP_NORMAL;
}

template<typename T>
std::string classify_as_string(T value)
{
	switch (classify(value))
	{
	case FP_ZERO:
		return "zero";
	case FP_SUBNORMAL:
		return "subnormal";
	case FP_INFINITE:
		return "inf";
	case FP_NAN:
		return "nan";
	default:
		return "normal";
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////

template<typename T>
void test(T value)
{
	std::cout << value << '\n'
			  << "Value    = " << Binary_prn{bit_cast(value)} << '\n'
			  << "Sign     = " << (sign(value) ? 1 : 0) << '\n'
			  << "Exponent =  " << Binary_prn{exponent(value), exponent_size<T>} << '\n'
			  << "Mantissa = " << std::string(exponent_size<T> + 1, ' ')
			  << Binary_prn{mantissa(value), mantissa_size<T>} << '\n'
			  << "Type     = " << classify_as_string(value) << '\n'
			  << std::endl;

	if (sign(value) != std::signbit(value))
		throw std::logic_error("Bad sign bit");

	if (classify(value) != std::fpclassify(value))
		throw std::logic_error("Bad floating-point class");
}

template<typename T>
void test()
{
	test(static_cast<T>(0));
	test(-static_cast<T>(0));
	test(static_cast<T>(0.1));
	test(-static_cast<T>(0.1));
	test(static_cast<T>(1));
	test(-static_cast<T>(1));

	test(std::numeric_limits<T>::min());
	test(std::numeric_limits<T>::lowest());
	test(std::numeric_limits<T>::max());
	test(std::numeric_limits<T>::denorm_min());
	test(std::numeric_limits<T>::epsilon());
	test(std::numeric_limits<T>::infinity());
	test(-std::numeric_limits<T>::infinity());
	test(std::numeric_limits<T>::quiet_NaN());
	test(std::numeric_limits<T>::signaling_NaN());
}

int main()
{
	try
	{
		test<float>();
		test<double>();
	}
	catch (std::exception& e)
	{
		std::cout << "Exception: " << e.what() << std::endl;
		return 1;
	}

	std::cout << "OK." << std::endl;
	return 0;
}
