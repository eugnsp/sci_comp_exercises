/*********************************************************************
IEEE 754
--------

Write functions to extract information from IEEE floating-point
numbers.

This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

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

template<class T>
struct Traits;

template<>
struct Traits<float>
{
	static_assert(std::numeric_limits<float>::is_iec559);

	using Bits = std::uint32_t;
	static constexpr std::size_t exponent_size = 8;
	static constexpr std::size_t mantissa_size = 23;
};

template<>
struct Traits<double>
{
	static_assert(std::numeric_limits<double>::is_iec559);

	using Bits = std::uint64_t;
	static constexpr std::size_t exponent_size = 11;
	static constexpr std::size_t mantissa_size = 52;
};

template<typename T>
using Bits = typename Traits<T>::Bits;

template<typename T>
inline constexpr auto exponent_size = Traits<T>::exponent_size;

template<typename T>
inline constexpr auto mantissa_size = Traits<T>::mantissa_size;

template<typename T>
Bits<T> bit_cast(T value)
{
	Bits<T> bits;
	std::memcpy(&bits, &value, sizeof(T));
	return bits;
}

template<typename T>
bool sign(T value)
{
	return (bit_cast(value) & (Bits<T>{1} << (CHAR_BIT * sizeof(T) - 1))) != 0;
}

template<typename T>
Bits<T> mantissa(T value)
{
	return (bit_cast(value)) & ((Bits<T>{1} << mantissa_size<T>) - 1);
}

template<typename T>
Bits<T> exponent(T value)
{
	return (bit_cast(value) >> mantissa_size<T>) & ~(Bits<T>{1} << exponent_size<T>);
}

template<typename T>
int classify(T value)
{
	const auto m = mantissa(value);
	const auto e = exponent(value);
	constexpr auto max_e = (Bits<T>{1} << exponent_size<T>) - 1;

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
			  << Binary{mantissa(value), mantissa_size<T>} << '\n'
			  << "Type     = " << classify_as_string(value) << '\n'
			  << std::endl;

	assert(sign(value) == std::signbit(value));
	assert(classify(value) == std::fpclassify(value));
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
	test<float>();
	test<double>();

	return 0;
}
