// This file is covered by the LICENSE file in the root of this project.

#pragma once
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

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
inline constexpr auto sign_mask = Bits<T>{1} << (exponent_size<T> + mantissa_size<T>);

template<typename T>
inline constexpr auto exponent_mask = ((Bits<T>{1} << exponent_size<T>) - 1) << mantissa_size<T>;

template<typename T>
inline constexpr auto mantissa_mask = (Bits<T>{1} << mantissa_size<T>) - 1;

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
	return (bit_cast(value) & sign_mask<T>) != 0;
}

template<typename T>
Bits<T> mantissa(T value)
{
	return (bit_cast(value)) & mantissa_mask<T>;
}

template<typename T>
Bits<T> exponent(T value)
{
	return (bit_cast(value) & exponent_mask<T>) >> mantissa_size<T>;
}

template<typename T>
std::size_t n_normals()
{
	constexpr auto n_exponents = (exponent_mask<T> >> mantissa_size<T>) + 1;
	constexpr auto n_mantissas = mantissa_mask<T> + 1;
	return 2 * (n_exponents - 2) * n_mantissas;
}

template<typename T>
std::size_t n_subnormals()
{
	constexpr auto n_mantissas = mantissa_mask<T> + 1;
	return 2 * (n_mantissas - 1);
}

template<typename T>
std::size_t n_numbers()
{
	constexpr auto n_exponents = (exponent_mask<T> >> mantissa_size<T>) + 1;
	constexpr auto n_mantissas = mantissa_mask<T> + 1;
	return 2 * n_mantissas * (n_exponents - 1);
}
