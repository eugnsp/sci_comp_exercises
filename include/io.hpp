/*********************************************************************
Input/output
------------

This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#pragma once
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

// template<typename T>
// std::ostream& operator<<(std::ostream& os, const Matrix<T>& mat)
// {
// 	os << std::fixed << std::setprecision(4);
// 	for (std::size_t row = 0; row < mat.rows(); ++row)
// 	{
// 		for (std::size_t col = 0; col < mat.cols(); ++col)
// 			os << std::setw(7) << mat(row, col) << ' ';
// 		os << '\n';
// 	}

// 	return os;
// }

template<class Vec, class... Vecs>
void write_vec(std::string file_name, const Vec& vec, const Vecs&... vecs)
{
	assert(((vec.size() == vecs.size()) && ...));

	std::ofstream file;
	file.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	file.open(file_name);

	for (std::size_t i = 0; i < vec.size(); ++i)
	{
		file << i << '\t' << vec[i];
		((file << '\t' << vecs[i]), ...);
		file << '\n';
	}
}

template<class Matrix, class X_titles, class Y_titles>
void write_gnuplot(std::string file_name, const Matrix& mat, X_titles x_titles, Y_titles y_titles)
{
	std::ofstream file;
	file.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	file.open(file_name, std::ofstream::binary);

	const auto write = [&file](auto v)
	{
		const auto val = static_cast<float>(v);
		file.write(reinterpret_cast<const char*>(&val), sizeof(float));
	};

	write(mat.rows());
	for (std::size_t ix = 0; ix < mat.rows(); ++ix)
		write(x_titles(ix));

	for (std::size_t iy = 0; iy < mat.cols(); ++iy)
	{
		write(y_titles(iy));
		for (std::size_t ix = 0; ix < mat.rows(); ++ix)
			write(mat(ix, iy));
	}
}
