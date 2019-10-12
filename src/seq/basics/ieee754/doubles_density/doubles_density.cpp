/*********************************************************************
Doubles density
---------------

Compute the relative density of double floating-point numbers.

This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#include "io.hpp"
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

int main()
{
	constexpr unsigned int n = 1'000;

	std::vector<double> ds;
	for (unsigned int i = 0; i < n; ++i)
		ds.push_back(std::pow(2, (4. * i / n) - 2));

	std::vector<double> density;
	for (const auto d : ds)
	{
		const auto dn = std::nextafter(d, std::numeric_limits<double>::max());
		density.push_back((dn - d) / d);
	}

	write_vec("doubles_density.txt", ds, density);
	return 0;
}
