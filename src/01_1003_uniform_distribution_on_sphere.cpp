/*********************************************************************
Uniform distribution on sphere
------------------------------
Problems and solutions in scientific computing by W.-H. Steeb et al.
Chapter 10, problem 3

Generate uniform sampling points on a sphere.

This file is covered by the LICENSE file in the root of this project.
**********************************************************************/

#define _USE_MATH_DEFINES

#include <array>
#include <cmath>
#include <fstream>
#include <random>

template<typename T>
class Uniform_distr_on_sphere
{
public:
	Uniform_distr_on_sphere() : rand_dist_s_(-1, 1), rand_dist_phi_(0, 2 * M_PI)
	{}

	template<class Random_generator>
	std::array<T, 3> operator()(Random_generator& gen, double radius = 1)
	{
		const auto s = rand_dist_s_(gen);
		const auto t = std::sqrt(1 - s * s);
		const auto phi = rand_dist_phi_(gen);

		return {radius * t * std::cos(phi), radius * t * std::sin(phi), radius * s};
	}

private:
	std::uniform_real_distribution<T> rand_dist_s_;
	std::uniform_real_distribution<T> rand_dist_phi_;
};

int main()
{
	std::random_device rand_dev;
	std::mt19937 rand_gen{rand_dev()};

	std::ofstream file("01_1003_uniform_distribution_on_sphere.txt");

	Uniform_distr_on_sphere<double> uds;
	for (unsigned int i = 0; i < 500; ++i)
	{
		const auto [x, y, z] = uds(rand_gen);
		file << x << ' ' << y << ' ' << z << std::endl;
	}

	return 0;
}
