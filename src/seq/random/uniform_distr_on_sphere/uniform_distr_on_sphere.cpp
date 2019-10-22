// This file is covered by the LICENSE file in the root of this project.

#include <esc/math.hpp>

#include <array>
#include <cmath>
#include <fstream>
#include <random>

template<typename T>
class Uniform_distr_on_sphere
{
public:
	Uniform_distr_on_sphere()
	:	rand_dist_s_(-1, 1),
		rand_dist_phi_(0, 2 * esc::pi)
	{}

	template<class Random_generator>
	auto operator()(
		Random_generator& gen,
		const double      radius = 1)
	-> std::array<T, 3>
	{
		const auto s   = rand_dist_s_(gen);
		const auto t   = std::sqrt(1 - s * s);
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

	std::ofstream file;
	file.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	file.open("distribution.txt");

	Uniform_distr_on_sphere<double> uds;
	for (unsigned int i = 0; i < 500; ++i)
	{
		const auto [x, y, z] = uds(rand_gen);
		file << x << '\t' << y << '\t' << z << '\n';
	}

	return 0;
}
