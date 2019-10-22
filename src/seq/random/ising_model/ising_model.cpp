// This file is covered by the LICENSE file in the root of this project.

#include <esl/dense.hpp>

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

template<typename T>
T sq(T x)
{
	return x * x;
}

class Ising_lattice
{
public:
	struct Stat_params
	{
		double energy;
		double magnetization;
		double sp_heat;
		double susceptibility;
	};

public:
	Ising_lattice(
		const std::size_t n,
		const double temp,
		const double coupling,
		const double field,
		const bool rand_init = true)
	:	n_(n),
		s_(n, n),
		temp_(temp),
		coupling_(coupling),
		field_(field),
		rand_gen_(rand_dev_()),
		rand_n_distr_(0, n_ - 1)
	{
		init_spins(rand_init);
	}

	void set_temp(const double temp)
	{
		temp_ = temp;
	}

	void sweep(const unsigned int n)
	{
		sweep(n, [](auto) {});
	}

	Stat_params stat_params(const unsigned int n)
	{
		Stat_params p{};

		double en = total_energy();
		auto magn = total_magnetization();
		std::size_t i = 1;

		sweep(n, [&](auto d_en_magn)
		{
			en   += d_en_magn.first;
			magn += d_en_magn.second;

			p.energy         += (en - p.energy) / i;
			p.magnetization  += (std::abs(magn) - p.magnetization) / i;
			p.sp_heat        += (sq(en) - p.sp_heat) / i;
			p.susceptibility += (sq(magn) - p.susceptibility) / i;

			++i;
		});

		p.sp_heat        = (p.sp_heat - sq(p.energy)) / sq(temp_);
		p.susceptibility = (p.susceptibility - sq(p.magnetization)) / temp_;

		p.energy         /= sq(n_);
		p.magnetization  /= sq(n_);
		p.sp_heat        /= sq(n_);
		p.susceptibility /= sq(n_);
		return p;
	}

	void write(const std::string& file_name) const
	{
		std::ofstream file;
		file.exceptions(std::ofstream::failbit | std::ofstream::badbit);
		file.open(file_name);

		for (std::size_t col = 0; col < n_; ++col)
		{
			for (std::size_t row = 0; row < n_; ++row)
				file << s_(row, col) << ' ';
			file << '\n';
		}
		file << '\n';
	}

private:
	std::size_t prev(const std::size_t i) const
	{
		return (i == 0) ? n_ - 1 : i - 1;
	}

	std::size_t next(const std::size_t i) const
	{
		return (i == n_ - 1) ? 0 : i + 1;
	}

	void init_spins(const bool rand_init)
	{
		std::bernoulli_distribution rand_distr;
		for (std::size_t col = 0; col < n_; ++col)
			for (std::size_t row = 0; row < n_; ++row)
				s_(row, col) = (!rand_init || rand_distr(rand_gen_)) ? 1 : -1;
	}

	bool take_step(const double d_en)
	{
		if (d_en <= 0)
			return true;

		const auto w = std::exp(-d_en / temp_);
		return w > rand_unit_distr_(rand_gen_);
	}

	int neighbour_sum(const std::size_t row,
					  const std::size_t col) const
	{
		return s_(prev(row), col) + s_(next(row), col) + s_(row, prev(col)) + s_(row, next(col));
	}

	std::pair<double, int> step()
	{
		const auto row = rand_n_distr_(rand_gen_);
		const auto col = rand_n_distr_(rand_gen_);

		const auto d_en = 2 * s_(row, col) * (coupling_ * neighbour_sum(row, col) + field_);
		if (take_step(d_en))
		{
			s_(row, col) = -s_(row, col);
			return {d_en, 2 * s_(row, col)};
		}

		return {0, 0};
	}

	template<class Fn>
	std::size_t sweep(std::size_t n,
					  Fn&& 		  fn)
	{
		const auto n_steps = n * sq(n_);
		for (std::size_t i = 0; i < n_steps; ++i)
		{
			const auto d = step();
			fn(d);
		}

		return n_steps;
	}

	double total_energy() const
	{
		double en = 0;
		for (std::size_t col = 0; col < n_; ++col)
			for (std::size_t row = 0; row < n_; ++row)
				en -= s_(row, col) * (.5 * coupling_ * neighbour_sum(row, col) + field_);

		return en;
	}

	int total_magnetization() const
	{
		int magn = 0;
		for (std::size_t col = 0; col < n_; ++col)
			for (std::size_t row = 0; row < n_; ++row)
				magn += s_(row, col);
		return magn;
	}

private:
	const std::size_t  n_;
	esl::Matrix_x<int> s_;

	double 		 temp_;
	const double coupling_;
	const double field_;

	std::random_device rand_dev_;
	std::mt19937       rand_gen_;
	std::uniform_int_distribution<std::size_t> rand_n_distr_;
	std::uniform_real_distribution<double>     rand_unit_distr_;
};

void params_vs_temp(
	const double 	   field,
	const std::string& file_name)
{
	std::ofstream file;
	file.exceptions(std::ofstream::failbit | std::ofstream::badbit);
	file.open(file_name);

	const std::size_t lattice_size = 15;
	const double coupling = 1;
	const double t_max = 8 * coupling;
	const std::size_t n_points = 100;

	Ising_lattice lattice(lattice_size, t_max / n_points, coupling, field * coupling, false);
	for (std::size_t i = 1; i <= n_points; ++i)
	{
		const double temp = t_max * i / n_points;
		lattice.set_temp(temp);
		lattice.sweep(10'000);

		const auto p = lattice.stat_params(100'000);
		file << temp / coupling << '\t' << p.energy / coupling << '\t' << p.magnetization << '\t'
			 << p.sp_heat << '\t' << p.susceptibility * coupling << '\n';
	}
}

void lattice_after_sweep()
{
	const std::size_t lattice_size = 100;
	const double coupling = 1;
	const double temp = 2 * coupling;
	Ising_lattice lattice(lattice_size, temp, coupling, 0, true);

	lattice.sweep(1'000);
	lattice.write("lattice.txt");
}

int main()
{
	params_vs_temp(0,  "mt0.txt");
	params_vs_temp(.1, "mt1.txt");
	params_vs_temp(.5, "mt2.txt");

	lattice_after_sweep();

	return 0;
}
