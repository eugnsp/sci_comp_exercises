// This file is covered by the LICENSE file in the root of this project.

#include <esf/matrix_based.hpp>
#include <esf/mesh/io.hpp>
#include <esl/sparse.hpp>
#include <esl/sparse/solver/pardiso_solver.hpp>
#include <esu/algorithm.hpp>

#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>

constexpr std::size_t element_order = 4;

using Bc = esf::Uniform_boundary_cond<esf::Lagrange<element_order>>;
using Var = esf::Var<esf::Lagrange<element_order>, 1, Bc, Bc>;
using System = esf::System<esf::Var_list<Var>>;

using Sp_solver = esl::Pardiso_solver<esl::Csr_matrix<double, esl::Symmetric_upper>>;

class Solver : public esf::Matrix_based_solver<System, Sp_solver>
{
private:
	using Base = esf::Matrix_based_solver<System, Sp_solver>;

public:
	Solver(const esf::Mesh2& mesh) : Base(mesh)
	{
		const auto br = mesh.bounding_rect();
		const esf::Linestring bnd1{br.bottom_left(), br.top_left()};
		const esf::Linestring bnd2{br.bottom_right(), br.top_right()};

		system().variable().set_bnd_cond<0>(mesh, bnd1, 0);
		system().variable().set_bnd_cond<1>(mesh, bnd2, 0.25);

		init();
		compute_and_set_sparsity_pattern(system(), matrix_);
	}

private:
	virtual void assemble() override
	{
		esf::Seq_cell_assembler assembler;
		assembler.assemble(system(), [this](const auto& cell) { assemble(cell); });
	}

	void assemble(const esf::Mesh2::Cell_view& cell)
	{
		const auto rhs_fn = [&cell](auto quadr_point_index)
		{
			auto pt = esf::point(quadr_point_index, cell);
			return std::cos(2 * pt.x()) * std::sin(2 * pt.y());
		};

		const double area = esf::area(cell);
		const auto mat = esf::stiffness_matrix<Var::Element>(cell, area);
		const auto rhs = esf::load_vector<Var::Element>(rhs_fn, area);

		const auto dofs = esf::dofs(system(), cell);
		for (std::size_t i1 = 0; i1 < dofs.size(); ++i1)
			if (const auto d1 = dofs[i1]; d1.is_free)
			{
				rhs_[d1.index] += rhs[i1];
				for (std::size_t i2 = 0; i2 <= i1; ++i2)
					if (auto d2 = dofs[i2]; d2.is_free)
					{
						const auto [dr, dc] = esu::sorted(d1.index, d2.index);
						matrix_(dr, dc) += mat(i2, i1);
					}
			}
			else
				for (std::size_t i2 = 0; i2 < dofs.size(); ++i2)
					if (auto d2 = dofs[i2]; d2.is_free)
						rhs_[d2.index] -= mat(i2, i1) * solution_[d1.index];
	}
};

int main()
{
	try
	{
		const auto mesh = esf::read_gmsh_mesh("mesh2.msh");
		std::cout << mesh << std::endl;

		Solver solver{mesh};
		solver.solve();

		esf::write_gnuplot("std.dat", solver.solution_view<0>());
		esf::write_interp("std_interp.dat", solver.solution_view<0>(), 0.01);
	}
	catch (const std::exception& e)
	{
		std::cout << "Exception: " << e.what() << std::endl;
		return -1;
	}

	std::cout << "Done.\n";
	return 0;
}
