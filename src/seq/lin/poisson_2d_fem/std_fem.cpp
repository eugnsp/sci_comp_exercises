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

template<
	class System,
	class Sp_solver>
class Solver : public esf::Matrix_based_solver<System, Sp_solver>
{
private:
	using Base = esf::Matrix_based_solver<System, Sp_solver>;

public:
	Solver(const esf::Mesh2& mesh)
	:	Base(mesh)
	{
		const auto br = mesh.bounding_rect();
		const esf::Linestring bnd1{br.bottom_left(), br.top_left()};
		const esf::Linestring bnd2{br.bottom_right(), br.top_right()};

		system().variable().template set_bnd_cond<0>(mesh, bnd1, 0);
		system().variable().template set_bnd_cond<1>(mesh, bnd2, 0.25);

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
		const auto mat    = esf::stiffness_matrix<Element>(cell, area);
		const auto rhs    = esf::load_vector<Element>(rhs_fn, area);

		const auto dofs = esf::dofs(system(), cell);
		for (std::size_t c = 0; c < dofs.size(); ++c)
			if (const auto dc = dofs[c]; dc.is_free)
			{
				rhs_[dc.index] += rhs[c];
				for (std::size_t r = 0; r <= c; ++r)
					if (auto dr = dofs[r]; dr.is_free)
					{
						const auto [d1, d2] = esu::sorted(dc.index, dr.index);
						matrix_(d1, d2) += mat(r, c);
					}
			}
			else
				for (std::size_t r = 0; r < dofs.size(); ++r)
					if (auto dr = dofs[r]; dr.is_free)
						rhs_[dr.index] -= mat(r, c) * solution_[dc.index];
	}

private:
	using Element = typename Base::System::template Var<0>::Element;

	using Base::system;
	using Base::init;

	using Base::solution_;
	using Base::rhs_;
	using Base::matrix_;
};

template<std::size_t element_order>
class Solver_type
{
private:
	using Bnd_cond = esf::Uniform_boundary_cond<esf::Lagrange<element_order>>;
	using Var      = esf::Var<esf::Lagrange<element_order>, 1, Bnd_cond, Bnd_cond>;
	using System   = esf::System<esf::Var_list<Var>>;

public:
	using Type = Solver<System, esl::Pardiso_solver<esl::Csr_matrix<double, esl::Symmetric_upper>>>;
};

int main()
{
	try
	{
		const auto mesh = esf::read_gmsh_mesh("mesh2.msh");
		std::cout << mesh << std::endl;

		Solver_type<4>::Type solver{mesh};
		solver.solve();

		esf::write_gnuplot("std_u.dat", solver.solution_view<0>());
		esf::write_interp("std_u_interp.dat", solver.solution_view<0>(), 0.01);
	}
	catch (const std::exception& e)
	{
		std::cout << "Exception: " << e.what() << std::endl;
		return -1;
	}

	std::cout << "Done.\n";
	return 0;
}
