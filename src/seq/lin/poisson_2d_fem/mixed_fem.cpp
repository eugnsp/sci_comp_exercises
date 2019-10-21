// This file is covered by the LICENSE file in the root of this project.

#include <esf/dof/tools.hpp>
#include <esf/matrix_based.hpp>
#include <esf/mesh/io.hpp>
#include <esf/quadr/quadr.hpp>
#include <esl/sparse.hpp>
#include <esl/sparse/solver/pardiso_solver.hpp>

#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>

constexpr std::size_t element_order = 1;

using Bc1 = esf::Uniform_boundary_cond<esf::Lagrange<element_order>>;
using Var1 = esf::Var<esf::Lagrange<element_order>, 1, Bc1, Bc1>;
using Var2 = esf::Var<esf::Discontinuous_lagrange<element_order - 1>, 2>;
using System = esf::System<esf::Var_list<Var1, Var2>, esf::Dof_mapper>;

using Sp_solver = esl::Pardiso_solver<esl::Csr_matrix<double, esl::Structural_symmetric>>;

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

		system().variable<0>().set_bnd_cond<0>(mesh, bnd1, 0);
		system().variable<0>().set_bnd_cond<1>(mesh, bnd2, .25);

		init();
		compute_and_set_sparsity_pattern(system(), matrix_);
	}

private:
	static auto a_matrix(const esf::Mesh2::Cell_view& cell, double scale)
	{
		using Quadr = esf::Quadr<Var1::Element::order - 1 + Var2::Element::order, 2>;
		const auto grads = esf::gradients<Var1::Element, Quadr>(inv_transp_jacobian(cell));

		constexpr auto n_dofs1 = Var1::Element::total_face_dofs;
		constexpr auto n_dofs2 = Var2::Element::total_face_dofs;
		esl::Matrix<esl::Vector_2d, n_dofs1, n_dofs2> mat;
		for (std::size_t i = 0; i < n_dofs1; ++i)
			for (std::size_t j = 0; j < n_dofs2; ++j)
				mat(i, j) = Quadr::sum([i, j, &grads, scale](auto iq) {
					constexpr auto basis = esf::Element_quadr<Var2::Element, Quadr>::basis();
					return scale * basis(iq, j) * grads(iq, i);
				});

		return mat;
	}

	virtual void assemble() override
	{
		esf::Seq_cell_assembler assembler;
		assembler.assemble(system(), [this](auto& cell) { assemble(cell); });
	}

	void assemble(const esf::Mesh2::Cell_view& cell)
	{
		const auto rhs_fn = [&cell](auto quadr_point_index) {
			auto pt = esf::point(quadr_point_index, cell);
			return std::cos(2 * pt.x()) * std::sin(2 * pt.y());
		};

		const double area = esf::area(cell);
		const auto mat_a12 = a_matrix(cell, area);
		const auto rhs1 = esf::load_vector<Var1::Element>(rhs_fn, area);

		const auto dofs1 = esf::dofs<0>(system(), cell);
		const auto dofs2 = esf::dofs<1>(system(), cell);
		for (std::size_t i1 = 0; i1 < dofs1.size(); ++i1)
			if (const auto d1 = dofs1[i1]; d1.is_free)
			{
				rhs_[d1.index] += rhs1[i1];
				for (std::size_t i2 = 0; i2 < dofs2.size(); ++i2)
				{
					const auto d2 = dofs2[i2];
					matrix_(d1.index, d2.index) += mat_a12(i1, i2)[0];
					matrix_(d2.index, d1.index) -= mat_a12(i1, i2)[0];

					matrix_(d1.index, d2.index + 1) += mat_a12(i1, i2)[1];
					matrix_(d2.index + 1, d1.index) -= mat_a12(i1, i2)[1];
				}
			}
			else
				for (std::size_t i2 = 0; i2 < dofs2.size(); ++i2)
				{
					const auto d2 = dofs2[i2];
					rhs_[d2.index] += mat_a12(i1, i2)[0] * solution_[d1.index];
					rhs_[d2.index + 1] += mat_a12(i1, i2)[1] * solution_[d1.index];
				}

		const auto mat_m2 = esf::mass_matrix<Var2::Element>(area);
		for (std::size_t i1 = 0; i1 < dofs2.size(); ++i1)
		{
			const auto d1 = dofs2[i1];
			for (std::size_t i2 = 0; i2 < dofs2.size(); ++i2)
			{
				const auto d2 = dofs2[i2];
				matrix_(d1.index, d2.index) += mat_m2(i1, i2);
				matrix_(d1.index + 1, d2.index + 1) += mat_m2(i1, i2);
			}
		}
	}
};

int main()
{
	try
	{
		const auto mesh = esf::read_gmsh_mesh("mesh.msh");
		std::cout << mesh << std::endl;

		Solver solver{mesh};
		solver.solve();

		esf::write_gnuplot("mixed1.dat", solver.solution_view<0>());
		esf::write_scattered("mixed2.dat", solver.solution_view<1>());
	}
	catch (const std::exception& e)
	{
		std::cout << "Exception: " << e.what() << std::endl;
		return -1;
	}

	std::cout << "Done.\n";
	return 0;
}
