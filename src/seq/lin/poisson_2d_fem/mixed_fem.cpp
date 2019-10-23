// This file is covered by the LICENSE file in the root of this project.

#include <esf/matrix_based.hpp>
#include <esf/mesh/io.hpp>
#include <esl/sparse.hpp>
#include <esl/sparse/solver/pardiso_solver.hpp>

#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>

template<class System, class Sp_solver>
class Solver : public esf::Matrix_based_solver<System, Sp_solver>
{
private:
	using Base     = esf::Matrix_based_solver<System, Sp_solver>;
	using Element0 = typename esf::Var_type<System, 0>::Element;
	using Element1 = typename esf::Var_type<System, 1>::Element;

public:
	Solver(const esf::Mesh2& mesh) : Base(mesh)
	{
		const auto br = mesh.bounding_rect();
		const esf::Linestring bnd1{br.bottom_left(), br.top_left()};
		const esf::Linestring bnd2{br.bottom_right(), br.top_right()};

		system().template variable<0>().template set_bnd_cond<0>(mesh, bnd1, 0);
		system().template variable<0>().template set_bnd_cond<1>(mesh, bnd2, .25);

		init();
		compute_and_set_sparsity_pattern(system(), matrix_);
	}

private:
	static auto a_matrix(const esf::Mesh2::Cell_view& cell, const double scale)
	{
		using Quadr = esf::Quadr<Element0::order - 1 + Element1::order, esf::Dim2>;
		auto grads = esf::gradients<Element0, Quadr>(inv_transp_jacobian(cell));
		grads *= scale;

		return esl::make_matrix<Element0::total_face_dofs, Element1::total_face_dofs>(
			[&grads](std::size_t row, std::size_t col)
			{
				return Quadr::sum([row, col, &grads](auto iq)
				{
					constexpr auto basis = esf::Element_quadr<Element1, Quadr>::basis();
					return basis(iq, col) * grads(iq, row);
				}).eval();
			});
	}

	virtual void assemble() override
	{
		esf::Seq_cell_assembler assembler;
		assembler.assemble(system(), [this](auto& cell) { assemble(cell); });
	}

	void assemble(const esf::Mesh2::Cell_view& cell)
	{
		const auto rhs_fn = [&cell](auto quadr_point_index)
		{
			auto pt = esf::point(quadr_point_index, cell);
			return std::cos(2 * pt.x()) * std::sin(2 * pt.y());
		};

		const double area = esf::area(cell);
		const auto mat_a12 = a_matrix(cell, area);
		const auto rhs1 = esf::load_vector<Element0>(rhs_fn, area);

		const auto dofs0 = esf::dofs<0>(system(), cell);
		const auto dofs1 = esf::dofs<1>(system(), cell);

		for (std::size_t r = 0; r < dofs0.size(); ++r)
			if (const auto d0r = dofs0[r]; d0r.is_free)
			{
				rhs_[d0r.index] += rhs1[r];
				for (std::size_t c = 0; c < dofs1.size(); ++c)
				{
					const auto d1c = dofs1[c];
					const auto d2c = d1c + 1;

					matrix_(d0r.index, d1c.index) += mat_a12(r, c)[0];
					matrix_(d1c.index, d0r.index) -= mat_a12(r, c)[0];

					matrix_(d0r.index, d2c.index) += mat_a12(r, c)[1];
					matrix_(d2c.index, d0r.index) -= mat_a12(r, c)[1];
				}
			}
			else
				for (std::size_t c = 0; c < dofs1.size(); ++c)
				{
					const auto d1c = dofs1[c];
					const auto d2c = d1c + 1;

					rhs_[d1c.index] += mat_a12(r, c)[0] * solution_[d0r.index];
					rhs_[d2c.index] += mat_a12(r, c)[1] * solution_[d0r.index];
				}

		const auto mat_m2 = esf::mass_matrix<Element1>(area);
		for (std::size_t r = 0; r < dofs1.size(); ++r)
		{
			const auto d1 = dofs1[r];
			for (std::size_t c = 0; c < dofs1.size(); ++c)
			{
				const auto d2 = dofs1[c];
				const auto d3 = d2 + 1;

				matrix_(d1.index, d2.index) += mat_m2(r, c);
				matrix_(d3.index, d3.index) += mat_m2(r, c);
			}
		}
	}

private:
	using Base::init;
	using Base::system;

	using Base::matrix_;
	using Base::rhs_;
	using Base::solution_;
};

template<std::size_t element_order>
class Solver_type
{
private:
	using Bnd_cond = esf::Uniform_boundary_cond<esf::Lagrange<element_order>>;
	using Var1     = esf::Var<esf::Lagrange<element_order>, 1, Bnd_cond, Bnd_cond>;
	using Var2     = esf::Var<esf::Discontinuous_lagrange<element_order - 1>, 2>;
	using System   = esf::System<esf::Var_list<Var1, Var2>, esf::Dof_mapper>;

public:
	using Type =
		Solver<System, esl::Pardiso_solver<esl::Csr_matrix<double, esl::Structural_symmetric>>>;
};

int main()
{
	try
	{
		const auto mesh = esf::read_gmsh_mesh("mesh.msh");
		std::cout << mesh << std::endl;

		Solver_type<1>::Type solver1{mesh};
		solver1.solve();

		esf::write_gnuplot(  "mixed_u.dat",   solver1.solution_view<0>());
		esf::write_gnuplot(  "mixed_s.dat",   solver1.solution_view<1>());
		esf::write_scattered("mixed_sv.dat", solver1.solution_view<1>());
	}
	catch (const std::exception& e)
	{
		std::cout << "Exception: " << e.what() << std::endl;
		return -1;
	}

	std::cout << "Done.\n";
	return 0;
}
