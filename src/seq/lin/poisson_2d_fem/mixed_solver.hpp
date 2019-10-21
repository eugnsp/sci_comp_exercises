#pragma once
#include "mixed_system.hpp"

#include <esf/geometry/function.hpp>
#include <esf/math/jacobian.hpp>
#include <esf/matrix_based/solver.hpp>
#include <esf/mesh/mesh2.hpp>
#include <esl/io/matfile_writer.hpp>
#include <esl/sparse/solver/pardiso_solver.hpp>
#include <esu/numeric.hpp>

#include <esf/matrix_based/par_cell_assembler.hpp>
#include <esf/matrix_based/seq_cell_assembler.hpp>

//#include <esl/function.hpp>
#include <esf/dof/layered_dof_mapper.hpp>
#include <esf/dof/tools.hpp>
#include <esf/element/lagrange.hpp>
#include <esf/io/gnuplot_writer2.hpp>
#include <esf/io/matlab_writer2.hpp>
#include <esf/math/matrix.hpp>
#include <esf/math/quadr.hpp>
#include <esf/var_list.hpp>

#include <algorithm>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

//#include <es_util/timer.hpp>

namespace mixed_fem
{
using Sp_matrix = esl::Csr_matrix<double, esl::Structural_symmetric>;
using Sp_solver = esl::Pardiso_solver<Sp_matrix>;

class Solver : public esf::Matrix_based_solver<System, Sp_solver>
{
private:
	using Base = esf::Matrix_based_solver<System, Sp_solver>;
	using System = typename Base::System;

public:
	using Phi_element = Phi::Element;
	using Exy_element = Exy::Element;

private:
	static constexpr auto n_dofs_tau = Phi_element::n_total_face_dofs;
	static constexpr auto n_dofs_chi = Exy_element::n_total_face_dofs;

public:
	using Base::mesh;
	using Base::system;

public:
	Solver(const esf::Mesh2& mesh) : Base(mesh)
	{
		system().variable<0>().set_name("phi");
		system().variable<1>().set_name("exy");
	}

	void init()
	{
		//system().variable<0>().set_bnd_cond<0>(mesh(), 0);
		//system().variable<1>().set_bnd_cond<1>(mesh(), 0);

		Base::init();
		esf::compute_and_set_sparsity_pattern(system(), matrix_);
	}

private:
	virtual void set_bnd_values() override
	{}

public:
	template<class Quadr, class Grads>
	auto ax_matrix(const Grads& grads)
	{
		static_assert(Grads::rows() == Quadr::size);
		static_assert(Grads::cols() == n_dofs_chi);

		esl::Matrix_d<n_dofs_chi, n_dofs_tau> m;
		for (std::size_t i = 0; i < n_dofs_chi; ++i)
			for (std::size_t j = 0; j < n_dofs_tau; ++j)
				m(i, j) = Quadr::sum([i, j, &grads](auto iq) {
					auto g = grads(iq, i);
					auto basis = esf::Element_quadr<Phi_element, Quadr>::basis();
					return g[0] * basis(iq, j);
				});

		return m;
	}

	template<class Quadr, class Grads>
	auto ay_matrix(const Grads& grads)
	{
		static_assert(Grads::rows() == Quadr::size);
		static_assert(Grads::cols() == n_dofs_chi);

		esl::Matrix_d<n_dofs_chi, n_dofs_tau> m;
		for (std::size_t i = 0; i < n_dofs_chi; ++i)
			for (std::size_t j = 0; j < n_dofs_tau; ++j)
				m(i, j) = Quadr::sum([i, j, &grads](auto iq) {
					auto g = grads(iq, i);
					auto basis = esf::Element_quadr<Phi_element, Quadr>::basis();
					return g[1] * basis(iq, j);
				});

		return m;
	}

	template<class Quadr>
	auto b_matrix()
	{
		esl::Matrix_d<n_dofs_chi, n_dofs_chi> m;
		for (std::size_t i = 0; i < n_dofs_chi; ++i)
			for (std::size_t j = 0; j < n_dofs_chi; ++j)
				m(i, j) = Quadr::sum([i, j](auto iq) {
					auto basis = esf::Element_quadr<Exy_element, Quadr>::basis();
					return basis(iq, i) * basis(iq, j);
				});

		return m;
	}

	template<class Quadr>
	auto load_vector()
	{
		return esl::make_vector<n_dofs_tau>([](auto i) {
			return Quadr::sum([i](std::size_t q) {
				auto constexpr basis = esf::Element_quadr<Phi_element, Quadr>::basis();
				return basis(q, i);
			});
		});
	}

	virtual void assemble() override
	{
		for (const auto& cell : mesh().cells())
			assemble_on_cell(cell);

		// {
		// 	const auto& v = system().variable<0>();
		// 	v.for_each_non_ess_bnd_cond([this](const auto& bc) {
		// 		for (const auto& halfedge : bc.halfedges())
		// 			assemble_on_halfedge(halfedge, bc);
		// 	});
		// }

		// {
		// 	const auto& v = system().variable<1>();
		// 	v.for_each_non_ess_bnd_cond([this](const auto& bc) {
		// 		for (const auto& halfedge : bc.halfedges())
		// 			assemble_on_halfedge(halfedge, bc);
		// 	});
		// }
	}

	void assemble_on_cell(const esf::Mesh2::Cell_view& cell)
	{
		const double area = esf::area(cell);

		// Precompute gradients of the basis functions
		using Q1 = esf::Quadr<Phi_element::order + Exy_element::order - 1, 2>;
		const auto grads = esf::gradients<Exy_element, Q1>(esf::inv_transp_jacobian(cell));

		const auto phi_dofs = system().dof_mapper().dofs<0>(cell);
		const auto exy_dofs = system().dof_mapper().dofs<1>(cell);

		// Assemble the local matrices and the load vector
		const auto ax = ax_matrix<Q1>(grads);
		const auto ay = ay_matrix<Q1>(grads);

		for (std::size_t r = 0; r < n_dofs_tau; ++r)
		{
			if (!phi_dofs[r].is_free)
				continue;

			auto i1 = phi_dofs[r].index;
			for (std::size_t c = 0; c < n_dofs_chi; ++c)
			{
				auto i2 = exy_dofs[c].index;
				matrix_(i1, i2) += -area * ax(c, r);
				matrix_(i1, i2 + 1) += -area * ay(c, r);

				matrix_(i2, i1) += area * ax(c, r);
				matrix_(i2 + 1, i1) += area * ay(c, r);
			}
		}

		using Q2 = esf::Quadr<2 * Exy_element::order, 2>;
		const auto b = b_matrix<Q2>();

		for (std::size_t r = 0; r < n_dofs_chi; ++r)
			for (std::size_t c = 0; c < n_dofs_chi; ++c)
			{
				auto i1 = exy_dofs[r].index;
				auto i2 = exy_dofs[c].index;

				matrix_(i1, i2) += area * b(r, c);
				matrix_(i1 + 1, i2 + 1) += area * b(r, c);
			}

		//const auto stiffness_matrix = fe::stiffness_matrix<Element, Stiff_quadr>(grads);
		//const auto mass_matrix = fe::mass_matrix<Element, Mass_quadr>();
		//const auto mass_matrix = fe::mass_matrix<Element, Mass_quadr>();

		using Q3 = esf::Quadr<Phi_element::order, 2>;
		const auto rhs = load_vector<Q3>();

		for (std::size_t r = 0; r < n_dofs_tau; ++r)
			if (phi_dofs[r].is_free)
				rhs_[phi_dofs[r].index] += area * rhs[r];

		// Add values to the global stiffness matrix
		// 		add_matrix_to_global(dofs, -area * stiffness_matrix);
		//		add_vector_to_global(dofs, area * load_vector);
	}

	// template<class Bc>
	// void assemble_on_halfedge(esf::Halfedge_index halfedge, const Bc& bc)
	// {
	// 	auto halfedge_view = mesh().view(halfedge);
	// 	auto face = halfedge_view.face_view();

	// 	const auto phi_dofs = system().dof_mapper().dofs<0>(face, halfedge);

	// 	const double edge_length = length(halfedge_view);

	// 	const double j = 0.5;
	// 	const double alpha = bc.alpha();
	// 	const double beta = bc.beta();

	// 	using Q1 = esf::Quadr<2 * Phi_element::order, 1>;
	// 	auto m = esl::make_matrix<n_dofs_tau, n_dofs_tau>(
	// 		[](auto i, auto j)
	// 		{
	// 			return Q1::sum(
	// 				[i, j](auto q)
	// 				{
	// 					auto const basis = esf::Element_quadr<Phi_element, Q1>::basis();
	// 					return basis(q, i) * basis(q, j);
	// 				});
	// 		});

	// 	for (esf::Local_index r = 0; r < n_dofs_tau; ++r)
	// 		for (esf::Local_index c = 0; c < n_dofs_tau; ++c)
	// 		{
	// 			auto i1 = phi_dofs[r].index;
	// 			auto i2 = phi_dofs[c].index;

	// 			matrix_(i1, i2) += j * alpha * edge_length * m(r, c);
	// 		}

	// 	using Q2 = esf::Quadr<Phi_element::order, 1>;
	// 	auto v = esl::make_vector<n_dofs_tau>(
	// 		[](auto i)
	// 		{
	// 			return Q2::sum(
	// 				[i](std::size_t q)
	// 				{
	// 					auto const basis = esf::Element_quadr<Phi_element, Q2>::basis();
	// 					return basis(q, i);
	// 				});
	// 		});

	// 	for (esf::Local_index r = 0; r < n_dofs_tau; ++r)
	// 		rhs_[phi_dofs[r].index] += j * beta * edge_length * v[r];
	// }

	virtual void after_solve() override
	{
		esl::Vector_xd phi(*mesh().n_faces(), 0);

		for (esf::Face_index face{0}; face < mesh().n_faces(); ++face)
		{
			typename Base::System::template Var_cell_dofs<0> cell_dofs;
			system().dof_mapper().template cell_dofs<0>(face, cell_dofs);
			phi[*face] = solution_[cell_dofs[0].index];
		}

		esf::Gnuplot_writer2 m("p0+p1.dat", mesh());
		m.write_face_field(phi);
	}

private:
	using Base::matrix_;
	using Base::rhs_;
	using Base::solution_;
};

} // namespace mixed_fem
