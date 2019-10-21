#pragma once

#include <esf/boundary_cond.hpp>
#include <esf/dof/dof_mapper.hpp>
#include <esf/element/lagrange.hpp>
#include <esf/system.hpp>
#include <esf/var.hpp>
#include <esf/var_list.hpp>

#include <string>

namespace mixed_fem
{
//using Phi = fe::Var<fe::Lagrange<1>, 1>;
//using Exy = fe::Var<fe::Discontinuous_lagrange<0>, 2, Dirichlet_boundary_cond, Robin_boundary_cond>;
using Phi =
	esf::Var<esf::Discontinuous_lagrange<0>, 1>; //, esf::Uniform_boundary_cond<esf::Lagrange<2>>>;
//using Exy = esf::Var<esf::Discontinuous_lagrange<0>, 2, esf::Uniform_boundary_cond<esf::Discontinuous_lagrange<0>>>;
using Exy = esf::Var<esf::Lagrange<1>, 2>;

class System final : public esf::System<esf::Var_list<Phi, Exy>, esf::Dof_mapper>
{
private:
	using Base = esf::System<esf::Var_list<Phi, Exy>, esf::Dof_mapper>;

public:
	using Base::Base;
};
} // namespace mixed_fem
