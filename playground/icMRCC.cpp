#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/sequant.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/wick.hpp>

#include "Utils.hpp"

#include <iostream>

using namespace sequant;

ExprPtr c0(std::size_t nAct) {
  std::vector<Index> indices;
  for (std::size_t i = 0; i < nAct; ++i) {
    indices.push_back(create_index(active));
  }

  return ex<Constant>(rational{1, 2}) *
         make_op(Tensor(L"{C_0}", indices, std::vector<Index>{},
                        Symmetry::antisymm));
}

ExprPtr c0dagger(std::size_t nAct) {
  std::vector<Index> indices;
  for (std::size_t i = 0; i < nAct; ++i) {
    indices.push_back(create_index(active));
  }

  return ex<Constant>(rational{1, 2}) *
         make_op(Tensor(L"{C_0^\\dagger}", std::vector<Index>{}, indices,
                        Symmetry::antisymm));
}

ExprPtr t() {
  return make_op(Tensor(L"t", std::vector<Index>{create_index(active)},
                        std::vector<Index>{create_index(occ)},
                        Symmetry::antisymm));
}

int main() {
  set_default_context(SeQuant(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
                              BraKetSymmetry::conjugate, SPBasis::spinorbital));

  setConvention();

  constexpr std::size_t nAct = 2;

  auto expr = c0dagger(nAct) * H() * t() * c0(nAct);

  auto contracted = FWickTheorem{expr}
                        .set_external_indices(std::array<Index, 0>{})
                        .full_contractions(true)
                        .compute();

  std::wcout << to_latex(expr) << "\n=\n" << to_latex(contracted) << std::endl;
}
