#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/context.hpp>
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
  set_locale();
  set_default_context(Context(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
                              BraKetSymmetry::conjugate, SPBasis::spinorbital));

  setConvention();

  constexpr std::size_t nAct = 2;

  auto expr = c0dagger(nAct) * H() * t() * c0(nAct);

  auto contracted = FWickTheorem{expr}
                        .set_external_indices(std::array<Index, 0>{})
                        .full_contractions(true)
                        .compute();

  std::wcout << "icMRCC equations:\n"
             << to_latex_align(contracted) << std::endl;

  for (auto &currentSummand : contracted->as<Sum>().summands()) {
    auto &currentProduct = currentSummand->as<Product>();

    auto c0Iter = std::find_if(
        currentProduct.factors().begin(), currentProduct.factors().end(),
        [](const ExprPtr &expr) {
          return expr.is<Tensor>() && expr.as<Tensor>().label() == L"{C_0}";
        });
    auto c0daggerIter =
        std::find_if(currentProduct.factors().begin(),
                     currentProduct.factors().end(), [](const ExprPtr &expr) {
                       return expr.is<Tensor>() &&
                              expr.as<Tensor>().label() == L"{C_0^\\dagger}";
                     });

    if (c0Iter == currentProduct.factors().end() ||
        c0daggerIter == currentProduct.factors().end()) {
      std::wcout << "Found term without density" << std::endl;
      continue;
    }

    auto braIndices = c0Iter->as<Tensor>().bra();
    auto ketIndices = c0daggerIter->as<Tensor>().ket();

    // Remove reference coefficients
    // TODO: The order of C0 and C0^+ may be different in general and we have to
    // look into iterator invalidation issues
    currentProduct.factors().erase(c0daggerIter);
    currentProduct.factors().erase(c0Iter);

    // instead, insert a density
    currentProduct.factors().push_back(
        ex<Tensor>(L"γ", braIndices, ketIndices, Symmetry::antisymm));
  }

  std::wcout << "icMRCC equations with densities:\n"
             << to_latex_align(contracted) << std::endl;
}
