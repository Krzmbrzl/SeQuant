#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/sequant.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/wick.hpp>

#include "Utils.hpp"

#include <iostream>
#include <vector>

using namespace sequant;

ExprPtr T() {
  return ex<Constant>(1.0 / 2) *
         make_op(Tensor(
             L"t", std::vector<Index>{create_index(virt), create_index(virt)},
             std::vector<Index>{create_index(occ), create_index(occ)},
             Symmetry::antisymm));
}

ExprPtr Lambda() {
  Tensor TOp = (*T()->begin())->as<Tensor>();

  return ex<Constant>(1.0 / 2) *
         make_op(Tensor(L"Λ", TOp.ket(), TOp.bra(), Symmetry::antisymm));
}

int main() {
  set_default_context(SeQuant(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
                              BraKetSymmetry::conjugate, SPBasis::spinorbital));

  setConvention();

  ExprPtr energy = FWickTheorem{H() * T()}
                       .full_contractions(true)
                       .set_external_indices(std::vector<Index>{})
                       .compute();

  ExprPtr amplitudes =
      FWickTheorem{Lambda() * H() * T()}
          .full_contractions(true)
          .set_external_indices(std::vector<Index>{
              /* TODO: Can we use this to replace the Lambda operator? */})
          .compute();

  std::wcout << "E_CCD = " << to_latex(energy) << "\n\n";
  std::wcout << "Doubles residuum = " << to_latex(amplitudes) << "\n";
}
