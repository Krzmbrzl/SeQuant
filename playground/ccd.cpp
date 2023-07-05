#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/sequant.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/wick.hpp>

#include "Utils.hpp"

#include <chrono>
#include <iostream>
#include <vector>

using namespace sequant;

ExprPtr T() {
  return ex<Constant>(rational{1, 4}) *
         make_op(Tensor(
             L"t", std::vector<Index>{create_index(virt), create_index(virt)},
             std::vector<Index>{create_index(occ), create_index(occ)},
             Symmetry::antisymm));
}

ExprPtr Lambda() {
  Tensor TOp = (*T()->begin())->as<Tensor>();

  return ex<Constant>(rational{1, 4}) *
         make_op(Tensor(L"Λ", TOp.ket(), TOp.bra(), Symmetry::antisymm));
}

ExprPtr commutator(ExprPtr A, ExprPtr B) { return A * B - B * A; }

ExprPtr bch() {
  ExprPtr expr = H();

  ExprPtr comm = commutator(H(), T());
  simplify(comm);
  expr += comm;

  comm = commutator(commutator(H(), T()), T());
  simplify(comm);
  expr += ex<Constant>(rational{1, 2}) * comm;

  comm = commutator(commutator(commutator(H(), T()), T()), T());
  simplify(comm);
  expr += ex<Constant>(rational{1, 6}) * comm;

  comm = commutator(commutator(commutator(commutator(H(), T()), T()), T()), T());
  simplify(comm);
  expr += ex<Constant>(rational{1, 24}) * comm;

  return expr;
}

int main() {
  set_locale();
  set_default_context(SeQuant(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
                              BraKetSymmetry::conjugate, SPBasis::spinorbital));

  setConvention();

  ExprPtr pre_equations = (ex<Constant>(1) + Lambda()) * bch();
  expand(pre_equations);
  // std::wcout << pre_equations->size() << std::endl;
  // std::size_t counter = 1;
  // for (auto current : *pre_equations) {
  //  auto begin = std::chrono::steady_clock::now();
  //  std::wcout << "Term " << counter++ << ": " << to_latex(current) << "\n";
  //  std::wcout << "-> "
  //             << to_latex(FWickTheorem{current}
  //                             .full_contractions(true)
  //                             .set_external_indices(std::vector<Index>{})
  //                             .compute())
  //             << "\n";
  //  std::wcout << "  contracted in "
  //             << std::chrono::duration_cast<std::chrono::seconds>(
  //                    std::chrono::steady_clock::now() - begin)
  //                    .count()
  //             << "s\n\n";
  //}
  std::wcout << "Now in total...\n";
  auto begin = std::chrono::steady_clock::now();
  ExprPtr equations =
      FWickTheorem{pre_equations}
          .full_contractions(true)
          .set_external_indices(std::vector<Index>{
              /* TODO: Can we use this to replace the Lambda operator? */})
          .compute();

  std::wcout << "  Contracting together took "
             << std::chrono::duration_cast<std::chrono::seconds>(
                    std::chrono::steady_clock::now() - begin)
                    .count()
             << "s\n";

  simplify(equations);

  std::wcout << "CCD terms:\n" << to_latex_align(equations) << "\n";
}
