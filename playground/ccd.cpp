#include <SeQuant/core/context.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>

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

  // return ex<Constant>(rational{1, 4}) *
  //        make_op(Tensor(L"Λ", TOp.ket(), TOp.bra(), Symmetry::antisymm));

  // SeQuant treats a tensor with name "A" as an antisymmetrization operator
  return ex<Constant>(rational{1, 4}) *
         make_op(Tensor(L"A", TOp.ket(), TOp.bra(), Symmetry::antisymm));
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

  // comm = commutator(commutator(commutator(H(), T()), T()), T());
  // simplify(comm);
  // expr += ex<Constant>(rational{1, 6}) * comm;

  /*
        These don't contribute for CCD, but take quite a while to Wick in Debug
  builds comm = commutator(commutator(commutator(commutator(H(), T()), T()),
  T()), T()); simplify(comm); expr += ex<Constant>(rational{1, 24}) * comm;
  */

  return expr;
}

struct Idx2Size {
  std::size_t operator()(const Index &idx) const {
    if (idx.space() == occ) {
      return 10;
    } else if (idx.space() == virt) {
      return 100;
    } else if (idx.space() == active) {
      return 5;
    } else {
      throw std::runtime_error("Encountered unexpected space in Idx2Size");
    }
  }
};

int main() {
  set_locale();
  set_default_context(Context(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
                              BraKetSymmetry::conjugate, SPBasis::spinorbital));

  setConvention();

  //ExprPtr pre_equations = (ex<Constant>(1) + Lambda()) * bch();
  ExprPtr pre_equations = Lambda() * bch();
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
  std::wcout << "Contracting..." << std::endl;
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

  std::wcout << "CCD terms:\n" << to_latex_align(equations) << "\n\n\n";

  auto ext_indices = external_indices(equations);

  std::wcout << "Spintracing...\n";
  for (const ExprPtr &current : equations->as<Sum>()) {
    std::wcout << to_latex_align(current) /*<< "\ntraces out to\n" << to_latex_align(spintrace(current)) */
               << "\ntraces to\n" << to_latex_align(simplify(spintrace(current, ext_indices)))
               << "\n\n";
  }
  //ExprPtr spinTracedEqs = simplify(spintrace(equations));

  //std::wcout << "Spintraced CCD terms:\n"
  //           << to_latex_align(spinTracedEqs) << "\n";

  //ExprPtr optimizedEqs = simplify(optimize(spinTracedEqs, Idx2Size{}));

  //std::wcout << "Optimized CCD terms:\n"
  //           << to_latex_align(optimizedEqs) << "\n";
}
