#include <SeQuant/core/export/itf.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/utility/expr.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>
#include <SeQuant/domain/mbpt/sr.hpp>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <unordered_map>

using namespace sequant;
using namespace sequant::mbpt::sr;

struct Idx2Size {
  static const size_t nocc = 10;
  static const size_t nvirt = 100;
  size_t operator()(Index const& idx) const {
    if (idx.space() == IndexSpace::active_occupied)
      return nocc;
    else if (idx.space() == IndexSpace::active_unoccupied)
      return nvirt;
    else
      throw std::runtime_error("Unsupported IndexSpace type encountered");
  }
};

ExprPtr T(const std::vector<std::size_t>& projectionManifold) {
  auto T = std::make_shared<Sum>();
  for (std::size_t current : projectionManifold) {
    if (current == 0) {
      continue;
    }

    T->append(op::T_(current));
  }

  return T;
}

std::size_t getTermCount(const ExprPtr& expr) {
  if (!expr.is<Sum>()) {
    return 1;
  }

  return expr.as<Sum>().summands().size();
}

int main(int argc, const char** argv) {
  sequant::detail::OpIdRegistrar op_id_registrar;
  sequant::set_default_context(
      Context(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
              BraKetSymmetry::conjugate, SPBasis::spinorbital));
  mbpt::set_default_convention();
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  std::size_t maxExcitation = 2;
  bool includeSingles = true;
  if (argc == 2) {
    if (std::strcmp(argv[1], "ccsd") == 0) {
      maxExcitation = 2;
      includeSingles = true;
      std::wcout << "Generating equations for ccsd\n";
    } else if (std::strcmp(argv[1], "ccd") == 0) {
      maxExcitation = 2;
      includeSingles = false;
      std::wcout << "Generating equations for ccd\n";
    } else {
      std::wcout << "Unknown/Unsupported CC method: " << argv[1] << std::endl;
      return 1;
    }
  } else {
    maxExcitation = 2;
    includeSingles = true;
    std::wcout << "Generating equations for ccsd\n";
  }

  std::vector<std::size_t> projectionManifold(maxExcitation + 1);
  std::iota(projectionManifold.begin(), projectionManifold.end(), 0);
  if (!includeSingles && maxExcitation > 1) {
    projectionManifold.erase(projectionManifold.begin() + 1);
  }

  std::wcout << L"Chosen projection manifold is { ";
  for (std::size_t i = 0; i < projectionManifold.size(); ++i) {
    std::wcout << L"<" << std::to_wstring(projectionManifold[i]) << "|"
               << (i + 1 < projectionManifold.size() ? ", " : " ");
  }
  std::wcout << "}\n\n";

  // 1. Construct Hbar
  ExprPtr hbar = op::H();
  ExprPtr H_Tk = hbar;
  for (int64_t k = 1; k <= 4; ++k) {
    H_Tk =
        simplify(ex<Constant>(rational{1, k}) * H_Tk * T(projectionManifold));
    hbar += H_Tk;
  }

  std::wcout << L"Hbar:\n" << to_latex_align(hbar) << "\n\n";

  // 2. project onto doubles manifold, screen, lower to tensor form and wick it
  std::vector<ExprPtr> equations;

  for (std::size_t p : projectionManifold) {
    // 2.a. screen out terms that cannot give nonzero after projection onto
    // <p|
    std::shared_ptr<Sum> screendedTerms;

    for (const ExprPtr& term : *hbar) {
      assert(term->is<Product>() || term->is<op_t>());

      if (op::raises_vacuum_to_rank(term, p)) {
        if (!screendedTerms) {
          screendedTerms = std::make_shared<Sum>(ExprPtrList{term});
        } else {
          screendedTerms->append(term);
        }
      }
    }

    if (p > 0) {
      // 2.b project onto <p|, i.e. multiply by P(p)
      ExprPtr P_hbar = simplify(op::P(p) * screendedTerms);

      // 2.c compute vacuum expectation value (Wick theorem)
      equations.push_back(op::vac_av(P_hbar));
    } else {
      // Use equation as-is (no projection required)
      equations.push_back(op::vac_av(screendedTerms));
    }

    equations.back() = simplify(equations.back());
  }

  std::vector<itf::Result> results;

  for (std::size_t i = 0; i < projectionManifold.size(); ++i) {
    std::size_t currentProjection = projectionManifold.at(i);
    std::wcout << L"Equations for projection on <" << currentProjection
               << L"|:\n=============================\nRaw ("
               << getTermCount(equations.at(i)) << L"):\n"
               << to_latex_align(equations.at(i)) << "\n\n";

    // Spintrace
    equations[i] =
        simplify(closed_shell_CC_spintrace(equations[i], currentProjection));

    std::wcout << L"Spin-traced (" << getTermCount(equations[i]) << L"):\n"
               << to_latex_align(equations.at(i)) << "\n\n";

    // Remove symmetrization operator as this is not a tensor (but the optimize
    // function would treat it as such)
    // -> from here on the final symmetrization is implicit!
    equations[i] = remove_tensor(equations[i], L"S");

    // Optimize
    equations[i] = optimize(equations[i], Idx2Size{});

    std::wcout << L"Optimized (" << getTermCount(equations[i]) << L"):\n"
               << to_latex_align(equations[i]) << "\n\n";

    // Replace amplitudes with names as expected in ITF code
    equations[i]->visit(
        [](ExprPtr& expr) {
          if (!expr.is<Tensor>()) {
            return;
          }
          const Tensor& tensor = expr.as<Tensor>();
          if (tensor.label() != L"t") {
            return;
          }

          if (tensor.braket().size() == 2) {
            expr = ex<Tensor>(L"T1", tensor.bra(), tensor.ket());
          } else if (tensor.braket().size() == 4) {
            expr = ex<Tensor>(L"T2", tensor.bra(), tensor.ket());
          }
        },
        true);

    IndexGroups externals = non_repeated_indices(equations[i]);
    std::wstring resultName = [&]() -> std::wstring {
      if (currentProjection == 0) {
        return L"ECC";
      }
      return L"R" + std::to_wstring(currentProjection) +
             (currentProjection > 1 ? L"u" : L"");
    }();
    Tensor resultTensor(resultName, externals.bra, externals.ket);

    results.push_back(
        itf::Result(equations[i], resultTensor, currentProjection <= 1));

    if (resultName[resultName.size() - 1] == L'u') {
      // Generate symmetrization
      Tensor symmetrizedResult(resultName.substr(0, resultName.size() - 1),
                               resultTensor.bra(), resultTensor.ket());

      assert(externals.bra.size() == externals.ket.size());
      assert(externals.bra.size() == maxExcitation);

      ExprPtr symmetrization = ex<Sum>(ExprPtrList{});
      for (std::size_t i = 0; i < externals.bra.size(); ++i) {
        std::vector<Index> symBra;
        std::vector<Index> symKet;

        for (std::size_t j = 0; j < externals.bra.size(); ++j) {
          symBra.push_back(externals.bra[(i + j) % externals.bra.size()]);
          symKet.push_back(externals.ket[(i + j) % externals.ket.size()]);
        }

        symmetrization +=
            ex<Tensor>(resultName, std::move(symBra), std::move(symKet));
      }

      results.push_back(itf::Result(symmetrization, symmetrizedResult, true));
    }
  }

  // TODO: Factor out K4E

  std::wcout << "ITF code:\n\n"
             << to_itf(itf::CodeBlock(L"Residual", results)) << "\n";
}
