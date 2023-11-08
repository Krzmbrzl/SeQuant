#include <SeQuant/core/export/itf.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/parse_expr.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/tensor_canonicalizer.hpp>
#include <SeQuant/core/utility/indices.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>
#include <SeQuant/domain/mbpt/sr.hpp>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
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

void defaultSetup() {
  sequant::detail::OpIdRegistrar op_id_registrar;
  sequant::set_default_context(
      Context(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
              BraKetSymmetry::conjugate, SPBasis::spinorbital));
  mbpt::set_default_convention();
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());
}

struct Options {
  std::size_t maxExcitation;
  bool includeSingles;
  std::string name;
};

Options parseOptions(int argc, const char** argv) {
  Options options;
  if (argc == 2) {
    if (std::strcmp(argv[1], "ccsd") == 0) {
      options.maxExcitation = 2;
      options.includeSingles = true;
      options.name = "ccsd";
    } else if (std::strcmp(argv[1], "ccd") == 0) {
      options.maxExcitation = 2;
      options.includeSingles = false;
      options.name = "ccd";
    } else {
      throw std::runtime_error(std::string("Unknown/Unsupported CC method: ") +
                               argv[1]);
    }
  } else {
    options.maxExcitation = 2;
    options.includeSingles = true;
    options.name = "ccsd";
  }
  std::wcout << "Generating equations for " << options.name.c_str() << "\n";

  return options;
}

std::vector<std::size_t> createProjectionManifold(std::size_t maxExcitation,
                                                  bool includeSingles) {
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

  return projectionManifold;
}

ExprPtr Hbar(const std::vector<std::size_t>& projectionManifold) {
  ExprPtr hbar = op::H();
  ExprPtr H_Tk = hbar;
  for (std::size_t k = 1; k <= 4; ++k) {
    H_Tk =
        simplify(ex<Constant>(rational{1, k}) * H_Tk * T(projectionManifold));
    hbar += H_Tk;
  }

  return hbar;
}

std::vector<ExprPtr> generateWorkingEquations(
    const std::vector<std::size_t>& projectionManifold, const ExprPtr& hbar) {
  // Project onto all functions in the projections manifold (<0|, <1|, etc.),
  // perform a pre-screening to filter out terms of which we know that they must
  // be zero and then apply Wick's theorem in order to obtain tensor
  // expressions.
  std::vector<ExprPtr> equations;

  for (std::size_t p : projectionManifold) {
    // Screen out terms that cannot give nonzero after projection onto <p|
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
      // Project onto <p|, i.e. multiply by P(p)
      ExprPtr P_hbar = simplify(op::P(p) * screendedTerms);

      // Compute vacuum expectation value (Wick theorem)
      equations.push_back(op::vac_av(P_hbar));
    } else {
      // Use equation as-is (no projection required)
      equations.push_back(op::vac_av(screendedTerms));
    }

    equations.back() = simplify(equations.back());
  }

  return equations;
}

ExprPtr processExpression(const ExprPtr& expr) {
  ExprPtr result;
  if (expr.is<Sum>()) {
    result = expr;
  } else {
    // The spintracing routine expects a sum
    result = ex<Sum>(ExprPtrList{expr});
  }

  // Spintrace
  result = simplify(closed_shell_CC_spintrace(result));

  // Remove symmetrization operator as this is not a tensor (but the
  // optimize function would treat it as such)
  // -> from here on the final symmetrization is implicit!
  result = remove_tensor(result, L"S");

  // Optimize
  result = optimize(result, Idx2Size{});

  return result;
}

itf::Result toItfResult(ExprPtr& expr, std::size_t projectionLevel) {
  // Replace amplitudes with names as expected in ITF code
  expr->visit(
      [](ExprPtr& expr) {
        if (!expr.is<Tensor>()) {
          return;
        }
        const Tensor& tensor = expr.as<Tensor>();
        if (tensor.label() != L"t") {
          return;
        }

        if (tensor.braket().size() == 2) {
          expr = ex<Tensor>(L"T1", tensor.bra(), tensor.ket(), tensor.auxiliary());
        } else if (tensor.braket().size() == 4) {
          expr = ex<Tensor>(L"T2", tensor.bra(), tensor.ket(), tensor.auxiliary());
        }
      },
      true);

  // Assemble the result tensor name
  std::wstring resultName = [&]() -> std::wstring {
    if (projectionLevel == 0) {
      return L"ECC";
    }
    return L"R" + std::to_wstring(projectionLevel) +
           (projectionLevel > 1 ? L"u" : L"");
  }();

  // Assemble the result tensor itself
  IndexGroups externals = get_unique_indices(expr);
  Tensor resultTensor(resultName, externals.bra, externals.ket, externals.aux);

  return itf::Result(expr, resultTensor, projectionLevel <= 1);
}

std::optional<itf::Result> generateResultSymmetrization(
    const Tensor& resultTensor) {
  std::wstring_view resultName = resultTensor.label();
  if (resultName[resultName.size() - 1] != L'u') {
    return {};
  }

  Tensor symmetrizedResult(resultName.substr(0, resultName.size() - 1),
                           resultTensor.bra(), resultTensor.ket(),
                           resultTensor.auxiliary());

  assert(resultTensor.bra_rank() == resultTensor.ket_rank());

  // Note: we're only symmetrizing over bra-ket, not over auxiliary indices
  ExprPtr symmetrization = ex<Sum>(ExprPtrList{});
  for (std::size_t i = 0; i < resultTensor.bra_rank(); ++i) {
    std::vector<Index> symBra;
    std::vector<Index> symKet;

    for (std::size_t j = 0; j < resultTensor.ket_rank(); ++j) {
      symBra.push_back(resultTensor.bra()[(i + j) % resultTensor.bra_rank()]);
      symKet.push_back(resultTensor.ket()[(i + j) % resultTensor.ket_rank()]);
    }

    symmetrization += ex<Tensor>(resultName, std::move(symBra),
                                 std::move(symKet), resultTensor.auxiliary());
  }

  // We assume that the symmetrized result always is an ITF-internal tensor that
  // we have to import rather than create
  return itf::Result(symmetrization, symmetrizedResult, true);
}

int main(int argc, const char** argv) {
  try {
    defaultSetup();
    Options options = parseOptions(argc, argv);

    // Generate desired equations
    std::vector<std::size_t> projectionManifold =
        createProjectionManifold(options.maxExcitation, options.includeSingles);

    ExprPtr hbar = Hbar(projectionManifold);
    std::wcout << L"Hbar:\n" << to_latex_align(hbar) << "\n\n";

    std::vector<ExprPtr> equations =
        generateWorkingEquations(projectionManifold, hbar);

    // Process generated equations
    std::vector<itf::Result> results;

    for (std::size_t i = 0; i < projectionManifold.size(); ++i) {
      const std::size_t currentProjection = projectionManifold.at(i);
      ExprPtr processedExpr = processExpression(equations.at(i));

      std::wcout << L"Equations for projection on <" << currentProjection
                 << L"|:\n=============================\nRaw ("
                 << getTermCount(processedExpr) << L"):\n"
                 << to_latex_align(processedExpr) << "\n\n";

      // Transform the current equation into a form that can later be translated
      // to ITF
      itf::Result result = toItfResult(processedExpr, currentProjection);
      results.push_back(result);

      // If necessary, generate the symmetrization code for the result computed
      // by processedExpr.
      std::optional<itf::Result> resultSymmetrization =
          generateResultSymmetrization(result.resultTensor);
      if (resultSymmetrization.has_value()) {
        results.push_back(resultSymmetrization.value());
      }
    }

    // TODO: Factor out K4E

    // Translate the intermediary representation into actual ITF code and print
    // that
    std::wfstream stream(options.name + ".itfaa", std::fstream::out);
    stream << to_itf(itf::CodeBlock(L"Residual", results));
    stream.close();
    std::wcout << "ITF code written to file " << options.name.c_str()
               << ".itfaa\n";

  } catch (const std::exception& e) {
    std::wcout << "[ERROR]: " << e.what() << std::endl;
    return 1;
  }
}
