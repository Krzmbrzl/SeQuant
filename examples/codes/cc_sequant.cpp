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
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>

using namespace sequant;
using namespace sequant::mbpt::sr;

struct Idx2Size {
  static const size_t nocc = 10;
  static const size_t nvirt = 100;
  static const size_t naux = 100;
  static const size_t nact = 4;
  size_t operator()(Index const& idx) const {
    if (idx.space() == IndexSpace::active_occupied ||
        idx.space() == IndexSpace::occupied)
      return nocc;
    else if (idx.space() == IndexSpace::active_unoccupied ||
             idx.space() == IndexSpace::unoccupied)
      return nvirt;
    else if (idx.space() == IndexSpace::active)
		return nact;
    else if (idx.space() == IndexSpace::all_active)
      return naux;
    else
      throw std::runtime_error("Unsupported IndexSpace type encountered in Idx2Size");
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
  bool densityFitting;
  bool termByTerm;
  std::optional<std::string> inputFile;
};

Options parseOptions(int argc, const char** argv) {
  Options options;
  if (argc >= 2) {
    std::string name(argv[1]);

    options.densityFitting = name.size() > 3 && name.substr(0, 3) == "df-";
    if (options.densityFitting) {
      name = name.substr(3);
    }

    options.termByTerm =
        name.size() > 5 && name.substr(name.size() - 4, 4) == "_tbt";
    if (options.termByTerm) {
      name = name.substr(0, name.size() - 4);
    }

    if (name == "ccsd") {
      options.maxExcitation = 2;
      options.includeSingles = true;
      options.name = argv[1];
    } else if (name == "ccd") {
      options.maxExcitation = 2;
      options.includeSingles = false;
      options.name = argv[1];
    } else if (name == "read") {
      if (argc < 3) {
        throw std::runtime_error("Missing argument for read");
      }

      options.inputFile = argv[2];
      name = argv[2];

      if (auto idx = name.rfind("."); idx != std::string::npos && idx > 1) {
        name = name.substr(0, idx);
      }
      if (options.termByTerm) {
        name += "_tbt";
      }
      if (options.densityFitting) {
        name = "df-" + name;
      }
      options.name = name;
    } else {
      throw std::runtime_error(std::string("Unknown/Unsupported CC method: ") +
                               argv[1]);
    }
  } else {
    options.maxExcitation = 2;
    options.includeSingles = true;
    options.name = "ccsd";
    options.densityFitting = false;
    options.termByTerm = false;
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

IndexGroups<std::vector<Index>> external_indices(const ExprPtr& expr) {
  IndexGroups groups;

  if (!expr.is<Sum>() && !expr.is<Product>()) {
    return groups;
  }

  const Tensor& symmetrizer = [&]() {
    if (expr.is<Sum>()) {
      return expr.as<Sum>().summand(0).as<Product>().factor(0).as<Tensor>();
    } else {
      return expr.as<Product>().factor(0).as<Tensor>();
    }
  }();

  if (symmetrizer.label() == L"A" || symmetrizer.label() == L"S") {
    groups.bra.insert(groups.bra.end(), symmetrizer.bra().begin(),
                      symmetrizer.bra().end());
    groups.ket.insert(groups.ket.end(), symmetrizer.ket().begin(),
                      symmetrizer.ket().end());
  }

  return groups;
}

std::tuple<ExprPtr, IndexGroups<std::vector<Index>>> processExpression(
    const ExprPtr& expr) {
  ExprPtr result;
  if (expr.is<Sum>()) {
    result = expr;
  } else {
    // The spintracing routine expects a sum
    result = ex<Sum>(ExprPtrList{expr});
  }

  // Spintrace
  //result = simplify(closed_shell_CC_spintrace(result));
  result = simplify(spintrace(result));

  IndexGroups externals = external_indices(result);

  // Remove symmetrization operator as this is not a tensor (but the
  // optimize function would treat it as such)
  // -> from here on the final symmetrization is implicit!
  result = remove_tensor(result, L"S");

  // Optimize
  result = optimize(result, Idx2Size{});

  return {result, externals};
}

itf::Result toItfResult(ExprPtr& expr, std::size_t projectionLevel,
                        const IndexGroups<std::vector<Index>>& externals) {
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
          expr =
              ex<Tensor>(L"T1", tensor.bra(), tensor.ket(), tensor.auxiliary());
        } else if (tensor.braket().size() == 4) {
          expr =
              ex<Tensor>(L"T2", tensor.bra(), tensor.ket(), tensor.auxiliary());
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

ExprPtr insertDensityFitting(ExprPtr expr) {
  expr->visit(
      [](ExprPtr& current) {
        if (!current->is<Tensor>()) {
          return;
        }
        Tensor& tensor = current.as<Tensor>();
        if (tensor.label() != L"g" || tensor.bra_rank() != 2 ||
            tensor.ket_rank() != 2 || tensor.auxiliary_rank() != 0) {
          return;
        }

        bool allVirtual = true;
        for (const Index& idx : tensor.const_indices()) {
          if (idx.space() != IndexSpace::instance(L"a_1")) {
            allVirtual = false;
            break;
          }
        }

        // Don't replace 4ext integrals - we want to make use of KExt
        if (allVirtual && false) {
          return;
        }

        const auto& braIdx = tensor.bra();
        const auto& ketIdx = tensor.ket();

        // IndexFactory factory;
        // Index contractionIdx = factory.make(Index(L"x_1"));
        Index contractionIdx(L"x_1");

        // std::wcout << "Before: " << deparse_expr(current) << "\n";

        ExprPtr df1 = ex<Tensor>(L"DF", std::vector<Index>{braIdx[0]},
                                 std::vector<Index>{ketIdx[0]},
                                 std::vector<Index>{contractionIdx},
                                 Symmetry::nonsymm, BraKetSymmetry::symm);
        ExprPtr df2 = ex<Tensor>(L"DF", std::vector<Index>{braIdx[1]},
                                 std::vector<Index>{ketIdx[1]},
                                 std::vector<Index>{contractionIdx},
                                 Symmetry::nonsymm, BraKetSymmetry::symm);
        ExprPtr df3 = ex<Tensor>(L"DF", std::vector<Index>{braIdx[0]},
                                 std::vector<Index>{ketIdx[1]},
                                 std::vector<Index>{contractionIdx},
                                 Symmetry::nonsymm, BraKetSymmetry::symm);
        ExprPtr df4 = ex<Tensor>(L"DF", std::vector<Index>{braIdx[1]},
                                 std::vector<Index>{ketIdx[0]},
                                 std::vector<Index>{contractionIdx},
                                 Symmetry::nonsymm, BraKetSymmetry::symm);

        if (tensor.symmetry() == Symmetry::antisymm) {
          // Antisymmetric exchange
          current =
              ex<Sum>(ExprPtrList{ex<Product>(ExprPtrList{df1, df2}),
                                  ex<Product>(-1, ExprPtrList{df3, df4})});
        } else {
          // Symmetric exchange
          current = ex<Product>(ExprPtrList{df1, df2});
        }

        // std::wcout << "Inserted DF: " << deparse_expr(current) << "\n";
      },
      true);

  return expand(expr);
}

itf::Result processToItf(ExprPtr expr, std::size_t projectionLevel,
                         bool useDensityFitting) {
  if (useDensityFitting) {
    expr = insertDensityFitting(expr);
  }

  auto [processed, externals] = processExpression(expr);

  return toItfResult(processed, projectionLevel, externals);
}

int main(int argc, const char** argv) {
  try {
    defaultSetup();
    Options options = parseOptions(argc, argv);

    std::vector<ExprPtr> equations;
    std::vector<std::size_t> projectionManifold;

    if (options.inputFile.has_value()) {
      // Read in equations
      std::ifstream input(options.inputFile.value());
      std::string line;
      std::size_t currentProjection;
      std::string currentBlock;

      while (std::getline(input, line)) {
        if (line.find("level:") == 0) {
          if (!currentBlock.empty()) {
            equations.push_back(
                parse_expr(to_wstring(currentBlock), Symmetry::antisymm));
            projectionManifold.push_back(currentProjection);
            currentBlock.clear();
          }
          std::stringstream(line.substr(6)) >> currentProjection;
        } else {
          currentBlock += "\n" + line;
        }
      }

      if (!currentBlock.empty()) {
        equations.push_back(
            parse_expr(to_wstring(currentBlock), Symmetry::antisymm));
        projectionManifold.push_back(currentProjection);
        currentBlock.clear();
      }
    } else {
      // Generate desired equations
      projectionManifold = createProjectionManifold(options.maxExcitation,
                                                    options.includeSingles);

      ExprPtr hbar = Hbar(projectionManifold);
      std::wcout << L"Hbar:\n" << to_latex_align(hbar) << "\n\n";

      equations = generateWorkingEquations(projectionManifold, hbar);
    }

    assert(equations.size() <= projectionManifold.size());

    // Process generated equations
    std::vector<itf::Result> results;

    for (std::size_t i = 0; i < equations.size(); ++i) {
      const std::size_t currentProjection = projectionManifold.at(i);

      std::wcout << "Processing equations for projection on <"
                 << currentProjection << "|\n";

      if (options.termByTerm) {
        if (!equations[i].is<Sum>()) {
          equations[i] = ex<Sum>(ExprPtrList{equations[i]});
        }

        Sum& sum = equations[i].as<Sum>();
        for (std::size_t k = 0; k < sum.size(); ++k) {
          std::wcout << "Term #" << (k + 1) << ":\n  " << deparse_expr(sum.summand(k)) << "\n  processes to\n";

          itf::Result result = processToItf(sum.summand(k), currentProjection,
                                            options.densityFitting);

		  std::wcout << "  " << deparse_expr(ex<Tensor>(result.resultTensor))
                     << " += " << deparse_expr(result.expression) << "\n\n";

          results.push_back(std::move(result));
        }
      } else {
        itf::Result result = processToItf(equations[i], currentProjection,
                                          options.densityFitting);

        std::wcout << deparse_expr(ex<Tensor>(result.resultTensor)) << " = "
                   << deparse_expr(result.expression) << "\n\n";

        results.push_back(std::move(result));
      }

      std::optional<itf::Result> symmetrization =
          generateResultSymmetrization(results.back().resultTensor);
      if (symmetrization.has_value()) {
        results.push_back(std::move(symmetrization.value()));
      }

	  std::wcout << "\n\n";
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
