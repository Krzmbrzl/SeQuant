//
// Created by Robert Adam on 2023-09-04
//

#include "itf.hpp"

#include <SeQuant/core/utility/expr.hpp>

#include <array>
#include <cassert>
#include <cwchar>
#include <map>
#include <unordered_map>

namespace sequant {

std::wstring to_itf(const itf::CodeBlock &block) {
  itf::detail::ITFGenerator generator;
  generator.addBlock(block);
  return generator.generate();
}

namespace itf {

Result::Result(ExprPtr expression, Tensor resultTensor)
    : expression(std::move(expression)),
      resultTensor(std::move(resultTensor)) {}

Tensor generateResultTensor(ExprPtr expr) {
  BraKet externals = non_repeated_indices(expr);

  return Tensor(L"Result", std::move(externals.bra), std::move(externals.ket));
}

Result::Result(ExprPtr expression)
    : expression(std::move(expression)),
      resultTensor(generateResultTensor(expression)) {}

CodeBlock::CodeBlock(std::wstring blockName, Result result)
    : CodeBlock(std::move(blockName), std::vector<Result>{std::move(result)}) {}

CodeBlock::CodeBlock(std::wstring blockName, std::vector<Result> results)
    : name(std::move(blockName)), results(std::move(results)) {}

namespace detail {

std::vector<Contraction> to_contractions(const ExprPtr &expression,
                                         const Tensor &resultTensor);

std::vector<Contraction> to_contractions(const Product &product,
                                         const Tensor &resultTensor) {
  static std::size_t intermediateCounter = 1;

  if (product.factors().size() == 1) {
    assert(product.factor(0).is<Tensor>());
    assert(product.scalar().imag() == 0);

    return {Contraction{product.scalar().real(), resultTensor,
                        product.factor(0).as<Tensor>(),
                        Tensor(L"One", {}, {})}};
  }

  // We assume that we're dealing with a binary tree
  assert(product.factors().size() == 2);

  std::unordered_map<const ExprPtr::element_type *, Tensor> intermediates;
  std::vector<Contraction> contractions;

  // Handle intermediates
  for (const ExprPtr &factor : product.factors()) {
    if (factor.is<Product>()) {
      // Create intermediate that computes this nested product
      BraKet intermediateIndices = non_repeated_indices(factor);

      std::array<wchar_t, 64> intermediateName;
      swprintf(intermediateName.data(), intermediateName.size(), L"INTER%06u",
               intermediateCounter++);

      Tensor intermediate(intermediateName.data(),
                          std::move(intermediateIndices.bra),
                          std::move(intermediateIndices.ket));

      std::vector<Contraction> intermediateContractions =
          to_contractions(factor, intermediate);
      contractions.reserve(contractions.size() +
                           intermediateContractions.size());
      contractions.insert(
          contractions.end(),
          std::make_move_iterator(intermediateContractions.begin()),
          std::make_move_iterator(intermediateContractions.end()));

      intermediates.insert({factor.get(), std::move(intermediate)});
    } else if (factor.is<Sum>()) {
      // TODO: Handle on-the-fly antisymmetrization (K[abij] - K[baij])
      throw std::invalid_argument(
          "Products of sums can not yet be translated to ITF");
    }
  }

  // Now create the contraction for the two factors
  auto lhsIntermediate = intermediates.find(product.factor(0).get());
  auto rhsIntermediate = intermediates.find(product.factor(1).get());

  assert(product.scalar().imag() == 0);
  contractions.push_back(Contraction{
      product.scalar().real(), resultTensor,
      lhsIntermediate == intermediates.end() ? product.factor(0).as<Tensor>()
                                             : lhsIntermediate->second,
      rhsIntermediate == intermediates.end() ? product.factor(1).as<Tensor>()
                                             : rhsIntermediate->second});

  return contractions;
}

std::vector<Contraction> to_contractions(const ExprPtr &expression,
                                         const Tensor &resultTensor) {
  std::wstring itfCode;

  if (expression.is<Constant>()) {
    throw std::invalid_argument("Can't transform constants into contractions");
  } else if (expression.is<Tensor>()) {
    // Use the special One[] tensor that is will result in a no-op contraction
    return {Contraction{1, resultTensor, expression.as<Tensor>(),
                        Tensor(L"One", {}, {})}};
  } else if (expression.is<Product>()) {
    // Separate into binary contractions
    return to_contractions(expression.as<Product>(), resultTensor);
  } else if (expression.is<Sum>()) {
    // Process each summand
    std::vector<Contraction> contractions;

    for (const ExprPtr &summand : expression.as<Sum>().summands()) {
      std::vector<Contraction> currentContractions =
          to_contractions(summand, resultTensor);

      contractions.reserve(contractions.size() + currentContractions.size());
      contractions.insert(contractions.end(),
                          std::make_move_iterator(currentContractions.begin()),
                          std::make_move_iterator(currentContractions.end()));
    }

    return contractions;
  } else {
    throw std::invalid_argument(
        "Unhandled expression type in to_contractions function");
  }
}

ExprPtr replaceTwoElectronIntgrals(const ExprPtr expr) {
  ExprPtr expression = expr->clone();

  expression->visit(
      [](ExprPtr &expr) {
        // TODO: Get rid of hardcoded label
        if (expr.is<Tensor>() && expr.as<Tensor>().label() == L"g") {
          const Tensor &tensor = expr.as<Tensor>();
          assert(tensor.bra().size() == 2);
          assert(tensor.ket().size() == 2);

          bool firstParticleSameSpaces =
              tensor.bra()[0].space() == tensor.ket()[0].space();
          bool secondParticleSameSpaces =
              tensor.bra()[1].space() == tensor.ket()[1].space();

          if ((firstParticleSameSpaces && secondParticleSameSpaces &&
               tensor.bra()[0].space() == tensor.ket()[0].space()) ||
              (!firstParticleSameSpaces && !secondParticleSameSpaces)) {
            // Integrals are stored on the K tensor in Molpro, if
            // - all indices are of the same space
            // - the indices for particle one and those for particle two belong
            // to different spaces
            expr = ex<Tensor>(L"K", tensor.bra(), tensor.ket());
          } else if (firstParticleSameSpaces && secondParticleSameSpaces) {
            // Integrals are stored on the J tensor in Molpro, if
            // - the indices for particle one and those for particle two belong
            // to the same space

            // Note the exchange of the second and third index to store the
            // tensor in a form in which the first two indices belong to the
            // same space
            expr = ex<Tensor>(
                L"J", std::vector<Index>{tensor.bra()[0], tensor.ket()[0]},
                std::vector<Index>{tensor.bra()[1], tensor.ket()[1]});
          } else {
			  // TODO: Handle cases with index spaces 3-to-1 (e.g. K:eccc)
            throw std::runtime_error(
                "Encountered (yet) unsupported index configuration on "
                "g-tensor");
          }
        }
      },
      true);

  return expression;
}

void ITFGenerator::addBlock(const itf::CodeBlock &block) {
  m_codes.reserve(m_codes.size() + block.results.size());

  for (const Result &currentResult : block.results) {
    ExprPtr expression = replaceTwoElectronIntgrals(currentResult.expression);

    m_codes.push_back(ContractionBlock{
        block.name, to_contractions(expression, currentResult.resultTensor)});

    // Assume that all tensors that appear in the expressions as well as the
    // results that are calculated are tensors that are imported rather than
    // created from scratch.
    m_importedTensors.insert(currentResult.resultTensor);
    expression->visit(
        [this](const ExprPtr &expr) {
          if (expr.is<Tensor>()) {
            const Tensor &tensor = expr.as<Tensor>();
            m_importedTensors.insert(tensor);
            m_encounteredIndices.insert(tensor.braket().begin(),
                                        tensor.braket().end());
          }
        },
        true);

    // Now go through all result tensors of the contractions that we have
    // produced and add all new tensors to the set of created tensors
    for (const Contraction &currentContraction : m_codes.back().contractions) {
      if (m_importedTensors.find(currentContraction.result) ==
          m_importedTensors.end()) {
        m_createdTensors.insert(currentContraction.result);
      }
    }
  }
}

struct IndexComponents {
  IndexSpace space;
  std::size_t id;
};

IndexComponents decomposeIndex(const Index &idx) {
  // The labels are of the form <letter>_<number> and we want <number>
  // Note that SeQuant uses 1-based indexing, but we want 0-based
  int num =
      std::stoi(std::wstring(idx.label().substr(idx.label().find('_') + 1))) -
      1;
  assert(num >= 0 && num <= 6);

  return {idx.space(), static_cast<std::size_t>(num)};
}

std::map<IndexSpace, std::set<std::size_t>> indicesBySpace(
    const std::set<Index> &indices) {
  std::map<IndexSpace, std::set<std::size_t>> indexMap;

  for (const Index &current : indices) {
    IndexComponents components = decomposeIndex(current);

    indexMap[components.space].insert(components.id);
  }

  return indexMap;
}

std::wstring to_itf(const Tensor &tensor, bool includeIndexing = true) {
  std::wstring tags;
  std::wstring indices;

  for (const Index &current : tensor.braket()) {
    IndexComponents components = decomposeIndex(current);

    assert(components.id <= 7);

    if (components.space.type() == IndexSpace::active_occupied) {
      tags += L"c";
      indices += static_cast<wchar_t>(L'i' + components.id);
    } else if (components.space.type() == IndexSpace::active_unoccupied) {
      tags += L"e";
      indices += static_cast<wchar_t>(L'a' + components.id);
    } else {
      std::runtime_error("Encountered unhandled index space type");
    }
  }

  return std::wstring(tensor.label()) + (tags.empty() ? L"" : L":" + tags) +
         (includeIndexing ? L"[" + indices + L"]" : L"");
}

std::wstring ITFGenerator::generate() const {
  std::wstring itf =
      L"// This ITF algo file has been generated via SeQuant's ITF export\n\n";

  itf += L"----decl\n";

  // Index declarations
  std::map<IndexSpace, std::set<std::size_t>> indexGroups =
      indicesBySpace(m_encounteredIndices);
  for (auto iter = indexGroups.begin(); iter != indexGroups.end(); ++iter) {
    wchar_t baseLabel;
    std::wstring spaceLabel;
    std::wstring spaceTag;
    if (iter->first.type() == IndexSpace::active_occupied) {
      baseLabel = L'i';
      spaceLabel = L"Closed";
      spaceTag = L"c";
    } else if (iter->first.type() == IndexSpace::active_unoccupied) {
      baseLabel = L'a';
      spaceLabel = L"External";
      spaceTag = L"e";
    } else {
      std::runtime_error("Encountered unhandled index space type");
    }

    itf += L"index-space: ";
    for (std::size_t i : iter->second) {
      assert(i <= 7);
      itf += static_cast<wchar_t>(baseLabel + i);
    }
    itf += L", " + spaceLabel + L", " + spaceTag + L"\n";
  }

  itf += L"\n";

  // Tensor declarations
  for (const Tensor &current : m_importedTensors) {
    itf +=
        L"tensor: " + to_itf(current) + L", " + to_itf(current, false) + L"\n";
  }
  itf += L"\n";
  for (const Tensor &current : m_createdTensors) {
    itf += L"tensor: " + to_itf(current) + L", !Create{type:plain}\n";
  }
  itf += L"\n\n";

  // Actual code
  for (const ContractionBlock &currentBlock : m_codes) {
    itf += L"---- code(\"" + currentBlock.name + L"\")\n";

    std::set<Tensor, TensorBlockCompare> allocatedTensors;

    for (const Contraction &currentContraction : currentBlock.contractions) {
      // For now we'll do a really silly contribution-by-contribution
      // load-process-store strategy
      if (allocatedTensors.find(currentContraction.result) ==
          allocatedTensors.end()) {
        itf += L"alloc " + to_itf(currentContraction.result) + L"\n";
        allocatedTensors.insert(currentContraction.result);
      } else {
        itf += L"load " + to_itf(currentContraction.result) + L"\n";
      }
      itf += L"load " + to_itf(currentContraction.lhs) + L"\n";
      itf += L"load " + to_itf(currentContraction.rhs) + L"\n";

      itf += L"." + to_itf(currentContraction.result) + L" ";
      int sign = currentContraction.factor < 0 ? -1 : 1;

      itf += (sign < 0 ? L"-= " : L"+= ");
      if (currentContraction.factor * sign != 1) {
        itf += to_wstring(currentContraction.factor * sign) + L"*";
      }
      itf += to_itf(currentContraction.lhs) + L" " +
             to_itf(currentContraction.rhs) + L"\n";

      itf += L"drop " + to_itf(currentContraction.rhs) + L"\n";
      itf += L"drop " + to_itf(currentContraction.lhs) + L"\n";
      itf += L"store " + to_itf(currentContraction.result) + L"\n";
    }

    itf += L"\n";
  }

  return itf;
}

}  // namespace detail

}  // namespace itf

}  // namespace sequant
