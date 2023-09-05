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

Result::Result(ExprPtr expression, Tensor resultTensor, bool importResultTensor)
    : expression(std::move(expression)),
      resultTensor(std::move(resultTensor)),
      importResultTensor(importResultTensor) {}

Tensor generateResultTensor(ExprPtr expr) {
  BraKet externals = non_repeated_indices(expr);

  return Tensor(L"Result", std::move(externals.bra), std::move(externals.ket));
}

Result::Result(ExprPtr expression, bool importResultTensor)
    : expression(std::move(expression)),
      resultTensor(generateResultTensor(expression)),
      importResultTensor(importResultTensor) {}

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

    return {Contraction{product.scalar().real(),
                        resultTensor,
                        product.factor(0).as<Tensor>(),
                        {}}};
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
    return {Contraction{1, resultTensor, expression.as<Tensor>(), {}}};
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

template <typename IndexContainer, typename SpaceTypeContainer>
bool isSpacePattern(const IndexContainer &indices,
                    const SpaceTypeContainer &pattern) {
  static_assert(std::is_same_v<typename IndexContainer::value_type, Index>);
  static_assert(std::is_same_v<typename SpaceTypeContainer::value_type,
                               IndexSpace::Type>);
  assert(indices.size() == pattern.size());

  auto indexIter = indices.cbegin();
  auto typeIter = pattern.cbegin();

  while (indexIter != indices.cend() && typeIter != pattern.cend()) {
    if (indexIter->space().type() != *typeIter) {
      return false;
    }

    ++indexIter;
    ++typeIter;
  }

  return true;
}

struct IndexTypeComparer {
  constexpr bool operator()(const IndexSpace::Type &lhs,
                            const IndexSpace::Type &rhs) const {
    static_assert(IndexSpace::active_occupied < IndexSpace::active_unoccupied);
    return lhs < rhs;
  }
};

ExprPtr replaceTwoElectronIntgrals(const ExprPtr expr) {
  ExprPtr expression = expr->clone();

  // Within Molpro the two-electron integrals are stored in two separate
  // tensors: J and K The mapping to either of these two depends on the index
  // spaces of the associated indices It's K if one of the following applies
  // - all indices belong to the same space
  // - The two indices belonging to particle 1 and those of particle 2 belong to
  // different spaces
  // - There
  expression->visit(
      [](ExprPtr &expr) {
        // TODO: Get rid of hardcoded label
        if (expr.is<Tensor>() && expr.as<Tensor>().label() == L"g") {
          const Tensor &tensor = expr.as<Tensor>();
          assert(tensor.bra().size() == 2);
          assert(tensor.ket().size() == 2);

          // Copy indices as we might have to mutate them
          auto braIndices = tensor.bra();
          auto ketIndices = tensor.ket();

          constexpr IndexSpace::Type occ = IndexSpace::active_occupied;
          constexpr IndexSpace::Type virt = IndexSpace::active_unoccupied;

          // Assert occ < virt
          IndexTypeComparer cmp;
          static_assert(cmp(occ, virt));

          // Step 1: Use 8-fold permutational symmetry of spin-summed integrals
          // to bring indices into a canonical order in terms of the index
          // spaces they belong to. Note: This symmetry is generated by the two
          // individual bra-ket symmetries for indices for particle one and two
          // as well as the particle-1,2-symmetry (column-symmetry)
          //
          // The final goal is to order the indices in descending index space
          // size, where the assumed relative sizes are
          // occ < virt

          // Step 1a: Particle-intern bra-ket symmetry
          for (std::size_t i = 0; i < braIndices.size(); ++i) {
            if (cmp(braIndices.at(i).space().type(),
                    ketIndices.at(i).space().type())) {
              // This bra index belongs to a smaller space than the ket index ->
              // swap them
              std::swap(braIndices[i], ketIndices[i]);
            }
          }

          // Step 1b: Particle-1,2-symmetry
          if (braIndices[0].space().type() != braIndices[1].space().type()) {
            if (cmp(braIndices[0].space().type(),
                    braIndices[1].space().type())) {
              std::swap(braIndices[0], braIndices[1]);
              std::swap(ketIndices[0], ketIndices[1]);
            }
          } else if (cmp(ketIndices[0].space().type(),
                         ketIndices[1].space().type())) {
            std::swap(braIndices[0], braIndices[1]);
            std::swap(ketIndices[0], ketIndices[1]);
          }

          // Step 2: Look at the index space patterns to figure out whether
          // this is a K or a J integral. If the previously attempted sorting
          // of index spaces can be improved by switching the second and third
          // index, do that and thereby produce a J tensor. Otherwise, we retain
          // the index sequence as-is and thereby produce a K tensor.
          if (cmp(braIndices[1].space().type(), ketIndices[0].space().type())) {
            std::swap(braIndices[1], ketIndices[0]);

            expr =
                ex<Tensor>(L"J", std::move(braIndices), std::move(ketIndices));
          } else {
            expr =
                ex<Tensor>(L"K", std::move(braIndices), std::move(ketIndices));
          }
        }
      },
      true);

  return expression;
}

void ITFGenerator::addBlock(const itf::CodeBlock &block) {
  m_codes.reserve(m_codes.size() + block.results.size());

  std::vector<std::vector<Contraction>> contractionBlocks;

  for (const Result &currentResult : block.results) {
    ExprPtr expression = replaceTwoElectronIntgrals(currentResult.expression);

    contractionBlocks.push_back(
        to_contractions(expression, currentResult.resultTensor));

    if (currentResult.importResultTensor) {
      m_importedTensors.insert(currentResult.resultTensor);
    } else {
      m_createdTensors.insert(currentResult.resultTensor);
    }

    // If we encounter a tensor in an expression that we have not yet seen
    // before, it must be an imported tensor (otherwise the expression would be
    // invalid)
    expression->visit(
        [this](const ExprPtr &expr) {
          if (expr.is<Tensor>()) {
            const Tensor &tensor = expr.as<Tensor>();
            if (m_createdTensors.find(tensor) == m_createdTensors.end()) {
              m_importedTensors.insert(tensor);
            }
            m_encounteredIndices.insert(tensor.braket().begin(),
                                        tensor.braket().end());
          }
        },
        true);

    // Now go through all result tensors of the contractions that we have
    // produced and add all new tensors to the set of created tensors
    for (const Contraction &currentContraction : contractionBlocks.back()) {
      if (m_importedTensors.find(currentContraction.result) ==
          m_importedTensors.end()) {
        m_createdTensors.insert(currentContraction.result);
      }
    }
  }

  m_codes.push_back(CodeSection{block.name, std::move(contractionBlocks)});
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
    itf += L"tensor: " + to_itf(current) + L", !Create{type:disk}\n";
  }
  itf += L"\n\n";

  // Actual code
  for (const CodeSection &currentSection : m_codes) {
    itf += L"---- code(\"" + currentSection.name + L"\")\n";

    std::set<Tensor, TensorBlockCompare> allocatedTensors;

    for (const std::vector<Contraction> &currentBlock :
         currentSection.contractionBlocks) {
      for (const Contraction &currentContraction : currentBlock) {
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
        if (currentContraction.rhs.has_value()) {
          itf += L"load " + to_itf(currentContraction.rhs.value()) + L"\n";
        }

        itf += L"." + to_itf(currentContraction.result) + L" ";
        int sign = currentContraction.factor < 0 ? -1 : 1;

        itf += (sign < 0 ? L"-= " : L"+= ");
        if (currentContraction.factor * sign != 1) {
          itf += to_wstring(currentContraction.factor * sign) + L"*";
        }
        itf += to_itf(currentContraction.lhs) +
               (currentContraction.rhs.has_value()
                    ? L" " + to_itf(currentContraction.rhs.value())
                    : L"") +
               L"\n";

        if (currentContraction.rhs.has_value()) {
          itf += L"drop " + to_itf(currentContraction.rhs.value()) + L"\n";
        }
        itf += L"drop " + to_itf(currentContraction.lhs) + L"\n";
        itf += L"store " + to_itf(currentContraction.result) + L"\n";
      }

      itf += L"\n";
    }

    itf += L"\n---- end\n";
  }

  return itf;
}

}  // namespace detail

}  // namespace itf

}  // namespace sequant
