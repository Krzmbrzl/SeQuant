//
// Created by Robert Adam on 2023-09-04
//

#ifndef SEQUANT_CORE_EXPORT_ITF_HPP
#define SEQUANT_CORE_EXPORT_ITF_HPP

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/tensor.hpp>

#include <optional>
#include <set>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <vector>
#include <optional>

namespace sequant {

namespace itf {

struct Result {
  ExprPtr expression;
  Tensor resultTensor;
  bool importResultTensor;

  Result(ExprPtr expression, Tensor resultTensor, bool importResultTensor = true);
  Result(ExprPtr expression, bool importResultTensor = false);
};

struct CodeBlock {
  std::wstring name;
  std::vector<Result> results;

  CodeBlock(std::wstring blockName, Result result);
  CodeBlock(std::wstring blockName, std::vector<Result> results);
};

namespace detail {

/// Comparator that identifies Tensors only by their "block", which is defined
/// by its name, the amount of its indices as well as the space these indices
/// belong to. Note that it explicitly does not depend on the explicit index
/// labelling.
struct TensorBlockCompare {
  bool operator()(const Tensor &lhs, const Tensor &rhs) const {
    if (lhs.label() != rhs.label()) {
      return lhs.label() < rhs.label();
    }
    if (lhs.braket().size() != rhs.braket().size()) {
      return lhs.braket().size() < rhs.braket().size();
    }
    auto lhsBraket = lhs.braket();
    auto rhsBraket = rhs.braket();

    for (std::size_t i = 0; i < lhsBraket.size(); ++i) {
      if (lhsBraket.at(i).space() != rhsBraket.at(i).space()) {
        return lhsBraket.at(i).space() < rhsBraket.at(i).space();
      }
    }

    return false;
  }
};

struct Contraction {
  rational factor;

  Tensor result;
  Tensor lhs;
  std::optional<Tensor> rhs;
};

struct CodeSection {
  std::wstring name;
  std::vector<std::vector<Contraction>> contractionBlocks;
};

class ITFGenerator {
 public:
  ITFGenerator() = default;

  void addBlock(const itf::CodeBlock &block);

  std::wstring generate() const;

 private:
  std::set<Index> m_encounteredIndices;
  std::set<Tensor, TensorBlockCompare> m_importedTensors;
  std::set<Tensor, TensorBlockCompare> m_createdTensors;
  std::vector<CodeSection> m_codes;
};

}  // namespace detail

}  // namespace itf

std::wstring to_itf(const itf::CodeBlock &block);

template <typename Container,
          typename = std::enable_if_t<!std::is_same_v<
              std::remove_const_t<std::remove_reference_t<Container>>,
              itf::CodeBlock>>>
std::wstring to_itf(Container &&container) {
  static_assert(
      std::is_same_v<typename Container::value_type, itf::CodeBlock> ||
      std::is_same_v<std::remove_const_t<typename Container::value_type>,
                     ExprPtr>);
  itf::detail::ITFGenerator generator;

  if constexpr (std::is_same_v<typename Container::value_type,
                               itf::CodeBlock>) {
    for (const itf::CodeBlock &current : container) {
      generator.addBlock(current);
    }
  } else {
    static_assert(
        std::is_same_v<typename Container::value_type, itf::Result>,
        "Container::value_type must either be itf::CodeBlock or itf::Result");

    itf::CodeBlock block(L"Generate_Results",
                         std::forward<Container>(container));

    generator.addBlock(block);
  }

  return generator.generate();
}

}  // namespace sequant

#endif  // SEQUANT_CORE_EXPORT_ITF_HPP
