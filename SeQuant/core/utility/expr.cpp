#include "SeQuant/core/utility/expr.hpp"

#include <SeQuant/core/tensor.hpp>

namespace sequant {

template <typename Container, typename Element>
void remove_one(Container& container, const Element& e) {
  auto iter = std::find(container.begin(), container.end(), e);

  if (iter != container.end()) {
    container.erase(iter);
  }
}

IndexGroups non_repeated_indices(const ExprPtr& expr) {
  if (expr.is<Constant>()) {
    return {};
  } else if (expr.is<Tensor>()) {
    const Tensor& tensor = expr.as<Tensor>();
    return {{tensor.bra().begin(), tensor.bra().end()},
            {tensor.ket().begin(), tensor.ket().end()}};
  } else if (expr.is<Sum>()) {
    // In order for the sum to be valid, all summands must have the same
    // external indices, so it suffices to look only at the first one
    return non_repeated_indices(expr.as<Sum>().summand(0));
  } else if (expr.is<Product>()) {
    std::set<Index> encounteredIndices;
    IndexGroups groups;

    for (const ExprPtr& current : expr.as<Product>()) {
      IndexGroups currentGroups = non_repeated_indices(current);

      for (Index& current : currentGroups.bra) {
        if (encounteredIndices.find(current) == encounteredIndices.end()) {
          encounteredIndices.insert(current);
          groups.bra.push_back(std::move(current));
        } else {
          remove_one(groups.bra, current);
          remove_one(groups.ket, current);
        }
      }
      // Same for ket indices
      for (Index& current : currentGroups.ket) {
        if (encounteredIndices.find(current) == encounteredIndices.end()) {
          encounteredIndices.insert(current);
          groups.ket.push_back(std::move(current));
        } else {
          remove_one(groups.bra, current);
          remove_one(groups.ket, current);
        }
      }
    }

    return groups;
  } else {
    throw std::runtime_error("Weird expression type encountered");
  }
}

}  // namespace sequant
