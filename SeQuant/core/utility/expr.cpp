#include "SeQuant/core/utility/expr.hpp"

#include <SeQuant/core/tensor.hpp>

namespace sequant {

BraKet non_repeated_indices(const ExprPtr& expr) {
	using Bra = decltype(BraKet::bra);
	using Ket = decltype(BraKet::ket);

  if (expr.is<Constant>()) {
    return {};
  } else if (expr.is<Tensor>()) {
    const Tensor& tensor = expr.as<Tensor>();
    return {Bra(tensor.bra().begin(), tensor.bra().end()),
            Ket(tensor.ket().begin(), tensor.ket().end())};
  } else if (expr.is<Sum>()) {
    // In order for the sum to be valid, all summands must have the same
    // external indices, so it suffices to look only at the first one
    return non_repeated_indices(expr.as<Sum>().summand(0));
  } else if (expr.is<Product>()) {
    std::set<Index> braIndices;
    std::set<Index> ketIndices;

    for (const ExprPtr& current : expr.as<Product>()) {
      BraKet indices = non_repeated_indices(current);

      braIndices.insert(indices.bra.begin(), indices.bra.end());
      ketIndices.insert(indices.ket.begin(), indices.ket.end());
    }

    BraKet externals;
    std::set_difference(braIndices.begin(), braIndices.end(),
                        ketIndices.begin(), ketIndices.end(),
                        std::back_inserter(externals.bra));
    std::set_difference(ketIndices.begin(), ketIndices.end(),
                        braIndices.begin(), braIndices.end(),
                        std::back_inserter(externals.ket));

    return externals;
  } else {
    throw std::runtime_error("Weird expression type encountered");
  }
}

}  // namespace sequant
