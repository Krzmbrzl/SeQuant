#ifndef SEQUANT_CORE_UTILITY_EXPR_HPP
#define SEQUANT_CORE_UTILITY_EXPR_HPP

#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>

#include <vector>

namespace sequant {

struct BraKet {
  std::vector<Index> bra;
  std::vector<Index> ket;
};

BraKet non_repeated_indices(const ExprPtr &expr);

}  // namespace sequant

#endif  // SEQUANT_CORE_UTILITY_EXPR_HPP
