#ifndef SEQUANT_CALCULUS_DIFFERENTIATION_HPP
#define SEQUANT_CALCULUS_DIFFERENTIATION_HPP

#include <SeQuant/core/expr_fwd.hpp>

namespace sequant {

ExprPtr differentiate(const Expr &expr, const Variable &var);
ExprPtr differentiate(const Expr &expr, const Tensor &var);

}  // namespace sequant

#endif  // SEQUANT_CALCULUS_DIFFERENTIATION_HPP
