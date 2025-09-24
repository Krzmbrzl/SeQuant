#include <SeQuant/core/calculus/differentiation.hpp>
#include <SeQuant/core/expr.hpp>

#include <algorithm>
#include <cassert>

namespace sequant {

template <typename T, typename EqualityCmp>
ExprPtr differentiate(const Expr &expr, const T &var, EqualityCmp cmp = {}) {
  if (expr.is_atom()) {
    if (cmp(expr, var)) {
      return ex<Constant>(1);
    }

    return ex<Constant>(0);
  }

  ExprPtr result = ex<Constant>(0);

  if (expr.is<Sum>()) {
    for (const ExprPtr &current : expr.as<Sum>().summands()) {
      result += differentiate(*current, var, cmp);
    }

    return result;
  }

  if (!expr.is<Product>()) {
    assert(false);
    throw std::runtime_error(
        "expr is neither atom, sum or product -> unhandled case");
  }

  const Product &product = expr.as<Product>();

  auto pos = product.begin();
  auto it = std::find_if(pos, product.end(),
                         [&](const ExprPtr &e) { return contains(e, var); });
  while (it != product.end()) {
    container::svector<ExprPtr> factors;
    auto idx = std::distance(product.begin(), it);
    for (std::size_t i = 0; i < product.factors().size(); ++i) {
      if (i == idx) {
        continue;
      }

      factors.emplace_back(product.factor(i)->clone());
    }

    it = std::find_if(pos, expr.end(),
                      [&](const ExprPtr &e) { return contains(e, var); });
  }

  return result;
}

}  // namespace sequant
