#ifndef SEQUANT_EXPRESSIONS_EXPR_CONTAINER_HPP
#define SEQUANT_EXPRESSIONS_EXPR_CONTAINER_HPP

#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/traits.hpp>
#include <SeQuant/core/meta.hpp>

#include <concepts>
#include <memory>

namespace sequant {

class ExprContainer {
 public:
  template <typename E>
    requires(is_an_expr_v<E> && std::copyable<E>)
  explicit ExprContainer(E &&expr)
      : expr_(std::make_unique<std::remove_cvref<E>>(std::forward<E>(expr))) {}

  template <typename E>
    requires(is_an_expr_v<E> && std::copyable<E>)
  ExprContainer &operator=(E &&expr) {
    expr_ = std::make_unique<std::remove_cvref_t<E>>(std::forward<E>(expr));

    return *this;
  }

  explicit ExprContainer(const ExprContainer &expr);
  ExprContainer(ExprContainer &&expr) = default;

  ExprContainer &operator=(const ExprContainer &expr);
  ExprContainer &operator=(ExprContainer &&expr) = default;

  ~ExprContainer() = default;

  ExprContainer copy() const;

  ExprContainer operator+(const Expr &expr) const;
  ExprContainer operator-(const Expr &expr) const;
  ExprContainer operator*(const Expr &expr) const;

  ExprContainer &operator+=(const Expr &expr);
  ExprContainer &operator-=(const Expr &expr);
  ExprContainer &operator*=(const Expr &expr);

  operator const Expr &() const;
  operator Expr &();

 private:
  std::unique_ptr<Expr> expr_;
};

}  // namespace sequant

namespace std {

void swap(sequant::ExprContainer &lhs, sequant::ExprContainer &rhs);

}

#endif  // SEQUANT_EXPRESSIONS_EXPR_CONTAINER_HPP
