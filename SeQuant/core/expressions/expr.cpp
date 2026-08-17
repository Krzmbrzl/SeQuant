//
// Created by Eduard Valeyev on 2019-02-06.
//

#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/product.hpp>

namespace sequant {

std::wstring Expr::to_latex() const {
  throw Exception("to_latex not implemented for " + type_name());
}

bool proportional_to::operator()(const ExprPtr &expr1,
                                 const ExprPtr &expr2) const {
  if (expr1->type_id() !=
      expr2->type_id()) {  // if expr1 is a Product with single factor == expr2,
                           // or vice versa
    if (expr1.is<Product>()) {
      return expr1.as<Product>().factors().size() == 1 &&
             expr1.as<Product>().factors().front() == expr2;
    } else if (expr2.is<Product>()) {
      return expr2.as<Product>().factors().size() == 1 &&
             expr2.as<Product>().factors().front() == expr1;
    } else
      return false;
  }

  // expr1 and expr2 are same type

  if (expr1.is<Constant>()) {
    return true;
  }
  if (expr1.is<Product>()) {
    return expr1->hash_value() == expr2->hash_value() &&
           expr1.as<Product>().factors() == expr2.as<Product>().factors();
  }
  return expr1 == expr2;
}

}  // namespace sequant
