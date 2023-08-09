//
// Created by Bimal Gaudel on 5/24/21.
//

#ifndef SEQUANT_EVAL_NODE_HPP
#define SEQUANT_EVAL_NODE_HPP

#include "asy_cost.hpp"
#include "binary_node.hpp"
#include "eval_expr.hpp"

#include <SeQuant/core/math.hpp>

namespace sequant {

namespace {

template <typename, typename = void>
constexpr bool IsEvalExpr{};

template <typename T>
constexpr bool
    IsEvalExpr<T, std::enable_if_t<std::is_convertible_v<T, EvalExpr>>>{true};

template <typename, typename = void>
constexpr bool IsEvalNode{};

template <typename T>
constexpr bool IsEvalNode<FullBinaryNode<T>, std::enable_if_t<IsEvalExpr<T>>>{
    true};

template <typename T>
constexpr bool
    IsEvalNode<const FullBinaryNode<T>, std::enable_if_t<IsEvalExpr<T>>>{true};

}  // namespace

template <typename T,
          typename = std::enable_if_t<std::is_convertible_v<T, EvalExpr>>>
using EvalNode = FullBinaryNode<T>;

///
/// \brief Creates an evaluation tree from @c ExprPtr.
///
/// \tparam ExprT Can be @c EvalExpr or derived from it.
///
/// \param expr The expression to be converted to an evaluation tree.
///
/// \return A full-binary tree whose nodes are @c EvalExpr (or derived from it)
///         types.
///
template <typename ExprT,
          std::enable_if_t<IsEvalExpr<std::decay_t<ExprT>>, bool> = true>
EvalNode<ExprT> eval_node(ExprPtr const& expr) {
  using ranges::accumulate;
  using ranges::views::tail;
  using ranges::views::transform;

  if (expr->is<Tensor>()) return EvalNode<ExprT>{ExprT{expr->as<Tensor>()}};
  if (expr->is<Constant>()) return EvalNode<ExprT>{ExprT{expr->as<Constant>()}};
  if (expr->is<Variable>()) return EvalNode<ExprT>{ExprT{expr->as<Variable>()}};
  assert(expr->is<Sum>() || expr->is<Product>());

  auto subxprs = *expr | ranges::views::transform([](auto const& x) {
    return eval_node<ExprT>(x);
  }) | ranges::to_vector;

  if (expr->is<Product>()) {
    auto const& prod = expr->as<Product>();
    if (prod.scalar() != 1)
      subxprs.emplace_back(eval_node<ExprT>(ex<Constant>(prod.scalar())));
  }

  auto const op = expr->is<Sum>() ? EvalOp::Sum : EvalOp::Prod;

  auto bnode = ranges::accumulate(
      ranges::views::tail(subxprs), std::move(*subxprs.begin()),
      [op](auto& lnode, auto& rnode) {
        auto pxpr = ExprT{*lnode, *rnode, op};
        return EvalNode<ExprT>(std::move(pxpr), std::move(lnode),
                               std::move(rnode));
      });

  return bnode;
}

///
/// \brief Creates an evaluation tree from Expr.
///
/// \tparam EvalNodeT is @c EvalNode<ExprT>, where @c ExprT is @c EvalExpr or
/// derived.
///
/// \param expr The expression to be converted to an evaluation tree.
///
/// \return A full-binary tree whose type is @c EvalNodeT.
///
template <typename EvalNodeT,
          std::enable_if_t<IsEvalNode<std::decay_t<EvalNodeT>>, bool> = true>
EvalNodeT eval_node(ExprPtr const& expr) {
  return eval_node<typename EvalNodeT::value_type>(expr);
}

template <typename ExprT>
ExprPtr to_expr(EvalNode<ExprT> const& node) {
  auto const op = node->op_type();
  auto const& evxpr = *node;

  if (node.leaf()) return evxpr.expr();

  if (op == EvalOp::Prod) {
    auto prod = Product{};

    ExprPtr lexpr = to_expr(node.left());
    ExprPtr rexpr = to_expr(node.right());

    prod.append(1, lexpr, Product::Flatten::No);
    prod.append(1, rexpr, Product::Flatten::No);

    assert(!prod.empty());

    if (prod.size() == 1 && !prod.factor(0)->is<Tensor>()) {
      return ex<Product>(Product{prod.scalar(), prod.factor(0)->begin(),
                                 prod.factor(0)->end(), Product::Flatten::No});
    } else {
      return ex<Product>(std::move(prod));
    }

  } else {
    assert(op == EvalOp::Sum && "unsupported operation type");
    return ex<Sum>(Sum{to_expr(node.left()), to_expr(node.right())});
  }
}

template <typename ExprT>
ExprPtr linearize_eval_node(EvalNode<ExprT> const& node) {
  if (node.leaf()) return to_expr(node);

  ExprPtr lres = linearize_eval_node(node.left());
  ExprPtr rres = linearize_eval_node(node.right());

  assert(lres);
  assert(rres);

  if (node->op_type() == EvalOp::Sum) return ex<Sum>(ExprPtrList{lres, rres});
  assert(node->op_type() == EvalOp::Prod);
  return ex<Product>(
      Product{1, ExprPtrList{lres, rres}, Product::Flatten::Yes});
}

namespace {

enum NodePos { Left = 0, Right, This };

class ContractedIndexCount {
 public:
  template <typename NodeT, typename = std::enable_if_t<IsEvalNode<NodeT>>>
  explicit ContractedIndexCount(NodeT const& n) {
    auto const L = NodePos::Left;
    auto const R = NodePos::Right;
    auto const T = NodePos::This;

    assert(n->is_tensor() && n.left()->is_tensor() && n.right()->is_tensor());

    for (auto p : {L, R, T}) {
      auto const& t = (p == L ? n.left() : p == R ? n.right() : n)->as_tensor();
      std::tie(occs_[p], virts_[p]) = occ_virt(t);
      ranks_[p] = occs_[p] + virts_[p];
    }

    // no. of contractions in occupied index space (always a whole number)
    occ_ = (occs_[L] + occs_[R] - occs_[T]) / 2;

    // no. of contractions in virtual index space (always a whole number)
    virt_ = (virts_[L] + virts_[R] - virts_[T]) / 2;

    is_outerprod_ = ranks_[L] + ranks_[R] == ranks_[T];
  }

  [[nodiscard]] size_t occ(NodePos p) const noexcept { return occs_[p]; }

  [[nodiscard]] size_t virt(NodePos p) const noexcept { return virts_[p]; }

  [[nodiscard]] size_t rank(NodePos p) const noexcept { return ranks_[p]; }

  [[nodiscard]] size_t occ() const noexcept { return occ_; }

  [[nodiscard]] size_t virt() const noexcept { return virt_; }

  [[nodiscard]] bool is_outerpod() const noexcept { return is_outerprod_; }

  [[nodiscard]] size_t unique_occs() const noexcept {
    return occ(NodePos::Left) + occ(NodePos::Right) - occ();
  }

  [[nodiscard]] size_t unique_virts() const noexcept {
    return virt(NodePos::Left) + virt(NodePos::Right) - virt();
  }

  static std::pair<size_t, size_t> occ_virt(sequant::Tensor const& t) {
    auto bk_rank = t.bra_rank() + t.ket_rank();
    auto nocc =
        ranges::count_if(t.const_braket(), [](sequant::Index const& idx) {
          return idx.space() == sequant::IndexSpace::active_occupied;
        });
    auto nvirt = bk_rank - nocc;
    return {nocc, nvirt};
  }

 private:
  std::array<size_t, 3> occs_{0, 0, 0};
  std::array<size_t, 3> virts_{0, 0, 0};
  std::array<size_t, 3> ranks_{0, 0, 0};
  size_t occ_ = 0;
  size_t virt_ = 0;
  bool is_outerprod_ = false;
};
}  // namespace

struct Flops {
  template <typename NodeT, typename = std::enable_if_t<IsEvalNode<NodeT>>>
  [[nodiscard]] AsyCost operator()(NodeT const& n) const {
    if (n.leaf()) return AsyCost::zero();
    if (n->op_type() == EvalOp::Prod  //
        && n.left()->is_tensor()      //
        && n.right()->is_tensor()) {
      auto const idx_count = ContractedIndexCount{n};
      auto c = AsyCost{idx_count.unique_occs(), idx_count.unique_virts()};
      return idx_count.is_outerpod() ? c : 2 * c;
    } else if (n->is_tensor()) {
      // scalar times a tensor
      // or a tensor plus a tensor
      auto [o, v] = ContractedIndexCount::occ_virt(n->as_tensor());
      return AsyCost{o, v};
    } else /* scalar (+|*) scalar */
      return AsyCost::zero();
  }
};

struct FlopsWithSymm {
  template <typename NodeT, typename = std::enable_if_t<IsEvalNode<NodeT>>>
  [[nodiscard]] AsyCost operator()(NodeT const& n) const {
    auto cost = Flops{}(n);
    if (n.leaf() || !(n->is_tensor()            //
                      && n.left()->is_tensor()  //
                      && n.right()->is_tensor()))
      return cost;

    // confirmed:
    // left, right and this node
    // all have tensor expression
    auto const& t = n->as_tensor();
    auto const tsymm = t.symmetry();
    //

    // ------
    // the rules of cost reduction are taken from
    //   doi:10.1016/j.procs.2012.04.044
    // ------
    if (tsymm == Symmetry::symm || tsymm == Symmetry::antisymm) {
      auto const op = n->op_type();
      auto const tbrank = t.bra_rank();
      auto const tkrank = t.ket_rank();
      if (op == EvalOp::Sum)
        cost = tsymm == Symmetry::symm
                   ? cost / (factorial(tbrank) * factorial(tkrank))
                   : cost / factorial(tbrank);
      else if (op == EvalOp::Prod) {
        auto const lsymm = n.left()->as_tensor().symmetry();
        auto const rsymm = n.right()->as_tensor().symmetry();
        cost = (lsymm == rsymm && lsymm == Symmetry::nonsymm)
                   ? cost / factorial(tbrank)
                   : cost / (factorial(tbrank) * factorial(tkrank));
      } else
        assert(false &&
               "Unsupported evaluation operation for asymptotic cost "
               "computation.");
    }
    return cost;
  }
};

template <typename NodeT, typename F = Flops,
          typename = std::enable_if_t<IsEvalNode<NodeT>>,
          typename = std::enable_if_t<std::is_invocable_r_v<AsyCost, F, NodeT>>>
AsyCost asy_cost(NodeT const& node, F const& cost_fn = {}) {
  return node.leaf() ? cost_fn(node)
                     : asy_cost(node.left(), cost_fn) +
                           asy_cost(node.right(), cost_fn) + cost_fn(node);
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_NODE_HPP
