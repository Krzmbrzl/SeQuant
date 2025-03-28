#ifndef SEQUANT_EVAL_EVAL_HPP
#define SEQUANT_EVAL_EVAL_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/logger.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/parse.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/eval/cache_manager.hpp>
#include <SeQuant/domain/eval/eval_result.hpp>

#include <btas/btas.h>
#include <tiledarray.h>

#include <chrono>
#include <range/v3/numeric.hpp>
#include <range/v3/view.hpp>

#include <any>
#include <iostream>
#include <stdexcept>
#include <type_traits>

namespace sequant {

namespace {

///
/// Invokes @c fun on @c args and returns a pair of the result, and
/// the time duration as @c std::chrono::duration<double>.
///
template <typename F, typename... Args>
auto timed_eval(F&& fun, Args&&... args) {
  using Clock = std::chrono::high_resolution_clock;
  auto tstart = Clock::now();
  auto&& res = std::forward<F>(fun)(std::forward<Args>(args)...);
  auto tend = Clock::now();
  return std::pair{std::move(res),
                   std::chrono::duration<double>(tend - tstart)};
}

///
/// Invokes @c fun that returns void on the arguments @c args and returns the
/// time duration as @c std::chrono::duration<double>.
template <
    typename F, typename... Args,
    std::enable_if_t<std::is_invocable_r_v<void, F, Args...>, bool> = true>
auto timed_eval_inplace(F&& fun, Args&&... args) {
  using Clock = std::chrono::high_resolution_clock;
  auto tstart = Clock::now();
  std::forward<F>(fun)(std::forward<Args>(args)...);
  auto tend = Clock::now();
  return std::chrono::duration<double>(tend - tstart);
}

template <typename... Args>
void log_eval(Args const&... args) noexcept {
  auto& l = Logger::instance();
  if (l.eval.level > 0) write_log(l, "[EVAL] ", args...);
}

[[maybe_unused]] void log_cache_access(size_t key, CacheManager const& cm) {
  auto& l = Logger::instance();
  if (l.eval.level > 0) {
    assert(cm.exists(key));
    auto max_l = cm.max_life(key);
    auto cur_l = cm.life(key);
    write_log(l,                                    //
              "[CACHE] Accessed key: ", key, ". ",  //
              cur_l, "/", max_l, " lives remain.\n");
    if (cur_l == 0) {
      write_log(l,  //
                "[CACHE] Released key: ", key, ".\n");
    }
  }
}

[[maybe_unused]] void log_cache_store(size_t key, CacheManager const& cm) {
  auto& l = Logger::instance();
  if (l.eval.level > 0) {
    assert(cm.exists(key));
    write_log(l,  //
              "[CACHE] Stored key: ", key, ".\n");
    // because storing automatically implies immediately accessing it
    log_cache_access(key, cm);
  }
}

[[maybe_unused]] std::string perm_groups_string(
    container::svector<std::array<size_t, 3>> const& perm_groups) {
  std::string result;
  for (auto const& g : perm_groups)
    result += "(" + std::to_string(g[0]) + "," + std::to_string(g[1]) + "," +
              std::to_string(g[2]) + ") ";
  result.pop_back();  // remove last space
  return result;
}

template <typename, typename = void>
constexpr bool HasAnnotMethod{};

template <typename T>
constexpr bool HasAnnotMethod<
    T, std::void_t<decltype(std::declval<meta::remove_cvref_t<T>>().annot())>> =
    true;

template <typename, typename = void>
constexpr bool IsEvaluable{};

template <typename T>
constexpr bool IsEvaluable<
    FullBinaryNode<T>,
    std::enable_if_t<std::is_convertible_v<T, EvalExpr> && HasAnnotMethod<T>>> =
    true;

template <typename T>
constexpr bool IsEvaluable<
    const FullBinaryNode<T>,
    std::enable_if_t<std::is_convertible_v<T, EvalExpr> && HasAnnotMethod<T>>> =
    true;

template <typename, typename = void>
constexpr bool IsIterableOfEvaluableNodes{};

template <typename Iterable>
constexpr bool IsIterableOfEvaluableNodes<
    Iterable, std::enable_if_t<IsEvaluable<meta::range_value_t<Iterable>>>> =
    true;

}  // namespace

template <typename, typename, typename = void>
constexpr bool IsLeafEvaluator{};

template <typename NodeT>
constexpr bool IsLeafEvaluator<NodeT, CacheManager, void>{};

template <typename NodeT, typename Le>
constexpr bool IsLeafEvaluator<
    NodeT, Le,
    std::enable_if_t<
        IsEvaluable<NodeT> &&
        std::is_same_v<
            ERPtr, std::remove_reference_t<std::invoke_result_t<Le, NodeT>>>>> =
    true;

///
/// \brief This class extends the EvalExpr class by adding a annot() method so
///        that it can be used to evaluate using TiledArray.
///
class EvalExprTA final : public EvalExpr {
 public:
  template <typename... Args, typename = std::enable_if_t<
                                  std::is_constructible_v<EvalExpr, Args...>>>
  EvalExprTA(Args&&... args) : EvalExpr{std::forward<Args>(args)...} {
    annot_ = indices_annot();
  }

  [[nodiscard]] inline auto const& annot() const noexcept { return annot_; }

 private:
  std::string annot_;
};

///
/// \brief This class extends the EvalExpr class by adding a annot() method so
///        that it can be used to evaluate using BTAS.
///
class EvalExprBTAS final : public EvalExpr {
 public:
  using annot_t = container::svector<long>;

  template <typename... Args, typename = std::enable_if_t<
                                  std::is_constructible_v<EvalExpr, Args...>>>
  EvalExprBTAS(Args&&... args) : EvalExpr{std::forward<Args>(args)...} {
    annot_ = index_hash(canon_indices()) | ranges::to<annot_t>;
  }

  ///
  /// \return Annotation (container::svector<long>) for BTAS::Tensor.
  ///
  [[nodiscard]] inline annot_t const& annot() const noexcept { return annot_; }

 private:
  annot_t annot_;
};  // EvalExprBTAS

template <typename NodeT, typename Le,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool> = true>
ERPtr evaluate_crust(NodeT const&, Le const&);

template <typename NodeT, typename Le,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool> = true>
ERPtr evaluate_crust(NodeT const&, Le const&, CacheManager&);

template <typename NodeT, typename Le, typename... Args,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool> = true>
ERPtr evaluate_core(NodeT const& node, Le const& le, Args&&... args) {
  if (node.leaf()) {
    auto&& [res, time] = timed_eval(le, node);
    log_eval(node->is_constant()   ? "[CONSTANT] "
             : node->is_variable() ? "[VARIABLE] "
                                   : "[TENSOR] ",
             node->label(), "  ", time.count(), "\n");
    return res;
  } else {
    ERPtr const left =
        evaluate_crust(node.left(), le, std::forward<Args>(args)...);
    ERPtr const right =
        evaluate_crust(node.right(), le, std::forward<Args>(args)...);

    assert(left);
    assert(right);

    std::array<std::any, 3> const ann{node.left()->annot(),
                                      node.right()->annot(), node->annot()};

    if (node->op_type() == EvalOp::Sum) {
      auto&& [res, time] = timed_eval([&]() { return left->sum(*right, ann); });
      log_eval("[SUM] ", node.left()->label(), " + ", node.right()->label(),
               " = ", node->label(), "  ", time.count(), "\n");
      return res;
    } else {
      assert(node->op_type() == EvalOp::Prod);
      auto const de_nest =
          node.left()->tot() && node.right()->tot() && !node->tot();

      auto&& [res, time] = timed_eval([&]() {
        return left->prod(*right, ann,
                          de_nest ? TA::DeNest::True : TA::DeNest::False);
      });

      log_eval("[PRODUCT] ", node.left()->label(), " * ", node.right()->label(),
               " = ", node->label(), "  ", time.count(), "\n");
      return res;
    }
  }
}

template <typename NodeT, typename Le,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool>>
ERPtr evaluate_crust(NodeT const& node, Le const& le) {
  return evaluate_core(node, le);
}

template <typename NodeT>
inline ERPtr mult_by_phase(NodeT const& node, ERPtr res) {
  return node->canon_phase() == 1 ? res
                                  : res->mult_by_phase(node->canon_phase());
}

template <typename NodeT, typename Le,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool>>
ERPtr evaluate_crust(NodeT const& node, Le const& le, CacheManager& cache) {
  auto const h = hash::value(*node);
  if (auto ptr = cache.access(h); ptr) {
    log_cache_access(h, cache);
    return mult_by_phase(node, ptr);
  } else if (cache.exists(h)) {
    auto ptr =
        cache.store(h, mult_by_phase(node, evaluate_core(node, le, cache)));
    log_cache_store(h, cache);
    return mult_by_phase(node, ptr);
  } else {
    return evaluate_core(node, le, cache);
  }
}

///
/// \param node An EvalNode to be evaluated.
///
/// \param le A leaf evaluator that takes an EvalNode and returns a tensor
///           (TA::TArrayD, btas::Tensor<double>, etc.) or a constant (double,
///           complex<double>, etc.).
///
/// \param args Optional CacheManager object passed by reference.
///
/// \return ERPtr to the resulting EvalResult.
///
/// \see EvalResult to know more about the return type.
///
template <typename NodeT, typename Le, typename... Args,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool> = true>
auto evaluate(NodeT const& node, Le&& le, Args&&... args) {
  return evaluate_crust(node, le, std::forward<Args>(args)...);
}

///
/// \param nodes An iterable of EvalNode objects that will be evaluated turn by
///              turn and summed up.
///
/// \param le A leaf evaluator that takes an EvalNode and returns a tensor
///           (TA::TArrayD, btas::Tensor<double>, etc.) or a constant (double,
///           complex<double>, etc.).
///
/// \param args Optional CacheManager object passed by reference.
///
/// \return ERPtr to the resulting EvalResult.
///
/// \see EvalResult to know more about the return type.
///
template <typename NodesT, typename Le, typename... Args,
          std::enable_if_t<IsIterableOfEvaluableNodes<NodesT>, bool> = true,
          std::enable_if_t<IsLeafEvaluator<meta::range_value_t<NodesT>, Le>,
                           bool> = true>
auto evaluate(NodesT const& nodes, Le const& le, Args&&... args) {
  auto iter = std::begin(nodes);
  auto end = std::end(nodes);
  assert(iter != end);

  log_eval("[TERM] ", " ", to_string(deparse(to_expr(*iter))), "\n");

  auto result = evaluate(*iter, le, std::forward<Args>(args)...);
  auto const pnode_label = (*iter)->label();

  for (++iter; iter != end; ++iter) {
    log_eval("[TERM] ", " ", to_string(deparse(to_expr(*iter))), "\n");
    auto right = evaluate(*iter, le, std::forward<Args>(args)...);
    auto&& time = timed_eval_inplace([&]() { result->add_inplace(*right); });
    log_eval("[ADD_INPLACE] ",  //
             pnode_label,       //
             " += ",            //
             (*iter)->label(),  //
             "  ",              //
             time.count(),      //
             "\n");
  }

  return result;
}

///
/// \param node An EvalNode to be evaluated into a tensor.
/// \param layout The layout of the resulting tensor. It is a permutation of the
///               result of node->annot().
/// \param le A leaf evaluator that takes an EvalNode and returns a tensor
///           (TA::TArrayD, btas::Tensor<double>, etc.) or a constant (double,
///           complex<double>, etc.).
///
/// \param args Optional CacheManager object passed by reference.
///
/// \return ERPtr to the resulting tensor.
///
/// \see EvalResult to know more about the return type.
///
template <typename NodeT, typename Annot, typename Le, typename... Args,
          std::enable_if_t<IsEvaluable<NodeT>, bool> = true,
          std::enable_if_t<IsLeafEvaluator<NodeT, Le>, bool> = true>
auto evaluate(NodeT const& node,    //
              Annot const& layout,  //
              Le const& le, Args&&... args) {
  log_eval("[TERM] ", " ", to_string(deparse(to_expr(node))), "\n");
  auto result = evaluate_crust(node, le, std::forward<Args>(args)...);

  auto&& [res, time] = timed_eval([&]() {
    return result->permute(std::array<std::any, 2>{node->annot(), layout});
  });
  log_eval("[PERMUTE] ", node->label(), "  ", time.count(), "\n");
  return res;
}

///
/// \param nodes An iterable of EvalNode objects that will be evaluated turn by
///              turn and summed up into a tensor.
///
/// \param layout The layout of the resulting tensor. It is a permutation of the
///               result of node->annot().
///
/// \param le A leaf evaluator that takes an EvalNode and returns a tensor
///           (TA::TArrayD, btas::Tensor<double>, etc.) or a constant (double,
///           complex<double>, etc.).
///
/// \param args Optional CacheManager object passed by reference.
///
/// \return ERPtr to the resulting tensor.
///
/// \see EvalResult to know more about the return type.
///
template <typename NodesT, typename Annot, typename Le, typename... Args,
          std::enable_if_t<IsIterableOfEvaluableNodes<NodesT>, bool> = true,
          std::enable_if_t<IsLeafEvaluator<meta::range_value_t<NodesT>, Le>,
                           bool> = true>
auto evaluate(NodesT const& nodes,  //
              Annot const& layout,  //
              Le const& le, Args&&... args) {
  auto iter = std::begin(nodes);
  auto end = std::end(nodes);
  assert(iter != end);
  auto const pnode_label = (*iter)->label();

  log_eval("[TERM] ", " ", to_string(deparse(to_expr(*iter))), "\n");

  auto result = evaluate(*iter, layout, le, std::forward<Args>(args)...);
  for (++iter; iter != end; ++iter) {
    log_eval("[TERM] ", " ", to_string(deparse(to_expr(*iter))), "\n");
    auto right = evaluate(*iter, layout, le, std::forward<Args>(args)...);
    auto&& time = timed_eval_inplace([&]() { result->add_inplace(*right); });
    log_eval("[ADD_INPLACE] ",  //
             pnode_label,       //
             " += ",            //
             (*iter)->label(),  //
             "  ",              //
             time.count(),      //
             "\n");
  }
  return result;
}

///
/// \param node An EvalNode or an iterable of such nodes to be evaluated into a
///             tensor.
///
/// \param layout The layout of the resulting tensor. It is a permutation of the
///               result of node->annot().
///
/// \param le A leaf evaluator that takes an EvalNode and returns a tensor
///           (TA::TArrayD, btas::Tensor<double>, etc.) or a constant (double,
///           complex<double>, etc.).
///
/// \param args Optional CacheManager object passed by reference.
///
/// \return ERPtr to the resulting tensor.
///
/// \see EvalResult to know more about the return type.
///
template <typename NodeT, typename Annot, typename Le, typename... Args>
auto evaluate_symm(NodeT const& node, Annot const& layout, Le const& le,
                   Args&&... args) {
  auto result = evaluate(node, layout, le, std::forward<Args>(args)...);

  auto&& [res, time] = timed_eval([&]() { return result->symmetrize(); });
  log_eval("[SYMMETRIZE] (layout) ",  //
           "(", layout, ") ",         //
           time.count(),              //
           "\n");
  return res;
}

///
/// \param node An EvalNode or an iterable of such nodes to be evaluated into a
///             tensor.
///
/// \param layout The layout of the resulting tensor. It is a permutation of the
///               result of node->annot().
///
/// \param le A leaf evaluator that takes an EvalNode and returns a tensor
///           (TA::TArrayD, btas::Tensor<double>, etc.) or a constant (double,
///           complex<double>, etc.).
///
/// \param args Optional CacheManager object passed by reference.
///
/// \return ERPtr to the resulting tensor.
///
/// \see EvalResult to know more about the return type.
///
template <typename NodeT, typename Annot, typename Le,
          typename... Args>
auto evaluate_antisymm(NodeT const& node,    //
                       Annot const& layout,  //
                       Le const& le,         //
                       Args&&... args) {
  size_t bra_rank;
  {
    ExprPtr expr_ptr{};
    if constexpr (IsIterableOfEvaluableNodes<NodeT>) {
      expr_ptr = (*std::begin(node))->expr();
    } else {
      expr_ptr = node->expr();
    }
    assert(expr_ptr->is<Tensor>());
    auto const& t = expr_ptr->as<Tensor>();
    bra_rank = t.bra_rank();
  }

  auto result = evaluate(node, layout, le, std::forward<Args>(args)...);
  auto&& [res, time] =
      timed_eval([&]() { return result->antisymmetrize(bra_rank); });
  log_eval("[ANTISYMMETRIZE] (bra rank, layout) ",  //
           "(", bra_rank, ", ", layout, ") ",       //
           time.count(),                            //
           "\n");
  return res;
}

}  // namespace sequant

#endif  // SEQUANT_EVAL_EVAL_HPP
