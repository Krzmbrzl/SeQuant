#ifndef SEQUANT_COMMON_SUBEXPRESSION_ELIMINATION_HPP
#define SEQUANT_COMMON_SUBEXPRESSION_ELIMINATION_HPP

#include <SeQuant/core/binary_node.hpp>
#include <SeQuant/core/container.hpp>
#include <SeQuant/core/eval_expr.hpp>
#include <SeQuant/core/eval_node.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/hash.hpp>

#include <cassert>
#include <concepts>
#include <ranges>
#include <string>
#include <unordered_map>

namespace sequant {

namespace {

template <typename TreeNode>
struct TreeNodeHasher {
  using is_transparent = void;

  std::size_t operator()(const TreeNode *node) const { return (*this)(*node); }

  std::size_t operator()(const TreeNode &node) const {
    return hash::value(*node);
  }
};

template <typename TreeNode>
struct TreeNodeEqualityComparator {
  using is_transparent = void;

  bool operator()(const TreeNode *lhs, const TreeNode *rhs) const {
    return (*this)(*lhs, *rhs);
  }

  bool operator()(const TreeNode &lhs, const TreeNode *rhs) const {
    return (*this)(lhs, *rhs);
  }

  bool operator()(const TreeNode *lhs, const TreeNode &rhs) const {
    return (*this)(*lhs, rhs);
  }

  bool operator()(const TreeNode &lhs, const TreeNode &rhs) const {
    // TODO: do something to guard against hash collisions
    return hash::value(*lhs) == hash::value(*rhs);
  }
};

template <typename TreeNode>
using SubexpressionHashCollector =
    std::unordered_map<const TreeNode *, std::size_t, TreeNodeHasher<TreeNode>,
                       TreeNodeEqualityComparator<TreeNode>>;

template <typename TreeNode>
using SubexpressionUsageCounts =
    std::unordered_map<TreeNode, std::size_t, TreeNodeHasher<TreeNode>,
                       TreeNodeEqualityComparator<TreeNode>>;

template <typename TreeNode>
using SubexpressionNames =
    std::unordered_map<TreeNode, std::wstring, TreeNodeHasher<TreeNode>,
                       TreeNodeEqualityComparator<TreeNode>>;

template <typename TreeNode>
class SubexpressionIdentifier {
 public:
  SubexpressionIdentifier() = default;

  bool operator()(const TreeNode &tree) {
    using std::ranges::end;

    if (auto it = intermediate_hashs.find(&tree);
        it != end(intermediate_hashs)) {
      // The expression identified by tree has been used before -> stop visiting
      // of subtree
      ++it->second;
      return false;
    }

    intermediate_hashs.emplace(&tree, 1);

    return true;
  }

  SubexpressionUsageCounts<TreeNode> take_subexpression_map() {
    SubexpressionUsageCounts<TreeNode> usages;
    for (const auto &[node_ptr, usage_count] : intermediate_hashs) {
      if (usage_count < 2) {
        // Everything that is used less than 2 times is not a common
        // subexpression
        continue;
      }

      usages.emplace(*node_ptr, usage_count);
    }

    intermediate_hashs.clear();

    return usages;
  }

 private:
  // Note: we are using the hash collector rather than the
  // SubexpressionUsageCounts because we will collect many, many entries which
  // will only appear once and the hash collector only stores node pointers, the
  // latter stores full node objects which are more than an order of magnitude
  // larger, which can easily lead to memory bottlenecks for large expressions.
  SubexpressionHashCollector<TreeNode> intermediate_hashs;
};

template <typename TreeNode, typename Transformer>
class SubexpressionReplacer {
 public:
  SubexpressionReplacer(const SubexpressionUsageCounts<TreeNode> &map,
                        const Transformer &transformer)
      : subexpressions(map), expr_to_tree(transformer) {}

  bool operator()(TreeNode &tree) {
    using std::ranges::begin;
    using std::ranges::end;

    auto it = subexpressions.find(tree);

    if (it == end(subexpressions)) {
      return true;
    }

    auto label_it = cse_names.find(tree);

    std::wstring label = label_it == end(cse_names)
                             ? L"CSE" + std::to_wstring(name_counter++)
                             : label_it->second;

    ExprPtr expr = [&]() {
      if (tree->is_tensor()) {
        return ex<Tensor>(label, bra(), ket(), aux(tree->canon_indices()));
      }

      return ex<Variable>(label);
    }();

	ExprPtr original = std::move(tree->expr());

    if (label_it == end(cse_names)) {
      cse_names.emplace(tree, label);
      // Store tree as the way to compute the current CSE (under the designated
      // name and with canonical indexing)
	  // TODO: make it work for EvalExpr as well
	  tree->set_expr(expr);

      cse_definitions.emplace_back(tree);
    }

    if (tree->canon_phase() != 1) {
      expr *= ex<Constant>(tree->canon_phase());
    }

    if (tree.root()) {
      assert(original->is<Tensor>() ||
             original->is<Variable>());

      ResultExpr result = [&]() {
        if (original->is<Tensor>()) {
          return ResultExpr(original->as<Tensor>(),
                            std::move(expr));
        }

        return ResultExpr(original->as<Variable>(),
                          std::move(expr));
      }();

      tree = expr_to_tree(result);
    } else {
      tree = expr_to_tree(expr);
    }

    return false;
  }

  container::svector<TreeNode> take_cse_definitions() {
    auto copy = std::move(cse_definitions);

    cse_definitions.clear();

    return copy;
  }

 private:
  const Transformer &expr_to_tree;
  const SubexpressionUsageCounts<TreeNode> &subexpressions;
  SubexpressionNames<TreeNode> cse_names;
  std::size_t name_counter = 1;
  container::svector<TreeNode> cse_definitions;
};

// TODO: provide predicate that computes the FLOPs required to compute the CSE
// and compares that with the expected speed for reading the CSE from disk. Use
// typical FLOPs/s and I/O rates from e.g. https://ssd.userbenchmark.com/
// https://openbenchmarking.org/test/pts/mt-dgemm
// to determine whether or not pre-computing the CSE is worth it.
// Note: of course the number of reusages of the CSE plays into this.
struct AcceptAllPredicate {
  template <typename... Ts>
  bool operator()(const Ts &...) const {
    return true;
  }
};

}  // namespace

template <std::ranges::range TreeRange,
          std::regular_invocable<ExprPtr> Transformer,
          std::predicate<std::ranges::range_value_t<TreeRange>, std::size_t>
              Filter = AcceptAllPredicate>
  requires std::regular_invocable<Transformer, ResultExpr> &&
           meta::eval_node<std::ranges::range_value_t<TreeRange>>
container::svector<std::ranges::range_value_t<TreeRange>>
eliminate_common_subexpressions(TreeRange &&range,
                                const Transformer &expr_to_tree,
                                const Filter &filter_predicate = {}) {
  using std::ranges::begin;
  using std::ranges::end;

  using TreeNode = std::ranges::range_value_t<TreeRange>;

  SubexpressionIdentifier<TreeNode> identifier;

  for (const TreeNode &current : range) {
    current.visit_internal(identifier);
  }

  SubexpressionUsageCounts<TreeNode> subexpression_usages =
      identifier.take_subexpression_map();

  // Apply filter
  std::erase_if(subexpression_usages, [&filter_predicate](const auto &pair) {
    return !filter_predicate(pair.first, pair.second);
  });

  SubexpressionReplacer<TreeNode, Transformer> replacer(subexpression_usages,
                                                        expr_to_tree);

  for (TreeNode &current : range) {
    current.visit_internal(replacer);
  }

  return replacer.take_cse_definitions();
}

}  // namespace sequant

#endif  // SEQUANT_COMMON_SUBEXPRESSION_ELIMINATION_HPP
