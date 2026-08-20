#ifndef SEQUANT_EXPRESSIONS_SUM_HPP
#define SEQUANT_EXPRESSIONS_SUM_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expressions/constant.hpp>
#include <SeQuant/core/expressions/expr.hpp>
#include <SeQuant/core/expressions/expr_container.hpp>
#include <SeQuant/core/expressions/expr_iterator.hpp>
#include <SeQuant/core/expressions/expr_ptr.hpp>
#include <SeQuant/core/expressions/product.hpp>
#include <SeQuant/core/meta.hpp>
#include <SeQuant/core/runtime.hpp>
#include <SeQuant/core/utility/aggregate.hpp>
#include <SeQuant/core/utility/macros.hpp>

#include <range/v3/range/access.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>

#include <concepts>
#include <initializer_list>
#include <memory>
#include <optional>
#include <ranges>
#include <type_traits>

namespace sequant {

/// @brief sum of zero or more summands

/// Sum is associative and is flattened automatically.
class Sum : public Expr {
 public:
  using summands_type = container::svector<ExprContainer, 2>;

  Sum() = default;
  virtual ~Sum() = default;
  Sum(const Sum &) = default;
  Sum(Sum &&) = default;
  Sum &operator=(const Sum &) = default;
  Sum &operator=(Sum &&) = default;

  /// construct a Sum out of zero or more summands
  /// @param summands an initializer list of summands
  template <expr_holder E>
  Sum(std::initializer_list<E> summands) {
    // use append to flatten out Sum summands
    for (auto &&summand : summands) {
      if constexpr (std::same_as<ExprPtr, std::remove_cvref_t<E>>) {
        append(std::forward<decltype(summand)>(summand));
      } else {
        append(ExprContainer(std::forward<decltype(summand)>(summand)));
      }
    }
  }

  /// construct a Sum out of a range of summands
  /// @param begin the begin iterator
  /// @param end the end iterator
  template <typename Iterator>
  Sum(Iterator begin, Iterator end) : Sum(std::ranges::subrange(begin, end)) {}

  /// construct a Sum out of a range of summands
  /// @param rng a range
  template <std::ranges::range Range>
    requires(!std::same_as<std::remove_cvref_t<Range>, ExprPtrList> &&
             !std::same_as<std::remove_cvref_t<Range>, ExprContainerList>)
  explicit Sum(Range &&rng) {
    // N.B. use append to flatten out Sum summands
    constexpr auto is_expr = is_an_expr_v<std::remove_cvref_t<Range>>;
    constexpr auto is_expr_holder = expr_holder<Range>;

    if constexpr (is_expr || is_expr_holder) {
      ExprContainer cont(std::forward<decltype(rng)>(rng));
      this->append(std::move(cont));
    } else {
      for (auto &&v : rng) {
        append(std::forward<decltype(v)>(v));
      }
    }
  }

  /// tags ctor to move the summands directly
  struct move_only_tag {};

  /// construct a Sum by moving in the summands, no flattening is performed,
  /// but zeros will be omitted and constants added up
  /// @param summands the summands to move in
  explicit Sum(summands_type &&summands, move_only_tag);

  /// append a summand to the sum
  /// @param summand the summand
  Sum &append(ExprContainer summand);

  /// prepend a summand to the sum
  /// @param summand the summand
  Sum &prepend(ExprContainer summand);

  /// Summands accessor
  const summands_type &summands() const;

  /// Summand accessor
  /// @param i summand index
  /// @return ith summand
  const ExprContainer &summand(size_t i) const;

  /// Takes the first @c count elements of the sum
  ExprContainer take_n(size_t count) const;

  /// Takes the first @c count elements of the sum starting with element @c
  /// offset
  ExprContainer take_n(size_t offset, size_t count) const;

  /// @param f Boolean predicate
  /// @returns A sum containing only the summands for which f was true.
  template <std::predicate<const Expr &> Filter>
  ExprContainer filter(Filter &&f) const {
    return ex<Sum>(summands_ | ranges::views::filter(f));
  }

  /// @return true if the number of factors is zero
  bool empty() const;

  /// @return the number of summands in a Sum
  std::size_t size() const;

  std::wstring to_latex() const override;

  Expr::type_id_type type_id() const override;

  /// @brief adjoint of a Sum is a sum of adjoints of its factors
  virtual void adjoint() override;

  Sum &operator+=(const Expr &that);

  Sum &operator-=(const Expr &that);

  ExprIterator begin_subexpr() override;

  ExprIterator end_subexpr() override;

  ConstExprIterator begin_subexpr() const override;

  ConstExprIterator end_subexpr() const override;

 protected:
  std::unique_ptr<Expr> unique_copy() const override;

 private:
  summands_type summands_{};
  std::optional<size_t>
      constant_summand_idx_{};  // points to the constant summand, if any; used
                                // to sum up constants in append/prepend

  /// @return the hash of this object
  /// @note this ensures that hash of a Sum of a single summand is
  /// identical to the hash of the summand itself.
  hash_type memoizing_hash() const override;

  /// @param multipass if true, will do a multipass canonicalization, with extra
  /// cleanup pass after the deep canonization pass
  ExprPtr canonicalize_impl(bool multipass, CanonicalizeOptions opt);

  ExprPtr canonicalize(CanonicalizeOptions opt =
                           CanonicalizeOptions::default_options()) override;

  ExprPtr rapid_canonicalize(
      CanonicalizeOptions opts =
          CanonicalizeOptions::default_options().copy_and_set(
              CanonicalizationMethod::Rapid)) override;

  bool static_equal(const Expr &that) const override;
};  // class Sum

/// @brief utility for eagerly accumulating summands in a hash table
/// Intended to accumulate summands that are already in canonical form and
/// produces canonical Sum
class HashingAccumulator {
 public:
  /// @p summand expr to append to the sum
  /// @p flatten if true, and @p summand is a Sum, will flatten the sum
  HashingAccumulator &append(ExprContainer summand, bool flatten = true);

  ExprContainer make_sum();

  ExprContainer make_canonicalized_sum();

  /// @param canonicalize if true, will sort the summands to canonical order
  /// defined by ExprPtr::operator<
  /// @return summands as a Sum (if have more than 1 summand), Constant (if have
  /// zero summands), or the lone summand itself
  ExprContainer make_expr(bool canonicalize = true);

  bool empty() const;

 private:
  /// @brief Common implementation for make_sum and make_canonicalized_sum
  /// @param canonicalize if true, sort the summands by hash value
  SumPtr make_sum_impl(bool canonicalize);

  container::unordered_set<ExprContainer, sequant::hash::_<ExprContainer>,
                           proportional_to>
      summands_;
};

struct TransformSumExprOptions {
  SEQUANT_DESIGNATED_INIT_ONLY;
  bool canonicalize = true;
  bool flatten = true;
};

/// variant of sequant::transform_reduce for eager sum reduction of Expr's
/// @sa HashingAccumulator
template <typename SizedRange, typename UnaryMapOp>
  requires(meta::is_range_v<std::remove_cvref_t<SizedRange>>)
ExprPtr transform_sum_expr(SizedRange &&rng, const UnaryMapOp &map,
                           const TransformSumExprOptions &options = {}) {
  HashingAccumulator result_acc;
  std::mutex result_mtx;  // serializes updates of result

  auto task = [&result_acc, &result_mtx, &map,
               canonicalize = options.canonicalize,
               flatten = options.flatten](const ExprPtr &input) {
    auto task_result = map(input);
    if (task_result) {
      if (canonicalize) {
        auto bp = task_result->canonicalize();
        if (bp) {
          task_result = bp * task_result;
        }
      }

      std::scoped_lock<std::mutex> lock(result_mtx);
      result_acc.append(task_result, flatten);
    }
  };
  sequant::for_each(std::forward<SizedRange>(rng), task);
  return result_acc.make_expr(options.canonicalize);
}

}  // namespace sequant

#endif  // SEQUANT_EXPRESSIONS_SUM_HPP
