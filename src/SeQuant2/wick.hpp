//
// Created by Eduard Valeyev on 3/23/18.
//

#ifndef SEQUANT2_WICK_HPP
#define SEQUANT2_WICK_HPP

#include <bitset>
#include <mutex>
#include <utility>

#include "op.hpp"
#include "ranges.hpp"
#include "tensor.hpp"

namespace sequant2 {

/// Applies Wick's theorem to a sequence of normal-ordered operators.
///
/// @tparam S particle statistics
template<Statistics S>
class WickTheorem {
 public:
  static constexpr const Statistics statistics = S;
  static_assert(S == Statistics::FermiDirac,
                "WickTheorem not yet implemented for Bose-Einstein");

  explicit WickTheorem(const NormalOperatorSequence<S> &input) : input_(input) {
    assert(input.size() <= max_input_size);
    assert(input.empty() || input.vacuum() != Vacuum::Invalid);
    assert(input.empty() || input.vacuum() != Vacuum::Invalid);
  }

  /// Controls whether next call to compute() will full contractions only or all
  /// (including partial) contractions. By default compute() generates all
  /// contractions.
  /// @param sf if true, will complete full contractions only.
  /// @return reference to @c *this , for daisy-chaining
  WickTheorem &full_contractions(bool fc) {
    full_contractions_ = fc;
    return *this;
  }
  /// Controls whether next call to compute() will assume spin-free or
  /// spin-orbital normal-ordered operators By default compute() assumes
  /// spin-orbital operators.
  /// @param sf if true, will complete full contractions only.
  WickTheorem &spinfree(bool sf) {
    spinfree_ = sf;
    return *this;
  }
  /// Controls whether next call to compute() will reduce the result
  /// By default compute() will not perform reduction.
  /// @param r if true, compute() will reduce the result.
  WickTheorem &reduce(bool r) {
    reduce_ = r;
    return *this;
  }

  /// Specifies the external indices; by default assume all indices are summed
  /// over
  /// @param ext_inds external (nonsummed) indices
  template <typename IndexContainer>
  WickTheorem &set_external_indices(IndexContainer &&external_indices) {
    if constexpr (std::is_convertible_v<IndexContainer,
                                        decltype(external_indices_)>)
      external_indices_ = std::forward<IndexContainer>(external_indices);
    else {
      ranges::for_each(std::forward<IndexContainer>(external_indices),
                       [this](auto &&v) {
                         auto result = this->external_indices_.emplace(v);
                         assert(result.second);
                       });
    }
    return *this;
  }

  /// Ensures that the given pairs of normal operators are connected; by default
  /// will not constrain connectivity
  /// @param op_index_pairs the list of pairs of op indices to be connected in
  /// the result
  template<typename IndexPairContainer>
  WickTheorem &set_op_connections(IndexPairContainer &&op_index_pairs) {
    for (const auto &opidx_pair : op_index_pairs) {
      if (opidx_pair.first < 0 || opidx_pair.first >= input_.size()) {
        throw std::invalid_argument(
            "WickTheorem::set_op_connections: op index out of range");
      }
      if (opidx_pair.second < 0 || opidx_pair.second >= input_.size()) {
        throw std::invalid_argument(
            "WickTheorem::set_op_connections: op index out of range");
      }
    }
    if (!op_index_pairs.empty()) {
      op_connections_.resize(input_.size());
      for (auto &v : op_connections_) {
        v.set();
      }
      for (const auto &opidx_pair : op_index_pairs) {
        op_connections_[opidx_pair.first].reset(opidx_pair.second);
        op_connections_[opidx_pair.second].reset(opidx_pair.first);
      }
    }

    return *this;
  }

  /// Computes and returns the result
  /// @param count_only if true, will return a vector of default-initialized
  /// values, useful if only interested in the total count
  /// @return the result of applying Wick's theorem, i.e. a sum of {prefactor,
  /// normal operator} pairs
  ExprPtr compute(const bool count_only = false) const {
    if (!full_contractions_)
      throw std::logic_error(
          "WickTheorem::compute: full_contractions=false not yet supported");
    if (spinfree_)
      throw std::logic_error(
          "WickTheorem::compute: spinfree=true not yet supported");
    auto result = compute_nontensor_wick(count_only);
    if (reduce_ && !count_only) {
      reduce(result);
      canonicalize(result);
    }
    return std::move(result);
  }

 private:
  static constexpr size_t max_input_size =
      32;  // max # of operators in the input sequence
  const NormalOperatorSequence<S> &input_;
  bool full_contractions_ = false;
  bool spinfree_ = false;
  bool reduce_ = false;
  container::vector<Index> external_indices_;
  // for each operator specifies the reverse bitmask of connections (0 = must
  // connect)
  container::svector<std::bitset<max_input_size>> op_connections_;

  /// carries state down the stack of recursive calls
  struct NontensorWickState {
    NontensorWickState(const NormalOperatorSequence<S> &opseq)
        : opseq(opseq),
          level(0),
          count_only(false),
          op_connections(opseq.size()),
          adjacency_matrix(opseq.size() * (opseq.size() - 1) / 2, 0) {
      compute_size();
    }
    NormalOperatorSequence<S> opseq;  //!< current state of operator sequence
    std::size_t opseq_size;           //!< current size of opseq
    Product sp;                       //!< current prefactor
    int level;                        //!< level in recursive wick call stack
    bool count_only;                  //!< if true, only update result size
    container::svector<std::bitset<max_input_size>>
        op_connections;  //!< bitmask of connections for each op (1 = connected)
    container::svector<size_t>
        adjacency_matrix;  // number of connections between each normop, only
    // lower triangle is kept

    void compute_size() {
      opseq_size = 0;
      for (const auto &op : opseq) opseq_size += op.size();
    }
    void reset(const NormalOperatorSequence<S>& o) {
      sp = Product{};
      opseq = o;
      op_connections = decltype(op_connections)(opseq.size());
      adjacency_matrix =
          decltype(adjacency_matrix)(opseq.size() * (opseq.size() - 1) / 2, 0);
      compute_size();
    }

    template<typename T>
    static auto lowtri_idx(T i, T j) {
      assert(i != j);
      auto ii = std::max(i, j);
      auto jj = std::min(i, j);
      return ii * (ii - 1) / 2 + jj;
    }

    /// @brief Updates connectivity if contraction satisfies target connectivity

    /// If the target connectivity will be violated by this contraction, keep
    /// the state unchanged and return false
    template<typename Cursor>
    inline bool connect(const container::svector<std::bitset<max_input_size>> &
    target_op_connections,
                        const Cursor &op1_cursor, const Cursor &op2_cursor) {
      if (target_op_connections.empty())  // if no constraints, all is fair game
        return true;

      auto result = true;

      // local vars
      const auto op1_idx = op1_cursor.range_ordinal();
      const auto op2_idx = op2_cursor.range_ordinal();
      const auto op1_op2_connected = op_connections[op1_idx].test(op2_idx);

      // update the connectivity
      if (!op1_op2_connected) {
        op_connections[op1_idx].set(op2_idx);
        op_connections[op2_idx].set(op1_idx);
      }

      // test if op1 has enough remaining indices to satisfy
      const auto nidx_op1_remain = op1_cursor.range_iter()->size() - 1;  // how many indices op1 has minus this index
      const auto nidx_op1_needs =
          (op_connections[op1_idx] | target_op_connections[op1_idx])
              .flip()
              .count();
      if (nidx_op1_needs > nidx_op1_remain) {
        if (!op1_op2_connected) {
          op_connections[op1_idx].reset(op2_idx);
          op_connections[op2_idx].reset(op1_idx);
        }
        return false;
      }

      // test if op2 has enough remaining indices to satisfy
      const auto nidx_op2_remain = op2_cursor.range_iter()->size() - 1;  // how many indices op2 has minus this index
      const auto nidx_op2_needs =
          (op_connections[op2_idx] | target_op_connections[op2_idx])
              .flip()
              .count();
      if (nidx_op2_needs > nidx_op2_remain) {
        if (!op1_op2_connected) {
          op_connections[op1_idx].reset(op2_idx);
          op_connections[op2_idx].reset(op1_idx);
        }
        return false;
      }

      adjacency_matrix[lowtri_idx(op1_idx, op2_idx)] += 1;

      return true;  // not yet implemented
    }
    /// @brief Updates connectivity when contraction is reversed
    template<typename Cursor>
    inline void disconnect(
        const container::svector<std::bitset<max_input_size>> &
        target_op_connections,
        const Cursor &op1_cursor, const Cursor &op2_cursor) {
      if (target_op_connections.empty())  // if no constraints, all is fair game
        return;

      // local vars
      const auto op1_idx = op1_cursor.range_ordinal();
      const auto op2_idx = op2_cursor.range_ordinal();
      assert(op_connections[op1_idx].test(op2_idx));

      auto &adjval = adjacency_matrix[lowtri_idx(op1_idx, op2_idx)];
      assert(adjval > 0);
      adjval -= 1;
      if (adjval == 0) {
        op_connections[op1_idx].reset(op2_idx);
        op_connections[op2_idx].reset(op1_idx);
      }
    }
  };
  /// Applies most naive version of Wick's theorem, where sign rule involves
  /// counting Ops
  ExprPtr compute_nontensor_wick(const bool count_only) const {
    std::vector<std::pair<Product, NormalOperator<S>>>
        result;      //!< current value of the result
    std::mutex mtx;  // used in critical sections updating the result
    auto result_plus_mutex = std::make_pair(&result, &mtx);
    NontensorWickState state(input_);
    state.count_only = count_only;

    recursive_nontensor_wick(result_plus_mutex, state);

    // convert result to an Expr
    // if result.size() == 0, return null ptr
    // TODO revise if decide to use Constant(0)
    ExprPtr result_expr;
    if (result.size() == 1) {  // if result.size() == 1, return Product
      result_expr = ex<Product>(std::move(result[0].first));
    } else if (result.size() > 1) {
      auto sum = std::make_shared<Sum>();
      for (auto &term : result) {
        sum->append(ex<Product>(std::move(term.first)));
      }
      result_expr = sum;
    }
    return result_expr;
  }

  void recursive_nontensor_wick(
      std::pair<std::vector<std::pair<Product, NormalOperator<S>>> *,
                std::mutex *> &result,
      NontensorWickState &state) const {
    // if full contractions needed, make contractions involving first index with
    // another index, else contract any index i with index j (i<j)
    if (full_contractions_) {
      using opseq_view_type = flattened_rangenest<NormalOperatorSequence<S>>;
      auto opseq_view = opseq_view_type(&state.opseq);
      using std::begin;
      using std::end;
      auto opseq_view_begin = begin(opseq_view);

      // optimization: can't contract fully if first op is not a qp annihilator
      if (!is_qpannihilator(*opseq_view_begin, input_.vacuum())) return;

      auto op_iter = opseq_view_begin;
      ++op_iter;
      for (; op_iter != end(opseq_view);) {
        if (op_iter != opseq_view_begin &&
            ranges::get_cursor(op_iter).range_iter() !=
                ranges::get_cursor(opseq_view_begin)
                    .range_iter()  // can't contract within same normop
            ) {
          if (can_contract(*opseq_view_begin, *op_iter, input_.vacuum()) &&
              state.connect(op_connections_, ranges::get_cursor(op_iter),
                            ranges::get_cursor(
                                opseq_view_begin))  // check connectivity
            // constraints if needed
              ) {
            //            std::wcout << "level " << state.level << ":
            //            contracting " << to_latex(*opseq_view_begin) << " with
            //            " << to_latex(*op_iter) << std::endl; std::wcout << "
            //            current opseq = " << to_latex(state.opseq) <<
            //            std::endl;

            // update the phase, if needed
            double phase = 1;
            if (statistics == Statistics::FermiDirac) {
              const auto distance =
                  ranges::get_cursor(op_iter).ordinal() -
                      ranges::get_cursor(opseq_view_begin).ordinal() - 1;
              if (distance % 2) {
                phase *= -1;
              }
            }

            // update the prefactor and opseq
            Product sp_copy = state.sp;
            state.sp.append(
                phase, contract(*opseq_view_begin, *op_iter, input_.vacuum()));
            // remove from back to front
            Op<S> right = *op_iter;
            ranges::get_cursor(op_iter).erase();
            --state.opseq_size;
            Op<S> left = *opseq_view_begin;
            ranges::get_cursor(opseq_view_begin).erase();
            --state.opseq_size;

            //            std::wcout << "  opseq after contraction = " <<
            //            to_latex(state.opseq) << std::endl;

            // update the result if nothing left to contract and have a nonzero
            // result
            if (state.opseq_size == 0 && !state.sp.empty()) {
              result.second->lock();
              //              std::wcout << "got " << to_latex(state.sp) <<
              //              std::endl;
              if (!state.count_only)
                result.first->push_back(
                    std::make_pair(std::move(state.sp.deep_copy()), NormalOperator<S>{}));
              else
                result.first->resize(result.first->size() + 1);
              //              std::wcout << "now up to " << result.first->size()
              //              << " terms" << std::endl;
              result.second->unlock();
            }

            if (state.opseq_size != 0) {
              ++state.level;
              recursive_nontensor_wick(result, state);
              --state.level;
            }

            // restore the prefactor and opseq
            state.sp = std::move(sp_copy);
            // restore from front to back
            ranges::get_cursor(opseq_view_begin).insert(std::move(left));
            ++state.opseq_size;
            ranges::get_cursor(op_iter).insert(std::move(right));
            ++state.opseq_size;
            state.disconnect(op_connections_, ranges::get_cursor(op_iter),
                             ranges::get_cursor(opseq_view_begin));
            //            std::wcout << "  restored opseq = " <<
            //            to_latex(state.opseq) << std::endl;
          }
          ++op_iter;
        } else
          ++op_iter;
      }
    } else
      assert(false);  // full_contraction_=false not implemented yet, should
    // result in error earlier
  }

 public:
  static bool can_contract(const Op<S> &left, const Op<S> &right,
                           Vacuum vacuum = get_default_context().vacuum()) {
    if (is_qpannihilator<S>(left, vacuum) && is_qpcreator<S>(right, vacuum)) {
      const auto qpspace_left = qpannihilator_space<S>(left, vacuum);
      const auto qpspace_right = qpcreator_space<S>(right, vacuum);
      const auto qpspace_common = intersection(qpspace_left, qpspace_right);
      if (qpspace_common != IndexSpace::null_instance()) return true;
    }
    return false;
  }

  static std::shared_ptr<Expr> contract(
      const Op<S> &left, const Op<S> &right,
      Vacuum vacuum = get_default_context().vacuum()) {
    assert(can_contract(left, right, vacuum));
//    assert(
//        !left.index().has_proto_indices() &&
//            !right.index().has_proto_indices());  // I don't think the logic is
    // correct for dependent indices
    if (is_pure_qpannihilator<S>(left, vacuum) &&
        is_pure_qpcreator<S>(right, vacuum))
      return overlap(left.index(), right.index());
    else {
      const auto qpspace_left = qpannihilator_space<S>(left, vacuum);
      const auto qpspace_right = qpcreator_space<S>(right, vacuum);
      const auto qpspace_common = intersection(qpspace_left, qpspace_right);
      const auto index_common = Index::make_tmp_index(qpspace_common);

      // preserve bra/ket positions of left & right
      const auto left_is_ann = left.action() == Action::annihilate;
      assert(left_is_ann || right.action() == Action::annihilate);

      if (qpspace_common != left.index().space() &&
          qpspace_common !=
              right.index().space()) {  // may need 2 overlaps if neither space
        // is pure qp creator/annihilator
        auto result = std::make_shared<Product>();
        result->append(1, left_is_ann ? overlap(left.index(), index_common) : overlap(index_common, left.index()));
        result->append(1, left_is_ann ? overlap(index_common, right.index()) : overlap(right.index(), index_common));
        return result;
      } else {
        return left_is_ann ? overlap(left.index(), right.index()) : overlap(right.index(), left.index());
      }
    }
  }

 public:  // TODO make these members private once WickTheorem can work on full
  // expressions (not on sequences of normal operators) directly
  /// @param[in,out] on input, Wick theorem result, on output the result of
  /// reducing the overlaps
  void reduce(ExprPtr& expr) const;
};

using BWickTheorem = WickTheorem<Statistics::BoseEinstein>;
using FWickTheorem = WickTheorem<Statistics::FermiDirac>;

}  // namespace sequant2

#include "wick.impl.hpp"

#endif  // SEQUANT2_WICK_HPP
