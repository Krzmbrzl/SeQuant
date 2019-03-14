//
// Created by Eduard Valeyev on 3/31/18.
//

#ifndef SEQUANT_WICK_IMPL_HPP
#define SEQUANT_WICK_IMPL_HPP

#if __has_include(<execution>)
#include <execution>
#define SEQUANT_HAS_EXECUTION_HEADER
#endif

namespace sequant {

namespace detail {

struct zero_result : public std::exception {};

/// @brief computes index replacement rules

/// If using orthonormal representation, overlaps are Kronecker deltas, hence
/// summations can be reduced by index replacements. Reducing sums over dummy
/// (internal) indices uses 2 rules:
/// - if a Kronecker delta binds 2 internal indices I and J, replace them with a
/// new internal index
///   representing intersection of spaces of I and J, !!remove delta!!
/// - if a Kronecker delta binds an internal index J and an external index I:
///   - if space of J includes space of I, replace J with I, !!remove delta!!
///   - if space of J is a subset of space of I, replace J with a new internal
///     index representing intersection of spaces of I and J, !!keep the delta!!
/// @throw zero_result if @c product is zero
inline container::map<Index, Index> compute_index_replacement_rules(
    std::shared_ptr<Product> &product,
    const container::set<Index> &external_indices,
    const std::set<Index, Index::LabelCompare> &all_indices) {
  expr_range exrng(product);

  /// this ensures that all temporary indices have unique *labels* (not just
  /// unique *full labels*)
  auto index_validator = [&all_indices](const Index &idx) {
    return all_indices.find(idx) == all_indices.end();
  };
  IndexFactory idxfac(index_validator);
  container::map<Index /* src */, Index /* dst */> result;  // src->dst

  // computes an index in intersection of space1 and space2
  auto make_intersection_index = [&idxfac](const IndexSpace &space1,
                                           const IndexSpace &space2) {
    const auto intersection_space = intersection(space1, space2);
    if (intersection_space == IndexSpace::null_instance()) throw zero_result{};
    return idxfac.make(intersection_space);
  };

  // transfers proto indices from idx (if any) to img
  auto proto = [](const Index &img, const Index &idx) {
    if (idx.has_proto_indices()) {
      if (img.has_proto_indices()) {
        assert(img.proto_indices() == idx.proto_indices());
        return img;
      } else
        return Index(img, idx.proto_indices());
    } else {
      assert(!img.has_proto_indices());
      return img;
    }
  };

  // adds src->dst or src->intersection(dst,current_dst)
  auto add_rule = [&result, &proto, &make_intersection_index](const Index &src,
                                                              const Index &dst) {
    auto src_it = result.find(src);
    if (src_it == result.end())  // if brand new, add the rule
      result[src] = proto(dst, src);
    else {  // else modify the destination of the existing rule to the
      // intersection
      const auto &old_dst = src_it->second;
      assert(old_dst.proto_indices() == src.proto_indices());
      if (dst.space() != old_dst.space()) {
        src_it->second =
            proto(make_intersection_index(old_dst.space(), dst.space()), src);
      }
    }
  };

  // adds src1->dst and src2->dst; if src1->dst1 and/or src2->dst2 already
  // exist the existing rules are updated to map to the intersection of dst1,
  // dst2 and dst
  auto add_rules = [&result, &idxfac, &proto, &make_intersection_index](
                       const Index &src1, const Index &src2, const Index &dst) {
    // are there replacement rules already for src{1,2}?
    auto src1_it = result.find(src1);
    auto src2_it = result.find(src2);
    const auto has_src1_rule = src1_it != result.end();
    const auto has_src2_rule = src2_it != result.end();

    // which proto-indices should dst1 and dst2 inherit? a source index without
    // proto indices will inherit its source counterpart's indices, unless it
    // already has its own protoindices: <a_ij|p> = <a_ij|a_ij> (hence replace p
    // with a_ij), but <a_ij|p_kl> = <a_ij|a_kl> != <a_ij|a_ij> (hence replace
    // p_kl with a_kl)
    const auto &dst1_proto =
        !src1.has_proto_indices() && src2.has_proto_indices() ? src2 : src1;
    const auto &dst2_proto =
        !src2.has_proto_indices() && src1.has_proto_indices() ? src1 : src2;

    if (!has_src1_rule && !has_src2_rule) {  // if brand new, add the rules
      result.insert(std::make_pair(src1, proto(dst, dst1_proto)));
      assert(!result.empty());
      result.insert(std::make_pair(src2, proto(dst, dst2_proto)));
      assert(!result.empty());
    } else if (has_src1_rule &&
               !has_src2_rule) {  // update the existing rule for src1
      const auto &old_dst1 = src1_it->second;
      assert(old_dst1.proto_indices() == dst1_proto.proto_indices());
      if (dst.space() != old_dst1.space()) {
        src1_it->second = proto(
            make_intersection_index(old_dst1.space(), dst.space()), dst1_proto);
      }
      result[src2] = src1_it->second;
    } else if (!has_src1_rule &&
               has_src2_rule) {  // update the existing rule for src2
      const auto &old_dst2 = src2_it->second;
      assert(old_dst2.proto_indices() == dst2_proto.proto_indices());
      if (dst.space() != old_dst2.space()) {
        src2_it->second = proto(
            make_intersection_index(old_dst2.space(), dst.space()), dst2_proto);
      }
      result[src1] = src2_it->second;
    } else {  // update both of the existing rules
      const auto &old_dst1 = src1_it->second;
      const auto &old_dst2 = src2_it->second;
      const auto new_dst_space =
          (dst.space() != old_dst1.space() || dst.space() != old_dst2.space())
          ? intersection(old_dst1.space(), old_dst2.space(), dst.space())
          : dst.space();
      if (new_dst_space == IndexSpace::null_instance()) throw zero_result{};
      Index new_dst;
      if (new_dst_space == old_dst1.space()) {
        new_dst = old_dst1;
        if (new_dst_space == old_dst2.space() && old_dst2 < new_dst) {
          new_dst = old_dst2;
        }
        if (new_dst_space == dst.space() && dst < new_dst) {
          new_dst = dst;
        }
      } else if (new_dst_space == old_dst2.space()) {
        new_dst = old_dst2;
        if (new_dst_space == dst.space() && dst < new_dst) {
          new_dst = dst;
        }
      } else if (new_dst_space == dst.space()) {
        new_dst = dst;
      } else
        new_dst = idxfac.make(new_dst_space);
      result[src1] = proto(new_dst,dst1_proto);
      result[src2] = proto(new_dst,dst2_proto);
    }
  };

  /// this makes the list of replacements ... we do not mutate the expressions
  /// to keep the information about which indices are related
  for (auto it = ranges::begin(exrng); it != ranges::end(exrng);
       ++it) {
    const auto &factor = *it;
    if (factor->type_id() == Expr::get_type_id<Tensor>()) {
      const auto &tensor = static_cast<const Tensor &>(*factor);
      if (tensor.label() == L"S") {
        assert(tensor.bra().size() == 1);
        assert(tensor.ket().size() == 1);
        const auto &bra = tensor.bra().at(0);
        const auto &ket = tensor.ket().at(0);
        assert(bra != ket);

        const auto bra_is_ext = ranges::find(external_indices, bra) !=
            ranges::end(external_indices);
        const auto ket_is_ext = ranges::find(external_indices, ket) !=
            ranges::end(external_indices);

        const auto intersection_space = intersection(bra.space(), ket.space());
        assert(intersection_space != IndexSpace::null_instance());

        if (!bra_is_ext && !ket_is_ext) {  // int + int
          const auto new_dummy = idxfac.make(intersection_space);
          add_rules(bra, ket, new_dummy);
        } else if (bra_is_ext && !ket_is_ext) {  // ext + int
          if (includes(bra.space(), ket.space())) {
            add_rule(ket, bra);
          } else {
            add_rule(ket, idxfac.make(intersection_space));
          }
        } else if (!bra_is_ext && ket_is_ext) {  // int + ext
          if (includes(ket.space(), bra.space())) {
            add_rule(bra, ket);
          } else {
            add_rule(bra, idxfac.make(intersection_space));
          }
        }
      }
    }
  }

  return result;
}

/// @return true if made any changes
inline bool apply_index_replacement_rules(
    std::shared_ptr<Product> &product,
    const container::map<Index, Index> &const_replrules,
    const container::set<Index> &external_indices,
    std::set<Index, Index::LabelCompare> &all_indices) {
  // to be able to use map[]
  auto &replrules = const_cast<container::map<Index, Index> &>(const_replrules);

  expr_range exrng(product);

  /// this recursively applies replacement rules until result does not
  /// change removes the deltas that are no longer needed
  const bool tag_transformed_indices =
      true;  // to replace indices when maps image and domain overlap, tag
             // transformed indices
#ifndef NDEBUG
  // assert that tensors_ indices are not tagged if going to tag indices
  if (tag_transformed_indices) {
    for (auto it = ranges::begin(exrng); it != ranges::end(exrng); ++it) {
      const auto &factor = *it;
      if (factor->is<Tensor>()) {
        auto &tensor = factor->as<Tensor>();
        assert(ranges::none_of(tensor.const_braket(), [](const Index &idx) {
          return idx.tag().has_value();
        }));
      }
    }
  }
#endif
  bool mutated = false;
  bool pass_mutated = false;
  do {
    pass_mutated = false;

    for (auto it = ranges::begin(exrng); it != ranges::end(exrng);) {
      const auto &factor = *it;
      if (factor->is<Tensor>()) {
        bool erase_it = false;
        auto &tensor = factor->as<Tensor>();

        /// replace indices
        pass_mutated &=
            tensor.transform_indices(const_replrules, tag_transformed_indices);

        if (tensor.label() == L"S") {
          const auto &bra = tensor.bra().at(0);
          const auto &ket = tensor.ket().at(0);

          if (bra.proto_indices() == ket.proto_indices()) {
            const auto bra_is_ext = ranges::find(external_indices, bra) !=
                                    ranges::end(external_indices);
            const auto ket_is_ext = ranges::find(external_indices, ket) !=
                                    ranges::end(external_indices);

#ifndef NDEBUG
            const auto intersection_space =
                intersection(bra.space(), ket.space());
#endif

            if (!bra_is_ext && !ket_is_ext) {  // int + int
#ifndef NDEBUG
              if (replrules.find(bra) != replrules.end() &&
                  replrules.find(ket) != replrules.end())
                assert(replrules[bra].space() == replrules[ket].space());
#endif
              erase_it = true;
            } else if (bra_is_ext && !ket_is_ext) {  // ext + int
              if (includes(bra.space(), ket.space())) {
#ifndef NDEBUG
                if (replrules.find(ket) != replrules.end())
                  assert(replrules[ket].space() == bra.space());
#endif
                erase_it = true;
              } else {
#ifndef NDEBUG
                if (replrules.find(ket) != replrules.end())
                  assert(replrules[ket].space() == intersection_space);
#endif
              }
            } else if (!bra_is_ext && ket_is_ext) {  // int + ext
              if (includes(ket.space(), bra.space())) {
#ifndef NDEBUG
                if (replrules.find(bra) != replrules.end())
                  assert(replrules[bra].space() == ket.space());
#endif
                erase_it = true;
              } else {
#ifndef NDEBUG
                if (replrules.find(bra) != replrules.end())
                  assert(replrules[bra].space() == intersection_space);
#endif
              }
            }

            if (erase_it) {
              pass_mutated = true;
              *it = ex<Constant>(1);
            }
          }  // matching proto indices
        }    // Kronecker delta
      }
      ++it;
    }
    mutated |= pass_mutated;
  } while (pass_mutated);  // keep replacing til fixed point

  // assert that tensors_ indices are not tagged if going to tag indices
  if (tag_transformed_indices) {
    for (auto it = ranges::begin(exrng); it != ranges::end(exrng); ++it) {
      const auto &factor = *it;
      if (factor->is<Tensor>()) {
        factor->as<Tensor>().reset_tags();
      }
    }
  }

  // update all_indices
  std::set<Index, Index::LabelCompare> all_indices_new;
  ranges::for_each(
      all_indices, [&const_replrules, &all_indices_new](const Index &idx) {
        auto dst_it = const_replrules.find(idx);
        all_indices_new.insert(dst_it != const_replrules.end() ? dst_it->second
                                                               : idx);
      });
  std::swap(all_indices_new, all_indices);

  return mutated;
}

/// If using orthonormal representation, resolves Kronecker deltas (=overlaps
/// between indices in orthonormal spaces) in summations
/// @throw zero_result if @c expr is zero
inline void reduce_wick_impl(std::shared_ptr<Product> &expr,
                             const container::set<Index> &external_indices) {
  if (get_default_context().metric() == IndexSpaceMetric::Unit) {
    bool pass_mutated = false;
    do {
      pass_mutated = false;

      // extract current indices
      std::set<Index, Index::LabelCompare> all_indices;
      ranges::for_each(*expr, [&all_indices](const auto &factor) {
        if (factor->template is<Tensor>()) {
          ranges::for_each(factor->template as<const Tensor>().braket(),
                           [&all_indices](const Index &idx) {
                             auto result = all_indices.insert(idx);
                           });
        }
      });

      const auto replacement_rules =
          compute_index_replacement_rules(expr, external_indices, all_indices);

      if (Logger::get_instance().wick_reduce){
        std::wcout << "reduce_wick_impl(expr, external_indices):\n  expr = "
                   << expr->to_latex() << "\n  external_indices = ";
        ranges::for_each(external_indices, [](auto &index) {
          std::wcout << index.label() << " ";
        });
        std::wcout << "\n  replrules = ";
        ranges::for_each(replacement_rules, [](auto &index) {
          std::wcout << to_latex(index.first) << "\\to"
                     << to_latex(index.second) << "\\,";
        });
        std::wcout.flush();
      }

      if (!replacement_rules.empty()) {
        pass_mutated = apply_index_replacement_rules(
            expr, replacement_rules, external_indices, all_indices);
      }

      if (Logger::get_instance().wick_reduce) {
        std::wcout << "\n  result = " << expr->to_latex() << std::endl;
      }

    } while (pass_mutated);  // keep reducing until stop changing
  } else
    abort();  // programming error?
}

}  // namespace detail

template <Statistics S>
ExprPtr WickTheorem<S>::compute(const bool count_only) {
  // have an Expr as input? Apply recursively ...
  if (expr_input_) {
    /// expand, then apply recursively to products
    expand(expr_input_);
    // if sum, canonicalize and apply to each summand ...
    if (expr_input_->is<Sum>()) {
      canonicalize(expr_input_);
      assert(!expr_input_->as<Sum>().empty());

      // parallelize over summands
      auto result = std::make_shared<Sum>();
      std::mutex result_mtx;  // serializes updates of result
      auto summands = expr_input_->as<Sum>().summands();

#if SEQUANT_HAS_EXECUTION_HEADER
      auto wick_task = [&result, &result_mtx, this,
                        &count_only](const ExprPtr &input) {
        WickTheorem wt(input, *this);
        auto task_result = wt.compute(count_only);
        if (task_result) {
          std::scoped_lock<std::mutex> lock(result_mtx);
          result->append(task_result);
        }
      };
      std::for_each(std::execution::par_unseq, begin(summands), end(summands),
                    wick_task);
#else
      auto wick_task = [&summands, &result, &result_mtx, this,
                        &count_only](size_t task_id) {
        auto &summand = summands[task_id];
        WickTheorem wt(summand, *this);
        auto task_result = wt.compute(count_only);
        if (task_result) {
          std::scoped_lock<std::mutex> lock(result_mtx);
          result->append(task_result);
        }
      };
      parallel_for_each(wick_task, summands.size());
#endif

      // if the sum is empty return zero
      // if the sum has 1 summand, return it directly
      ExprPtr result_expr = result;
      if (result->summands().size() == 0) {
        result_expr = ex<Constant>(0);
      }
      if (result->summands().size() == 1)
        result_expr = std::move(result->summands()[0]);
      return result_expr;
    }
    // ... else if a product, find NormalOperatorSequence, if any, and compute
    // ...
    else if (expr_input_->is<Product>()) {
      canonicalize(expr_input_);
      // NormalOperators should be all at the end
      auto first_nop_it = ranges::find_if(
          *expr_input_,
          [](const ExprPtr &expr) { return expr->is<NormalOperator<S>>(); });
      // if have ops, split into prefactor and op sequence
      if (first_nop_it != ranges::end(*expr_input_)) {

        // TODO compute tensor-nop connections and topological partitions
        if (use_topology_) {
          //        abort(); // not yet implemented

          // hack for CC only ... assume that third and higher nops come from Ts, to determine their equivalence just sort them by rank
          container::multimap<size_t, size_t> rank_2_nop_idx; // maps rank to list of normal operator indices
          auto nops_view = ranges::view::filter(
              *expr_input_,
              [](const ExprPtr &expr) { return expr->is<NormalOperator<S>>(); });
//          assert(ranges::size(nops_view) >= 2);  // proj + hamiltonian
          size_t max_nop_rank = 0;
          size_t nop_idx = 0;
          ranges::for_each(nops_view, [&nop_idx, &max_nop_rank, &rank_2_nop_idx](const ExprPtr& expr) {
            // should actually use ranges::view::drop_exactly
            if (nop_idx>=2) {
              const auto& nop = expr->as<NormalOperator<S>>();
              const auto nop_rank = nop.rank();
              max_nop_rank = std::max(max_nop_rank, nop_rank);
              rank_2_nop_idx.insert(std::make_pair(nop_rank, nop_idx));
            }
            ++nop_idx;
          });

          // convert rank_2_nop_idx to nop partitions
          if (!rank_2_nop_idx.empty()) {
            container::vector<container::vector<size_t>> nop_partitions;
            assert(rank_2_nop_idx.upper_bound(8) == rank_2_nop_idx.end());
            assert(max_nop_rank >= 1);
            for(size_t rank=1; rank<=max_nop_rank; ++rank) {
              auto rank_nop_range = rank_2_nop_idx.equal_range(rank);
              const auto rank_nop_range_size = rank_nop_range.second - rank_nop_range.first;
              if (rank_nop_range_size > 1) {
                nop_partitions.push_back({});
                auto& partition = nop_partitions.back();
                partition.reserve(rank_nop_range_size);
                for(auto it=rank_nop_range.first; it != rank_nop_range.second; ++it) {
                  partition.push_back(it->second);
                }
              }
            }
            if (!nop_partitions.empty())
              this->set_op_partitions(nop_partitions);
          }

        }

        ExprPtr prefactor =
            ex<CProduct>(expr_input_->as<Product>().scalar(), ExprPtrList{});
        bool found_op = false;
        ranges::for_each(
            *expr_input_, [this, &found_op, &prefactor](const ExprPtr &factor) {
              if (factor->is<NormalOperator<S>>()) {
                input_.push_back(factor->as<NormalOperator<S>>());
                found_op = true;
              } else {
                assert(factor->is_cnumber());
                assert(!found_op);  // make sure that ops are at the end
                *prefactor *= *factor;
              }
            });
        if (!input_.empty()) {
          auto result = compute_nopseq(count_only);
          if (result) {  // simplify if obtained nonzero ...
            result = prefactor * result;
            expand(result);
            this->reduce(result);
            rapid_simplify(result);
            canonicalize(result);
            rapid_simplify(
                result);  // rapid_simplify again since canonization may produce
                          // new opportunities (e.g. terms cancel, etc.)
          } else
            result = ex<Constant>(0);
          return result;
        }
      } else {  // product does not include ops
        return expr_input_;
      }
    }
    // ... else if NormalOperatorSequence already, compute ...
    else if (expr_input_->is<NormalOperatorSequence<S>>()) {
      input_ = expr_input_->as<NormalOperatorSequence<S>>();
      // NB no simplification possible for a bare product w/ full contractions
      // ... partial contractions will need simplification
      return compute_nopseq(count_only);
    } else  // ... else do nothing
      return expr_input_;
  } else  // given a NormalOperatorSequence instead of an expression
    return compute_nopseq(count_only);
  abort();
}

template <Statistics S>
void WickTheorem<S>::reduce(ExprPtr &expr) const {
  // there are 2 possibilities: expr is a single Product, or it's a Sum of
  // Products
  if (expr->type_id() == Expr::get_type_id<Product>()) {
    auto expr_cast = std::static_pointer_cast<Product>(expr);
    try {
      detail::reduce_wick_impl(expr_cast, external_indices_);
      expr = expr_cast;
    } catch (detail::zero_result &) {
      expr = std::make_shared<Constant>(0);
    }
  } else {
    assert(expr->type_id() == Expr::get_type_id<Sum>());
    for (auto &&subexpr : *expr) {
      assert(subexpr->is<Product>());
      auto subexpr_cast = std::static_pointer_cast<Product>(subexpr);
      try {
        detail::reduce_wick_impl(subexpr_cast, external_indices_);
        subexpr = subexpr_cast;
      }
      catch (detail::zero_result &) {
        subexpr = std::make_shared<Constant>(0);
      }
    }
  }
}

}  // namespace sequant

#endif  // SEQUANT_WICK_IMPL_HPP
