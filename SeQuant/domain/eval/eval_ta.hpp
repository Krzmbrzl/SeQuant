#ifndef SEQUANT_EVAL_TA_HPP
#define SEQUANT_EVAL_TA_HPP

#include "eval.hpp"

#include <tiledarray.h>
#include <range/v3/all.hpp>

namespace sequant::eval {

auto const braket_to_annot = [](auto const& bk) {
  using ranges::views::join;
  using ranges::views::transform;
  using ranges::views::intersperse;
  return join(bk | transform([](auto const& idx) { return idx.label(); }) |
              intersperse(L",")) |
         ranges::to<std::string>;
};  // braket_to_annot

auto const ords_to_annot = [](auto const& ords) {
  using ranges::accumulate;
  using ranges::views::intersperse;
  using ranges::views::transform;
  auto to_str = [](auto x) { return std::to_string(x); };
  return ranges::accumulate(
      ords | transform(to_str) | intersperse(std::string{","}), std::string{},
      std::plus{});
};  // ords_to_annot

template <typename Tensor_t>
struct yield_leaf {
  size_t const no, nv;
  Tensor_t const &G, &F, &t_vo, &t_vvoo;
  yield_leaf(size_t nocc, size_t nvirt, Tensor_t const& fock,
             Tensor_t const& eri, Tensor_t const& ampl_vo,
             Tensor_t const& ampl_vvoo)
      : no{nocc},
        nv{nvirt},
        G{eri},
        F{fock},
        t_vo{ampl_vo},
        t_vvoo{ampl_vvoo}

  {}

  auto range1_limits(sequant::Tensor const& tensor) {
    return tensor.const_braket() |
           ranges::views::transform([this](auto const& idx) {
             auto ao = sequant::IndexSpace::active_occupied;
             auto au = sequant::IndexSpace::active_unoccupied;
             auto sp = idx.space();
             assert(sp == ao || sp == au);

             return sp == ao ? no : nv;
           });
  }

  Tensor_t operator()(sequant::Tensor const& tensor) {
    if (tensor.label() == L"t") {
      auto rank = tensor.rank();
      assert(rank == 1 || rank == 2);
      return rank == 1 ? t_vo : t_vvoo;
    }

    auto r1_limits = range1_limits(tensor);

    auto trange_vec = r1_limits | ranges::views::transform([](auto x) {
                        return TA::TiledRange1{0, x};
                      }) |
                      ranges::to_vector;

    auto iter_limits =
        r1_limits | ranges::views::transform([this](auto x) {
          return x == no ? std::pair{size_t{0}, no} : std::pair{no, no + nv};
        });

    auto tlabel = tensor.label();
    assert(tlabel == L"g" || tlabel == L"f");

    auto const& big_tensor = tlabel == L"g" ? G : F;

    auto slice =
        TA::TArrayD{big_tensor.world(),
                    TA::TiledRange{trange_vec.begin(), trange_vec.end()}};
    slice.fill(0);
    auto tile_orig = big_tensor.find(0).get();
    auto tile_dest = slice.find(0).get();

    assert(iter_limits.size() == 2 || iter_limits.size() == 4);
    if (iter_limits.size() == 2) {
      for (auto ii = iter_limits[0].first; ii < iter_limits[0].second; ++ii)
        for (auto jj = iter_limits[1].first; jj < iter_limits[1].second; ++jj) {
          tile_dest(ii - iter_limits[0].first,  //
                    jj - iter_limits[1].first) = tile_orig(ii, jj);
        }
    } else {  // 4 iterations
      for (auto ii = iter_limits[0].first; ii < iter_limits[0].second; ++ii)
        for (auto jj = iter_limits[1].first; jj < iter_limits[1].second; ++jj)
          for (auto kk = iter_limits[2].first; kk < iter_limits[2].second; ++kk)
            for (auto ll = iter_limits[3].first; ll < iter_limits[3].second;
                 ++ll) {
              tile_dest(ii - iter_limits[0].first, jj - iter_limits[1].first,
                        kk - iter_limits[2].first, ll - iter_limits[3].first) =
                  tile_orig(ii, jj, kk, ll);
            }
    }

    // return cache.store(hash, slice);
    return slice;
  }
};

template <typename Tensor_t>
Tensor_t inode_evaluate_ta(
    sequant::utils::binary_node<sequant::utils::eval_expr> const& node,
    Tensor_t const& leval, Tensor_t const& reval) {
  assert((node->op() == sequant::utils::eval_expr::eval_op::Sum ||
          node->op() == sequant::utils::eval_expr::eval_op::Prod) &&
         "unsupported intermediate operation");

  auto assert_imaginary_zero = [](sequant::Constant const& c) {
    assert(c.value().imag() == 0 &&
           "complex scalar unsupported for real tensor");
  };

  assert_imaginary_zero(node.left()->scalar());
  assert_imaginary_zero(node.right()->scalar());

  auto this_annot = braket_to_annot(node->tensor().const_braket());
  auto lannot = braket_to_annot(node.left()->tensor().const_braket());
  auto rannot = braket_to_annot(node.right()->tensor().const_braket());

  auto lscal = node.left()->scalar().value().real();
  auto rscal = node.right()->scalar().value().real();

  auto result = Tensor_t{};
  if (node->op() == sequant::utils::eval_expr::eval_op::Prod) {
    // prod
    result(this_annot) = (lscal * rscal) * leval(lannot) * reval(rannot);
  } else {
    // sum
    result(this_annot) = lscal * leval(lannot) + rscal * reval(rannot);
  }

  return result;
}

template <typename Tensor_t, typename Yielder>
Tensor_t evaluate_ta(
    sequant::utils::binary_node<sequant::utils::eval_expr> const& node,
    Yielder& yielder, sequant::utils::cache_manager<Tensor_t>& cman) {
  static_assert(
      std::is_invocable_r_v<Tensor_t, Yielder, sequant::Tensor const&>);

  auto const key = node->hash();

  if (auto&& exists = cman.access(key); exists && exists.value())
    return *exists.value();

  return node.leaf()
             ? *cman.store(key, yielder(node->tensor()))
             : *cman.store(key,
                           inode_evaluate_ta(
                               node, evaluate_ta(node.left(), yielder, cman),
                               evaluate_ta(node.right(), yielder, cman)));
}

struct eval_instance {
  sequant::utils::binary_node<sequant::utils::eval_expr> const& node;

  template <typename Tensor_t, typename Fetcher>
  auto evaluate(Fetcher& f, sequant::utils::cache_manager<Tensor_t>& man) {
    static_assert(
        std::is_invocable_r_v<Tensor_t, Fetcher, sequant::Tensor const&>);

    auto result = evaluate_ta(node, f, man);
    auto const annot = braket_to_annot(node->tensor().const_braket());
    auto scaled = decltype(result){};
    scaled(annot) = node->scalar().value().real() * result(annot);
    return scaled;
  }

  template <typename Tensor_t, typename Fetcher>
  auto evaluate_asymm(Fetcher& f,
                      sequant::utils::cache_manager<Tensor_t>& man) {
    auto result = evaluate(f, man);

    auto asymm_result = decltype(result){result.world(), result.trange()};
    asymm_result.fill(0);

    auto const lannot =
        ords_to_annot(ranges::views::iota(size_t{0}, result.trange().rank()) |
                      ranges::to_vector);

    auto asym_impl = [&result, &asymm_result,
                      &lannot](auto const& pwp) {  // pwp = perm with phase
      asymm_result(lannot) += pwp.phase * result(ords_to_annot(pwp.perm));
    };

    sequant::eval::antisymmetrize_tensor(result.trange().rank(), asym_impl);
    return asymm_result;
  }

  template <typename Tensor_t, typename Fetcher>
  auto evaluate_symm(Fetcher& f, sequant::utils::cache_manager<Tensor_t>& man) {
    auto result = evaluate(f, man);

    auto symm_result = decltype(result){result.world(), result.trange()};
    symm_result.fill(0);

    auto const lannot =
        ords_to_annot(ranges::views::iota(size_t{0}, result.trange().rank()) |
                      ranges::to_vector);

    auto sym_impl = [&result, &symm_result, &lannot](auto const& perm) {
      symm_result(lannot) += result(ords_to_annot(perm));
    };

    sequant::eval::symmetrize_tensor(result.trange().rank(), sym_impl);
    return symm_result;
  }
};  // eval_instance

}  // namespace sequant::eval

#endif  // SEQUANT_EVAL_TA_HPP
