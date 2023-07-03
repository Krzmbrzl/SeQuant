#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/sequant.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/rational.hpp>

#include "Utils.hpp"

using namespace sequant;

static IndexFactory idxFactory(nullptr, 1);

void setConvention() {
  IndexSpace::register_instance(L"f", IndexSpace::frozen_occupied,
                                IndexSpace::nullqns, true);
  IndexSpace::register_instance(L"z", IndexSpace::inactive_occupied,
                                IndexSpace::nullqns, true);
  IndexSpace::register_instance(L"y", IndexSpace::active_occupied,
                                IndexSpace::nullqns, true);
  IndexSpace::register_instance(L"o", IndexSpace::occupied, IndexSpace::nullqns,
                                true);
  IndexSpace::register_instance(L"a", IndexSpace::active_unoccupied,
                                IndexSpace::nullqns, true);
  IndexSpace::register_instance(L"v", IndexSpace::inactive_unoccupied,
                                IndexSpace::nullqns, true);
  IndexSpace::register_instance(L"x", IndexSpace::unoccupied,
                                IndexSpace::nullqns, true);
  IndexSpace::register_instance(L"w", IndexSpace::all_active,
                                IndexSpace::nullqns, true);
  IndexSpace::register_instance(L"v", IndexSpace::all, IndexSpace::nullqns,
                                true);
  IndexSpace::register_instance(L"u", IndexSpace::other_unoccupied,
                                IndexSpace::nullqns, true);
  IndexSpace::register_instance(L"t", IndexSpace::complete_unoccupied,
                                IndexSpace::nullqns, true);
  IndexSpace::register_instance(L"p", IndexSpace::complete, IndexSpace::nullqns,
                                true);

  TensorCanonicalizer::set_cardinal_tensor_labels(
      {L"t", L"f", L"g", L"{C_0}", L"{C_0^\\dagger}"});
}

ExprPtr make_op(Tensor tensor) {
  auto bra = tensor.bra();
  auto ket = tensor.ket();
  return ex<Tensor>(std::move(tensor)) *
         ex<FNOperator>(std::move(bra), std::move(ket),
                        get_default_context().vacuum());
}

Index create_index(const IndexSpace::Type &type) {
  return idxFactory.make(IndexSpace::instance(type));
}

ExprPtr f() {
  return make_op(Tensor(L"f", std::vector<Index>{create_index(general)},
                        std::vector<Index>{create_index(general)},
                        Symmetry::antisymm));
}

ExprPtr g() {
  return ex<Constant>(rational{1,4}) *
         make_op(Tensor(
             L"g",
             std::vector<Index>{create_index(general), create_index(general)},
             std::vector<Index>{create_index(general), create_index(general)},
             Symmetry::antisymm));
}

ExprPtr H() { return f() + g(); }
