#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/rational.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>

#include "Utils.hpp"

using namespace sequant;

static IndexFactory idxFactory(nullptr, 1);

void registerSpace(IndexSpace::Type type, std::wstring label) {
  IndexSpace::register_instance(label, type, IndexSpace::nullqns, true);
  IndexSpace::register_instance(label + L"↑", type, IndexSpace::alpha, true);
  IndexSpace::register_instance(label + L"↓", type, IndexSpace::beta, true);
}

void setConvention() {
  // Base index spaces (occ., act. & virt.)
  registerSpace(IndexSpace::active_occupied, L"o");
  registerSpace(IndexSpace::active_unoccupied, L"v");
  registerSpace(IndexSpace::active, L"a");

  // Internal (occ. + act.) and external (act. virt.) index spaces
  registerSpace(IndexSpace::maybe_occupied, L"I");
  registerSpace(IndexSpace::maybe_unoccupied, L"A");

  // General indices
  registerSpace(IndexSpace::complete, L"p");

  // Unused index spaces (we have to define them, just in case though)
  registerSpace(IndexSpace::frozen_occupied, L"l");
  registerSpace(IndexSpace::inactive_occupied, L"m");
  registerSpace(IndexSpace::occupied, L"n");
  registerSpace(IndexSpace::active_maybe_occupied, L"q");
  registerSpace(IndexSpace::active_maybe_unoccupied, L"r");
  registerSpace(IndexSpace::inactive_unoccupied, L"s");
  registerSpace(IndexSpace::unoccupied, L"t");
  registerSpace(IndexSpace::all_active, L"u");
  registerSpace(IndexSpace::all, L"w");
  registerSpace(IndexSpace::other_unoccupied, L"x");
  registerSpace(IndexSpace::complete_unoccupied, L"y");
  registerSpace(IndexSpace::complete_maybe_unoccupied, L"z");

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
  return ex<Constant>(rational{1, 4}) *
         make_op(Tensor(
             L"g",
             std::vector<Index>{create_index(general), create_index(general)},
             std::vector<Index>{create_index(general), create_index(general)},
             Symmetry::antisymm));
}

ExprPtr H() { return f() + g(); }
