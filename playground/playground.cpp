#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/context.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

int main() {
  using namespace sequant;

  set_default_context(Context(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
                              BraKetSymmetry::conjugate, SPBasis::spinorbital));
  //TensorCanonicalizer::set_cardinal_tensor_labels(std::vector<std::wstring>{L"f", L"t", L"λ"});

  mbpt::set_default_convention();

  // Index spaces:
  // i: active_unoccupied
  // m: occupied
  // a: active_unoccupied
  // e: unoccupied
  // x: all_active
  // p: all
  // α': other_unoccupied
  // α: complete_unoccupied
  // κ: complete

  // Create tensor object
  auto f = ex<Tensor>(L"f", std::vector<Index>{Index(L"κ_1")},
                      std::vector<Index>{Index(L"κ_2")}, Symmetry::antisymm);
  auto t1 = ex<Tensor>(L"t", std::vector<Index>{Index(L"a_1")},
                      std::vector<Index>{Index(L"m_1")}, Symmetry::antisymm);
  auto t2 = ex<Tensor>(L"t", std::vector<Index>{Index(L"a_3")},
                      std::vector<Index>{Index(L"m_3")}, Symmetry::antisymm);
  auto l = ex<Tensor>(L"λ", std::vector<Index>{Index(L"m_2")},
                      std::vector<Index>{Index(L"e_2")}, Symmetry::antisymm);

  // Create (Normal)Operator object that contains necessary indices
  auto fop = ex<FNOperator>(f->as<Tensor>().bra(), f->as<Tensor>().ket(),
                            get_default_context().vacuum());
  auto t1op = ex<FNOperator>(t1->as<Tensor>().bra(), t1->as<Tensor>().ket(),
                            get_default_context().vacuum());
  auto t2op = ex<FNOperator>(t2->as<Tensor>().bra(), t2->as<Tensor>().ket(),
                            get_default_context().vacuum());
  auto lop = ex<FNOperator>(l->as<Tensor>().bra(), l->as<Tensor>().ket(),
                            get_default_context().vacuum());

  // Multiply expressions as needed
  //auto overallExpr = l * lop * f * fop * t * top;
  auto overallExpr = t1 * t1op * f * fop * t2 * t2op;

  std::wcout << to_latex(overallExpr) << " = "
             << to_latex(FWickTheorem{overallExpr}
                             .set_external_indices(std::array<Index, 0>{})
                             .full_contractions(true)
                             .compute())
             << std::endl;

  return 0;
}
