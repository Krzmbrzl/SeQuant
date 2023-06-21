#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/op.hpp>
#include <SeQuant/core/sequant.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/core/wick.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>

int main() {
  using namespace sequant;

  set_default_context(SeQuant(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
                              BraKetSymmetry::conjugate, SPBasis::spinorbital));
  mbpt::set_default_convention();

  IndexSpace occ(IndexSpace::occupied, IndexSpace::nullqns);
  IndexSpace virt(IndexSpace::unoccupied, IndexSpace::nullqns);
  IndexSpace general(IndexSpace::complete, IndexSpace::nullqns);
  Index i(L"i", occ), j(L"j", occ), a(L"a", virt), b(L"b", virt),
      p(L"p", general), q(L"q", general);

  // Create tensor object
  auto f = ex<Tensor>(L"f", std::vector<Index>{p}, std::vector<Index>{q},
                      Symmetry::antisymm);
  auto t = ex<Tensor>(L"t", std::vector<Index>{a}, std::vector<Index>{i},
                      Symmetry::antisymm);

  // Create (Normal)Operator object that contains necessary indices
  auto fop = ex<FNOperator>(std::vector<Index>{p}, std::vector<Index>{q},
                            get_default_context().vacuum());
  auto top = ex<FNOperator>(std::vector<Index>{a}, std::vector<Index>{i},
                            get_default_context().vacuum());
  // Multiply expressions as needed
  auto overallExpr = ex<Constant>(0.5) * f * fop * t * top;

  std::wcout << to_latex(overallExpr) << " = "
             << to_latex(FWickTheorem{overallExpr}
                             .set_external_indices(std::array<Index, 0>{})
                             .full_contractions(true)
                             .compute())
             << std::endl;

  return 0;
}
