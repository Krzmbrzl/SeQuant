#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/optimize.hpp>
#include <SeQuant/core/tensor.hpp>
#include <SeQuant/domain/mbpt/convention.hpp>
#include <SeQuant/domain/mbpt/spin.hpp>
#include <SeQuant/domain/mbpt/sr.hpp>

#include <iostream>

using namespace sequant;
using namespace sequant::mbpt::sr;

int main(int argc, const char** argv) {
  sequant::detail::OpIdRegistrar op_id_registrar;
  sequant::set_default_context(
      Context(Vacuum::SingleProduct, IndexSpaceMetric::Unit,
              BraKetSymmetry::conjugate, SPBasis::spinorbital));
  mbpt::set_default_convention();
  TensorCanonicalizer::register_instance(
      std::make_shared<DefaultTensorCanonicalizer>());

  // 1. Construct Hbar
  ExprPtr hbar = op::H();
  ExprPtr H_Tk = hbar;
  for (int64_t k = 1; k <= 4; ++k) {
    H_Tk = simplify(ex<Constant>(rational{1, k}) * H_Tk * op::T_(2));
    hbar += H_Tk;
  }

  std::wcout << L"Hbar:\n" << to_latex_align(hbar) << "\n\n";

  // 2. project onto doubles manifold, screen, lower to tensor form and wick it
  ExprPtr equations;

  const std::size_t projection = 2;

  // 2.a. screen out terms that cannot give nonzero after projection onto
  // <p|
  std::shared_ptr<Sum> screendedTerms;

  for (const ExprPtr& term : *hbar) {
    assert(term->is<Product>() || term->is<op_t>());

    if (op::raises_vacuum_to_rank(term, projection)) {
      if (!screendedTerms) {
        screendedTerms = std::make_shared<Sum>(ExprPtrList{term});
      } else {
        screendedTerms->append(term);
      }
    }
  }

  if (projection > 0) {
    // 2.b project onto <p|, i.e. multiply by P(p)
    ExprPtr P_hbar = simplify(op::P(projection) * screendedTerms);

    // 2.c compute vacuum expectation value (Wick theorem)
    equations = op::vac_av(P_hbar);
  } else {
    // Use equation as-is (no projection required)
    equations = op::vac_av(screendedTerms);
  }

  equations = simplify(equations);

  auto ext_indices = external_indices(equations);

  std::wcout << "External index groups:\n";
  for (const auto& currentGroup : ext_indices) {
    assert(currentGroup.size() == 2);
    std::wcout << "{ " << to_latex(currentGroup.front()) << ", "
               << to_latex(currentGroup.back()) << " }\n";
  }
  std::wcout << "\n\n";

  ExprPtr specialTraced =
      simplify(closed_shell_CC_spintrace(equations, projection));
  ExprPtr generalTraced = simplify(biorthogonal_transform(
      simplify(spintrace(equations, ext_indices)), projection, ext_indices));

  std::wcout << "Special:\n"
             << to_latex_align(specialTraced) << "\n\nGeneral:\n"
             << to_latex_align(generalTraced)
             << "\n\nAre equal (after expand S operator): " << std::boolalpha
             << (*simplify(S_maps(specialTraced)) == *generalTraced) << "\n";
}
