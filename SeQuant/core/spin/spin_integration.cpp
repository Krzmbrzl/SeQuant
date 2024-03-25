#include <SeQuant/core/spin/spin_integration.hpp>

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/spin/restricted_diophantine_solver.hpp>
#include <SeQuant/core/tensor.hpp>

namespace sequant {

using Solver = RestrictedDiophantineSolver<Index>;

template <typename Range>
std::size_t count_unique_indices(const Range &range) {
  container::set<Index> indices;

  for (const Index &idx : range) {
    indices.insert(idx);
  }

  return indices.size();
}

ExprPtr spin_integrate(const Constant &constant) { return constant.clone(); }

ExprPtr spin_integrate(const Variable &variable) { return variable.clone(); }

ExprPtr spin_integrate(const Tensor &tensor) {
  const std::size_t nIndices = count_unique_indices(tensor.braket());

  // For now, we assume all tensors ought to have singlet-symmetry
  const int targetMultiplicity = 0;

  Solver::variable_list indices;
  Solver::solution_vectors solutions;
  if (tensor.symmetry() == Symmetry::nonsymm) {
    if (tensor.bra_rank() != tensor.ket_rank()) {
      throw std::runtime_error(
          "Don't know how to spin-integrate particle non-conserving operators "
          "without index symmetries");
    }
    const std::size_t nEquations = tensor.rank();

    Solver solver(nEquations, nIndices);

    auto bra_indices = tensor.bra();
    auto ket_indices = tensor.ket();

    for (std::size_t i = 0; i < tensor.rank(); ++i) {
      solver.addTerm(bra_indices[i], 1);
      solver.addTerm(ket_indices[i], -1);
      solver.endEquation(0);
    }

    solutions = solver.solve();
    indices = solver.getVariables();
  } else {
    const std::size_t nEquations = 1;
    Solver solver(nEquations, nIndices);

    for (const Index &bra_idx : tensor.bra()) {
      solver.addTerm(bra_idx, 1);
    }

    for (const Index &ket_idx : tensor.ket()) {
      solver.addTerm(ket_idx, -1);
    }

    solver.endEquation(0);

    solutions = solver.solve();
    indices = solver.getVariables();
  }

  if (solutions.empty()) {
    return ex<Constant>(0);
  }

  SumPtr integrated = std::dynamic_pointer_cast<Sum>(ex<Sum>());

  container::map<Index, Index> idx_replacements;
  idx_replacements.reserve(indices.size());
  for (std::size_t i = 0; i < solutions.size(); ++i) {
    const auto &current_solution = solutions[i];
    assert(current_solution.size() == indices.size());

    idx_replacements.clear();

    for (std::size_t k = 0; k < indices.size(); ++k) {
      assert(current_solution[k] == 1 || current_solution[k] == -1);

      Index old = indices[i];

	  // TODO: Creation of index with spin label is not working yet
	  // -> cmp. make_index_with_spincase in spin.cpp
      idx_replacements.insert(
          {old,
           Index(indices[i].label(),
                 IndexSpace(IndexSpace::nonnulltype,
                            current_solution[k] == -1 ? IndexSpace::beta
                                                      : IndexSpace::alpha))});
    }

    ExprPtr copy = tensor.clone();
    bool success = copy.as<Tensor>().transform_indices(idx_replacements);
    assert(success);
    (void)success;

    integrated->append(std::move(copy));
  }

  assert(integrated->size() > 0);

  return integrated->size() == 1 ? integrated->summand(0) : integrated;
}

ExprPtr spin_integrate(const Sum &sum) {
  container::vector<ExprPtr> integrated;
  for (const ExprPtr &summand : sum) {
    // TODO: parallelize?
    integrated.push_back(spin_integrate(summand));
  }

  if (integrated.empty()) {
    return ex<Constant>(0);
  }

  return integrated.size() == 1 ? integrated[0]
                                : ex<Sum>(std::move(integrated));
}

ExprPtr spin_integrate(const Product &product) {
  // TODO
  // TODO: expand all terms / implement dealing with them properly
}

ExprPtr spin_integrate(const ExprPtr &expression) {
  if (expression.is<Constant>()) {
    return spin_integrate(expression.as<Constant>());
  } else if (expression.is<Variable>()) {
    return spin_integrate(expression.as<Variable>());
  } else if (expression.is<Tensor>()) {
    return spin_integrate(expression.as<Tensor>());
  } else if (expression.is<Sum>()) {
    return spin_integrate(expression.as<Sum>());
  } else if (expression.is<Product>()) {
    return spin_integrate(expression.as<Product>());
  }

  throw std::runtime_error("Unhandled expression type in spin_integrate");
}

}  // namespace sequant
