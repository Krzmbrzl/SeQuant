#ifndef SEQUANT_DIOPHANTINE_SOLVER_HPP
#define SEQUANT_DIOPHANTINE_SOLVER_HPP

#include <SeQuant/core/container.hpp>

#include <Eigen/Dense>

#include <cstddef>
#include <string>

namespace sequant {

/// This solver can solve systems of linear diophantine equations such that the
/// solution vector elements are restricted to be in {-1,+1}.
/// A system of linear diophantine equations is a system of linear equations
/// where the solution vector elements are restricted to integer values.
class RestrictedDiophantineSolver {
 public:
	 using solution_vectors = container::vector<Eigen::VectorXi>;

  RestrictedDiophantineSolver(std::size_t numEquations,
                              std::size_t numVariables);

  void addTerm(std::wstring_view variable, int coefficient);

  void endEquation(int result);

  [[nodiscard]] solution_vectors solve() const;

  [[nodiscard]] const container::vector<std::wstring_view> &getVariables() const;

  void reset();

 private:
  Eigen::MatrixXi m_coefficientMatrix;
  Eigen::VectorXi m_inhomogeneity;
  container::vector<std::wstring_view> m_variableNames;
  std::size_t m_currentEquation = 0;
};

}  // namespace sequant

#endif
