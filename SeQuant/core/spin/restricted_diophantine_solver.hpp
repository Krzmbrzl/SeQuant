#ifndef SEQUANT_DIOPHANTINE_SOLVER_HPP
#define SEQUANT_DIOPHANTINE_SOLVER_HPP

#include <SeQuant/core/container.hpp>
#include <SeQuant/core/math.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <bitset>
#include <cassert>
#include <cstddef>
#include <string>

namespace sequant {

/// This solver can solve systems of linear diophantine equations such that the
/// solution vector elements are restricted to be in {-1,+1}.
/// A system of linear diophantine equations is a system of linear equations
/// where the solution vector elements are restricted to integer values.
template <typename Variable>
class RestrictedDiophantineSolver {
 public:
  using solution_vectors = container::vector<Eigen::VectorXi>;
  using variable_list = container::vector<Variable>;

  RestrictedDiophantineSolver(std::size_t numEquations,
                              std::size_t numVariables)
      : m_coefficientMatrix(Eigen::MatrixXi(numEquations, numVariables)),
        m_inhomogeneity(Eigen::VectorXi(numEquations)) {
    m_coefficientMatrix.setZero();
    m_inhomogeneity.setZero();
    m_variableNames.reserve(numVariables);
  }

  void addTerm(const Variable& variable, int coefficient) {
    auto iter =
        std::find(m_variableNames.begin(), m_variableNames.end(), variable);

    assert(m_currentEquation < m_coefficientMatrix.rows());

    if (iter != m_variableNames.end()) {
      // Variable already known
      auto pos = std::distance(m_variableNames.begin(), iter);
      assert(pos >= 0);

      assert(pos < m_coefficientMatrix.cols());

      m_coefficientMatrix(m_currentEquation, pos) = coefficient;
    } else {
      // Variable not seen yet
      std::size_t pos = m_variableNames.size();
      m_variableNames.push_back(variable);

      assert(pos < m_coefficientMatrix.cols());

      m_coefficientMatrix(m_currentEquation, pos) = coefficient;
    }
  }

  void endEquation(int result) {
    assert(m_currentEquation < m_inhomogeneity.size());
    m_inhomogeneity(m_currentEquation) = result;
    m_currentEquation++;
  }

  [[nodiscard]] solution_vectors solve() const {
    // There are some standard algorithms for solving linear diophantine
    // systems, namely via the Smith normal form or via the Hermite normal form.
    // Additionally, the field of linear programming also provides ways of
    // solving such systems. These algorithms would have to be adapted to
    // restrict the solution vector entries to the desired range.
    //
    // However, for now we follow a simple brute-force solving strategy,
    // assuming that the dimension of the solution vector space 2^N is small (N
    // being the amount of unknowns).
    assert(m_variableNames.size() == m_coefficientMatrix.cols());

    const std::size_t nVariables = m_variableNames.size();

    using counter_t = std::uint32_t;
    // -1 required to avoid endless loop below
    assert(m_inhomogeneity.size() < sizeof(counter_t) - 1);

    Eigen::VectorXi trialVector(nVariables);

    solution_vectors solutions;

    for (counter_t i = 0; i < pow2(nVariables); ++i) {
      const std::bitset<sizeof(counter_t)> bitset(i);

      // Prepare trial vector to reflect current trial solution (encoded in i's
      // bit representation)
      for (std::size_t k = 0; k < nVariables; ++k) {
        // 0 -> -1; 1 -> +1
        trialVector(k) = bitset.test(k) ? +1 : -1;
      }

      if (m_coefficientMatrix * trialVector != m_inhomogeneity) {
        // The trial vector does not represent a valid solution
        continue;
      }

      // Valid solution found
      solutions.push_back(trialVector);
    }

    return solutions;
  }

  [[nodiscard]] const variable_list& getVariables() const {
    return m_variableNames;
  }

  void reset() {
    m_coefficientMatrix.setZero();
    m_inhomogeneity.setZero();
    m_variableNames.clear();
    m_currentEquation = 0;
  }

 private:
  Eigen::MatrixXi m_coefficientMatrix;
  Eigen::VectorXi m_inhomogeneity;
  variable_list m_variableNames;
  std::size_t m_currentEquation = 0;
};

}  // namespace sequant

#endif
