#include <SeQuant/core/spin/RestrictedDiophantineSolver.hpp>

#include <SeQuant/core/math.hpp>

#include <algorithm>
#include <bitset>
#include <cassert>

namespace sequant {

RestrictedDiophantineSolver::RestrictedDiophantineSolver(
    std::size_t numEquations, std::size_t numVariables)
    : m_coefficientMatrix(Eigen::MatrixXi(numEquations, numVariables)),
      m_inhomogeneity(Eigen::VectorXi(numEquations)) {
  m_coefficientMatrix.setZero();
  m_inhomogeneity.setZero();
  m_variableNames.reserve(numVariables);
}

void RestrictedDiophantineSolver::addTerm(std::wstring_view variable,
                                          int coefficient) {
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
    m_variableNames.push_back(std::move(variable));

    assert(pos < m_coefficientMatrix.cols());

    m_coefficientMatrix(m_currentEquation, pos) = coefficient;
  }
}

void RestrictedDiophantineSolver::endEquation(int result) {
  assert(m_currentEquation < m_inhomogeneity.size());
  m_inhomogeneity(m_currentEquation) = result;
  m_currentEquation++;
}

RestrictedDiophantineSolver::solution_vectors
RestrictedDiophantineSolver::solve() const {
  // There are some standard algorithms for solving linear diophantine systems,
  // namely via the Smith normal form or via the Hermite normal form.
  // Additionally, the field of linear programming also provides ways of solving
  // such systems.
  // These algorithms would have to be adapted to restrict the solution vector
  // entries to the desired range.
  //
  // However, for now we follow a simple brute-force solving strategy, assuming
  // that the dimension of the solution vector space 2^N is small (N being the
  // amount of unknowns).
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

const container::vector<std::wstring_view>
    &RestrictedDiophantineSolver::getVariables() const {
  return m_variableNames;
}

void RestrictedDiophantineSolver::reset() {
  m_coefficientMatrix.setZero();
  m_inhomogeneity.setZero();
  m_variableNames.clear();
  m_currentEquation = 0;
}

}  // namespace sequant
