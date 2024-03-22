#include <SeQuant/core/spin/RestrictedDiophantineSolver.hpp>

#include <Eigen/Dense>

#include "catch.hpp"
#include "test_config.hpp"

#include <algorithm>
#include <iostream>
#include <sstream>

using namespace sequant;

struct VectorCmp {
  template <typename L, typename R>
  bool operator()(const L &lhs, const R &rhs) const {
    if (lhs.size() != rhs.size()) {
      return false;
    }

    return std::equal(lhs.begin(), lhs.end(), rhs.begin());
  }
};

template <typename T>
void write_vector(std::stringstream &stream, const T &vec) {
  stream << "(";
  for (const auto &current : vec) {
    stream << current << ",";
  }
  stream << ") ";
}

#define REQUIRE_SOLUTION_SET(solver, solutions)                               \
  {                                                                           \
    auto actual = solver.solve();                                             \
    REQUIRE(actual.size() == solutions.size());                               \
    if (!std::is_permutation(actual.begin(), actual.end(), solutions.begin(), \
                             VectorCmp{})) {                                  \
      std::stringstream sstream;                                              \
      sstream << "Expected ";                                                 \
      for (const auto &current : solutions) write_vector(sstream, current);   \
      sstream << "\nbut got ";                                                \
      for (const auto &current : actual) write_vector(sstream, current);      \
      FAIL(sstream.str());                                                    \
    }                                                                         \
  }

TEST_CASE("DiophantineSolver", "[spin]") {
    SECTION("single equation"){
      // a * x1 + b * x2 = b1
      RestrictedDiophantineSolver solver(1, 2);

      // x1 - x2 = 0
      solver.addTerm(L"a", 1);
      solver.addTerm(L"b", -1);
      solver.endEquation(0);
      std::vector<std::vector<int>> expectedSolutions = {{1, 1}, {-1, -1}};
      REQUIRE_SOLUTION_SET(solver, expectedSolutions);

      // x1 + x2 = 0
      solver.reset();
      solver.addTerm(L"a", 1);
      solver.addTerm(L"b", 1);
      solver.endEquation(0);
      expectedSolutions = {{1, -1}, {-1, 1}};
      REQUIRE_SOLUTION_SET(solver, expectedSolutions);

      // 4 * x1 + 3 * x2 = 0
      solver.reset();
      solver.addTerm(L"a", 4);
      solver.addTerm(L"b", 3);
      solver.endEquation(0);
      expectedSolutions = {};
      REQUIRE_SOLUTION_SET(solver, expectedSolutions);
    }

    SECTION("system of equations"){
      // a * x1 + b * x2 = b1
      // c * x1 + b * x2 = b2
      RestrictedDiophantineSolver solver(2, 3);

      // x1 - x2      = 0
      //      x2 - x3 = -2
      solver.addTerm(L"a", 1);
      solver.addTerm(L"b", -1);
      solver.endEquation(0);
      solver.addTerm(L"b", 1);
      solver.addTerm(L"c", -1);
      solver.endEquation(-2);

      std::vector<std::vector<int>> expectedSolutions = {{-1,-1,1}};
      REQUIRE_SOLUTION_SET(solver, expectedSolutions);

    }
}

#undef REQUIRE_SOLUTION_SET
