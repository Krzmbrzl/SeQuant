#include <SeQuant/core/expr.hpp>
#include <SeQuant/core/index.hpp>
#include <SeQuant/core/space.hpp>
#include <SeQuant/core/tensor.hpp>

// Define aliases for index spaces
static auto active = sequant::IndexSpace::active_unoccupied;
static auto occ = sequant::IndexSpace::occupied;
static auto virt = sequant::IndexSpace::inactive_unoccupied;
static auto general = sequant::IndexSpace::complete;

/**
 * Configures SeQuant to use our conventions
 */
void setConvention();

/**
 * Creates an expression containing the given Tensor as well as a corresponding
 * set of creator and annihilator operators
 */
sequant::ExprPtr make_op(sequant::Tensor tensor);

/**
 * Creates a unique index of the given type
 */
sequant::Index create_index(const sequant::IndexSpace::Type &type);

/**
 * @returns The fock operator
 */
sequant::ExprPtr f();

/**
 * @returns The two-electron interaction operator
 */
sequant::ExprPtr g();

/**
 * @returns The Hamiltonian (f + g)
 */
sequant::ExprPtr H();
