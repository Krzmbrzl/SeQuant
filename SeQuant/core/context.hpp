#ifndef SEQUANT_CORE_CONTEXT_HPP
#define SEQUANT_CORE_CONTEXT_HPP

#include "attr.hpp"
#include "index_space_registry.hpp"

#include "utility/context.hpp"

namespace sequant {

/// Specifies second quantization context, such as vacuum choice, whether index
/// spaces are orthonormal, sizes of index spaces, etc.
class Context {
 public:
  struct Defaults {
    constexpr static auto vacuum = Vacuum::Physical;
    constexpr static auto metric = IndexSpaceMetric::Unit;
    constexpr static auto braket_symmetry = BraKetSymmetry::conjugate;
    constexpr static auto spbasis = sequant::SPBasis::spinorbital;
    constexpr static auto first_dummy_index_ordinal = 100;
  };

  /// standard full-form constructor

  /// @param isr_sptr a shared_ptr to an IndexSpaceRegistry object
  /// @param vac a Vacuum object
  /// @param m an IndexSpaceMetric object
  /// @param bks a BraKetSymmetry object
  /// @param spb single-particle basis (spin-free or spin-dependent)
  /// @param fdio first dummy index ordinal
  explicit Context(std::shared_ptr<IndexSpaceRegistry> isr_sptr, Vacuum vac,
                   IndexSpaceMetric m = Defaults::metric,
                   BraKetSymmetry bks = Defaults::braket_symmetry,
                   SPBasis spb = Defaults::spbasis,
                   std::size_t fdio = Defaults::first_dummy_index_ordinal);

  /// @brief same as the standard ctor, using IndexSpaceRegistry passed by value

  /// @param isr an IndexSpaceRegistry object
  /// @param vac a Vacuum object
  /// @param m an IndexSpaceMetric object
  /// @param bks a BraKetSymmetry object
  /// @param spb single-particle basis (spin-free or spin-dependent)
  /// @param fdio first dummy index ordinal
  explicit Context(IndexSpaceRegistry isr, Vacuum vac,
                   IndexSpaceMetric m = Defaults::metric,
                   BraKetSymmetry bks = Defaults::braket_symmetry,
                   SPBasis spb = Defaults::spbasis,
                   std::size_t fdio = Defaults::first_dummy_index_ordinal);

  /// @brief same as the standard ctor, using default-constructed
  /// IndexSpaceRegistry

  /// @param vac a Vacuum object
  /// @param m an IndexSpaceMetric object
  /// @param bks a BraKetSymmetry object
  /// @param spb single-particle basis (spin-free or spin-dependent)
  /// @param fdio first dummy index ordinal
  explicit Context(Vacuum vac, IndexSpaceMetric m = Defaults::metric,
                   BraKetSymmetry bks = Defaults::braket_symmetry,
                   SPBasis spb = Defaults::spbasis,
                   std::size_t fdio = Defaults::first_dummy_index_ordinal);

  /// default constructor, equivalent to Context({}, Vacuum::Physical,
  /// IndexSpaceMetric::Unit, BraKetSymmetry::conjugate,
  /// sequant::SPBasis::spinorbital, 100)
  Context() = default;

  ~Context() = default;

  Context(const Context&) = default;
  Context& operator=(const Context&) = default;

  /// \return Vacuum of this context
  Vacuum vacuum() const;
  /// @return a constant pointer to the IndexSpaceRegistry for this context.
  std::shared_ptr<const IndexSpaceRegistry> index_space_registry() const;
  /// @return a pointer to the IndexSpaceRegistry for this context.
  std::shared_ptr<IndexSpaceRegistry> mutable_index_space_registry() const;
  /// \return IndexSpaceMetric of this context
  IndexSpaceMetric metric() const;
  /// \return BraKetSymmetry of this context
  BraKetSymmetry braket_symmetry() const;
  /// \return SPBasis of this context
  SPBasis spbasis() const;
  /// \return first ordinal of the dummy indices generated by calls to
  /// Index::next_tmp_index when this context is active
  std::size_t first_dummy_index_ordinal() const;

  /// Sets the Vacuum for this context, convenient for chaining
  /// \param vacuum Vacuum
  /// \return ref to `*this`, for chaining
  Context& set(Vacuum vacuum);
  /// sets the IndexSpaceRegistry for this context
  /// \param IndexSpaceRegistry ISR
  /// \return ref to '*this' for chaining
  Context& set(IndexSpaceRegistry ISR);
  /// Sets the IndexSpaceMetric for this context, convenient for chaining
  /// \param metric IndexSpaceMetric
  /// \return ref to `*this`, for chaining
  Context& set(IndexSpaceMetric metric);
  /// Sets the BraKetSymmetry for this context, convenient for chaining
  /// \param braket_symmetry BraKetSymmetry
  /// \return ref to `*this`, for chaining
  Context& set(BraKetSymmetry braket_symmetry);
  /// Sets the SPBasis for this context, convenient for chaining
  /// \param spbasis SPBasis
  /// \return ref to `*this`, for chaining
  Context& set(SPBasis spbasis);
  /// Sets the first dummy index ordinal for this context, convenient for
  /// chaining \param first_dummy_index_ordinal the first dummy index ordinal
  /// \return ref to `*this`, for chaining
  Context& set_first_dummy_index_ordinal(std::size_t first_dummy_index_ordinal);

  /// @return the IndexRegistry object

 private:
  std::shared_ptr<IndexSpaceRegistry> idx_space_reg_ =
      std::make_shared<IndexSpaceRegistry>();
  Vacuum vacuum_ = Defaults::vacuum;
  IndexSpaceMetric metric_ = Defaults::metric;
  BraKetSymmetry braket_symmetry_ = Defaults::braket_symmetry;
  SPBasis spbasis_ = Defaults::spbasis;
  std::size_t first_dummy_index_ordinal_ = Defaults::first_dummy_index_ordinal;
};

/// old name of Context is a deprecated alias
using SeQuant
    [[deprecated("use sequant::Context instead of sequant::SeQuant")]] =
        Context;

/// Context object equality comparison
/// \param ctx1
/// \param ctx2
/// \return true if \p ctx1 and \p ctx2 are equal
/// \warning does not compare index registries
bool operator==(const Context& ctx1, const Context& ctx2);

/// Context object inequality comparison
/// \param ctx1
/// \param ctx2
/// \return true if \p ctx1 and \p ctx2 are not equal
/// \warning does not compare index registries
bool operator!=(const Context& ctx1, const Context& ctx2);

/// \name manipulation of implicit context for SeQuant
/// \warning all of these are re-entrant

/// @{

/// @brief access default Context for the given Statistics
/// @param s Statistics
/// @return the default context used for Statistics @p s
const Context& get_default_context(Statistics s = Statistics::Arbitrary);

/// @brief sets default Context for the given Statistics
/// @param s Statistics
void set_default_context(const Context& ctx,
                         Statistics s = Statistics::Arbitrary);

/// @brief sets default Context for several given Statistics
/// @param ctxs a Statistics->Context map
void set_default_context(const std::map<Statistics, Context>& ctxs);

/// @brief resets default Contexts for all statistics to their initial values
void reset_default_context();

/// @brief changes default contexts
/// @param ctx Context objects for one or more statistics
/// @return a move-only ContextResetter object whose destruction will reset the
/// default context to the previous value
/// @example
/// ```cpp
/// {
///   auto resetter = set_scoped_default_context({{Statistics::Arbitrary,
///   ctx}});
///   // ctx is now the default context for all statistics
/// } // leaving scope, resetter is destroyed, default context is reset back to
/// the old value
/// ```
[[nodiscard]] detail::ImplicitContextResetter<std::map<Statistics, Context>>
set_scoped_default_context(const std::map<Statistics, Context>& ctx);

/// @brief changes default context for arbitrary statistics
/// @note equivalent to `set_scoped_default_context({{Statistics::Arbitrary,
/// ctx}})`
[[nodiscard]] detail::ImplicitContextResetter<std::map<Statistics, Context>>
set_scoped_default_context(const Context& ctx);

///@}

}  // namespace sequant

#endif
