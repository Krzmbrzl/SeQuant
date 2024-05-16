//
// Created by Conner Masteran on 4/16/24.
//

#ifndef SEQUANT_INDEX_SPACE_REGISTRY_HPP
#define SEQUANT_INDEX_SPACE_REGISTRY_HPP

#include "space.hpp"

#include <range/v3/algorithm/sort.hpp>
#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/unique.hpp>

#include <boost/hana.hpp>
#include <boost/hana/ext/std/integral_constant.hpp>

namespace sequant {

inline namespace space_tags {
struct IsVacuumOccupied {};
struct IsReferenceOccupied {};
struct IsComplete {};
struct IsHole {};
struct IsParticle {};

constexpr auto is_vacuum_occupied = IsVacuumOccupied{};
constexpr auto is_reference_occupied = IsReferenceOccupied{};
constexpr auto is_complete = IsComplete{};
constexpr auto is_hole = IsHole{};
constexpr auto is_particle = IsParticle{};

}  // namespace space_tags

/// @brief set of known IndexSpace objects

/// Each IndexSpace object has hardwired base key (label) that gives
/// indexed expressions appropriate semantics; e.g., spaces referred to by
/// indices in \f$ t_{p_1}^{i_1} \f$ are defined if IndexSpace objects with
/// base keys \f$ p \f$ and \f$ i \f$ are registered.
/// Since index spaces have set-theoretic semantics, the user must
/// provide complete set of unions/intersects of the base spaces to
/// cover all possible IndexSpace objects that can be generated in their
/// program.
///
/// Registry contains 2 parts: set of IndexSpace objects (managed by a
/// `std::shared_ptr`, see IndexSpaceRegistry::spaces()) and specification of
/// various spaces (vacuum, reference, complete, etc.). Copy semantics is thus
/// partially shallow, with spaces shared between copies. This allows to have
/// multiple registries share same set of spaces but have different
/// specifications of vacuum, reference, etc.; this is useful for providing
/// different contexts for fermions and bosons, for example.
///
/// @note  @c IndexSpaces with @c typeattr corresponding to occupied indices
/// will always be lower than typeattr corresponding to unoccupied orbitals.
class IndexSpaceRegistry {
 public:
  IndexSpaceRegistry()
      : spaces_(std::make_shared<
                container::set<IndexSpace, IndexSpace::KeyCompare>>()) {
    // register nullspace
    this->add(IndexSpace::null);
  }

  /// constructs an IndexSpaceRegistry from an existing set of IndexSpace
  /// objects
  IndexSpaceRegistry(
      std::shared_ptr<container::set<IndexSpace, IndexSpace::KeyCompare>>
          spaces)
      : spaces_(std::move(spaces)) {}

  IndexSpaceRegistry(const IndexSpaceRegistry& other) = default;
  IndexSpaceRegistry(IndexSpaceRegistry&& other) = default;
  IndexSpaceRegistry& operator=(const IndexSpaceRegistry& other) = default;
  IndexSpaceRegistry& operator=(IndexSpaceRegistry&& other) = default;

  const auto& spaces() const { return spaces_; }

  decltype(auto) begin() const { return spaces_->cbegin(); }
  decltype(auto) end() const { return spaces_->cend(); }

  /// @brief retrieve an IndexSpace from the registry by the label
  /// @param label can be numbered or the @c base_key
  /// @return IndexSpace associated with that key
  /// @throw IndexSpace::bad_key if matching space is not found
  template <typename S, typename = meta::EnableIfAnyStringConvertible<S>>
  const IndexSpace& retrieve(S&& label) const {
    auto it =
        spaces_->find(IndexSpace::reduce_key(to_basic_string_view(label)));
    if (it == spaces_->end()) {
      throw IndexSpace::bad_key(label);
    }
    return *it;
  }

  /// @brief retrieve an IndexSpace from the registry by its type and quantum
  /// numbers
  /// @param type IndexSpace::Type
  /// @param qns IndexSpace::QuantumNumbers
  /// @return IndexSpace associated with that key.
  /// @throw std::invalid_argument if matching space is not found
  const IndexSpace& retrieve(const IndexSpace::Type& type,
                             const IndexSpace::QuantumNumbers& qns) const {
    auto it = std::find_if(
        spaces_->begin(), spaces_->end(),
        [&](const auto& is) { return is.type() == type && is.qns() == qns; });
    if (it == spaces_->end()) {
      throw std::invalid_argument(
          "IndexSpaceRegistry::retrieve(type,qn): missing { IndexSpace::Type=" +
          std::to_string(type.to_int32()) + " , IndexSpace::QuantumNumbers=" +
          std::to_string(qns.to_int32()) + " } combination");
    }
    return *it;
  }

  /// @brief retrieve an IndexSpace from the registry by the IndexSpace::Attr
  /// @param space_attr an IndexSpace::Attr
  /// @return IndexSpace associated with that key.
  /// @throw std::invalid_argument if matching space is not found
  const IndexSpace& retrieve(const IndexSpace::Attr& space_attr) const {
    auto it = std::find_if(
        spaces_->begin(), spaces_->end(),
        [&space_attr](const IndexSpace& s) { return s.attr() == space_attr; });
    using std::cend;
    if (it == cend(*spaces_)) {
      throw std::invalid_argument(
          "IndexSpaceRegistry::retrieve(attr): missing { IndexSpace::Type=" +
          std::to_string(space_attr.type().to_int32()) +
          " , IndexSpace::QuantumNumbers=" +
          std::to_string(space_attr.qns().to_int32()) + " } combination");
    }
    return *it;
  }

  /// @name adding IndexSpace objects to the registry
  /// @{

  /// @brief add an IndexSpace to this registry.
  /// @param IS an IndexSpace
  /// @return reference to `this`
  /// @throw std::invalid_argument if `IS.base_key()` or `IS.attr()` matches
  /// an already registered IndexSpace
  IndexSpaceRegistry& add(const IndexSpace& IS) {
    auto it = spaces_->find(IS.base_key());
    if (it != spaces_->end()) {
      throw std::invalid_argument(
          "IndexSpaceRegistry::add(is): already have an IndexSpace associated "
          "with is.base_key(); if you are trying to replace the IndexSpace use "
          "IndexSpaceRegistry::replace(is)");
    } else {
      // make sure there are no duplicate IndexSpaces whose attribute is
      // IS.attr()
      if (ranges::any_of(*spaces_,
                         [&IS](auto&& is) { return IS.attr() == is.attr(); })) {
        throw std::invalid_argument(
            "IndexSpaceRegistry::add(is): already have an IndexSpace "
            "associated with is.attr(); if you are trying to replace the "
            "IndexSpace use IndexSpaceRegistry::replace(is)");
      }
      spaces_->emplace(IS);
    }

    return clear_memoized_data_and_return_this();
  }

  /// @brief add an IndexSpace to this registry.
  /// @param type_label a label that will denote the space type,
  ///                   must be convertible to a std::string
  /// @param type an IndexSpace::Type
  /// @param args optional arguments consisting of a mix of zero or more of
  /// the following:
  ///   - IndexSpace::QuantumNumbers
  ///   - approximate size of the space (unsigned long)
  ///   - any of { is_vacuum_occupied , is_reference_occupied , is_complete ,
  ///   is_hole , is_particle }
  /// @return reference to `this`
  /// @throw std::invalid_argument if `type_label` or `type` matches
  /// an already registered IndexSpace
  template <typename S, typename... OptionalArgs,
            typename = meta::EnableIfAnyStringConvertible<S>>
  IndexSpaceRegistry& add(S&& type_label, IndexSpace::Type type,
                          OptionalArgs&&... args) {
    auto h_args = boost::hana::make_tuple(args...);

    // process IndexSpace::QuantumNumbers, set to default is not given
    auto h_qns = boost::hana::filter(h_args, [](auto arg) {
      return boost::hana::type_c<decltype(arg)> ==
             boost::hana::type_c<IndexSpace::QuantumNumbers>;
    });
    constexpr auto nqns = boost::hana::size(h_qns);
    static_assert(
        nqns == boost::hana::size_c<0> || nqns == boost::hana::size_c<1>,
        "IndexSpaceRegistry::add: only one IndexSpace::QuantumNumbers argument "
        "is allowed");
    constexpr auto have_qns = nqns == boost::hana::size_c<1>;
    IndexSpace::QuantumNumbersAttr qns;
    if constexpr (have_qns) {
      qns = boost::hana::at_c<0>(h_qns);
    }

    // process approximate_size, set to default is not given
    auto h_ints = boost::hana::filter(h_args, [](auto arg) {
      return boost::hana::traits::is_integral(boost::hana::decltype_(arg));
    });
    constexpr auto nints = boost::hana::size(h_ints);
    static_assert(
        nints == boost::hana::size_c<0> || nints == boost::hana::size_c<1>,
        "IndexSpaceRegistry::add: only one integral argument is allowed");
    constexpr auto have_approximate_size = nints == boost::hana::size_c<1>;
    unsigned long approximate_size = 10;
    if constexpr (have_approximate_size) {
      approximate_size = boost::hana::at_c<0>(h_ints);
    }

    // make space
    IndexSpace space(std::forward<S>(type_label), type, qns, approximate_size);
    this->add(space);

    // process attribute tags
    auto h_attributes = boost::hana::filter(h_args, [](auto arg) {
      return !boost::hana::traits::is_integral(
                 boost::hana::type_c<decltype(arg)>) &&
             boost::hana::type_c<decltype(arg)> !=
                 boost::hana::type_c<IndexSpace::QuantumNumbers>;
    });
    process_attribute_tags(h_attributes, type);

    return clear_memoized_data_and_return_this();
  }

  /// @brief add a union of IndexSpace objects to this registry.
  /// @param type_label a label that will denote the space type,
  ///                   must be convertible to a std::string
  /// @param components sequence of IndexSpace objects or labels (known to this)
  /// whose union will be known by @p type_label
  /// @param args optional arguments consisting of a mix of zero or more of
  /// { is_vacuum_occupied , is_reference_occupied , is_complete , is_hole ,
  /// is_particle }
  /// @return reference to `this`
  template <typename S, typename IndexSpaceOrLabel, typename... OptionalArgs,
            typename = meta::EnableIfAnyStringConvertible<S>,
            typename = std::enable_if_t<
                (std::is_same_v<std::decay_t<IndexSpaceOrLabel>, IndexSpace> ||
                 meta::is_basic_string_convertible_v<
                     std::decay_t<IndexSpaceOrLabel>>)>>
  IndexSpaceRegistry& add_unIon(
      S&& type_label, std::initializer_list<IndexSpaceOrLabel> components,
      OptionalArgs&&... args) {
    assert(components.size() > 1);

    auto h_args = boost::hana::make_tuple(args...);

    // make space
    IndexSpace::Attr space_attr;
    long count = 0;
    if (components.size() <= 1) {
      throw std::invalid_argument(
          "IndexSpaceRegistry::add_union: must have at least two components");
    }
    for (auto&& component : components) {
      const IndexSpace* component_ptr;
      if constexpr (std::is_same_v<std::decay_t<IndexSpaceOrLabel>,
                                   IndexSpace>) {
        component_ptr = &component;
      } else {
        component_ptr = &(this->retrieve(component));
      }
      if (count == 0)
        space_attr = component_ptr->attr();
      else
        space_attr = space_attr.unIon(component_ptr->attr());
      ++count;
    }
    const auto approximate_size = compute_approximate_size(space_attr);

    IndexSpace space(std::forward<S>(type_label), space_attr.type(),
                     space_attr.qns(), approximate_size);
    this->add(space);
    auto type = space.type();

    // process attribute tags
    auto h_attributes = boost::hana::filter(h_args, [](auto arg) {
      return !boost::hana::traits::is_integral(
                 boost::hana::type_c<decltype(arg)>) &&
             boost::hana::type_c<decltype(arg)> !=
                 boost::hana::type_c<IndexSpace::QuantumNumbers>;
    });
    process_attribute_tags(h_attributes, type);

    return clear_memoized_data_and_return_this();
  }

  /// alias to add_unIon
  template <typename S, typename IndexSpaceOrLabel, typename... OptionalArgs,
            typename = meta::EnableIfAnyStringConvertible<S>,
            typename = std::enable_if_t<
                (std::is_same_v<std::decay_t<IndexSpaceOrLabel>, IndexSpace> ||
                 meta::is_basic_string_convertible_v<
                     std::decay_t<IndexSpaceOrLabel>>)>>
  IndexSpaceRegistry& add_union(
      S&& type_label, std::initializer_list<IndexSpaceOrLabel> components,
      OptionalArgs&&... args) {
    return this->add_unIon(std::forward<S>(type_label), components,
                           std::forward<OptionalArgs>(args)...);
  }

  /// @brief add a union of IndexSpace objects to this registry.
  /// @param type_label a label that will denote the space type,
  ///                   must be convertible to a std::string
  /// @param components sequence of IndexSpace objects or labels (known to this)
  /// whose intersection will be known by @p type_label
  /// @param args optional arguments consisting of a mix of zero or more of
  /// { is_vacuum_occupied , is_reference_occupied , is_complete , is_hole ,
  /// is_particle }
  /// @return reference to `this`
  template <typename S, typename IndexSpaceOrLabel, typename... OptionalArgs,
            typename = meta::EnableIfAnyStringConvertible<S>,
            typename = std::enable_if_t<
                (std::is_same_v<std::decay_t<IndexSpaceOrLabel>, IndexSpace> ||
                 meta::is_basic_string_convertible_v<
                     std::decay_t<IndexSpaceOrLabel>>)>>
  IndexSpaceRegistry& add_intersection(
      S&& type_label, std::initializer_list<IndexSpaceOrLabel> components,
      OptionalArgs&&... args) {
    assert(components.size() > 1);

    auto h_args = boost::hana::make_tuple(args...);

    // make space
    IndexSpace::Attr space_attr;
    long count = 0;
    if (components.size() <= 1) {
      throw std::invalid_argument(
          "IndexSpaceRegistry::add_intersection: must have at least two "
          "components");
    }
    for (auto&& component : components) {
      const IndexSpace* component_ptr;
      if constexpr (std::is_same_v<std::decay_t<IndexSpaceOrLabel>,
                                   IndexSpace>) {
        component_ptr = &component;
      } else {
        component_ptr = &(this->retrieve(component));
      }
      if (count == 0)
        space_attr = component_ptr->attr();
      else
        space_attr = space_attr.intersection(component_ptr->attr());
      ++count;
    }
    const auto approximate_size = compute_approximate_size(space_attr);

    IndexSpace space(std::forward<S>(type_label), space_attr.type(),
                     space_attr.qns(), approximate_size);
    this->add(space);
    auto type = space.type();

    // process attribute tags
    auto h_attributes = boost::hana::filter(h_args, [](auto arg) {
      return !boost::hana::traits::is_integral(
                 boost::hana::type_c<decltype(arg)>) &&
             boost::hana::type_c<decltype(arg)> !=
                 boost::hana::type_c<IndexSpace::QuantumNumbers>;
    });
    process_attribute_tags(h_attributes, type);

    return clear_memoized_data_and_return_this();
  }

  /// @}

  /// @brief removes an IndexSpace associated with `IS.base_key()` from this
  /// @param IS an IndexSpace
  /// @return reference to `this`
  IndexSpaceRegistry& remove(const IndexSpace& IS) {
    auto it = spaces_->find(IS.base_key());
    if (it != spaces_->end()) {
      spaces_->erase(IS);
    }
    return clear_memoized_data_and_return_this();
  }

  /// @brief equivalent to `remove(this->retrieve(label))`
  /// @param label space label
  /// @return reference to `this`
  template <typename S, typename = meta::EnableIfAnyStringConvertible<S>>
  IndexSpaceRegistry& remove(S&& label) {
    auto&& IS = this->retrieve(std::forward<S>(label));
    return this->remove(IS);
  }

  /// @brief replaces an IndexSpace registered in the registry under
  /// IS.base_key()
  ///        with @p IS
  /// @param IS an IndexSpace
  /// @return reference to `this`
  IndexSpaceRegistry& replace(const IndexSpace& IS) {
    this->remove(IS);
    return this->add(IS);
  }

  /// @brief returns the list of _basis_ IndexSpace::Type objects

  /// A base IndexSpace::Type object has 1 bit in its bitstring.
  /// @sa IndexSpaceRegistry::is_base
  /// @return (memoized) set of base IndexSpace::Type objects, sorted in
  /// increasing order
  const std::vector<IndexSpace::Type>& base_space_types() const {
    if (!base_space_types_) {
      auto types = *spaces_ | ranges::views::transform([](const auto& s) {
        return s.type();
      }) | ranges::views::filter([](const auto& t) { return is_base(t); }) |
                   ranges::views::unique | ranges::to_vector;
      ranges::sort(types, [](auto t1, auto t2) { return t1 < t2; });
      base_space_types_ = std::move(types);
    }
    return *base_space_types_;
  }

  /// @brief returns the list of _basis_ IndexSpace objects

  /// A base IndexSpace object has 1 bit in its type() bitstring.
  /// @sa IndexSpaceRegistry::is_base
  /// @return (memoized) set of base IndexSpace objects, sorted in the order of
  /// increasing type()
  const std::vector<IndexSpace>& base_spaces() const {
    if (!base_spaces_) {
      auto spaces =
          *spaces_ |
          ranges::views::filter([](const auto& s) { return is_base(s); }) |
          ranges::views::unique | ranges::to_vector;
      ranges::sort(spaces,
                   [](auto s1, auto s2) { return s1.type() < s2.type(); });
      base_spaces_ = std::move(spaces);
    }
    return *base_spaces_;
  }

  /// @brief checks if an IndexSpace is in the basis
  /// @param IS IndexSpace
  /// @return true if @p IS is in the basis
  /// @sa base_spaces
  static bool is_base(const IndexSpace& IS) {
    return has_single_bit(IS.type().to_int32());
  }

  /// @brief checks if an IndexSpace::Type is in the basis
  /// @param t IndexSpace::Type
  /// @return true if @p t is in the basis
  /// @sa space_type_basis
  static bool is_base(const IndexSpace::Type& t) {
    return has_single_bit(t.to_int32());
  }

  /// @brief clear the IndexSpaceRegistry map
  /// @return reference to `this`
  IndexSpaceRegistry& clear_registry() {
    spaces_->clear();
    return this->add(IndexSpace::null);
  }

  /// @brief Is the result of a binary operation null or not registered return
  /// false.
  /// a user may wish to know if an operation returns a space they have
  /// registered.
  /// @param i1 IndexSpace
  /// @param i2 IndexSpace
  /// @param op a function which takes two int32 as arguments and returns an
  /// int32
  /// @note IndexSpaces must have the same @c QuantumNumberAttr to be a valid
  /// bitop
  bool valid_bitop(const IndexSpace& i1, const IndexSpace& i2,
                   const std::function<int32_t(int32_t, int32_t)>& op) {
    auto bitop_int = op(i1.type().to_int32(), i2.type().to_int32());
    bool same_qn = i1.qns() == i2.qns();
    if (!same_qn) return false;
    auto& temp_space = find_by_attr({bitop_int, i1.qns()});
    return temp_space == IndexSpace::null ? false : true;
  }

  /// @brief return the resulting space corresponding to a bitwise intersection
  /// between two spaces.
  /// @param space1
  /// @param space2
  /// @return the resulting space after intesection
  /// @note can return nullspace
  /// @note throw invalid_argument if the bitwise result is not registered
  const IndexSpace& intersection(const IndexSpace& space1,
                                 const IndexSpace& space2) const {
    if (space1 == space2) {
      return space1;
    } else {
      const auto target_qns = space1.qns().intersection(space2.qns());
      bool same_qns = space1.qns() == space2.qns();
      if (!target_qns && !same_qns) {  // spaces with different quantum numbers
                                       // do not intersect.
        return IndexSpace::null;
      }
      auto intersection_attr = space1.type().intersection(space2.type());
      const IndexSpace& intersection_space =
          find_by_attr({intersection_attr, space1.qns()});
      // the nullspace is a reasonable return value for intersection
      if (intersection_space == IndexSpace::null && intersection_attr) {
        throw std::invalid_argument(
            "The resulting space is not registered in this context. Add this "
            "space to the registry with a label to use it.");
      } else {
        return intersection_space;
      }
    }
  }

  /// @param space1
  /// @param space2
  /// @return the union of two spaces.
  /// @note can only return registered spaces
  /// @note never returns nullspace
  const IndexSpace& unIon(const IndexSpace& space1,
                          const IndexSpace& space2) const {
    if (space1 == space2) {
      return space1;
    } else {
      bool same_qns = space1.qns() == space2.qns();
      if (!same_qns) {
        throw std::invalid_argument(
            "asking for the intersection of spaces with incompatible quantum "
            "number attributes.");
      }
      auto unIontype = space1.type().unIon(space2.type());
      const IndexSpace& unIonSpace = find_by_attr({unIontype, space1.qns()});
      if (unIonSpace == IndexSpace::null) {
        throw std::invalid_argument(
            "The resulting space is not registered in this context. Add this "
            "space to the registry with a label to use it.");
      } else {
        return unIonSpace;
      }
    }
  }
  /// @brief which spaces result from XOR bitwise operation of two spaces. Only
  /// keep connected spaces.
  /// @param space1
  /// @param space2
  /// @return a list of spaces
  /// @note nullspace is not a valid return
  /// @note finding unregistered @c typeattr will throw
  std::vector<IndexSpace> non_overlapping_spaces(
      const IndexSpace& space1, const IndexSpace& space2) const {
    auto attributes = space1.attr().excluded_spaces(space2.attr());
    std::vector<IndexSpace> result;
    for (int i = 0; i < attributes.size(); i++) {
      auto excluded_space = find_by_attr(attributes[i]);
      if (excluded_space == IndexSpace::null) {
        throw std::invalid_argument(
            "The resulting space is not registered in this context. Add this "
            "space to the registry with a label to use it.");
      }
      result.push_back(excluded_space);
    }
    return result;
  }

  ///@brief do two spaces have non-overlapping bitsets.
  /// @note does not probe the registry for these spaces
  bool has_non_overlapping_spaces(const IndexSpace& space1,
                                  const IndexSpace& space2) const {
    return space1.type().xOr(space2.type()).to_int32() != 0;
  }

  ///@brief an @c IndexSpace is occupied with respect to the fermi vacuum or a
  /// subset of that space
  /// @note only makes sense to ask this if in a SingleProduct vacuum context.
  bool is_pure_occupied(const IndexSpace& IS) const {
    if (!IS) {
      return false;
    }
    if (IS.type().to_int32() <=
        vacuum_occupied_space(IS.qns()).type().to_int32()) {
      return true;
    } else {
      return false;
    }
  }

  ///@brief all states are unoccupied in the fermi vacuum
  ///@note again, this only makes sense to ask if in a SingleProduct vacuum
  /// context.
  bool is_pure_unoccupied(const IndexSpace& IS) const {
    if (!IS) {
      return false;
    } else {
      return !IS.type().intersection(vacuum_occupied_space(IS.qns()).type());
    }
  }

  ///@brief some states are fermi vacuum occupied
  bool contains_occupied(const IndexSpace& IS) const {
    return IS.type().intersection(vacuum_occupied_space(IS.qns()).type()) !=
           IndexSpace::Type::null;
  }

  ///@brief some states are fermi vacuum unoccupied
  bool contains_unoccupied(const IndexSpace& IS) const {
    if (IS == nullspace) {
      return false;
    } else {
      return vacuum_occupied_space(IS.qns()).type() < IS.type();
    }
  }

  /// @name  specifies which spaces have nonzero occupancy in the vacuum wave
  ///        function
  /// @note needed for applying Wick theorem with Fermi vacuum
  /// @{

  /// @param t an IndexSpace::Type specifying which base spaces have nonzero
  ///          occupancy in
  ///          the vacuum wave function by default (i.e. for any quantum number
  ///          choice); to specify occupied space per specific QN set use the
  ///          other overload
  /// @return reference to `this`
  IndexSpaceRegistry& vacuum_occupied_space(const IndexSpace::Type& t) {
    throw_if_missing(t, "vacuum_occupied_space");
    std::get<0>(vacocc_) = t;
    return *this;
  }

  /// @param qn2type for each quantum number specifies which base spaces have
  ///                nonzero occupancy in the reference wave function
  /// @return reference to `this`
  IndexSpaceRegistry& vacuum_occupied_space(
      std::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
    throw_if_missing_any(qn2type, "vacuum_occupied_space");
    std::get<1>(vacocc_) = std::move(qn2type);
    return *this;
  }

  /// equivalent to `vacuum_occupied_space(s.type())`
  /// @note QuantumNumbers attribute of `s` ignored
  /// @param s an IndexSpace
  /// @return reference to `this`
  IndexSpaceRegistry& vacuum_occupied_space(const IndexSpace& s) {
    return vacuum_occupied_space(s.type());
  }

  /// equivalent to `vacuum_occupied_space(retrieve(l).type())`
  /// @param l label of a known IndexSpace
  /// @return reference to `this`
  template <typename S, typename = meta::EnableIfAnyStringConvertible<S>>
  IndexSpaceRegistry& vacuum_occupied_space(S&& l) {
    return vacuum_occupied_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return the space occupied in vacuum state for any set of quantum numbers
  /// @throw std::invalid_argument if @p nulltype_ok is false and
  /// vacuum_occupied_space had not been specified
  const IndexSpace::Type& vacuum_occupied_space(
      bool nulltype_ok = false) const {
    if (!std::get<0>(vacocc_)) {
      if (nulltype_ok) return IndexSpace::Type::null;
      throw std::invalid_argument(
          "vacuum occupied space has not been specified, invoke "
          "vacuum_occupied_space(IndexSpace::Type) or "
          "vacuum_occupied_space(std::map<IndexSpace::QuantumNumbers,"
          "IndexSpace::Type>)");
    } else
      return std::get<0>(vacocc_);
  }

  /// @param qn the quantum numbers of the space
  /// @return the space occupied in vacuum state for the given set of quantum
  /// numbers
  const IndexSpace& vacuum_occupied_space(
      const IndexSpace::QuantumNumbers& qn) const {
    auto it = std::get<1>(vacocc_).find(qn);
    if (it != std::get<1>(vacocc_).end()) {
      return retrieve(it->second, qn);
    } else {
      return retrieve(this->vacuum_occupied_space(), qn);
    }
  }

  /// @}

  /// @name  assign which spaces have nonzero occupancy in the reference wave
  ///        function (i.e., the wave function uses to compute reference
  ///        expectation value)
  /// @note needed for computing expectation values when the vacuum state does
  /// not match the wave function of interest.
  /// @{

  /// @param t an IndexSpace::Type specifying which base spaces have nonzero
  /// occupancy in
  ///          the reference wave function by default (i.e., for any choice of
  ///          quantum numbers); to specify occupied space per specific QN set
  ///          use the other overload
  /// @return reference to `this`
  IndexSpaceRegistry& reference_occupied_space(const IndexSpace::Type& t) {
    throw_if_missing(t, "reference_occupied_space");
    std::get<0>(refocc_) = t;
    return *this;
  }

  /// @param qn2type for each quantum number specifies which base spaces have
  /// nonzero occupancy in
  ///          the reference wave function
  /// @return reference to `this`
  IndexSpaceRegistry& reference_occupied_space(
      std::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
    throw_if_missing_any(qn2type, "reference_occupied_space");
    std::get<1>(refocc_) = std::move(qn2type);
    return *this;
  }

  /// equivalent to `reference_occupied_space(s.type())`
  /// @note QuantumNumbers attribute of `s` ignored
  /// @param s an IndexSpace
  /// @return reference to `this`
  IndexSpaceRegistry& reference_occupied_space(const IndexSpace& s) {
    return reference_occupied_space(s.type());
  }

  /// equivalent to `reference_occupied_space(retrieve(l).type())`
  /// @param l label of a known IndexSpace
  /// @return reference to `this`
  template <typename S, typename = meta::EnableIfAnyStringConvertible<S>>
  IndexSpaceRegistry& reference_occupied_space(S&& l) {
    return reference_occupied_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return the space occupied in reference state for any set of quantum
  /// numbers
  /// @throw std::invalid_argument if @p nulltype_ok is false and
  /// reference_occupied_space had not been specified
  const IndexSpace::Type& reference_occupied_space(
      bool nulltype_ok = false) const {
    if (!std::get<0>(refocc_)) {
      if (nulltype_ok) return IndexSpace::Type::null;
      throw std::invalid_argument(
          "reference occupied space has not been specified, invoke "
          "reference_occupied_space(IndexSpace::Type) or "
          "reference_occupied_space(std::map<IndexSpace::QuantumNumbers,"
          "IndexSpace::Type>)");
    } else
      return std::get<0>(refocc_);
  }

  /// @param qn the quantum numbers of the space
  /// @return the space occupied in vacuum state for the given set of quantum
  /// numbers
  const IndexSpace& reference_occupied_space(
      const IndexSpace::QuantumNumbers& qn) const {
    auto it = std::get<1>(refocc_).find(qn);
    if (it != std::get<1>(refocc_).end()) {
      return retrieve(it->second, qn);
    } else {
      return retrieve(this->reference_occupied_space(), qn);
    }
  }

  /// @}

  /// @name  specifies which spaces comprise the entirety of Hilbert space
  /// @note needed for creating general operators in mbpt/op
  /// @{

  /// @param t an IndexSpace::Type specifying the complete Hilbert space;
  ///          to specify occupied space per specific QN set use the other
  ///          overload
  IndexSpaceRegistry& complete_space(const IndexSpace::Type& s) {
    throw_if_missing(s, "complete_space");
    std::get<0>(complete_) = s;
    return *this;
  }

  /// @param qn2type for each quantum number specifies which base spaces have
  /// nonzero occupancy in
  ///          the reference wave function
  IndexSpaceRegistry& complete_space(
      std::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
    throw_if_missing_any(qn2type, "complete_space");
    std::get<1>(complete_) = std::move(qn2type);
    return *this;
  }

  /// equivalent to `complete_space(s.type())`
  /// @note QuantumNumbers attribute of `s` ignored
  /// @param s an IndexSpace
  /// @return reference to `this`
  IndexSpaceRegistry& complete_space(const IndexSpace& s) {
    return complete_space(s.type());
  }

  /// equivalent to `complete_space(retrieve(l).type())`
  /// @param l label of a known IndexSpace
  /// @return reference to `this`
  template <typename S, typename = meta::EnableIfAnyStringConvertible<S>>
  IndexSpaceRegistry& complete_space(S&& l) {
    return complete_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return the complete Hilbert space for any set of quantum numbers
  /// @throw std::invalid_argument if @p nulltype_ok is false and complete_space
  /// had not been specified
  const IndexSpace::Type& complete_space(bool nulltype_ok = false) const {
    if (!std::get<0>(complete_)) {
      if (nulltype_ok) return IndexSpace::Type::null;
      throw std::invalid_argument(
          "complete space has not been specified, call "
          "complete_space(IndexSpace::Type)");
    } else
      return std::get<0>(complete_);
  }

  /// @param qn the quantum numbers of the space
  /// @return the complete Hilbert space for the given set of quantum numbers
  const IndexSpace& complete_space(const IndexSpace::QuantumNumbers& qn) const {
    auto it = std::get<1>(complete_).find(qn);
    if (it != std::get<1>(complete_).end()) {
      return retrieve(it->second, qn);
    } else {
      return retrieve(this->complete_space(), qn);
    }
  }

  /// @}

  /// @name specifies in which space holes can be created successfully from the
  /// reference wave function
  /// @note convenience for making operators
  /// @{

  /// @param t an IndexSpace::Type specifying where holes can be created;
  ///          to specify hole space per specific QN set use the other
  ///          overload
  IndexSpaceRegistry& hole_space(const IndexSpace::Type& t) {
    throw_if_missing(t, "hole_space");
    std::get<0>(hole_space_) = t;
    return *this;
  }

  /// @param qn2type for each quantum number specifies the space in which holes
  /// can be created
  IndexSpaceRegistry& hole_space(
      std::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
    throw_if_missing_any(qn2type, "hole_space");
    std::get<1>(hole_space_) = std::move(qn2type);
    return *this;
  }

  /// equivalent to `hole_space(s.type())`
  /// @note QuantumNumbers attribute of `s` ignored
  /// @param s an IndexSpace
  /// @return reference to `this`
  IndexSpaceRegistry& hole_space(const IndexSpace& s) {
    return hole_space(s.type());
  }

  /// equivalent to `hole_space(retrieve(l).type())`
  /// @param l label of a known IndexSpace
  /// @return reference to `this`
  template <typename S, typename = meta::EnableIfAnyStringConvertible<S>>
  IndexSpaceRegistry& hole_space(S&& l) {
    return hole_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return default space in which holes can be created
  /// @throw std::invalid_argument if @p nulltype_ok is false and
  /// hole_space had not been specified
  const IndexSpace::Type& hole_space(bool nulltype_ok = false) const {
    if (!std::get<0>(hole_space_)) {
      if (nulltype_ok) return IndexSpace::Type::null;
      throw std::invalid_argument(
          "active hole space has not been specified, invoke "
          "hole_space(IndexSpace::Type) or "
          "hole_space(std::map<IndexSpace::QuantumNumbers,IndexSpace::"
          "Type>)");
    } else
      return std::get<0>(hole_space_);
  }

  /// @param qn the quantum numbers of the space
  /// @return the space in which holes can be created for the given set of
  /// quantum numbers
  const IndexSpace& hole_space(const IndexSpace::QuantumNumbers& qn) const {
    auto it = std::get<1>(hole_space_).find(qn);
    if (it != std::get<1>(hole_space_).end()) {
      return this->retrieve(it->second, qn);
    } else {
      return this->retrieve(this->hole_space(), qn);
    }
  }

  /// @}

  /// @name specifies in which space particles can be created successfully from
  /// the reference wave function
  /// @note convenience for making operators
  /// @{

  /// @param t an IndexSpace::Type specifying where particles can be created;
  ///          to specify particle space per specific QN set use the other
  ///          overload
  IndexSpaceRegistry& particle_space(const IndexSpace::Type& t) {
    throw_if_missing(t, "particle_space");
    std::get<0>(particle_space_) = t;
    return *this;
  }

  /// @param qn2type for each quantum number specifies the space in which
  /// particles can be created
  IndexSpaceRegistry& particle_space(
      std::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type) {
    throw_if_missing_any(qn2type, "particle_space");
    std::get<1>(particle_space_) = std::move(qn2type);
    return *this;
  }

  /// equivalent to `particle_space(s.type())`
  /// @note QuantumNumbers attribute of `s` ignored
  /// @param s an IndexSpace
  /// @return reference to `this`
  IndexSpaceRegistry& particle_space(const IndexSpace& s) {
    return particle_space(s.type());
  }

  /// equivalent to `particle_space(retrieve(l).type())`
  /// @param l label of a known IndexSpace
  /// @return reference to `this`
  template <typename S, typename = meta::EnableIfAnyStringConvertible<S>>
  IndexSpaceRegistry& particle_space(S&& l) {
    return particle_space(this->retrieve(std::forward<S>(l)).type());
  }

  /// @return default space in which particles can be created
  /// @throw std::invalid_argument if @p nulltype_ok is false and
  /// particle_space had not been specified
  const IndexSpace::Type& particle_space(bool nulltype_ok = false) const {
    if (!std::get<0>(particle_space_)) {
      if (nulltype_ok) return IndexSpace::Type::null;
      throw std::invalid_argument(
          "active particle space has not been specified, invoke "
          "particle_space(IndexSpace::Type) or "
          "particle_space(std::map<IndexSpace::QuantumNumbers,"
          "IndexSpace::Type>)");
    } else
      return std::get<0>(particle_space_);
  }

  /// @param qn the quantum numbers of the space
  /// @return the space in which particles can be created for the given set of
  /// quantum numbers
  const IndexSpace& particle_space(const IndexSpace::QuantumNumbers& qn) const {
    auto it = std::get<1>(particle_space_).find(qn);
    if (it != std::get<1>(particle_space_).end()) {
      return this->retrieve(it->second, qn);
    } else {
      return this->retrieve(this->particle_space(), qn);
    }
  }

  /// @}

 private:
  // N.B. need transparent comparator, see https://stackoverflow.com/a/35525806
  std::shared_ptr<container::set<IndexSpace, IndexSpace::KeyCompare>> spaces_;

  // memoized data
  mutable std::optional<std::vector<IndexSpace::Type>> base_space_types_;
  mutable std::optional<std::vector<IndexSpace>> base_spaces_;
  IndexSpaceRegistry& clear_memoized_data_and_return_this() {
    base_space_types_.reset();
    base_spaces_.reset();
    return *this;
  }

  /// TODO use c++20 std::has_single_bit() when we update to this version
  static bool has_single_bit(std::uint32_t bits) {
    return bits & (((bool)(bits & (bits - 1))) - 1);
  }

  ///@brief find an IndexSpace from its type. return nullspace if not present.
  ///@param find the IndexSpace via it's @c attr
  const IndexSpace& find_by_attr(const IndexSpace::Attr& attr) const {
    for (auto&& space : *spaces_) {
      if (space.attr() == attr) {
        return space;
      }
    }
    return IndexSpace::null;
  }

  void throw_if_missing(const IndexSpace::Type& t,
                        const IndexSpace::QuantumNumbers& qn,
                        std::string call_context = "") {
    for (auto&& space : *spaces_) {
      if (space.type() == t && space.qns() == qn) {
        return;
      }
    }
    throw std::invalid_argument(
        call_context +
        ": missing { IndexSpace::Type=" + std::to_string(t.to_int32()) +
        " , IndexSpace::QuantumNumbers=" + std::to_string(qn.to_int32()) +
        " } combination");
  }

  // same as above, but ignoring qn
  void throw_if_missing(const IndexSpace::Type& t,
                        std::string call_context = "") {
    for (auto&& space : *spaces_) {
      if (space.type() == t) {
        return;
      }
    }
    throw std::invalid_argument(call_context + ": missing { IndexSpace::Type=" +
                                std::to_string(t.to_int32()) +
                                " , any IndexSpace::QuantumNumbers } space");
  }

  void throw_if_missing_any(
      const std::map<IndexSpace::QuantumNumbers, IndexSpace::Type>& qn2type,
      std::string call_context = "") {
    std::map<IndexSpace::QuantumNumbers, IndexSpace::Type> qn2type_found;
    for (auto&& space : *spaces_) {
      for (auto&& [qn, t] : qn2type) {
        if (space.type() == t && space.qns() == qn) {
          auto [it, found] = qn2type_found.try_emplace(qn, t);
          assert(!found);
          // found all? return
          if (qn2type_found.size() == qn2type.size()) {
            return;
          }
        }
      }
    }

    std::string errmsg;
    for (auto&& [qn, t] : qn2type) {
      if (!qn2type_found.contains(qn)) {
        errmsg +=
            call_context +
            ": missing { IndexSpace::Type=" + std::to_string(t.to_int32()) +
            " , IndexSpace::QuantumNumbers=" + std::to_string(qn.to_int32()) +
            " } combination\n";
      }
    }
    throw std::invalid_argument(errmsg);
  }

  // Need to define defaults for various traits, like which spaces are occupied
  // in vacuum, etc. Makes sense to make these part of the registry to avoid
  // having to pass these around in every call N.B. default and QN-specific
  // space selections merged into single tuple

  // defines active bits in TypeAttr; used by general operators in mbpt/op
  std::tuple<IndexSpace::Type,
             std::map<IndexSpace::QuantumNumbers, IndexSpace::Type>>
      complete_ = {{}, {}};

  // used for fermi vacuum wick application
  std::tuple<IndexSpace::Type,
             std::map<IndexSpace::QuantumNumbers, IndexSpace::Type>>
      vacocc_ = {{}, {}};

  // used for MR MBPT to take average over multiconfiguration reference
  std::tuple<IndexSpace::Type,
             std::map<IndexSpace::QuantumNumbers, IndexSpace::Type>>
      refocc_ = {{}, {}};

  // both needed to make excitation and de-excitation operators. not
  // necessarily equivalent in the case of multi-reference context.
  std::tuple<IndexSpace::Type,
             std::map<IndexSpace::QuantumNumbers, IndexSpace::Type>>
      particle_space_ = {{}, {}};
  std::tuple<IndexSpace::Type,
             std::map<IndexSpace::QuantumNumbers, IndexSpace::Type>>
      hole_space_ = {{}, {}};

  // Boost.Hana snippet to process attribute tag arguments
  template <typename ArgsHanaTuple>
  void process_attribute_tags(ArgsHanaTuple h_tuple,
                              const IndexSpace::Type& type) {
    boost::hana::for_each(h_tuple, [this, &type](auto arg) {
      if constexpr (boost::hana::type_c<decltype(arg)> ==
                    boost::hana::type_c<space_tags::IsVacuumOccupied>) {
        this->vacuum_occupied_space(type);
      } else if constexpr (boost::hana::type_c<decltype(arg)> ==
                           boost::hana::type_c<
                               space_tags::IsReferenceOccupied>) {
        this->reference_occupied_space(type);
      } else if constexpr (boost::hana::type_c<decltype(arg)> ==
                           boost::hana::type_c<space_tags::IsComplete>) {
        this->complete_space(type);
      } else if constexpr (boost::hana::type_c<decltype(arg)> ==
                           boost::hana::type_c<space_tags::IsHole>) {
        this->hole_space(type);
      } else if constexpr (boost::hana::type_c<decltype(arg)> ==
                           boost::hana::type_c<space_tags::IsParticle>) {
        this->particle_space(type);
      } else {
        static_assert(meta::always_false<decltype(arg)>::value,
                      "IndexSpaceRegistry::add{,_union,_intersect}: unknown "
                      "attribute tag");
      }
    });
  }

  /// @brief computes the approximate size of the space

  /// for a base space return its extent, for a composite space compute as a sum
  /// of extents of base subspaces
  /// @param s an IndexSpace object
  /// @return the approximate size of the space
  unsigned long compute_approximate_size(
      const IndexSpace::Attr& space_attr) const {
    if (is_base(space_attr.type())) {
      return this->retrieve(space_attr).approximate_size();
    } else {
      // compute_approximate_size is used when populating the registry
      // so don't use base_spaces() here
      unsigned long size = ranges::accumulate(
          *spaces_ | ranges::views::filter([&space_attr](auto& s) {
            return s.qns() == space_attr.qns() && is_base(s) &&
                   space_attr.type().intersection(s.type());
          }),
          0ul, [](unsigned long size, const IndexSpace& s) {
            return size + s.approximate_size();
          });
      return size;
    }
  }

  friend bool operator==(const IndexSpaceRegistry& isr1,
                         const IndexSpaceRegistry& isr2) {
    return *isr1.spaces_ == *isr2.spaces_;
  }
};

}  // namespace sequant
#endif  // SEQUANT_INDEX_SPACE_REGISTRY_HPP
