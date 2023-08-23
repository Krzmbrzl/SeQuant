//
// Created by Robert Adam on 8/23/23.
//

#ifndef SEQUANT_TYPETRAITS_HPP
#define SEQUANT_TYPETRAITS_HPP

#include <type_traits>

namespace sequant {

// Make remove_cvref available also in pre-C++20 code
#if __cplusplus >= 202002L
template <typename T>
using remove_cvref = std::remove_cvref<T>;

using std::remove_cvref_t;
#else
template <class T>
struct remove_cvref {
  using type = std::remove_cv_t<::std::remove_reference_t<T>>;
};

template <typename T>
using remove_cvref_t = typename remove_cvref<T>::type;
#endif

/// Checks whether \c T is a \c Base (is either the same class or a sub-class
/// ignoring CV and reference qualifiers
template <typename Base, typename T>
using is_a = std::is_base_of<remove_cvref_t<Base>, remove_cvref_t<T>>;
template <typename Base, typename T>
constexpr bool is_a_v = is_a<Base, T>::value;

/// Checks whether \c T and \c U are the same type, ignoring any CV and
/// reference qualifiers
template <typename T, typename U>
struct is
    : std::bool_constant<std::is_same_v<remove_cvref_t<T>, remove_cvref_t<U>>> {
};

template <typename T, typename U>
constexpr bool is_v = is<T, U>::value;

}  // namespace sequant

#endif  // SEQUANT_TYPETRAITS_HPP
