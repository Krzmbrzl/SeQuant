#ifndef SEQUANT_EXTERNAL_BOOST_INTERVAL_H
#define SEQUANT_EXTERNAL_BOOST_INTERVAL_H

// boost/numeric/interval does not know about arm rounding .. on arm64/macos use
// c99 rounding
#if defined(__arm64__) && defined(__APPLE__) && !defined(__USE_ISOC99)
#define __USE_ISOC99 1
#include <boost/numeric/interval.hpp>
#undef __USE_ISOC99
#else
#include <boost/numeric/interval.hpp>
#endif

namespace boost::numeric {

template <typename T>
inline auto hash_value(const interval<T>& i) {
  using boost::hash_value;
  auto val = hash_value(i.lower());
  using boost::hash_combine;
  hash_combine(val, i.upper());
  return val;
}

}  // namespace boost::numeric

#endif  // SEQUANT_EXTERNAL_BOOST_INTERVAL_H
