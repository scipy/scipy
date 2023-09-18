#ifndef PYTHONIC_INCLUDE_UTILS_INT_HPP
#define PYTHONIC_INCLUDE_UTILS_INT_HPP

PYTHONIC_NS_BEGIN

namespace utils
{
  template <size_t>
  struct int_ {
  }; // compile-time counter
}
PYTHONIC_NS_END

#endif
