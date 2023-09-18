#ifndef PYTHONIC_INCLUDE_UTILS_TAGS_HPP
#define PYTHONIC_INCLUDE_UTILS_TAGS_HPP

#include "pythonic/include/types/traits.hpp"

PYTHONIC_NS_BEGIN

namespace purity
{
  struct unknown_tag {
  };

  struct pure_tag {
  };
}

template <class T>
struct purity_of {
  using type =
      typename std::conditional<types::is_pure<T>::value, purity::pure_tag,
                                purity::unknown_tag>::type;
};
PYTHONIC_NS_END

#endif
