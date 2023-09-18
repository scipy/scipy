#ifndef PYTHONIC_INCLUDE_NUMPY_ISNAN_HPP
#define PYTHONIC_INCLUDE_NUMPY_ISNAN_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class T>
    bool isnan(std::complex<T> const &v);
    template <class T>
    auto isnan(T const &v) -> typename std::enable_if<
        std::is_floating_point<typename std::decay<T>::type>::value,
        bool>::type;
    template <class T>
    auto isnan(T const &v) -> typename std::enable_if<
        !std::is_floating_point<typename std::decay<T>::type>::value,
        bool>::type;
  }

#define NUMPY_NARY_FUNC_NAME isnan
#define NUMPY_NARY_FUNC_SYM wrapper::isnan
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif
