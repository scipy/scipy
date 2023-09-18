#ifndef PYTHONIC_NUMPY_NDARRAY_FILL_HPP
#define PYTHONIC_NUMPY_NDARRAY_FILL_HPP

#include "pythonic/include/numpy/ndarray/fill.hpp"

#include "pythonic/utils/functor.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/builtins/None.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace ndarray
  {
    template <class E, class F>
    types::none_type fill(E &&e, F f)
    {
      std::fill(e.begin(), e.end(), f);
      return builtins::None;
    }

    template <class T, class pS, class F>
    types::none_type fill(types::ndarray<T, pS> &e, F f)
    {
      std::fill(e.fbegin(), e.fend(), f);
      return builtins::None;
    }
  }
}
PYTHONIC_NS_END

#endif
