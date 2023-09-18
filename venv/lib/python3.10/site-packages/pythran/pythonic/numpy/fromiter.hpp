#ifndef PYTHONIC_NUMPY_FROMITER_HPP
#define PYTHONIC_NUMPY_FROMITER_HPP

#include "pythonic/include/numpy/fromiter.hpp"

#include "pythonic/utils/functor.hpp"
#include "pythonic/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class Iterable, class dtype>
  types::ndarray<typename std::remove_cv<typename std::remove_reference<
                     Iterable>::type>::type::value_type,
                 types::pshape<long>>
  fromiter(Iterable &&iterable, dtype d, long count)
  {
    using T = typename std::remove_cv<
        typename std::remove_reference<Iterable>::type>::type::value_type;
    if (count < 0) {
      types::list<T> buffer(0);
      std::copy(iterable.begin(), iterable.end(), std::back_inserter(buffer));
      return {buffer};
    } else {
      utils::shared_ref<types::raw_array<T>> buffer(count);
      std::copy_n(iterable.begin(), count, buffer->data);
      types::array<long, 1> shape = {count};
      return {buffer, shape};
    }
  }
}
PYTHONIC_NS_END

#endif
