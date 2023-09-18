#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_SHUFFLE_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_SHUFFLE_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/NoneType.hpp"
#include "pythonic/include/numpy/random/generator.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace random
  {
    template <class T>
    types::none_type shuffle(T &seq);

    DEFINE_FUNCTOR(pythonic::numpy::random, shuffle);
  }
}

PYTHONIC_NS_END

#endif
