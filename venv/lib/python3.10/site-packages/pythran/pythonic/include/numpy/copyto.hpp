#ifndef PYTHONIC_INCLUDE_NUMPY_COPYTO_HPP
#define PYTHONIC_INCLUDE_NUMPY_COPYTO_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  template <class T, class pS, class E>
  types::none_type copyto(types::ndarray<T, pS> &out, E const &expr);

  template <class T, class pS, class E>
  types::none_type copyto(types::ndarray<T, pS> &&out, E const &expr);

  template <class T, class pS, class E>
  types::none_type copyto(types::numpy_texpr<types::ndarray<T, pS>> &out, E const &expr);

  template <class T, class pS, class E>
  types::none_type copyto(types::numpy_texpr<types::ndarray<T, pS>> &&out, E const &expr);

  // pythran extensions
  template <class E, class F>
  types::none_type copyto(E &out, F const &expr) {
    out[types::fast_contiguous_slice(0, types::none_type{})] = expr;
    return {};
  }

  DEFINE_FUNCTOR(pythonic::numpy, copyto);
}
PYTHONIC_NS_END

#endif
