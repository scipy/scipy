#ifndef PYTHONIC_NUMPY_FROMFUNCTION_HPP
#define PYTHONIC_NUMPY_FROMFUNCTION_HPP

#include "pythonic/include/numpy/fromfunction.hpp"

#include "pythonic/utils/functor.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/builtins/None.hpp"
#include "pythonic/utils/tags.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class F, size_t N, class dtype, class Tags>
  struct fromfunction_helper;

  template <class F, class dtype, class purity_tag>
  template <class pS>
  types::ndarray<typename std::remove_cv<typename std::remove_reference<
                     typename std::result_of<F(dtype)>::type>::type>::type,
                 pS> fromfunction_helper<F, 1, dtype, purity_tag>::
  operator()(F &&f, pS const &shape, dtype d)
  {
    types::ndarray<typename std::remove_cv<typename std::remove_reference<
                       typename std::result_of<F(dtype)>::type>::type>::type,
                   pS> out(shape, builtins::None);
    long n = out.template shape<0>();
    for (long i = 0; i < n; ++i)
      out[i] = f(i);
    return out;
  }

  template <class F, class dtype, class purity_tag>
  template <class pS>
  types::ndarray<
      typename std::remove_cv<typename std::remove_reference<
          typename std::result_of<F(dtype, dtype)>::type>::type>::type,
      pS> fromfunction_helper<F, 2, dtype, purity_tag>::
  operator()(F &&f, pS const &shape, dtype d)
  {
    types::ndarray<
        typename std::remove_cv<typename std::remove_reference<
            typename std::result_of<F(dtype, dtype)>::type>::type>::type,
        pS> out(shape, builtins::None);
    long n = out.template shape<0>();
    long m = out.template shape<1>();
    for (long i = 0; i < n; ++i)
      for (long j = 0; j < m; ++j)
        out[i][j] = f(i, j);
    return out;
  }

  template <class F, class pS, class dtype>
  auto fromfunction(F &&f, pS const &shape, dtype d)
      -> decltype(fromfunction_helper<F, std::tuple_size<pS>::value, dtype,
                                      typename pythonic::purity_of<F>::type>()(
          std::forward<F>(f), shape))
  {
    return fromfunction_helper<F, std::tuple_size<pS>::value, dtype,
                               typename pythonic::purity_of<F>::type>()(
        std::forward<F>(f), shape);
  }

  /* TODO: must specialize for higher order */
}
PYTHONIC_NS_END

#endif
