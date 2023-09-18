#ifndef PYTHONIC_INCLUDE_NUMPY_FROMFUNCTION_HPP
#define PYTHONIC_INCLUDE_NUMPY_FROMFUNCTION_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/builtins/None.hpp"
#include "pythonic/include/utils/tags.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class F, size_t N, class dtype, class Tags>
  struct fromfunction_helper;

  template <class F, class dtype, class purity_tag>
  struct fromfunction_helper<F, 1, dtype, purity_tag> {
    template <class pS>
    types::ndarray<typename std::remove_cv<typename std::remove_reference<
                       typename std::result_of<F(dtype)>::type>::type>::type,
                   pS>
    operator()(F &&f, pS const &shape, dtype d = dtype());
  };

  template <class F, class dtype, class purity_tag>
  struct fromfunction_helper<F, 2, dtype, purity_tag> {
    template <class pS>
    types::ndarray<
        typename std::remove_cv<typename std::remove_reference<
            typename std::result_of<F(dtype, dtype)>::type>::type>::type,
        pS>
    operator()(F &&f, pS const &shape, dtype d = dtype());
  };

  template <class F, class pS, class dtype = double>
  auto fromfunction(F &&f, pS const &shape, dtype d = dtype())
      -> decltype(fromfunction_helper<F, std::tuple_size<pS>::value, dtype,
                                      typename pythonic::purity_of<F>::type>()(
          std::forward<F>(f), shape));

  /* TODO: must specialize for higher order */
  DEFINE_FUNCTOR(pythonic::numpy, fromfunction);
}
PYTHONIC_NS_END

#endif
