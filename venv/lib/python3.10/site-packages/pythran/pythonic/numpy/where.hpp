#ifndef PYTHONIC_NUMPY_WHERE_HPP
#define PYTHONIC_NUMPY_WHERE_HPP

#include "pythonic/include/numpy/where.hpp"

#include "pythonic/numpy/asarray.hpp"
#include "pythonic/numpy/nonzero.hpp"
#include "pythonic/numpy/copy.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace impl
  {
    template <class E, class F, class G>
    typename __combined<F, G>::type where(E const &cond, F const &true_,
                                          G const &false_)
    {
      if (cond)
        return true_;
      else
        return false_;
    }
  }

#define NUMPY_NARY_FUNC_NAME where
#define NUMPY_NARY_FUNC_SYM impl::where
#define NUMPY_NARY_RESHAPE_MODE reshape_type
#include "pythonic/types/numpy_nary_expr.hpp"
}

namespace types
{

  template <>
  struct Dereferencer<numpy::functor::where> {

    template <class Ts>
    auto operator()(Ts const &iters, utils::index_sequence<0, 1, 2>) ->
        typename std::enable_if<
            types::is_dtype<
                typename std::remove_cv<typename std::remove_reference<
                    decltype(*std::get<0>(iters))>::type>::type>::value &&
                types::is_dtype<
                    typename std::remove_cv<typename std::remove_reference<
                        decltype(*std::get<1>(iters))>::type>::type>::value &&
                types::is_dtype<
                    typename std::remove_cv<typename std::remove_reference<
                        decltype(*std::get<2>(iters))>::type>::type>::value,
            decltype(numpy::impl::where(*std::get<0>(iters),
                                        *std::get<1>(iters),
                                        *std::get<2>(iters)))>::type
    {
      if (*std::get<0>(iters))
        return *std::get<1>(iters);
      else
        return *std::get<2>(iters);
    }

    template <class Ts, size_t... I>
    auto operator()(Ts const &iters, utils::index_sequence<I...>, ...)
        -> decltype(numpy::functor::where{}(*std::get<I>(iters)...))
    {
      return numpy::functor::where{}(*std::get<I>(iters)...);
    }
  };
}
PYTHONIC_NS_END

#endif
