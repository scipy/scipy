#ifndef PYTHONIC_INCLUDE_NUMPY_FLOORDIVIDE_HPP
#define PYTHONIC_INCLUDE_NUMPY_FLOORDIVIDE_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/numpy_broadcast.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include//numpy/floor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  namespace wrapper
  {
    template <class Arg0, class Arg1>
    std::complex<typename std::common_type<Arg0, Arg1>::type>
    divfloor(std::complex<Arg0> const &arg0, std::complex<Arg1> const &arg1)
    {
      return {functor::floor{}(std::real(arg0 / arg1)), 0};
    }

    template <class Arg0, class Arg1>
    auto divfloor(Arg0 const &arg0, Arg1 const &arg1) ->
        typename std::enable_if<(std::is_integral<Arg0>::value &&
                                 std::is_integral<Arg1>::value),
                                decltype(arg0 / arg1)>::type
    {
      bool opposite_sign = (arg0 >= 0 && arg1 < 0) || (arg0 < 0 && arg1 >= 0);
      return (arg0 + opposite_sign * (-arg1 + 1)) / arg1;
    }

    template <class Arg0, class Arg1>
    auto divfloor(Arg0 const &arg0, Arg1 const &arg1) ->
        typename std::enable_if<!std::is_integral<Arg0>::value ||
                                    !std::is_integral<Arg1>::value,
                                decltype(functor::floor{}(arg0 / arg1))>::type
    {
      return functor::floor{}(arg0 / arg1);
    }
  }
#define NUMPY_NARY_FUNC_NAME floor_divide
#define NUMPY_NARY_FUNC_SYM wrapper::divfloor
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif
