#ifndef PYTHONIC_INCLUDE_NUMPY_LOGICALNOT_HPP
#define PYTHONIC_INCLUDE_NUMPY_LOGICALNOT_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/numpy_traits.hpp"
#include "pythonic/include/operator_/not_.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

#define NUMPY_NARY_FUNC_NAME logical_not
#define NUMPY_NARY_FUNC_SYM pythonic::operator_::not_
#include "pythonic/include/types/numpy_nary_expr.hpp"
}
PYTHONIC_NS_END

#endif
