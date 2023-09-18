#ifndef PYTHONIC_BUILTIN_FLOAT_HPP
#define PYTHONIC_BUILTIN_FLOAT_HPP

#include "pythonic/include/builtins/float_.hpp"

#include "pythonic/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace functor
  {
    template <class T>
    float_::type float_::operator()(T &&t) const
    {
      return static_cast<float_::type>(t);
    }

    inline float_::type float_::operator()() const
    {
      return 0.;
    }
  } // namespace functor
} // namespace builtins
PYTHONIC_NS_END

#endif
