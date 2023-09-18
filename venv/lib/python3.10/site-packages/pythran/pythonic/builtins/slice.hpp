#ifndef PYTHONIC_BUILTIN_SLICE_HPP
#define PYTHONIC_BUILTIN_SLICE_HPP

#include "pythonic/include/builtins/slice.hpp"
#include "pythonic/types/slice.hpp"

PYTHONIC_NS_BEGIN

namespace builtins
{

  namespace anonymous
  {
    inline types::cstride_slice<1> slice(types::none<long> stop)
    {
      return {types::none<long>(), stop};
    }
    inline types::cstride_slice<1> slice(types::none<long> start,
                                         types::none<long> stop)
    {
      return {start, stop};
    }
    inline types::slice slice(types::none<long> start, types::none<long> stop,
                              types::none<long> step)
    {
      return {start, stop, step};
    }
  } // namespace anonymous
} // namespace builtins
PYTHONIC_NS_END

#endif
