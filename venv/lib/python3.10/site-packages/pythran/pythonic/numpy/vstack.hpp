#ifndef PYTHONIC_NUMPY_VSTACK_HPP
#define PYTHONIC_NUMPY_VSTACK_HPP

#include <pythonic/include/numpy/vstack.hpp>
#include <pythonic/numpy/concatenate.hpp>

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class ArraySequence>
  auto vstack(ArraySequence &&seq) ->
      typename std::enable_if<(impl::vstack_helper<ArraySequence>::value > 1),
                              impl::vstack_helper<ArraySequence>>::type
  {

    return concatenate(std::forward<ArraySequence>(seq), 0);
  }

  template <class ArraySequence>
  auto vstack(ArraySequence &&seq) -> typename std::enable_if<
      (impl::vstack_helper<ArraySequence>::value == 1),
      decltype(std::declval<impl::vstack_helper<ArraySequence>>().reshape(
          std::declval<types::array<long, 2>>()))>::type
  {
    auto &&temp = concatenate(std::forward<ArraySequence>(seq), 0);
    long const seq_size = seq.size(), temp_size = temp.size();
    types::array<long, 2> new_shape{{seq_size, temp_size / seq_size}};
    return temp.reshape(new_shape);
  }
}
PYTHONIC_NS_END

#endif
