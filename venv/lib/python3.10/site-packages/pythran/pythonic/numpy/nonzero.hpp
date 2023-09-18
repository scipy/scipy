#ifndef PYTHONIC_NUMPY_NONZERO_HPP
#define PYTHONIC_NUMPY_NONZERO_HPP

#include "pythonic/include/numpy/nonzero.hpp"

#include "pythonic/utils/functor.hpp"
#include "pythonic/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace
  {
    template <class I, class O, size_t M>
    void _nonzero(I begin, I end, O &out, types::array<long, M> &curr,
                  utils::int_<1>)
    {
      I start = begin;
      for (; begin != end; ++begin) {
        curr[M - 1] = begin - start;
        if (*begin)
          for (size_t i = 0; i < M; ++i) {
            *(out[i]) = curr[i];
            ++out[i];
          }
      }
    }

    template <class I, class O, size_t M, size_t N>
    void _nonzero(I begin, I end, O &out, types::array<long, M> &curr,
                  utils::int_<N>)
    {
      I start = begin;
      for (; begin != end; ++begin) {
        curr[M - N] = begin - start;
        _nonzero((*begin).begin(), (*begin).end(), out, curr,
                 utils::int_<N - 1>());
      }
    }
  }

  template <size_t... Is>
  types::array<utils::shared_ref<types::raw_array<long>>, sizeof...(Is)>
  init_buffers(long sz, utils::index_sequence<Is...>)
  {
    auto fwd = [](long ret, long) { return ret; }; // just to avoid a warning
    return {{fwd(sz, Is)...}};                     // too much memory used
  }

  template <class E>
  auto nonzero(E const &expr)
      -> types::array<types::ndarray<long, types::array<long, 1>>, E::value>
  {
    constexpr long N = E::value;
    typedef types::array<types::ndarray<long, types::array<long, 1>>, E::value>
        out_type;
    long sz = expr.flat_size();

    types::array<utils::shared_ref<types::raw_array<long>>, N> out_buffers =
        init_buffers(sz, utils::make_index_sequence<N>());
    types::array<long *, N> out_iters;
    for (size_t i = 0; i < N; ++i)
      out_iters[i] = out_buffers[i]->data;

    types::array<long, N> indices;
    _nonzero(expr.begin(), expr.end(), out_iters, indices, utils::int_<N>());

    types::array<long, 1> shape = {
        {(long)(out_iters[0] - out_buffers[0]->data)}};

    out_type out;
    for (size_t i = 0; i < N; ++i)
      out[i] = types::ndarray<long, types::array<long, 1>>(
          std::move(out_buffers[i]), shape);

    return out;
  }
}
PYTHONIC_NS_END

#endif
