#ifndef PYTHONIC_INCLUDE_NUMPY_CONCATENATE_HPP
#define PYTHONIC_INCLUDE_NUMPY_CONCATENATE_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, size_t M, class V>
  types::ndarray<typename E::dtype, types::array<long, E::value>>
  concatenate(types::array_base<E, M, V> const &args, long axis = 0);

  template <class... Types>
  auto concatenate(std::tuple<Types...> const &args, long axis = 0)
      -> types::ndarray<
          typename __combined<typename std::decay<Types>::type::dtype...>::type,
          types::array<
              long, std::tuple_element<0, std::tuple<Types...>>::type::value>>;

  template <class E>
  types::ndarray<typename E::dtype, types::array<long, E::value>>
  concatenate(types::list<E> const &args, long axis = 0);

  DEFINE_FUNCTOR(pythonic::numpy, concatenate);
}
PYTHONIC_NS_END

#endif
