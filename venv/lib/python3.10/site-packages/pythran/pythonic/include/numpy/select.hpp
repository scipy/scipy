#ifndef PYTHONIC_INCLUDE_NUMPY_SELECT_HPP
#define PYTHONIC_INCLUDE_NUMPY_SELECT_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/int_.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class C, class L>
  types::ndarray<typename L::dtype, types::array<long, L::value - 1>>
  select(C const &condlist, L const &choicelist,
         typename L::dtype _default = 0);

  template <class T, class TpS, class U, class UpS>
  typename std::enable_if<
      std::tuple_size<TpS>::value == std::tuple_size<UpS>::value,
      types::ndarray<T, types::array<long, std::tuple_size<TpS>::value>>>::type
  select(types::list<types::ndarray<U, UpS>> const &condlist,
         types::list<types::ndarray<T, TpS>> const &choicelist, T _default = 0);

  template <class T, class TpS, class U, class UpS, size_t M>
  typename std::enable_if<std::tuple_size<TpS>::value ==
                              std::tuple_size<UpS>::value,
                          types::ndarray<T, TpS>>::type
  select(types::static_list<types::ndarray<U, UpS>, M> const &condlist,
         types::static_list<types::ndarray<T, TpS>, M> const &choicelist,
         T _default = 0);

  template <class T, class TpS, class U, class UpS, size_t M>
  typename std::enable_if<std::tuple_size<TpS>::value ==
                              std::tuple_size<UpS>::value,
                          types::ndarray<T, TpS>>::type
  select(types::static_list<types::ndarray<U, UpS>, M> const &condlist,
         types::list<types::ndarray<T, TpS>> const &choicelist, T _default = 0);

  template <class T, class TpS, class U, class UpS, size_t M>
  typename std::enable_if<std::tuple_size<TpS>::value ==
                              std::tuple_size<UpS>::value,
                          types::ndarray<T, TpS>>::type
  select(types::list<types::ndarray<U, UpS>> const &condlist,
         types::static_list<types::ndarray<T, TpS>, M> const &choicelist,
         T _default = 0);

  DEFINE_FUNCTOR(pythonic::numpy, select);
}
PYTHONIC_NS_END

#endif
