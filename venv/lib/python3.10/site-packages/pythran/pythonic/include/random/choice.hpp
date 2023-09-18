#ifndef PYTHONIC_INCLUDE_RANDOM_CHOICE_HPP
#define PYTHONIC_INCLUDE_RANDOM_CHOICE_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/random/random.hpp"

PYTHONIC_NS_BEGIN

namespace random
{

  template <class Seq>
  typename Seq::value_type choice(Seq const &seq);

  DEFINE_FUNCTOR(pythonic::random, choice);
}
PYTHONIC_NS_END

#endif
