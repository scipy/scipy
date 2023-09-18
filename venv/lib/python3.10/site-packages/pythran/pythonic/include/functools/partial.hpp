#ifndef PYTHONIC_INCLUDE_FUNCTOOLS_PARTIAL_HPP
#define PYTHONIC_INCLUDE_FUNCTOOLS_PARTIAL_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/seq.hpp"

#include <utility>
#include <tuple>

PYTHONIC_NS_BEGIN

namespace functools
{

  namespace details
  {

    /* a task that captures its environment for later call */
    template <typename... ClosureTypes>
    struct task {

      using callable = void;
      friend std::ostream &operator<<(std::ostream &os, task)
      {
        return os << "partial_function_wrapper";
      }

      mutable std::tuple<ClosureTypes...> closure; // closure associated to
                                                   // the task, mutable
                                                   // because pythran assumes
                                                   // all function calls are
                                                   // const

      task();
      task(task const &) = default;
      task(ClosureTypes const &... types);

      template <std::size_t... S, typename... Types>
      auto call(utils::index_sequence<S...>, Types &&... types) const
          -> decltype(std::get<0>(closure)(std::get<S + 1>(closure)...,
                                           std::forward<Types>(types)...))
      {
        return std::get<0>(closure)(std::get<S + 1>(closure)...,
                                    std::forward<Types>(types)...);
      }

      template <typename... Types>
      auto operator()(Types &&... types) const -> decltype(
          this->call(utils::make_index_sequence<sizeof...(ClosureTypes)-1>(),
                     std::forward<Types>(types)...));
    };
  }

  template <typename... Types>
  // remove references as closure capture the env by copy
  details::task<typename std::remove_cv<
      typename std::remove_reference<Types>::type>::type...>
  partial(Types &&... types);

  DEFINE_FUNCTOR(pythonic::functools, partial);
}
PYTHONIC_NS_END

#endif
