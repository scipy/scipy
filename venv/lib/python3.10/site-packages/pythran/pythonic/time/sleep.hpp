#ifndef PYTHONIC_TIME_SLEEP_HPP
#define PYTHONIC_TIME_SLEEP_HPP

#include "pythonic/include/time/sleep.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/builtins/None.hpp"

#include <thread>
#include <chrono>

PYTHONIC_NS_BEGIN

namespace time
{

  types::none_type sleep(double const value)
  {
    std::this_thread::sleep_for(std::chrono::duration<double>(value));
    return builtins::None;
  }
}
PYTHONIC_NS_END

#endif
