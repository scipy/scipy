#ifndef PYTHONIC_INCLUDE_TYPES_IMMEDIATE_HPP
#define PYTHONIC_INCLUDE_TYPES_IMMEDIATE_HPP

PYTHONIC_NS_BEGIN

namespace types
{

  template <class T, T Val>
  struct immediate {
    immediate() = default;
    immediate(immediate const &) = default;
    immediate(immediate &&) = default;

    operator T() const
    {
      return Val;
    }

    template <class U, U Wal,
              class _ = typename std::enable_if<Val == (T)Wal, void>::type>
    immediate(std::integral_constant<U, Wal>)
    {
    }
  };

  using true_immediate = immediate<bool, true>;
  using false_immediate = immediate<bool, false>;
}

PYTHONIC_NS_END

#endif
