#ifndef PYTHONIC_INCLUDE_TYPES_NONE_HPP
#define PYTHONIC_INCLUDE_TYPES_NONE_HPP

#include "pythonic/include/operator_/mod.hpp"
#include "pythonic/include/types/assignable.hpp"
#include <ostream>

PYTHONIC_NS_BEGIN

namespace types
{

  static const intptr_t NONE_ID = 0x1331;

  struct none_type {
    none_type();
    intptr_t id() const;
  };

  inline std::ostream &operator<<(std::ostream &os, none_type const &)
  {
    return os << "None";
  }

  template <class T, bool is_fundamental = std::is_fundamental<T>::value>
  struct none;

  /* Type adapator to simulate an option type
   *
   * see http://en.wikipedia.org/wiki/Option_type
   */
  template <class T>
  struct none<T, false> : T {

    bool is_none; // set to true if the type is none

    none(none_type const &);

    none() : T(), is_none{true}
    {
    }
    none(none const &other) = default;
    none(T const &arg) : T(arg), is_none(false)
    {
    }
    template <class OT>
    none(OT const &arg) : none(T(arg))
    {
    }

    bool operator==(none_type const &) const;

    template <class O>
    bool operator==(O const &t) const;

    bool operator!=(none_type const &) const;

    template <class O>
    bool operator!=(O const &t) const;

    explicit operator bool() const;

    intptr_t id() const;

    template <class T0>
    friend std::ostream &operator<<(std::ostream &os, none<T0, false> const &);
  };

  /* specialization of none for integral types we cannot derive from
   */
  template <class P, class T>
  struct none_data {
    explicit operator bool() const
    {
      return !static_cast<P const *>(this)->is_none &&
             static_cast<P const *>(this)->data;
    }
    operator T() const
    {
      return static_cast<P const *>(this)->data;
    }
  };
  template <class P>
  struct none_data<P, bool> {
    operator bool() const
    {
      return !static_cast<P const *>(this)->is_none &&
             static_cast<P const *>(this)->data;
    }
  };
  template <class T>
  struct none<T, true> : none_data<none<T, true>, T> {
    T data;
    template <class T1>
    friend std::ostream &operator<<(std::ostream &, none<T1, true> const &);
    template <class T1>
    friend T1 operator+(none<T1, true> const &t0, T1 const &t1);
    template <class T1>
    friend T1 operator+(T1 const &t0, none<T1, true> const &t1);
    template <class T1>
    friend none<T1, true> operator+(none<T1, true> const &t0,
                                    none<T1, true> const &t1);
    template <class T1>
    friend bool operator>(none<T1, true> const &t0, T1 const &t1);
    template <class T1>
    friend bool operator>(T1 const &t0, none<T1, true> const &t1);
    template <class T1>
    friend none<bool> operator>(none<T1, true> const &t0,
                                none<T1, true> const &t1);
    template <class T1>
    friend bool operator>=(none<T1, true> const &t0, T1 const &t1);
    template <class T1>
    friend bool operator>=(T1 const &t0, none<T1, true> const &t1);
    template <class T1>
    friend none<bool> operator>=(none<T1, true> const &t0,
                                 none<T1, true> const &t1);
    template <class T1>
    friend bool operator<(none<T1, true> const &t0, T1 const &t1);
    template <class T1>
    friend bool operator<(T1 const &t0, none<T1, true> const &t1);
    template <class T1>
    friend none<bool> operator<(none<T1, true> const &t0,
                                none<T1, true> const &t1);
    template <class T1>
    friend bool operator<=(none<T1, true> const &t0, T1 const &t1);
    template <class T1>
    friend bool operator<=(T1 const &t0, none<T1, true> const &t1);
    template <class T1>
    friend none<bool> operator<=(none<T1, true> const &t0,
                                 none<T1, true> const &t1);
    template <class T1>
    friend T1 operator-(none<T1, true> const &t0, T1 const &t1);
    template <class T1>
    friend T1 operator-(T1 const &t0, none<T1, true> const &t1);
    template <class T1>
    friend none<T1, true> operator-(none<T1, true> const &t0,
                                    none<T1, true> const &t1);
    template <class T1>
    friend T1 operator*(none<T1, true> const &t0, T1 const &t1);
    template <class T1>
    friend T1 operator*(T1 const &t0, none<T1, true> const &t1);
    template <class T1>
    friend none<T1, true> operator*(none<T1, true> const &t0,
                                    none<T1, true> const &t1);
    template <class T1>
    friend T1 operator/(none<T1, true> const &t0, T1 const &t1);
    template <class T1>
    friend T1 operator/(T1 const &t0, none<T1, true> const &t1);
    template <class T1>
    friend none<T1, true> operator/(none<T1, true> const &t0,
                                    none<T1, true> const &t1);

    template <class T1>
    none &operator+=(T1 other);
    template <class T1>
    none &operator-=(T1 other);
    template <class T1>
    none &operator*=(T1 other);
    template <class T1>
    none &operator/=(T1 other);

  public:
    bool is_none;
    none();
    none(none_type const &);
    none(T const &data);
    bool operator==(none_type const &) const;
    template <class O>
    bool operator==(O const &t) const;
    bool operator!=(none_type const &) const;
    template <class O>
    bool operator!=(O const &t) const;
    T &operator=(T const &t);
    intptr_t id() const;
    template <class T1>
    operator none<T1, true>()
    {
      if (is_none)
        return {none_type{}};
      else
        return {static_cast<T1>(data)};
    }
  };
  template <class T>
  T operator+(none<T, true> const &t0, T const &t1);
  template <class T>
  T operator+(T const &t0, none<T, true> const &t1);
  template <class T>
  none<T, true> operator+(none<T, true> const &t0, none<T, true> const &t1);
  template <class T>
  bool operator>(none<T, true> const &t0, T const &t1);
  template <class T>
  bool operator>(T const &t0, none<T, true> const &t1);
  template <class T>
  none<bool> operator>(none<T, true> const &t0, none<T, true> const &t1);
  template <class T>
  bool operator>=(none<T, true> const &t0, T const &t1);
  template <class T>
  bool operator>=(T const &t0, none<T, true> const &t1);
  template <class T>
  none<bool> operator>=(none<T, true> const &t0, none<T, true> const &t1);
  template <class T>
  bool operator<(none<T, true> const &t0, T const &t1);
  template <class T>
  bool operator<(T const &t0, none<T, true> const &t1);
  template <class T>
  none<bool> operator<(none<T, true> const &t0, none<T, true> const &t1);
  template <class T>
  bool operator<=(none<T, true> const &t0, T const &t1);
  template <class T>
  bool operator<=(T const &t0, none<T, true> const &t1);
  template <class T>
  none<bool> operator<=(none<T, true> const &t0, none<T, true> const &t1);
  template <class T>
  T operator-(none<T, true> const &t0, T const &t1);
  template <class T>
  T operator-(T const &t0, none<T, true> const &t1);
  template <class T>
  none<T, true> operator-(none<T, true> const &t0, none<T, true> const &t1);
  template <class T>
  T operator*(none<T, true> const &t0, T const &t1);
  template <class T>
  T operator*(T const &t0, none<T, true> const &t1);
  template <class T>
  none<T, true> operator*(none<T, true> const &t0, none<T, true> const &t1);
  template <class T>
  T operator/(none<T, true> const &t0, T const &t1);
  template <class T>
  T operator/(T const &t0, none<T, true> const &t1);
  template <class T>
  none<T, true> operator/(none<T, true> const &t0, none<T, true> const &t1);
  template <class T0, class T1>
  decltype(operator_::mod(std::declval<T0>(), std::declval<T1>()))
  operator%(none<T0, true> const &t0, T1 const &t1);
  template <class T0, class T1>
  decltype(operator_::mod(std::declval<T0>(), std::declval<T1>()))
  operator%(T0 const &t0, none<T1, true> const &t1);
  template <class T0, class T1>
  none<decltype(operator_::mod(std::declval<T0>(), std::declval<T1>())), true>
  operator%(none<T0, true> const &t0, none<T1, true> const &t1);
  template <class T>
  std::ostream &operator<<(std::ostream &os, none<T, true> const &v);

  template <class T>
  struct is_none {
    static const bool value = false;
  };

  template <class T>
  struct is_none<none<T>> {
    static const bool value = true;
  };
} // namespace types

template <class T>
struct assignable<types::none<T>> {
  using type = types::none<typename assignable<T>::type>;
};
PYTHONIC_NS_END

namespace std
{
  /* std::get overload */
  template <size_t I, class T0>
  auto get(pythonic::types::none<T0> const &t)
      -> decltype(std::get<I>((T0 const &)t));

  template <size_t I, class T0>
  struct tuple_element<I, pythonic::types::none<T0>> {
    using type = typename std::tuple_element<I, T0>::type;
  };
} // namespace std

/* type inference stuff { */
#include "pythonic/include/types/combined.hpp"

template <class T0, class T1>
struct __combined<pythonic::types::none<T0>, T1> {
  static_assert(!pythonic::types::is_none<T1>::value,
                "none of none should'nt exist");
  using type = pythonic::types::none<typename __combined<T0, T1>::type>;
};

template <class T0, class T1>
struct __combined<T1, pythonic::types::none<T0>> {
  static_assert(!pythonic::types::is_none<T0>::value,
                "none of none should'nt exist");
  using type = pythonic::types::none<typename __combined<T0, T1>::type>;
};

template <class T0, class T1>
struct __combined<pythonic::types::none<T1>, pythonic::types::none<T0>> {
  static_assert(!pythonic::types::is_none<T0>::value,
                "none of none shouldn't exist");
  static_assert(!pythonic::types::is_none<T1>::value,
                "none of none shouldn't exist");
  using type = pythonic::types::none<typename __combined<T0, T1>::type>;
};

template <class T>
struct __combined<pythonic::types::none_type, T> {
  static_assert(!pythonic::types::is_none<T>::value,
                "none of none shouldn't exist");
  using type = pythonic::types::none<T>;
};

template <class T>
struct __combined<pythonic::types::none_type, pythonic::types::none<T>> {
  static_assert(!pythonic::types::is_none<T>::value,
                "none of none shouldn't exist");
  using type = pythonic::types::none<T>;
};

template <class T>
struct __combined<T, pythonic::types::none_type> {
  static_assert(!pythonic::types::is_none<T>::value,
                "none of none shouldn't exist");
  using type = pythonic::types::none<T>;
};
template <class T>
struct __combined<pythonic::types::none<T>, pythonic::types::none_type> {
  static_assert(!pythonic::types::is_none<T>::value,
                "none of none shouldn't exist");
  using type = pythonic::types::none<T>;
};

template <>
struct __combined<pythonic::types::none_type, pythonic::types::none_type> {
  using type = pythonic::types::none_type;
};

/* } */

#ifdef ENABLE_PYTHON_MODULE

#include "pythonic/python/core.hpp"

PYTHONIC_NS_BEGIN
template <>
struct to_python<types::none_type> {
  static PyObject *convert(types::none_type);
};

template <class T>
struct to_python<types::none<T>> {
  static PyObject *convert(types::none<T> const &n);
};

template <>
struct from_python<types::none_type> {

  static bool is_convertible(PyObject *obj);

  static types::none_type convert(PyObject *obj);
};
PYTHONIC_NS_END

#endif
#endif
