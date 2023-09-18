#ifndef PYTHONIC_INCLUDE_TYPES_NDITERATOR_HPP
#define PYTHONIC_INCLUDE_TYPES_NDITERATOR_HPP

#include <iterator>

#ifdef USE_XSIMD
#include <xsimd/xsimd.hpp>
#endif

PYTHONIC_NS_BEGIN

namespace types
{
  struct fast {
  };

  template <class T>
  auto fast_begin(T const &e) ->
      typename std::enable_if<has_fast_iterator<T>::value,
                              decltype(e.begin(fast{}))>::type
  {
    return e.begin(fast{});
  }
  template <class T>
  auto fast_begin(T const &e) ->
      typename std::enable_if<!has_fast_iterator<T>::value,
                              decltype(e.begin())>::type
  {
    return e.begin();
  }
  template <class T>
  auto fast_end(T const &e) ->
      typename std::enable_if<has_fast_iterator<T>::value,
                              decltype(e.end(fast{}))>::type
  {
    return e.end(fast{});
  }
  template <class T>
  auto fast_end(T const &e) ->
      typename std::enable_if<!has_fast_iterator<T>::value,
                              decltype(e.end())>::type
  {
    return e.end();
  }

  /* Iterator over whatever provides a fast(long) method to access its element
   */
  template <class E>
  struct nditerator
      : public std::iterator<std::random_access_iterator_tag,
                             typename std::remove_reference<decltype(
                                 std::declval<E &>().fast(0))>::type> {
    E &data;
    long index;
    nditerator(E &data, long index);

    auto operator*() -> decltype(data.fast(index));
    auto operator*() const -> decltype(data.fast(index));
    nditerator<E> &operator++();
    nditerator<E> &operator--();
    nditerator<E> &operator+=(long i);
    nditerator<E> &operator-=(long i);
    nditerator<E> operator+(long i) const;
    nditerator<E> operator-(long i) const;
    long operator-(nditerator<E> const &other) const;
    bool operator!=(nditerator<E> const &other) const;
    bool operator==(nditerator<E> const &other) const;
    bool operator<(nditerator<E> const &other) const;
    nditerator &operator=(nditerator const &other);
  };

  /* Const iterator over whatever provides a fast(long) method to access its
   * element
   */
  template <class E>
  struct const_nditerator
      : public std::iterator<std::random_access_iterator_tag,
                             typename std::remove_reference<decltype(
                                 std::declval<E &>().fast(0))>::type> {
    E const &data;
    long index;
    const_nditerator(E const &data, long index);

    auto operator*() const -> decltype(data.fast(index));
    const_nditerator<E> &operator++();
    const_nditerator<E> &operator--();
    const_nditerator<E> &operator+=(long i);
    const_nditerator<E> &operator-=(long i);
    const_nditerator<E> operator+(long i) const;
    const_nditerator<E> operator-(long i) const;
    long operator-(const_nditerator<E> const &other) const;
    bool operator!=(const_nditerator<E> const &other) const;
    bool operator==(const_nditerator<E> const &other) const;
    bool operator<(const_nditerator<E> const &other) const;
    const_nditerator &operator=(const_nditerator const &other);
  };
#ifdef USE_XSIMD
  template <class E>
  struct const_simd_nditerator
      : public std::iterator<std::random_access_iterator_tag,
                             xsimd::batch<typename E::dtype>> {

    using vector_type = typename xsimd::batch<typename E::dtype>;
    typename E::dtype const *data;
    static const std::size_t vector_size = vector_type::size;

    const_simd_nditerator(typename E::dtype const *data);

    auto operator*() const -> decltype(xsimd::load_unaligned(data));
    const_simd_nditerator &operator++();
    const_simd_nditerator &operator+=(long);
    const_simd_nditerator operator+(long);
    const_simd_nditerator &operator--();
    long operator-(const_simd_nditerator const &other) const;
    bool operator!=(const_simd_nditerator const &other) const;
    bool operator==(const_simd_nditerator const &other) const;
    bool operator<(const_simd_nditerator const &other) const;
    const_simd_nditerator &operator=(const_simd_nditerator const &other);
    void store(xsimd::batch<typename E::dtype> const &);
  };
  template <class E>
  struct const_simd_nditerator_nostep : const_simd_nditerator<E> {
    const_simd_nditerator_nostep &operator++()
    {
      return *this;
    }
    const_simd_nditerator_nostep &operator+=(long)
    {
      return *this;
    }
    const_simd_nditerator_nostep &operator--()
    {
      return *this;
    }
    const_simd_nditerator_nostep &
    operator=(const_simd_nditerator_nostep const &other) = default;
  };
#endif

  // build an iterator over T, selecting a raw pointer if possible
  template <bool is_strided>
  struct make_nditerator {
    template <class T>
    auto operator()(T &self, long i) -> decltype(nditerator<T>(self, i)) const;
  };

  template <>
  struct make_nditerator<false> {
    template <class T>
    typename T::dtype *operator()(T &self, long i) const;
  };

  template <bool is_strided>
  struct make_const_nditerator {
    template <class T>
    auto operator()(T const &self, long i)
        -> decltype(const_nditerator<T>(self, i)) const;
  };

  template <>
  struct make_const_nditerator<false> {
    template <class T>
    typename T::dtype const *operator()(T const &self, long i) const;
  };
}
PYTHONIC_NS_END

#endif
