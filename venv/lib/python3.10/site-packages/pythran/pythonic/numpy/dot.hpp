#ifndef PYTHONIC_NUMPY_DOT_HPP
#define PYTHONIC_NUMPY_DOT_HPP

#include "pythonic/include/numpy/dot.hpp"

#include "pythonic/numpy/multiply.hpp"
#include "pythonic/numpy/sum.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/traits.hpp"

#ifdef PYTHRAN_BLAS_NONE
#error pythran configured without BLAS but BLAS seem needed
#endif

#if defined(PYTHRAN_BLAS_ATLAS) || defined(PYTHRAN_BLAS_SATLAS)
extern "C" {
#endif
#include <cblas.h>
#if defined(PYTHRAN_BLAS_ATLAS) || defined(PYTHRAN_BLAS_SATLAS)
}
#endif

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, class F>
  typename std::enable_if<types::is_dtype<E>::value &&
                              types::is_dtype<F>::value,
                          decltype(std::declval<E>() * std::declval<F>())>::type
  dot(E const &e, F const &f)
  {
    return e * f;
  }

  template <class E>
  struct blas_buffer_t {
    typename E::dtype const *operator()(E const &e) const
    {
      return e.buffer;
    }
  };
  template <class T>
  struct blas_buffer_t<types::list<T>> {
    T const *operator()(types::list<T> const &e) const
    {
      return &e.fast(0);
    }
  };
  template <class T, size_t N>
  struct blas_buffer_t<types::array<T, N>> {
    T const *operator()(types::array<T, N> const &e) const
    {
      return e.data();
    }
  };

  template <class E>
  auto blas_buffer(E const &e) -> decltype(blas_buffer_t<E>{}(e))
  {
    return blas_buffer_t<E>{}(e);
  }

  template <class E, class F>
  typename std::enable_if<
      types::is_numexpr_arg<E>::value &&
          types::is_numexpr_arg<F>::value   // Arguments are array_like
          && E::value == 1 && F::value == 1 // It is a two vectors.
          && (!is_blas_array<E>::value || !is_blas_array<F>::value ||
              !std::is_same<typename E::dtype, typename F::dtype>::value),
      typename __combined<typename E::dtype, typename F::dtype>::type>::type
  dot(E const &e, F const &f)
  {
    return sum(functor::multiply{}(e, f));
  }

  template <class E, class F>
  typename std::enable_if<E::value == 1 && F::value == 1 &&
                              std::is_same<typename E::dtype, float>::value &&
                              std::is_same<typename F::dtype, float>::value &&
                              is_blas_array<E>::value &&
                              is_blas_array<F>::value,
                          float>::type
  dot(E const &e, F const &f)
  {
    return cblas_sdot(e.size(), blas_buffer(e), 1, blas_buffer(f), 1);
  }

  template <class E, class F>
  typename std::enable_if<E::value == 1 && F::value == 1 &&
                              std::is_same<typename E::dtype, double>::value &&
                              std::is_same<typename F::dtype, double>::value &&
                              is_blas_array<E>::value &&
                              is_blas_array<F>::value,
                          double>::type
  dot(E const &e, F const &f)
  {
    return cblas_ddot(e.size(), blas_buffer(e), 1, blas_buffer(f), 1);
  }

  template <class E, class F>
  typename std::enable_if<
      E::value == 1 && F::value == 1 &&
          std::is_same<typename E::dtype, std::complex<float>>::value &&
          std::is_same<typename F::dtype, std::complex<float>>::value &&
          is_blas_array<E>::value && is_blas_array<F>::value,
      std::complex<float>>::type
  dot(E const &e, F const &f)
  {
    std::complex<float> out;
    cblas_cdotu_sub(e.size(), blas_buffer(e), 1, blas_buffer(f), 1, &out);
    return out;
  }

  template <class E, class F>
  typename std::enable_if<
      E::value == 1 && F::value == 1 &&
          std::is_same<typename E::dtype, std::complex<double>>::value &&
          std::is_same<typename F::dtype, std::complex<double>>::value &&
          is_blas_array<E>::value && is_blas_array<F>::value,
      std::complex<double>>::type
  dot(E const &e, F const &f)
  {
    std::complex<double> out;
    cblas_zdotu_sub(e.size(), blas_buffer(e), 1, blas_buffer(f), 1, &out);
    return out;
  }

  /// Matrice / Vector multiplication

#define MV_DEF(T, L)                                                           \
  inline void mv(int m, int n, T *A, T *B, T *C)                               \
  {                                                                            \
    cblas_##L##gemv(CblasRowMajor, CblasNoTrans, n, m, 1, A, m, B, 1, 0, C,    \
                    1);                                                        \
  }

  MV_DEF(double, d)
  MV_DEF(float, s)

#undef MV_DEF

#define TV_DEF(T, L)                                                           \
  inline void tv(int m, int n, T *A, T *B, T *C)                               \
  {                                                                            \
    cblas_##L##gemv(CblasRowMajor, CblasTrans, m, n, 1, A, n, B, 1, 0, C, 1);  \
  }

  TV_DEF(double, d)
  TV_DEF(float, s)

#undef TV_DEF

#define MV_DEF(T, K, L)                                                        \
  inline void mv(int m, int n, T *A, T *B, T *C)                               \
  {                                                                            \
    T alpha = 1, beta = 0;                                                     \
    cblas_##L##gemv(CblasRowMajor, CblasNoTrans, n, m, (K *)&alpha, (K *)A, m, \
                    (K *)B, 1, (K *)&beta, (K *)C, 1);                         \
  }
  MV_DEF(std::complex<float>, float, c)
  MV_DEF(std::complex<double>, double, z)
#undef MV_DEF

  template <class E, class pS0, class pS1>
  typename std::enable_if<is_blas_type<E>::value &&
                              std::tuple_size<pS0>::value == 2 &&
                              std::tuple_size<pS1>::value == 1,
                          types::ndarray<E, types::pshape<long>>>::type
  dot(types::ndarray<E, pS0> const &f, types::ndarray<E, pS1> const &e)
  {
    types::ndarray<E, types::pshape<long>> out(
        types::pshape<long>{f.template shape<0>()}, builtins::None);
    const int m = f.template shape<1>(), n = f.template shape<0>();
    mv(m, n, f.buffer, e.buffer, out.buffer);
    return out;
  }

  template <class E, class pS0, class pS1>
  typename std::enable_if<is_blas_type<E>::value &&
                              std::tuple_size<pS0>::value == 2 &&
                              std::tuple_size<pS1>::value == 1,
                          types::ndarray<E, types::pshape<long>>>::type
  dot(types::numpy_texpr<types::ndarray<E, pS0>> const &f,
      types::ndarray<E, pS1> const &e)
  {
    types::ndarray<E, types::pshape<long>> out(
        types::pshape<long>{f.template shape<0>()}, builtins::None);
    const int m = f.template shape<1>(), n = f.template shape<0>();
    tv(m, n, f.arg.buffer, e.buffer, out.buffer);
    return out;
  }

// The trick is to not transpose the matrix so that MV become VM
#define VM_DEF(T, L)                                                           \
  inline void vm(int m, int n, T *A, T *B, T *C)                               \
  {                                                                            \
    cblas_##L##gemv(CblasRowMajor, CblasTrans, n, m, 1, A, m, B, 1, 0, C, 1);  \
  }

  VM_DEF(double, d)
  VM_DEF(float, s)

#undef VM_DEF
#define VT_DEF(T, L)                                                           \
  inline void vt(int m, int n, T *A, T *B, T *C)                               \
  {                                                                            \
    cblas_##L##gemv(CblasRowMajor, CblasNoTrans, m, n, 1, A, n, B, 1, 0, C,    \
                    1);                                                        \
  }

  VT_DEF(double, d)
  VT_DEF(float, s)

#undef VM_DEF
#define VM_DEF(T, K, L)                                                        \
  inline void vm(int m, int n, T *A, T *B, T *C)                               \
  {                                                                            \
    T alpha = 1, beta = 0;                                                     \
    cblas_##L##gemv(CblasRowMajor, CblasTrans, n, m, (K *)&alpha, (K *)A, m,   \
                    (K *)B, 1, (K *)&beta, (K *)C, 1);                         \
  }
  VM_DEF(std::complex<float>, float, c)
  VM_DEF(std::complex<double>, double, z)
#undef VM_DEF

  template <class E, class pS0, class pS1>
  typename std::enable_if<is_blas_type<E>::value &&
                              std::tuple_size<pS0>::value == 1 &&
                              std::tuple_size<pS1>::value == 2,
                          types::ndarray<E, types::pshape<long>>>::type
  dot(types::ndarray<E, pS0> const &e, types::ndarray<E, pS1> const &f)
  {
    types::ndarray<E, types::pshape<long>> out(
        types::pshape<long>{f.template shape<1>()}, builtins::None);
    const int m = f.template shape<1>(), n = f.template shape<0>();
    vm(m, n, f.buffer, e.buffer, out.buffer);
    return out;
  }

  template <class E, class pS0, class pS1>
  typename std::enable_if<is_blas_type<E>::value &&
                              std::tuple_size<pS0>::value == 1 &&
                              std::tuple_size<pS1>::value == 2,
                          types::ndarray<E, types::pshape<long>>>::type
  dot(types::ndarray<E, pS0> const &e,
      types::numpy_texpr<types::ndarray<E, pS1>> const &f)
  {
    types::ndarray<E, types::pshape<long>> out(
        types::pshape<long>{f.template shape<1>()}, builtins::None);
    const int m = f.template shape<1>(), n = f.template shape<0>();
    vt(m, n, f.arg.buffer, e.buffer, out.buffer);
    return out;
  }

  // If arguments could be use with blas, we evaluate them as we need pointer
  // on array for blas
  template <class E, class F>
  typename std::enable_if<
      types::is_numexpr_arg<E>::value &&
          types::is_numexpr_arg<F>::value // It is an array_like
          && (!(types::is_ndarray<E>::value && types::is_ndarray<F>::value) ||
              !std::is_same<typename E::dtype, typename F::dtype>::value) &&
          is_blas_type<typename E::dtype>::value &&
          is_blas_type<typename F::dtype>::value // With dtype compatible with
                                                 // blas
          && E::value == 2 && F::value == 1,     // And it is matrix / vect
      types::ndarray<
          typename __combined<typename E::dtype, typename F::dtype>::type,
          types::pshape<long>>>::type
  dot(E const &e, F const &f)
  {
    types::ndarray<
        typename __combined<typename E::dtype, typename F::dtype>::type,
        typename E::shape_t>
        e_ = e;
    types::ndarray<
        typename __combined<typename E::dtype, typename F::dtype>::type,
        typename F::shape_t>
        f_ = f;
    return dot(e_, f_);
  }

  // If arguments could be use with blas, we evaluate them as we need pointer
  // on array for blas
  template <class E, class F>
  typename std::enable_if<
      types::is_numexpr_arg<E>::value &&
          types::is_numexpr_arg<F>::value // It is an array_like
          && (!(types::is_ndarray<E>::value && types::is_ndarray<F>::value) ||
              !std::is_same<typename E::dtype, typename F::dtype>::value) &&
          is_blas_type<typename E::dtype>::value &&
          is_blas_type<typename F::dtype>::value // With dtype compatible with
                                                 // blas
          && E::value == 1 && F::value == 2,     // And it is vect / matrix
      types::ndarray<
          typename __combined<typename E::dtype, typename F::dtype>::type,
          types::pshape<long>>>::type
  dot(E const &e, F const &f)
  {
    types::ndarray<
        typename __combined<typename E::dtype, typename F::dtype>::type,
        typename E::shape_t>
        e_ = e;
    types::ndarray<
        typename __combined<typename E::dtype, typename F::dtype>::type,
        typename F::shape_t>
        f_ = f;
    return dot(e_, f_);
  }

  // If one of the arg doesn't have a "blas compatible type", we use a slow
  // matrix vector multiplication.
  template <class E, class F>
  typename std::enable_if<
      (!is_blas_type<typename E::dtype>::value ||
       !is_blas_type<typename F::dtype>::value) &&
          E::value == 1 && F::value == 2, // And it is vect / matrix
      types::ndarray<
          typename __combined<typename E::dtype, typename F::dtype>::type,
          types::pshape<long>>>::type
  dot(E const &e, F const &f)
  {
    types::ndarray<
        typename __combined<typename E::dtype, typename F::dtype>::type,
        types::pshape<long>>
        out(types::pshape<long>{f.template shape<1>()}, 0);
    for (long i = 0; i < out.template shape<0>(); i++)
      for (long j = 0; j < f.template shape<0>(); j++)
        out[i] += e[j] * f[types::array<long, 2>{{j, i}}];
    return out;
  }

  // If one of the arg doesn't have a "blas compatible type", we use a slow
  // matrix vector multiplication.
  template <class E, class F>
  typename std::enable_if<
      (!is_blas_type<typename E::dtype>::value ||
       !is_blas_type<typename F::dtype>::value) &&
          E::value == 2 && F::value == 1, // And it is vect / matrix
      types::ndarray<
          typename __combined<typename E::dtype, typename F::dtype>::type,
          types::pshape<long>>>::type
  dot(E const &e, F const &f)
  {
    types::ndarray<
        typename __combined<typename E::dtype, typename F::dtype>::type,
        types::pshape<long>>
        out(types::pshape<long>{e.template shape<0>()}, 0);
    for (long i = 0; i < out.template shape<0>(); i++)
      for (long j = 0; j < f.template shape<0>(); j++)
        out[i] += e[types::array<long, 2>{{i, j}}] * f[j];
    return out;
  }

  /// Matrix / Matrix multiplication

#define MM_DEF(T, L)                                                           \
  inline void mm(int m, int n, int k, T *A, T *B, T *C)                        \
  {                                                                            \
    cblas_##L##gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, 1, A,  \
                    k, B, n, 0, C, n);                                         \
  }
  MM_DEF(double, d)
  MM_DEF(float, s)
#undef MM_DEF
#define MM_DEF(T, K, L)                                                        \
  inline void mm(int m, int n, int k, T *A, T *B, T *C)                        \
  {                                                                            \
    T alpha = 1, beta = 0;                                                     \
    cblas_##L##gemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k,        \
                    (K *)&alpha, (K *)A, k, (K *)B, n, (K *)&beta, (K *)C, n); \
  }
  MM_DEF(std::complex<float>, float, c)
  MM_DEF(std::complex<double>, double, z)
#undef MM_DEF

  template <class E, class pS0, class pS1>
  typename std::enable_if<is_blas_type<E>::value &&
                              std::tuple_size<pS0>::value == 2 &&
                              std::tuple_size<pS1>::value == 2,
                          types::ndarray<E, types::array<long, 2>>>::type
  dot(types::ndarray<E, pS0> const &a, types::ndarray<E, pS1> const &b)
  {
    int n = b.template shape<1>(), m = a.template shape<0>(),
        k = b.template shape<0>();

    types::ndarray<E, types::array<long, 2>> out(types::array<long, 2>{{m, n}},
                                                 builtins::None);
    mm(m, n, k, a.buffer, b.buffer, out.buffer);
    return out;
  }

  template <class E, class pS0, class pS1, class pS2>
  typename std::enable_if<
      is_blas_type<E>::value && std::tuple_size<pS0>::value == 2 &&
          std::tuple_size<pS1>::value == 2 && std::tuple_size<pS2>::value == 2,
      types::ndarray<E, pS2>>::type &
  dot(types::ndarray<E, pS0> const &a, types::ndarray<E, pS1> const &b,
      types::ndarray<E, pS2> &c)
  {
    int n = b.template shape<1>(), m = a.template shape<0>(),
        k = b.template shape<0>();

    mm(m, n, k, a.buffer, b.buffer, c.buffer);
    return c;
  }

#define TM_DEF(T, L)                                                           \
  inline void tm(int m, int n, int k, T *A, T *B, T *C)                        \
  {                                                                            \
    cblas_##L##gemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k, 1, A, m, \
                    B, n, 0, C, n);                                            \
  }
  TM_DEF(double, d)
  TM_DEF(float, s)
#undef TM_DEF
#define TM_DEF(T, K, L)                                                        \
  inline void tm(int m, int n, int k, T *A, T *B, T *C)                        \
  {                                                                            \
    T alpha = 1, beta = 0;                                                     \
    cblas_##L##gemm(CblasRowMajor, CblasTrans, CblasNoTrans, m, n, k,          \
                    (K *)&alpha, (K *)A, m, (K *)B, n, (K *)&beta, (K *)C, n); \
  }
  TM_DEF(std::complex<float>, float, c)
  TM_DEF(std::complex<double>, double, z)
#undef TM_DEF

  template <class E, class pS0, class pS1>
  typename std::enable_if<is_blas_type<E>::value &&
                              std::tuple_size<pS0>::value == 2 &&
                              std::tuple_size<pS1>::value == 2,
                          types::ndarray<E, types::array<long, 2>>>::type
  dot(types::numpy_texpr<types::ndarray<E, pS0>> const &a,
      types::ndarray<E, pS1> const &b)
  {
    int n = b.template shape<1>(), m = a.template shape<0>(),
        k = b.template shape<0>();

    types::ndarray<E, types::array<long, 2>> out(types::array<long, 2>{{m, n}},
                                                 builtins::None);
    tm(m, n, k, a.arg.buffer, b.buffer, out.buffer);
    return out;
  }

#define MT_DEF(T, L)                                                           \
  inline void mt(int m, int n, int k, T *A, T *B, T *C)                        \
  {                                                                            \
    cblas_##L##gemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k, 1, A, k, \
                    B, k, 0, C, n);                                            \
  }
  MT_DEF(double, d)
  MT_DEF(float, s)
#undef MT_DEF
#define MT_DEF(T, K, L)                                                        \
  inline void mt(int m, int n, int k, T *A, T *B, T *C)                        \
  {                                                                            \
    T alpha = 1, beta = 0;                                                     \
    cblas_##L##gemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, n, k,          \
                    (K *)&alpha, (K *)A, k, (K *)B, k, (K *)&beta, (K *)C, n); \
  }
  MT_DEF(std::complex<float>, float, c)
  MT_DEF(std::complex<double>, double, z)
#undef MT_DEF

  template <class E, class pS0, class pS1>
  typename std::enable_if<is_blas_type<E>::value &&
                              std::tuple_size<pS0>::value == 2 &&
                              std::tuple_size<pS1>::value == 2,
                          types::ndarray<E, types::array<long, 2>>>::type
  dot(types::ndarray<E, pS0> const &a,
      types::numpy_texpr<types::ndarray<E, pS1>> const &b)
  {
    int n = b.template shape<1>(), m = a.template shape<0>(),
        k = b.template shape<0>();

    types::ndarray<E, types::array<long, 2>> out(types::array<long, 2>{{m, n}},
                                                 builtins::None);
    mt(m, n, k, a.buffer, b.arg.buffer, out.buffer);
    return out;
  }

#define TT_DEF(T, L)                                                           \
  inline void tt(int m, int n, int k, T *A, T *B, T *C)                        \
  {                                                                            \
    cblas_##L##gemm(CblasRowMajor, CblasTrans, CblasTrans, m, n, k, 1, A, m,   \
                    B, k, 0, C, n);                                            \
  }
  TT_DEF(double, d)
  TT_DEF(float, s)
#undef TT_DEF
#define TT_DEF(T, K, L)                                                        \
  inline void tt(int m, int n, int k, T *A, T *B, T *C)                        \
  {                                                                            \
    T alpha = 1, beta = 0;                                                     \
    cblas_##L##gemm(CblasRowMajor, CblasTrans, CblasTrans, m, n, k,            \
                    (K *)&alpha, (K *)A, m, (K *)B, k, (K *)&beta, (K *)C, n); \
  }
  TT_DEF(std::complex<float>, float, c)
  TT_DEF(std::complex<double>, double, z)
#undef TT_DEF

  template <class E, class pS0, class pS1>
  typename std::enable_if<is_blas_type<E>::value &&
                              std::tuple_size<pS0>::value == 2 &&
                              std::tuple_size<pS1>::value == 2,
                          types::ndarray<E, types::array<long, 2>>>::type
  dot(types::numpy_texpr<types::ndarray<E, pS0>> const &a,
      types::numpy_texpr<types::ndarray<E, pS1>> const &b)
  {
    int n = b.template shape<1>(), m = a.template shape<0>(),
        k = b.template shape<0>();

    types::ndarray<E, types::array<long, 2>> out(types::array<long, 2>{{m, n}},
                                                 builtins::None);
    tt(m, n, k, a.arg.buffer, b.arg.buffer, out.buffer);
    return out;
  }

  // If arguments could be use with blas, we evaluate them as we need pointer
  // on array for blas
  template <class E, class F>
  typename std::enable_if<
      types::is_numexpr_arg<E>::value &&
          types::is_numexpr_arg<F>::value // It is an array_like
          && (!(types::is_ndarray<E>::value && types::is_ndarray<F>::value) ||
              !std::is_same<typename E::dtype, typename F::dtype>::value) &&
          is_blas_type<typename E::dtype>::value &&
          is_blas_type<typename F::dtype>::value // With dtype compatible with
                                                 // blas
          && E::value == 2 && F::value == 2,     // And both are matrix
      types::ndarray<
          typename __combined<typename E::dtype, typename F::dtype>::type,
          types::array<long, 2>>>::type
  dot(E const &e, F const &f)
  {
    types::ndarray<
        typename __combined<typename E::dtype, typename F::dtype>::type,
        typename E::shape_t>
        e_ = e;
    types::ndarray<
        typename __combined<typename E::dtype, typename F::dtype>::type,
        typename F::shape_t>
        f_ = f;
    return dot(e_, f_);
  }

  // If one of the arg doesn't have a "blas compatible type", we use a slow
  // matrix multiplication.
  template <class E, class F>
  typename std::enable_if<
      (!is_blas_type<typename E::dtype>::value ||
       !is_blas_type<typename F::dtype>::value) &&
          E::value == 2 && F::value == 2, // And it is matrix / matrix
      types::ndarray<
          typename __combined<typename E::dtype, typename F::dtype>::type,
          types::array<long, 2>>>::type
  dot(E const &e, F const &f)
  {
    types::ndarray<
        typename __combined<typename E::dtype, typename F::dtype>::type,
        types::array<long, 2>>
        out(types::array<long, 2>{{e.template shape<0>(),
                                   f.template shape<1>()}},
            0);
    for (long i = 0; i < out.template shape<0>(); i++)
      for (long j = 0; j < out.template shape<1>(); j++)
        for (long k = 0; k < e.template shape<1>(); k++)
          out[types::array<long, 2>{{i, j}}] +=
              e[types::array<long, 2>{{i, k}}] *
              f[types::array<long, 2>{{k, j}}];
    return out;
  }

  template <class E, class F>
  typename std::enable_if<
      (E::value >= 3 && F::value == 1), // And it is matrix / matrix
      types::ndarray<
          typename __combined<typename E::dtype, typename F::dtype>::type,
          types::array<long, E::value - 1>>>::type
  dot(E const &e, F const &f)
  {
    auto out = dot(
        e.reshape(types::array<long, 2>{{sutils::prod_head(e), f.size()}}), f);
    types::array<long, E::value - 1> out_shape;
    auto tmp = sutils::getshape(e);
    std::copy(tmp.begin(), tmp.end() - 1, out_shape.begin());
    return out.reshape(out_shape);
  }

  template <class E, class F>
  typename std::enable_if<
      (E::value >= 3 && F::value >= 2),
      types::ndarray<
          typename __combined<typename E::dtype, typename F::dtype>::type,
          types::array<long, E::value - 1>>>::type
  dot(E const &e, F const &f)
  {
    static_assert(E::value == 0, "not implemented yet");
  }
} // namespace numpy
PYTHONIC_NS_END

#endif
