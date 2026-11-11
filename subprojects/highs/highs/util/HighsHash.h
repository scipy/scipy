/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_UTIL_HASH_H_
#define HIGHS_UTIL_HASH_H_

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iterator>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include "util/HighsType.h"

#ifdef HIGHS_HAVE_BITSCAN_REVERSE
#include <intrin.h>
#pragma intrinsic(_BitScanReverse)
#ifdef _WIN64
#pragma intrinsic(_BitScanReverse64)
#pragma intrinsic(__popcnt64)
#else
#pragma intrinsic(__popcnt)
#endif
#endif

#if __GNUG__ && __GNUC__ < 5 && !defined(__clang__)
#define IS_TRIVIALLY_COPYABLE(T) __has_trivial_copy(T)
#else
#define IS_TRIVIALLY_COPYABLE(T) std::is_trivially_copyable<T>::value
#endif

template <typename T>
struct HighsHashable : std::integral_constant<bool, IS_TRIVIALLY_COPYABLE(T)> {
};

template <typename U, typename V>
struct HighsHashable<std::pair<U, V>>
    : public std::integral_constant<bool, HighsHashable<U>::value &&
                                              HighsHashable<V>::value> {};

template <typename U, typename V>
struct HighsHashable<std::tuple<U, V>> : public HighsHashable<std::pair<U, V>> {
};
template <typename U, typename V, typename W, typename... Args>
struct HighsHashable<std::tuple<U, V, W, Args...>>
    : public std::integral_constant<
          bool, HighsHashable<U>::value &&
                    HighsHashable<std::tuple<V, W, Args...>>::value> {};

struct HighsHashHelpers {
  using u8 = std::uint8_t;
  using i8 = std::int8_t;

  using u16 = std::uint16_t;
  using i16 = std::int16_t;

  using u32 = std::uint32_t;
  using i32 = std::int32_t;

  using u64 = std::uint64_t;
  using i64 = std::uint64_t;

  static constexpr u64 c[] = {
      u64{0xc8497d2a400d9551}, u64{0x80c8963be3e4c2f3}, u64{0x042d8680e260ae5b},
      u64{0x8a183895eeac1536}, u64{0xa94e9c75f80ad6de}, u64{0x7e92251dec62835e},
      u64{0x07294165cb671455}, u64{0x89b0f6212b0a4292}, u64{0x31900011b96bf554},
      u64{0xa44540f8eee2094f}, u64{0xce7ffd372e4c64fc}, u64{0x51c9d471bfe6a10f},
      u64{0x758c2a674483826f}, u64{0xf91a20abe63f8b02}, u64{0xc2a069024a1fcc6f},
      u64{0xd5bb18b70c5dbd59}, u64{0xd510adac6d1ae289}, u64{0x571d069b23050a79},
      u64{0x60873b8872933e06}, u64{0x780481cc19670350}, u64{0x7a48551760216885},
      u64{0xb5d68b918231e6ca}, u64{0xa7e5571699aa5274}, u64{0x7b6d309b2cfdcf01},
      u64{0x04e77c3d474daeff}, u64{0x4dbf099fd7247031}, u64{0x5d70dca901130beb},
      u64{0x9f8b5f0df4182499}, u64{0x293a74c9686092da}, u64{0xd09bdab6840f52b3},
      u64{0xc05d47f3ab302263}, u64{0x6b79e62b884b65d6}, u64{0xa581106fc980c34d},
      u64{0xf081b7145ea2293e}, u64{0xfb27243dd7c3f5ad}, u64{0x5211bf8860ea667f},
      u64{0x9455e65cb2385e7f}, u64{0x0dfaf6731b449b33}, u64{0x4ec98b3c6f5e68c7},
      u64{0x007bfd4a42ae936b}, u64{0x65c93061f8674518}, u64{0x640816f17127c5d1},
      u64{0x6dd4bab17b7c3a74}, u64{0x34d9268c256fa1ba}, u64{0x0b4d0c6b5b50d7f4},
      u64{0x30aa965bc9fadaff}, u64{0xc0ac1d0c2771404d}, u64{0xc5e64509abb76ef2},
      u64{0xd606b11990624a36}, u64{0x0d3f05d242ce2fb7}, u64{0x469a803cb276fe32},
      u64{0xa4a44d177a3e23f4}, u64{0xb9d9a120dcc1ca03}, u64{0x2e15af8165234a2e},
      u64{0x10609ba2720573d4}, u64{0xaa4191b60368d1d5}, u64{0x333dd2300bc57762},
      u64{0xdf6ec48f79fb402f}, u64{0x5ed20fcef1b734fa}, u64{0x4c94924ec8be21ee},
      u64{0x5abe6ad9d131e631}, u64{0xbe10136a522e602d}, u64{0x53671115c340e779},
      u64{0x9f392fe43e2144da}};

  /// mersenne prime 2^61 - 1
  static constexpr u64 M61() { return u64{0x1fffffffffffffff}; };

#ifdef HIGHS_HAVE_BUILTIN_CLZ
  static int log2i(uint64_t n) { return 63 - __builtin_clzll(n); }

  static int log2i(uint32_t n) { return 31 - __builtin_clz(n); }

  static int popcnt(uint64_t x) { return __builtin_popcountll(x); }

#elif defined(HIGHS_HAVE_BITSCAN_REVERSE)
  static int log2i(uint64_t n) {
    unsigned long result;
#ifdef _WIN64
    _BitScanReverse64(&result, n);
#else
    if (_BitScanReverse(&result, (n >> 32)))
      result += 32;
    else
      _BitScanReverse(&result, (n & 0xffffffffu));
#endif
    return static_cast<int>(result);
  }

  static int log2i(uint32_t n) {
    unsigned long result;
    _BitScanReverse(&result, static_cast<unsigned long>(n));
    return static_cast<int>(result);
  }

  static int popcnt(uint64_t x) {
#ifdef _WIN64
    return static_cast<int>(__popcnt64(x));
#else
    return __popcnt(x & 0xffffffffu) + __popcnt(x >> 32);
#endif
  }
#else
  // integer log2 algorithm without floating point arithmetic. It uses an
  // unrolled loop and requires few instructions that can be well optimized.
  static int log2i(uint64_t n) {
    int x = 0;

    auto log2Iteration = [&](int p) {
      if (n >= uint64_t{1} << p) {
        x += p;
        n >>= p;
      }
    };

    log2Iteration(32);
    log2Iteration(16);
    log2Iteration(8);
    log2Iteration(4);
    log2Iteration(2);
    log2Iteration(1);

    return x;
  }

  static int log2i(uint32_t n) {
    int x = 0;

    auto log2Iteration = [&](int p) {
      if (n >= 1u << p) {
        x += p;
        n >>= p;
      }
    };

    log2Iteration(16);
    log2Iteration(8);
    log2Iteration(4);
    log2Iteration(2);
    log2Iteration(1);

    return x;
  }

  static int popcnt(uint64_t x) {
    constexpr uint64_t m1 = 0x5555555555555555ull;
    constexpr uint64_t m2 = 0x3333333333333333ull;
    constexpr uint64_t m4 = 0x0f0f0f0f0f0f0f0full;
    constexpr uint64_t h01 = 0x0101010101010101ull;

    x -= (x >> 1) & m1;
    x = (x & m2) + ((x >> 2) & m2);
    x = (x + (x >> 4)) & m4;

    return (x * h01) >> 56;
  }

#endif

  /// compute a * b mod 2^61-1
  static u64 multiply_modM61(u64 a, u64 b) {
    u64 ahi = a >> 32;
    u64 bhi = b >> 32;
    u64 alo = a & 0xffffffffu;
    u64 blo = b & 0xffffffffu;

    // compute the different order terms with adicities 2^64, 2^32, 2^0
    u64 term_64 = ahi * bhi;
    u64 term_32 = ahi * blo + bhi * alo;
    u64 term_0 = alo * blo;

    // Partially reduce term_0 and term_32 modulo M61() individually to not deal
    // with a possible carry bit (thanks @https://github.com/WTFHCN for catching
    // the bug with this). We do not need to completely reduce by an additional
    // check for the range of the resulting term as this is done in the end in
    // any case and the reduced sizes do not cause troubles with the available
    // 64 bits.
    term_0 = (term_0 & M61()) + (term_0 >> 61);
    term_0 += ((term_32 >> 29) + (term_32 << 32)) & M61();

    // The lower 61 bits of term_0 are now the lower 61 bits of the result that
    // we need. Now extract the upper 61 of the result so that we can compute
    // the result of the multiplication modulo M61()
    u64 ab61 = (term_64 << 3) | (term_0 >> 61);

    // finally take the result modulo M61 which is computed by exploiting
    // that M61 is a mersenne prime, particularly, if a * b = q * 2^61 + r
    // then a * b = (q + r) (mod 2^61 - 1)
    u64 result = (term_0 & M61()) + ab61;
    if (result >= M61()) result -= M61();
    return result;
  }

  static u64 modexp_M61(u64 a, u64 e) {
    // the exponent need to be greater than zero
    assert(e > 0);
    u64 result = a;

    while (e != 1) {
      // square
      result = multiply_modM61(result, result);

      // multiply with a if exponent is odd
      if (e & 1) result = multiply_modM61(result, a);

      // shift to next bit
      e = e >> 1;
    }

    return result;
  }

  /// mersenne prime 2^31 - 1
  static constexpr u64 M31() { return u32{0x7fffffff}; };

  /// compute a * b mod 2^31-1
  static u32 multiply_modM31(u32 a, u32 b) {
    u64 result = u64(a) * u64(b);
    result = (result >> 31) + (result & M31());
    if (result >= M31()) result -= M31();
    return static_cast<u32>(result);
  }

  static u32 modexp_M31(u32 a, u64 e) {
    // the exponent need to be greater than zero
    assert(e > 0);
    u32 result = a;

    while (e != 1) {
      // square
      result = multiply_modM31(result, result);

      // multiply with a if exponent is odd
      if (e & 1) result = multiply_modM31(result, a);

      // shift to next bit
      e = e >> 1;
    }

    return result;
  }

  template <HighsInt k>
  static u64 pair_hash(u32 a, u32 b) {
    return (static_cast<u64>(a) + c[2 * k]) *
           (static_cast<u64>(b) + c[2 * k + 1]);
  }

  static void sparse_combine(u64& hash, HighsInt index, u64 value) {
    // we take each value of the sparse hash as coefficient for a polynomial
    // of the finite field modulo the mersenne prime 2^61-1 where the monomial
    // for a sparse entry has the degree of its index. We evaluate the
    // polynomial at a random constant. This allows to compute the hashes of
    // sparse vectors independently of each others nonzero contribution and
    // therefore allows to use the order of best access patterns for cache
    // performance. E.g. we can compute a strong hash value for parallel row and
    // column detection and only need to loop over the nonzeros once in
    // arbitrary order. This comes at the expense of more expensive hash
    // calculations as it would be more efficient to evaluate the polynomial
    // with horners scheme, but allows for parallelization and arbitrary order.
    // Since we have 64 random constants available, we slightly improve
    // the scheme by using a lower degree polynomial with 64 variables
    // which we evaluate at the random vector of 64.

    // make sure input value is never zero and at most 61bits are used
    value = ((value << 1) & M61()) | 1;

    // make sure that the constant has at most 61 bits, as otherwise the modulo
    // algorithm for multiplication mod M61 might not work properly due to
    // overflow
    u64 a = c[index & 63] & M61();
    u64 degree = (static_cast<u64>(index) >> 6) + 1;

    hash += multiply_modM61(value, modexp_M61(a, degree));
    hash = (hash >> 61) + (hash & M61());
    if (hash >= M61()) hash -= M61();
    assert(hash < M61());
  }

  static void sparse_inverse_combine(u64& hash, HighsInt index, u64 value) {
    // same hash algorithm as sparse_combine(), but for updating a hash value to
    // the state before it was changed with a call to sparse_combine(). This is
    // easily possible as the hash value just uses finite field arithmetic. We
    // can simply add the additive inverse of the previous hash value. This is a
    // very useful routine for symmetry detection. During partition refinement
    // the hashes do not need to be recomputed but can be updated with this
    // procedure.

    // make sure input value is never zero and at most 61bits are used
    value = ((value << 1) & M61()) | 1;

    u64 a = c[index & 63] & M61();
    u64 degree = (static_cast<u64>(index) >> 6) + 1;
    // add the additive inverse (M61() - hashvalue) instead of the hash value
    // itself
    hash += M61() - multiply_modM61(value, modexp_M61(a, degree));
    hash = (hash >> 61) + (hash & M61());
    if (hash >= M61()) hash -= M61();
    assert(hash < M61());
  }

  /// overload that is not taking a value and saves one multiplication call
  /// useful for sparse hashing of bit vectors
  static void sparse_combine(u64& hash, HighsInt index) {
    u64 a = c[index & 63] & M61();
    u64 degree = (static_cast<u64>(index) >> 6) + 1;

    hash += modexp_M61(a, degree);
    hash = (hash >> 61) + (hash & M61());
    if (hash >= M61()) hash -= M61();
    assert(hash < M61());
  }

  /// overload that is not taking a value and saves one multiplication call
  /// useful for sparse hashing of bit vectors
  static void sparse_inverse_combine(u64& hash, HighsInt index) {
    // same hash algorithm as sparse_combine(), but for updating a hash value to
    // the state before it was changed with a call to sparse_combine(). This is
    // easily possible as the hash value just uses finite field arithmetic. We
    // can simply add the additive inverse of the previous hash value. This is a
    // very useful routine for symmetry detection. During partition refinement
    // the hashes do not need to be recomputed but can be updated with this
    // procedure.

    u64 a = c[index & 63] & M61();
    u64 degree = (static_cast<u64>(index) >> 6) + 1;
    // add the additive inverse (M61() - hashvalue) instead of the hash value
    // itself
    hash += M61() - modexp_M61(a, degree);
    hash = (hash >> 61) + (hash & M61());
    if (hash >= M61()) hash -= M61();
    assert(hash < M61());
  }

  static void sparse_combine32(u32& hash, HighsInt index, u64 value) {
    // we take each value of the sparse hash as coefficient for a polynomial
    // of the finite field modulo the mersenne prime 2^61-1 where the monomial
    // for a sparse entry has the degree of its index. We evaluate the
    // polynomial at a random constant. This allows to compute the hashes of
    // sparse vectors independently of each others nonzero contribution and
    // therefore allows to use the order of best access patterns for cache
    // performance. E.g. we can compute a strong hash value for parallel row and
    // column detection and only need to loop over the nonzeros once in
    // arbitrary order. This comes at the expense of more expensive hash
    // calculations as it would be more efficient to evaluate the polynomial
    // with horners scheme, but allows for parallelization and arbitrary order.
    // Since we have 16 random constants available, we slightly improve
    // the scheme by using a lower degree polynomial with 16 variables
    // which we evaluate at the random vector of 16.

    // make sure input value is never zero and at most 31bits are used
    value = (pair_hash<0>(static_cast<u32>(value), value >> 32) >> 33) | 1;

    // make sure that the constant has at most 31 bits, as otherwise the modulo
    // algorithm for multiplication mod M31 might not work properly due to
    // overflow
    u32 a = static_cast<u32>(c[index & 63] & M31());
    u64 degree = (static_cast<u64>(index) >> 6) + 1;

    u64 result = hash;
    result += multiply_modM31(static_cast<u32>(value), modexp_M31(a, degree));
    result = (result >> 31) + (result & M31());
    if (result >= M31()) result -= M31();
    assert(result < M31());
    hash = static_cast<u32>(result);
  }

  static void sparse_inverse_combine32(u32& hash, HighsInt index, u64 value) {
    // same hash algorithm as sparse_combine(), but for updating a hash value to
    // the state before it was changed with a call to sparse_combine(). This is
    // easily possible as the hash value just uses finite field arithmetic. We
    // can simply add the additive inverse of the previous hash value. This is a
    // very useful routine for symmetry detection. During partition refinement
    // the hashes do not need to be recomputed but can be updated with this
    // procedure.

    // make sure input value is never zero and at most 31bits are used
    value = (pair_hash<0>(static_cast<u32>(value), value >> 32) >> 33) | 1;

    u32 a = static_cast<u32>(c[index & 63] & M31());
    u64 degree = (static_cast<u64>(index) >> 6) + 1;
    // add the additive inverse (M31() - hashvalue) instead of the hash value
    // itself
    u64 result = hash;
    result +=
        M31() - multiply_modM31(static_cast<u32>(value), modexp_M31(a, degree));
    result = (result >> 31) + (result & M31());
    if (result >= M31()) result -= M31();
    assert(result < M31());
    hash = static_cast<u32>(result);
  }

  static constexpr u64 fibonacci_muliplier() { return u64{0x9e3779b97f4a7c15}; }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value, int>::type = 0>
  static u64 vector_hash(const T* vals, size_t numvals) {
    std::array<u32, 2> pair{};
    u64 hash = 0;
    HighsInt k = 0;

    const char* dataptr = (const char*)vals;
    const char* dataend = (const char*)(vals + numvals);

    while (dataptr != dataend) {
      using std::size_t;
      size_t numBytes = std::min(size_t(dataend - dataptr), size_t{256});
      size_t numPairs = (numBytes + 7) / 8;
      size_t lastPairBytes = numBytes - (numPairs - 1) * 8;
      u64 chunkhash[] = {u64{0}, u64{0}};

#define HIGHS_VECHASH_CASE_N(N, B)                         \
  std::memcpy(&pair[0], dataptr, B);                       \
  chunkhash[N & 1] += pair_hash<32 - N>(pair[0], pair[1]); \
  dataptr += B;

      switch (numPairs) {
        case 32:
          if (hash != 0) {
            // make sure hash is reduced mod M61() before multiplying with the
            // next random constant. For vectors at most 240 bytes we never
            // get here and only use the fast pair hashing scheme
            // for vectors with 240 bytes to 256 bytes we do have the one
            // additional check for hash != 0 above which will return false
            // and only for longer vectors we ever reduce modulo M61
            if (hash >= M61()) hash -= M61();
            hash = multiply_modM61(hash, c[(k++) & 63] & M61());
          }
          HIGHS_VECHASH_CASE_N(32, 8)
          // fall through
        case 31:
          HIGHS_VECHASH_CASE_N(31, 8)
          // fall through
        case 30:
          HIGHS_VECHASH_CASE_N(30, 8)
          // fall through
        case 29:
          HIGHS_VECHASH_CASE_N(29, 8)
          // fall through
        case 28:
          HIGHS_VECHASH_CASE_N(28, 8)
          // fall through
        case 27:
          HIGHS_VECHASH_CASE_N(27, 8)
          // fall through
        case 26:
          HIGHS_VECHASH_CASE_N(26, 8)
          // fall through
        case 25:
          HIGHS_VECHASH_CASE_N(25, 8)
          // fall through
        case 24:
          HIGHS_VECHASH_CASE_N(24, 8)
          // fall through
        case 23:
          HIGHS_VECHASH_CASE_N(23, 8)
          // fall through
        case 22:
          HIGHS_VECHASH_CASE_N(22, 8)
          // fall through
        case 21:
          HIGHS_VECHASH_CASE_N(21, 8)
          // fall through
        case 20:
          HIGHS_VECHASH_CASE_N(20, 8)
          // fall through
        case 19:
          HIGHS_VECHASH_CASE_N(19, 8)
          // fall through
        case 18:
          HIGHS_VECHASH_CASE_N(18, 8)
          // fall through
        case 17:
          HIGHS_VECHASH_CASE_N(17, 8)
          // fall through
        case 16:
          HIGHS_VECHASH_CASE_N(16, 8)
          // fall through
        case 15:
          HIGHS_VECHASH_CASE_N(15, 8)
          // fall through
        case 14:
          HIGHS_VECHASH_CASE_N(14, 8)
          // fall through
        case 13:
          HIGHS_VECHASH_CASE_N(13, 8)
          // fall through
        case 12:
          HIGHS_VECHASH_CASE_N(12, 8)
          // fall through
        case 11:
          HIGHS_VECHASH_CASE_N(11, 8)
          // fall through
        case 10:
          HIGHS_VECHASH_CASE_N(10, 8)
          // fall through
        case 9:
          HIGHS_VECHASH_CASE_N(9, 8)
          // fall through
        case 8:
          HIGHS_VECHASH_CASE_N(8, 8)
          // fall through
        case 7:
          HIGHS_VECHASH_CASE_N(7, 8)
          // fall through
        case 6:
          HIGHS_VECHASH_CASE_N(6, 8)
          // fall through
        case 5:
          HIGHS_VECHASH_CASE_N(5, 8)
          // fall through
        case 4:
          HIGHS_VECHASH_CASE_N(4, 8)
          // fall through
        case 3:
          HIGHS_VECHASH_CASE_N(3, 8)
          // fall through
        case 2:
          HIGHS_VECHASH_CASE_N(2, 8)
          // fall through
        case 1:
          HIGHS_VECHASH_CASE_N(1, lastPairBytes)
      }

      hash += (chunkhash[0] >> 3) ^ (chunkhash[1] >> 32);
    }

#undef HIGHS_VECHASH_CASE_N

    return hash * fibonacci_muliplier();
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value &&
                                        (sizeof(T) <= 8) && (sizeof(T) >= 1),
                                    int>::type = 0>
  static u64 hash(const T& val) {
    std::array<u32, 2> bytes;
    if (sizeof(T) < 4) bytes[0] = 0;
    if (sizeof(T) < 8) bytes[1] = 0;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return pair_hash<1>(bytes[0], bytes[1]) ^
           pair_hash<0>(bytes[0], bytes[1]) >> 32;
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value &&
                                        (sizeof(T) >= 9) && (sizeof(T) <= 16),
                                    int>::type = 0>
  static u64 hash(const T& val) {
    std::array<u32, 4> bytes;
    if (sizeof(T) < 12) bytes[2] = 0;
    if (sizeof(T) < 16) bytes[3] = 0;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return (pair_hash<0>(bytes[0], bytes[1]) ^
            (pair_hash<1>(bytes[2], bytes[3]) >> 32)) *
           fibonacci_muliplier();
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value &&
                                        (sizeof(T) >= 17) && (sizeof(T) <= 24),
                                    int>::type = 0>
  static u64 hash(const T& val) {
    std::array<u32, 6> bytes;
    if (sizeof(T) < 20) bytes[4] = 0;
    if (sizeof(T) < 24) bytes[5] = 0;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return (pair_hash<0>(bytes[0], bytes[1]) ^
            ((pair_hash<1>(bytes[2], bytes[3]) +
              pair_hash<2>(bytes[4], bytes[5])) >>
             32)) *
           fibonacci_muliplier();
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value &&
                                        (sizeof(T) >= 25) && (sizeof(T) <= 32),
                                    int>::type = 0>
  static u64 hash(const T& val) {
    std::array<u32, 8> bytes;
    if (sizeof(T) < 28) bytes[6] = 0;
    if (sizeof(T) < 32) bytes[7] = 0;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return ((pair_hash<0>(bytes[0], bytes[1]) +
             pair_hash<1>(bytes[2], bytes[3])) ^
            ((pair_hash<2>(bytes[4], bytes[5]) +
              pair_hash<3>(bytes[6], bytes[7])) >>
             32)) *
           fibonacci_muliplier();
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value &&
                                        (sizeof(T) >= 33) && (sizeof(T) <= 40),
                                    int>::type = 0>
  static u64 hash(const T& val) {
    std::array<u32, 10> bytes;
    if (sizeof(T) < 36) bytes[8] = 0;
    if (sizeof(T) < 40) bytes[9] = 0;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return ((pair_hash<0>(bytes[0], bytes[1]) +
             pair_hash<1>(bytes[2], bytes[3])) ^
            ((pair_hash<2>(bytes[4], bytes[5]) +
              pair_hash<3>(bytes[6], bytes[7]) +
              pair_hash<4>(bytes[8], bytes[9])) >>
             32)) *
           fibonacci_muliplier();
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value &&
                                        (sizeof(T) >= 41) && (sizeof(T) <= 48),
                                    int>::type = 0>
  static u64 hash(const T& val) {
    std::array<u32, 12> bytes;
    if (sizeof(T) < 44) bytes[10] = 0;
    if (sizeof(T) < 48) bytes[11] = 0;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return ((pair_hash<0>(bytes[0], bytes[1]) +
             pair_hash<1>(bytes[2], bytes[3]) +
             pair_hash<2>(bytes[4], bytes[5])) ^
            ((pair_hash<3>(bytes[6], bytes[7]) +
              pair_hash<4>(bytes[8], bytes[9]) +
              pair_hash<5>(bytes[10], bytes[11])) >>
             32)) *
           fibonacci_muliplier();
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value &&
                                        (sizeof(T) >= 49) && (sizeof(T) <= 56),
                                    int>::type = 0>
  static u64 hash(const T& val) {
    std::array<u32, 14> bytes;
    if (sizeof(T) < 52) bytes[12] = 0;
    if (sizeof(T) < 56) bytes[13] = 0;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return ((pair_hash<0>(bytes[0], bytes[1]) +
             pair_hash<1>(bytes[2], bytes[3]) +
             pair_hash<2>(bytes[4], bytes[5])) ^
            ((pair_hash<3>(bytes[6], bytes[7]) +
              pair_hash<4>(bytes[8], bytes[9]) +
              pair_hash<5>(bytes[10], bytes[11]) +
              pair_hash<6>(bytes[12], bytes[13])) >>
             32)) *
           fibonacci_muliplier();
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value &&
                                        (sizeof(T) >= 57) && (sizeof(T) <= 64),
                                    int>::type = 0>
  static u64 hash(const T& val) {
    std::array<u32, 16> bytes;
    if (sizeof(T) < 60) bytes[14] = 0;
    if (sizeof(T) < 64) bytes[15] = 0;
    std::memcpy(&bytes[0], &val, sizeof(T));
    return ((pair_hash<0>(bytes[0], bytes[1]) +
             pair_hash<1>(bytes[2], bytes[3]) +
             pair_hash<2>(bytes[4], bytes[5]) +
             pair_hash<3>(bytes[6], bytes[7])) ^
            ((pair_hash<4>(bytes[8], bytes[9]) +
              pair_hash<5>(bytes[10], bytes[11]) +
              pair_hash<6>(bytes[12], bytes[13]) +
              pair_hash<7>(bytes[14], bytes[15])) >>
             32)) *
           fibonacci_muliplier();
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value && (sizeof(T) > 64),
                                    int>::type = 0>
  static u64 hash(const T& val) {
    return vector_hash(&val, 1);
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value, int>::type = 0>
  static u64 hash(const std::vector<T>& val) {
    return vector_hash(val.data(), val.size());
  }

  template <typename T, typename std::enable_if<
                            std::is_same<decltype(*reinterpret_cast<T*>(0) ==
                                                  *reinterpret_cast<T*>(0)),
                                         bool>::value,
                            int>::type = 0>
  static bool equal(const T& a, const T& b) {
    return a == b;
  }

  template <typename T,
            typename std::enable_if<HighsHashable<T>::value, int>::type = 0>
  static bool equal(const std::vector<T>& a, const std::vector<T>& b) {
    if (a.size() != b.size()) return false;
    return std::memcmp(a.data(), b.data(), sizeof(T) * a.size()) == 0;
  }

  static constexpr double golden_ratio_reciprocal() {
    return 0.61803398874989484;
  }

  static u32 double_hash_code(double val) {
    // we multiply by some irrational number, so that the buckets in which we
    // put the real numbers do not break on a power of two pattern. E.g.
    // consider the use case for detecting parallel rows when we have two
    // parallel rows scaled to have their largest coefficient 1.0 and another
    // coefficient which is 0.5
    // +- epsilon. Clearly we want to detect those rows as parallel and give
    // them the same hash value for small enough epsilon. The exponent,
    // however will switch to -2 for the value just below 0.5 and the hashcodes
    // will differ. when multiplying with the reciprocal of the golden ratio the
    // exact 0.5 will yield 0.30901699437494742 and 0.5 - 1e-9 will yield
    // 0.3090169937569134 which has the same exponent and matches in the most
    // significant bits. Hence it yields the same hashcode. Obviously there will
    // now be different values which exhibit the same pattern as the 0.5 case,
    // but they do not have a small denominator like 1/2 in their rational
    // representation but are power of two multiples of the golden ratio and
    // therefore irrational, which we do not expect in non-artificial input
    // data.
    int exponent;
    double hashbits = std::frexp(val * golden_ratio_reciprocal(), &exponent);

    // some extra casts to be more verbose about what is happening.
    // We want the exponent to use only 16bits so that the remaining 16 bits
    // are used for the most significant bits of the mantissa and the sign bit.
    // casting to unsigned 16bits first ensures that the value after the cast is
    // defined to be UINT16_MAX - |exponent| when the exponent is negative.
    // casting the exponent to a uint32_t directly would give wrong promotion
    // of negative exponents as UINT32_MAX - |exponent| and take up to many bits
    // or possibly lose information after the 16 bit shift. For the mantissa we
    // take the 15 most significant bits, even though we could squeeze out a few
    // more of the exponent. We don't need more bits as this would make the
    // buckets very small and might miss more values that are equal within
    // epsilon. Therefore the most significant 15 bits of the mantissa and the
    // sign is encoded in the 16 lower bits of the hashcode and the upper 16bits
    // encode the sign and value of the exponent.
    u32 hashvalue = (u32)(u16)(i16)exponent;
    hashvalue = (hashvalue << 16) | (u32)(u16)(i16)std::ldexp(hashbits, 15);

    return hashvalue;
  }
};

struct HighsHasher {
  template <typename T>
  size_t operator()(const T& x) const {
    return HighsHashHelpers::hash(x);
  }
};

struct HighsVectorHasher {
  template <typename T>
  size_t operator()(const std::vector<T>& vec) const {
    return HighsHashHelpers::vector_hash(vec.data(), vec.size());
  }
};

struct HighsVectorEqual {
  template <typename T>
  bool operator()(const std::vector<T>& vec1,
                  const std::vector<T>& vec2) const {
    if (vec1.size() != vec2.size()) return false;
    return std::equal(vec1.begin(), vec1.end(), vec2.begin());
  }
};

template <typename K, typename V = void>
struct HighsHashTableEntry {
 private:
  K key_;
  V value_;

 public:
  HighsHashTableEntry(HighsHashTableEntry<K, V>&&) = default;
  HighsHashTableEntry(const HighsHashTableEntry<K, V>&) = default;
  ~HighsHashTableEntry() = default;
  HighsHashTableEntry() = default;
  HighsHashTableEntry<K, V>& operator=(HighsHashTableEntry<K, V>&&) = default;
  HighsHashTableEntry<K, V>& operator=(const HighsHashTableEntry<K, V>&) =
      default;

  // add a constructor to pass an argument to initialize the key with a value
  // and the value as default
  // the enable if statement makes sure this overload is never selected
  // when the type of the single argument is HighsHashTableEntry<K,V> so that
  // the default move and copy constructors are preferred when they match
  // and this is only used to initialize the key type from a single argument.
  template <
      typename K_,
      typename std::enable_if<
          !std::is_same<typename std::remove_cv<
                            typename std::remove_reference<K_>::type>::type,
                        HighsHashTableEntry<K, V>>::value,
          int>::type = 0>
  HighsHashTableEntry(K_&& k) : key_(std::forward<K_>(k)), value_() {}

  template <typename K_, typename V_>
  HighsHashTableEntry(K_&& k, V_&& v)
      : key_(std::forward<K_>(k)), value_(std::forward<V_>(v)) {}

  const K& key() const { return key_; }
  const V& value() const { return value_; }
  V& value() { return value_; }

  template <typename Func>
  auto forward(Func&& f) -> decltype(f(key_, value_)) {
    const K& keyRef = key_;
    return f(keyRef, value_);
  }

  template <typename Func>
  auto forward(Func&& f) const -> decltype(f(key_)) {
    const K& keyRef = key_;
    return f(keyRef);
  }

  template <typename Func>
  auto forward(Func&& f) const -> decltype(f(key_, value_)) {
    return f(key_, value_);
  }
};

template <typename T>
struct HighsHashTableEntry<T, void> {
 private:
  T value_;

 public:
  HighsHashTableEntry(HighsHashTableEntry<T, void>&&) = default;
  HighsHashTableEntry(const HighsHashTableEntry<T, void>&) = default;
  ~HighsHashTableEntry() = default;
  HighsHashTableEntry() = default;
  HighsHashTableEntry<T, void>& operator=(HighsHashTableEntry<T, void>&&) =
      default;
  HighsHashTableEntry<T, void>& operator=(const HighsHashTableEntry<T, void>&) =
      default;

  // Add a constructor to accept an arbitrary argument pack for initialize the
  // underlying value of type T. The enable if statement makes sure this
  // overload is never selected when the type of the single argument is
  // HighsHashTableEntry<T,void> so that the default move and copy constructors
  // are preferred when they match and this is only used to initialize the value
  // of type from a set of arguments which are properly forwarded.
  // The std::tuple usage in enable_if is a work-around to make the statement
  // legal when multiple arguments are passed in, since std::is_same expects a
  // single type. In that case is_same will obviously return false and the
  // overload is appropriate to initialize the value_ with multiple forwarded
  // arguments.
  template <typename... Args,
            typename std::enable_if<
                !std::is_same<
                    std::tuple<typename std::remove_cv<
                        typename std::remove_reference<Args>::type>::type...>,
                    std::tuple<HighsHashTableEntry<T, void>>>::value,
                int>::type = 0>
  HighsHashTableEntry(Args&&... args) : value_(std::forward<Args>(args)...) {}

  const T& key() const { return value_; }
  const T& value() const { return value_; }

  template <typename Func>
  auto forward(Func&& f) -> decltype(f(value_)) {
    return f(value_);
  }

  template <typename Func>
  auto forward(Func&& f) const -> decltype(f(value_)) {
    return f(value_);
  }
};

template <typename K, typename V = void>
class HighsHashTable {
  struct OpNewDeleter {
    void operator()(void* ptr) { ::operator delete(ptr); }
  };

 public:
  using u8 = std::uint8_t;
  using i8 = std::int8_t;

  using u16 = std::uint16_t;
  using i16 = std::int16_t;

  using u32 = std::uint32_t;
  using i32 = std::int32_t;

  using u64 = std::uint64_t;
  using i64 = std::int64_t;

  using Entry = HighsHashTableEntry<K, V>;
  using KeyType = K;
  using ValueType =
      typename std::remove_reference<decltype(Entry().value())>::type;

  std::unique_ptr<Entry, OpNewDeleter> entries;
  std::unique_ptr<u8[]> metadata;
  u64 tableSizeMask;
  u64 numHashShift;
  u64 numElements = 0;

  template <typename IterType>
  class HashTableIterator {
    u8* pos;
    u8* end;
    Entry* entryEnd;

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = IterType;
    using pointer = IterType*;
    using reference = IterType&;
    using iterator_category = std::forward_iterator_tag;
    HashTableIterator(u8* pos_, u8* end_, Entry* entryEnd_)
        : pos(pos_), end(end_), entryEnd(entryEnd_) {}
    HashTableIterator() = default;

    HashTableIterator<IterType> operator++(int) {
      // postfix
      HashTableIterator<IterType> oldpos = *this;
      for (++pos; pos != end; ++pos)
        if ((*pos) & 0x80u) break;

      return oldpos;
    }

    HashTableIterator<IterType>& operator++() {
      // prefix
      for (++pos; pos != end; ++pos)
        if ((*pos) & 0x80u) break;

      return *this;
    }

    reference operator*() const { return *(entryEnd - (end - pos)); }
    pointer operator->() const { return (entryEnd - (end - pos)); }
    HashTableIterator<IterType> operator+(difference_type v) const {
      for (difference_type k = 0; k != v; ++k) ++(*this);
    }

    bool operator==(const HashTableIterator<IterType>& rhs) const {
      return pos == rhs.pos;
    }
    bool operator!=(const HashTableIterator<IterType>& rhs) const {
      return pos != rhs.pos;
    }
  };

  using const_iterator = HashTableIterator<const Entry>;
  using iterator = HashTableIterator<Entry>;

  HighsHashTable() { makeEmptyTable(128); }
  HighsHashTable(u64 minCapacity) {
    u64 initCapacity = u64{1} << (u64)std::ceil(std::log2(std::max(
                           128.0, 8 * static_cast<double>(minCapacity) / 7)));
    makeEmptyTable(initCapacity);
  }

  HighsHashTable(const HighsHashTable& hashTable)
      : tableSizeMask(hashTable.tableSizeMask),
        numHashShift(hashTable.numHashShift),
        numElements(hashTable.numElements) {
    u64 capacity = tableSizeMask + 1;
    metadata = decltype(metadata)(new u8[capacity]);
    entries =
        decltype(entries)((Entry*)::operator new(sizeof(Entry) * capacity));

    std::copy(hashTable.metadata.get(), hashTable.metadata.get() + capacity,
              metadata.get());
    std::copy(hashTable.entries.get(), hashTable.entries.get() + capacity,
              entries.get());
  }

  iterator end() {
    u64 capacity = tableSizeMask + 1;
    return iterator{metadata.get() + capacity, metadata.get() + capacity,
                    entries.get() + capacity};
  };

  const_iterator end() const {
    u64 capacity = tableSizeMask + 1;
    return const_iterator{metadata.get() + capacity, metadata.get() + capacity,
                          entries.get() + capacity};
  };

  const_iterator begin() const {
    if (numElements == 0) return end();
    u64 capacity = tableSizeMask + 1;
    const_iterator iter{metadata.get(), metadata.get() + capacity,
                        entries.get() + capacity};
    if (!occupied(metadata[0])) ++iter;

    return iter;
  };

  iterator begin() {
    if (numElements == 0) return end();
    u64 capacity = tableSizeMask + 1;
    iterator iter{metadata.get(), metadata.get() + capacity,
                  entries.get() + capacity};
    if (!occupied(metadata[0])) ++iter;

    return iter;
  };

 private:
  u8 toMetadata(u64 hash) const { return (hash >> numHashShift) | 0x80u; }

  static constexpr u64 maxDistance() { return 127; }

  void makeEmptyTable(u64 capacity) {
    tableSizeMask = capacity - 1;
    numHashShift = 64 - HighsHashHelpers::log2i(capacity);
    assert(capacity == (u64{1} << (64 - numHashShift)));
    numElements = 0;

    metadata = decltype(metadata)(new u8[capacity]{});
    entries =
        decltype(entries)((Entry*)::operator new(sizeof(Entry) * capacity));
  }

  bool occupied(u8 meta) const { return meta & 0x80; }

  u64 distanceFromIdealSlot(u64 pos) const {
    // we store 7 bits of the hash in the metadata. Assuming a decent
    // hashfunction it is practically never happening that an item travels more
    // than 127 slots from its ideal position, therefore, we can compute the
    // distance from the ideal position just as it would normally be done
    // assuming there is at most one overflow. Consider using 3 bits which gives
    // values from 0 to 7. When an item is at a position with lower bits 7 and
    // is placed 3 positions after its ideal position, the lower bits of the
    // hash value will overflow and yield the value 2. With the assumption that
    // an item never cycles through one full cycle of the range 0 to 7, its
    // position would never be placed in a position with lower bits 7 other than
    // its ideal position. This allows us to compute the distance from its ideal
    // position by simply ignoring an overflow. In our case the correct answer
    // would be 3, but we get (2 - 7)=-5. This, however, is the correct result 3
    // when promoting to an unsigned value and looking at the lower 3 bits.

    return ((pos - metadata[pos])) & 0x7f;
  }

  void growTable() {
    decltype(entries) oldEntries = std::move(entries);
    decltype(metadata) oldMetadata = std::move(metadata);
    u64 oldCapactiy = tableSizeMask + 1;

    makeEmptyTable(2 * oldCapactiy);

    for (u64 i = 0; i != oldCapactiy; ++i)
      if (occupied(oldMetadata[i])) insert(std::move(oldEntries.get()[i]));
  }

  void shrinkTable() {
    decltype(entries) oldEntries = std::move(entries);
    decltype(metadata) oldMetadata = std::move(metadata);
    u64 oldCapactiy = tableSizeMask + 1;

    makeEmptyTable(oldCapactiy / 2);

    for (u64 i = 0; i != oldCapactiy; ++i)
      if (occupied(oldMetadata[i])) insert(std::move(oldEntries.get()[i]));
  }

  bool findPosition(const KeyType& key, u8& meta, u64& startPos, u64& maxPos,
                    u64& pos) const {
    u64 hash = HighsHashHelpers::hash(key);
    startPos = hash >> numHashShift;
    maxPos = (startPos + maxDistance()) & tableSizeMask;
    meta = toMetadata(hash);

    const Entry* entryArray = entries.get();
    pos = startPos;
    do {
      if (!occupied(metadata[pos])) return false;
      if (metadata[pos] == meta &&
          HighsHashHelpers::equal(key, entryArray[pos].key()))
        return true;

      u64 currentDistance = (pos - startPos) & tableSizeMask;
      if (currentDistance > distanceFromIdealSlot(pos)) return false;

      pos = (pos + 1) & tableSizeMask;
    } while (pos != maxPos);

    return false;
  }

 public:
  void clear() {
    if (numElements) {
      u64 capacity = tableSizeMask + 1;
      if (!std::is_trivially_destructible<Entry>::value) {
        for (u64 i = 0; i < capacity; ++i)
          if (occupied(metadata[i])) entries.get()[i].~Entry();
      }
      if (capacity == 128) {
        std::memset(metadata.get(), 0, 128);
        numElements = 0;
      } else
        makeEmptyTable(128);
    }
  }

  const ValueType* find(const KeyType& key) const {
    u64 pos, startPos, maxPos;
    u8 meta;
    if (findPosition(key, meta, startPos, maxPos, pos))
      return &(entries.get()[pos].value());

    return nullptr;
  }

  ValueType* find(const KeyType& key) {
    u64 pos, startPos, maxPos;
    u8 meta;
    if (findPosition(key, meta, startPos, maxPos, pos))
      return &(entries.get()[pos].value());

    return nullptr;
  }

  ValueType& operator[](const KeyType& key) {
    Entry* entryArray = entries.get();
    u64 pos, startPos, maxPos;
    u8 meta;
    if (findPosition(key, meta, startPos, maxPos, pos))
      return entryArray[pos].value();

    if (numElements == ((tableSizeMask + 1) * 7) / 8 || pos == maxPos) {
      growTable();
      return (*this)[key];
    }

    using std::swap;
    ValueType& insertLocation = entryArray[pos].value();
    Entry entry(key);
    ++numElements;

    do {
      if (!occupied(metadata[pos])) {
        metadata[pos] = meta;
        new (&entryArray[pos]) Entry{std::move(entry)};
        return insertLocation;
      }

      u64 currentDistance = (pos - startPos) & tableSizeMask;
      u64 distanceOfCurrentOccupant = distanceFromIdealSlot(pos);
      if (currentDistance > distanceOfCurrentOccupant) {
        // steal the position
        swap(entry, entryArray[pos]);
        swap(meta, metadata[pos]);

        startPos = (pos - distanceOfCurrentOccupant) & tableSizeMask;
        maxPos = (startPos + maxDistance()) & tableSizeMask;
      }
      pos = (pos + 1) & tableSizeMask;
    } while (pos != maxPos);

    growTable();
    insert(std::move(entry));
    return (*this)[key];
  }

  template <typename... Args>
  bool insert(Args&&... args) {
    Entry entry(std::forward<Args>(args)...);

    u64 pos, startPos, maxPos;
    u8 meta;
    if (findPosition(entry.key(), meta, startPos, maxPos, pos)) return false;

    if (numElements == ((tableSizeMask + 1) * 7) / 8 || pos == maxPos) {
      growTable();
      return insert(std::move(entry));
    }

    using std::swap;
    Entry* entryArray = entries.get();
    ++numElements;

    do {
      if (!occupied(metadata[pos])) {
        metadata[pos] = meta;
        new (&entryArray[pos]) Entry{std::move(entry)};
        return true;
      }

      u64 currentDistance = (pos - startPos) & tableSizeMask;
      u64 distanceOfCurrentOccupant = distanceFromIdealSlot(pos);
      if (currentDistance > distanceOfCurrentOccupant) {
        // steal the position
        swap(entry, entryArray[pos]);
        swap(meta, metadata[pos]);

        startPos = (pos - distanceOfCurrentOccupant) & tableSizeMask;
        maxPos = (startPos + maxDistance()) & tableSizeMask;
      }
      pos = (pos + 1) & tableSizeMask;
    } while (pos != maxPos);

    growTable();
    insert(std::move(entry));
    return true;
  }

  bool erase(const KeyType& key) {
    u64 pos, startPos, maxPos;
    u8 meta;
    if (!findPosition(key, meta, startPos, maxPos, pos)) return false;
    // delete element at position pos
    Entry* entryArray = entries.get();
    entryArray[pos].~Entry();
    metadata[pos] = 0;

    // retain at least a quarter of slots occupied, otherwise shrink the table
    // if its not at its minimum size already
    --numElements;
    u64 capacity = tableSizeMask + 1;
    if (capacity != 128 && numElements < capacity / 4) {
      shrinkTable();
      return true;
    }

    // shift elements after pos backwards
    while (true) {
      u64 shift = (pos + 1) & tableSizeMask;
      if (!occupied(metadata[shift])) return true;

      u64 dist = distanceFromIdealSlot(shift);
      if (dist == 0) return true;

      entryArray[pos] = std::move(entryArray[shift]);
      metadata[pos] = metadata[shift];
      metadata[shift] = 0;
      pos = shift;
    }
  }

  i64 size() const { return numElements; }

  HighsHashTable(HighsHashTable<K, V>&&) = default;
  HighsHashTable<K, V>& operator=(HighsHashTable<K, V>&&) = default;

  ~HighsHashTable() {
    if (!std::is_trivially_destructible<Entry>::value && metadata) {
      u64 capacity = tableSizeMask + 1;
      for (u64 i = 0; i < capacity; ++i) {
        if (occupied(metadata[i])) entries.get()[i].~Entry();
      }
    }
  }
};

#endif
