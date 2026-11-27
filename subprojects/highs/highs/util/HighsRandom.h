/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsRandom.h
 * @brief Random number generators for HiGHS
 */
#ifndef UTIL_HIGHSRANDOM_H_
#define UTIL_HIGHSRANDOM_H_

#include <cassert>

#include "util/HighsHash.h"

/**
 * @brief Class for HiGHS random number generators
 */
class HighsRandom {
 private:
  uint32_t drawUniform(uint32_t sup, int nbits) {
    // draw uniformly in interval [0,sup) where nbits is the maximal number of
    // bits the results can have we draw random numbers with nbits many bits
    // until we get one that is in the desired range we first use all available
    // output functions for the same state before we advance the state again. We
    // expect nbits to be at most 32 for this 32 bit version.
    assert(sup <= uint32_t{1} << nbits);
    while (true) {
      advance();
      uint32_t lo = static_cast<uint32_t>(state);
      uint32_t hi = state >> 32;

      uint32_t val;

#define HIGHS_RAND_TRY_OUTPUT(n)                                        \
  val = static_cast<uint32_t>(HighsHashHelpers::pair_hash<n>(lo, hi) >> \
                              (64 - nbits));                            \
  if (val < sup) return val;

      HIGHS_RAND_TRY_OUTPUT(0);
      HIGHS_RAND_TRY_OUTPUT(1);
      HIGHS_RAND_TRY_OUTPUT(2);
      HIGHS_RAND_TRY_OUTPUT(3);
      HIGHS_RAND_TRY_OUTPUT(4);
      HIGHS_RAND_TRY_OUTPUT(5);
      HIGHS_RAND_TRY_OUTPUT(6);
      HIGHS_RAND_TRY_OUTPUT(7);
      HIGHS_RAND_TRY_OUTPUT(9);
      HIGHS_RAND_TRY_OUTPUT(10);
      HIGHS_RAND_TRY_OUTPUT(11);
      HIGHS_RAND_TRY_OUTPUT(12);
      HIGHS_RAND_TRY_OUTPUT(13);
      HIGHS_RAND_TRY_OUTPUT(14);
      HIGHS_RAND_TRY_OUTPUT(15);
      HIGHS_RAND_TRY_OUTPUT(16);
      HIGHS_RAND_TRY_OUTPUT(17);
      HIGHS_RAND_TRY_OUTPUT(18);
      HIGHS_RAND_TRY_OUTPUT(19);
      HIGHS_RAND_TRY_OUTPUT(20);
      HIGHS_RAND_TRY_OUTPUT(21);
      HIGHS_RAND_TRY_OUTPUT(22);
      HIGHS_RAND_TRY_OUTPUT(23);
      HIGHS_RAND_TRY_OUTPUT(24);
      HIGHS_RAND_TRY_OUTPUT(25);
      HIGHS_RAND_TRY_OUTPUT(26);
      HIGHS_RAND_TRY_OUTPUT(27);
      HIGHS_RAND_TRY_OUTPUT(28);
      HIGHS_RAND_TRY_OUTPUT(29);
      HIGHS_RAND_TRY_OUTPUT(30);
      HIGHS_RAND_TRY_OUTPUT(31);

#undef HIGHS_RAND_TRY_OUTPUT
    }
  }

  uint64_t drawUniform(uint64_t sup, int nbits) {
    // 64 bit version for drawUniform. If result fits in 32 bits use 32 bit
    // version above so that it is only called when we actually need more than
    // 32 bits
    if (nbits <= 32) return drawUniform(uint32_t(sup), nbits);

    assert(sup <= uint64_t{1} << nbits);
    while (true) {
      advance();
      uint32_t lo = static_cast<uint32_t>(state);
      uint32_t hi = state >> 32;

      uint64_t val;

#define HIGHS_RAND_TRY_OUTPUT(n)                                       \
  val = (HighsHashHelpers::pair_hash<2 * n>(lo, hi) >> (64 - nbits)) ^ \
        (HighsHashHelpers::pair_hash<2 * n + 1>(lo, hi) >> 32);        \
  if (val < sup) return val;

      HIGHS_RAND_TRY_OUTPUT(0);
      HIGHS_RAND_TRY_OUTPUT(1);
      HIGHS_RAND_TRY_OUTPUT(2);
      HIGHS_RAND_TRY_OUTPUT(3);
      HIGHS_RAND_TRY_OUTPUT(4);
      HIGHS_RAND_TRY_OUTPUT(5);
      HIGHS_RAND_TRY_OUTPUT(6);
      HIGHS_RAND_TRY_OUTPUT(7);
      HIGHS_RAND_TRY_OUTPUT(9);
      HIGHS_RAND_TRY_OUTPUT(10);
      HIGHS_RAND_TRY_OUTPUT(11);
      HIGHS_RAND_TRY_OUTPUT(12);
      HIGHS_RAND_TRY_OUTPUT(13);
      HIGHS_RAND_TRY_OUTPUT(14);
      HIGHS_RAND_TRY_OUTPUT(15);

#undef HIGHS_RAND_TRY_OUTPUT
    }
  }

 public:
  /**
   * @brief Initialisations
   */
  HighsRandom(HighsUInt seed = 0) { initialise(seed); }

  /**
   * @brief (Re-)initialise the random number generator
   */
  void initialise(HighsUInt seed = 0) {
    state = seed;
    do {
      state = HighsHashHelpers::pair_hash<0>(static_cast<uint32_t>(state),
                                             state >> 32);
      state ^= (HighsHashHelpers::pair_hash<1>(state >> 32, seed) >> 32);
    } while (state == 0);
  }

  void advance() {
    // advance state with simple xorshift the outputs are produced by applying
    // strongly universal hash functions to the state so that the least
    // significant bits are as well distributed as the most significant bits.
    state ^= state >> 12;
    state ^= state << 25;
    state ^= state >> 27;
  }

  /**
   * @brief Return a random integer between 0 and 2147483647
   */
  HighsInt integer() {
    advance();
#ifdef HIGHSINT64
    // use 63 bits of the first hash result and use second hash for lower half
    // of bits
    return (HighsHashHelpers::pair_hash<0>(state, state >> 32) >> 1) ^
           (HighsHashHelpers::pair_hash<1>(state, state >> 32) >> 32);
#else
    // use 31 bits of the 64 bit result
    return HighsHashHelpers::pair_hash<0>(static_cast<uint32_t>(state),
                                          state >> 32) >>
           33;
#endif
  }

  /**
   * @brief Return a random integer between [0,sup)
   */
  HighsInt integer(HighsInt sup) {  // let overload resolution select the 32bit
                                    // or the 64bit version
    if (sup <= 1) return 0;
    int nbits = HighsHashHelpers::log2i(HighsUInt(sup - 1)) + 1;
    return static_cast<HighsInt>(drawUniform(HighsUInt(sup), nbits));
  }

  /**
   * @brief Return a random integer between [min,sup)
   */
  HighsInt integer(HighsInt min, HighsInt sup) {
    return min + integer(sup - min);
  }

  /**
   * @brief Return a random fraction - real in (0, 1)
   */
  double fraction() {
    advance();
    // 52 bit output is in interval [0,2^52-1]
    uint64_t output = (HighsHashHelpers::pair_hash<0>(
                           static_cast<uint32_t>(state), state >> 32) >>
                       (64 - 52)) ^
                      (HighsHashHelpers::pair_hash<1>(
                           static_cast<uint32_t>(state), state >> 32) >>
                       (64 - 26));
    // compute (1+output) / (2^52+1) which is strictly between 0 and 1
    return static_cast<double>(1 + output) * 2.2204460492503125e-16;
  }

  /**
   * @brief Return a random fraction - real in [0, 1]
   */
  double closedFraction() {
    advance();
    // 53 bit result is in interval [0,2^53-1]
    uint64_t output = (HighsHashHelpers::pair_hash<0>(
                           static_cast<uint32_t>(state), state >> 32) >>
                       (64 - 53)) ^
                      (HighsHashHelpers::pair_hash<1>(
                           static_cast<uint32_t>(state), state >> 32) >>
                       32);
    // compute output / (2^53-1) in double precision which is in the closed
    // interval [0,1]
    return static_cast<double>(output) * 1.1102230246251566e-16;
  }

  /**
   * @brief Return a random real value in the interval [a,b]
   */
  double real(double a, double b) { return a + (b - a) * closedFraction(); }

  /**
   * @brief Return a random bit
   */
  bool bit() {
    advance();
    return HighsHashHelpers::pair_hash<0>(static_cast<uint32_t>(state),
                                          state >> 32) >>
           63;
  }

  /**
   * @brief shuffle the given data array
   */
  template <typename T>
  void shuffle(T* data, HighsInt N) {
    for (HighsInt i = N; i > 1; --i) {
      HighsInt pos = integer(i);
      std::swap(data[pos], data[i - 1]);
    }
  }

 private:
  uint64_t state;
};

#endif
