/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Written and engineered 2008-2020 at the University of Edinburgh    */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/**@file util/HighsRandom.h
 * @brief Random number generators for HiGHS
 * @author Julian Hall, Ivet Galabova, Qi Huangfu and Michael Feldmeier
 */
#ifndef UTIL_HIGHSRANDOM_H_
#define UTIL_HIGHSRANDOM_H_

const int initial_random_mw = 1985;
const int initial_random_mz = 2012;

/**
 * @brief Class for HiGHS random number generators
 */
class HighsRandom {
 public:
  /**
   * @brief Initialisations
   */
  HighsRandom() {
    /**
     * @brief Initialise the two seeds to default values
     */
    random_mw = initial_random_mw;
    random_mz = initial_random_mz;
  }

  /**
   * @brief (Re-)initialise the random number generator
   */
  void initialise() {
    random_mw = initial_random_mw;
    random_mz = initial_random_mz;
  }

  /**
   * @brief Return a random integer between 0 and 2147483647
   */
  int integer() {
    random_mz = 36969 * (random_mz & 65535) + (random_mz >> 16);
    random_mw = 18000 * (random_mw & 65535) + (random_mw >> 16);
    unsigned result = (random_mz << 16) + random_mw;
    return result >> 1;
  }

  /**
   * @brief Return a random fraction - real in (0, 1)
   */
  double fraction() {
    random_mz = 36969 * (random_mz & 65535) + (random_mz >> 16);
    random_mw = 18000 * (random_mw & 65535) + (random_mw >> 16);
    unsigned result = (random_mz << 16) + random_mw;
    double returnValue = (result + 1.0) * 2.328306435454494e-10;
    return returnValue;
  }

 private:
  unsigned random_mw;
  unsigned random_mz;
};
#endif
