/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#ifndef HIGHS_UTIL_INTEGERS_H_
#define HIGHS_UTIL_INTEGERS_H_

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "util/HighsCDouble.h"
#include "util/HighsType.h"

class HighsIntegers {
 public:
  static int64_t mod(int64_t a, int64_t m) {
    int64_t r = a % m;
    return r + (r < 0) * m;
  }

  static double mod(double a, double m) {
    return std::trunc(std::fmod(a, m)) + (a < 0) * m;
  }

  static int64_t nearestInteger(double x) {
    return (int64_t)(x + std::copysign(0.5, x));
  }

  static bool isIntegral(double x, double eps) {
    double y = std::fabs(x - std::trunc(x));
    return std::min(y, 1.0 - y) <= eps;
  }

  static int64_t modularInverse(int64_t a, int64_t m) {
    int64_t y = 0;
    int64_t x = 1;

    if (m == 1) return 0;

    a = mod(a, m);

    while (a > 1) {
      // compute quotient q = a / m and remainder r = a % m
      int64_t q = a / m;
      int64_t r = a - q * m;

      // update (a,m) = (m,r)
      a = m;
      m = r;

      // update x and y of extended euclidean algorithm
      r = x - q * y;
      x = y;
      y = r;
    }

    return x;
  }

  static int64_t gcd(int64_t a, int64_t b) {
    assert(a != std::numeric_limits<int64_t>::min());
    assert(b != std::numeric_limits<int64_t>::min());

    int64_t h;
    if (a < 0) a = -a;
    if (b < 0) b = -b;

    if (a == 0) return b;
    if (b == 0) return a;

    do {
      h = a % b;
      a = b;
      b = h;
    } while (b != 0);

    return a;
  }

  // computes a rational approximation with given maximal denominator
  static int64_t denominator(double x, double eps, int64_t maxdenom) {
    int64_t ai = (int64_t)x;
    int64_t m[] = {ai, 1, 1, 0};

    HighsCDouble xi = x;
    HighsCDouble fraction = xi - double(ai);

    while (fraction > eps) {
      xi = 1.0 / fraction;
      if (double(xi) > double(int64_t{1} << 53)) break;

      ai = (int64_t)(double)xi;
      int64_t t = m[2] * ai + m[3];
      if (t > maxdenom) break;

      m[3] = m[2];
      m[2] = t;

      t = m[0] * ai + m[1];
      m[1] = m[0];
      m[0] = t;

      fraction = xi - ai;
    }

    ai = (maxdenom - m[3]) / m[2];
    m[1] += m[0] * ai;
    m[3] += m[2] * ai;

    double x0 = static_cast<double>(m[0]) / static_cast<double>(m[2]);
    double x1 = static_cast<double>(m[1]) / static_cast<double>(m[3]);
    x = std::abs(x);
    double err0 = std::abs(x - x0);
    double err1 = std::abs(x - x1);

    if (err0 < err1) return m[2];
    return m[3];
  }

  static double integralScale(const double* vals, HighsInt numVals,
                              double deltadown, double deltaup) {
    if (numVals == 0) return 0.0;

    auto minmax = std::minmax_element(
        vals, vals + numVals,
        [](double a, double b) { return std::abs(a) < std::abs(b); });
    const double minval = *minmax.first;
    const double maxval = *minmax.second;

    int expshift = 0;

    // to cover many small denominators at once use a denominator of 75 * 2^n
    // with n-3 being large enough so that the smallest value is not below 0.5
    // but ignore tiny values bew deltadown/deltaup.
    if (minval < -deltadown || minval > deltaup) std::frexp(minval, &expshift);
    expshift = std::max(-expshift, 0) + 3;

    // guard against making the largest value too big which may cause overflows
    // with intermediate gcd values
    int expMaxVal;
    std::frexp(maxval, &expMaxVal);
    expMaxVal = std::min(expMaxVal, 32);
    if (expMaxVal + expshift > 32) expshift = 32 - expMaxVal;

    uint64_t denom = uint64_t{75} << expshift;
    int64_t startdenom = denom;
    // now check if the values are integral and if not compute a common
    // denominator for their remaining fraction
    HighsCDouble val = startdenom * HighsCDouble(vals[0]);
    HighsCDouble downval = floor(val + deltaup);
    HighsCDouble fraction = val - downval;

    if (fraction > deltadown) {
      // use a continued fraction algorithm to compute small missing
      // denominators for the remaining fraction
      denom *= denominator(double(fraction), deltaup, 1000);
      val = denom * HighsCDouble(vals[0]);
      downval = floor(val + deltaup);
      fraction = val - downval;

      // if this is not sufficient for reaching integrality, we stop here
      if (fraction > deltadown) return 0.0;
    }

    uint64_t currgcd = (uint64_t)std::abs(double(downval));

    for (HighsInt i = 1; i != numVals; ++i) {
      val = denom * HighsCDouble(vals[i]);
      downval = floor(val + deltaup);
      fraction = val - downval;

      if (fraction > deltadown) {
        val = startdenom * HighsCDouble(vals[i]);
        fraction = val - floor(val);
        denom *= denominator(double(fraction), deltaup, 1000);
        val = denom * HighsCDouble(vals[i]);
        downval = floor(val + deltaup);
        fraction = val - downval;

        if (fraction > deltadown) return 0.0;
      }

      if (currgcd != 1) {
        currgcd = gcd(currgcd, (int64_t) double(downval));

        // if the denominator is large, divide by the current gcd to prevent
        // unnecessary overflows
        if (denom > std::numeric_limits<unsigned int>::max()) {
          denom /= currgcd;
          if (startdenom != 1) startdenom /= gcd(currgcd, startdenom);
          currgcd = 1;
        }
      }
    }

    return static_cast<double>(denom) / static_cast<double>(currgcd);
  }

  static double integralScale(const std::vector<double>& vals, double deltadown,
                              double deltaup) {
    return integralScale(vals.data(), vals.size(), deltadown, deltaup);
  }
};

#endif
