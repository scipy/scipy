/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                       */
/*    This file is part of the HiGHS linear optimization suite           */
/*                                                                       */
/*    Available as open-source under the MIT License                     */
/*                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file util/HighsCDouble.h
 * @brief Quad precision type implemented with two standard double precision
 *        representing the value and a compensation term
 */
#ifndef UTIL_HIGHSCDOUBLE_H_
#define UTIL_HIGHSCDOUBLE_H_

#include <cmath>
#include <cstdint>

/// A compensated double number achieving roughly quad precision on the
/// supported operations

class HighsCDouble {
 private:
  double hi;
  double lo;

  // The following functions are implemented as described in:
  // Rump, Siegfried M. "High precision evaluation of nonlinear functions."
  // Proceedings of. 2005.

  /// performs an exact transformation such that x + y = a + b
  /// and x = double(a + b). The operation uses 6 flops (addition/subtraction).
  static void two_sum(double& x, double& y, double a, double b) {
    x = a + b;
    double z = x - a;
    y = (a - (x - z)) + (b - z);
  }

  /// splits a 53 bit double precision number into two 26 bit parts
  /// such that x + y = a holds exactly
  static void split(double& x, double& y, double a) {
    constexpr double factor = double((1 << 27) + 1);
    double c = factor * a;
    x = c - (c - a);
    y = a - x;
  }

  /// performs an exact transformation such that x + y = a * b
  /// and x = double(a * b). The operation uses 10 flops for
  /// addition/subtraction and 7 flops for multiplication.
  static void two_product(double& x, double& y, double a, double b) {
    x = a * b;
    double a1, a2, b1, b2;
    split(a1, a2, a);
    split(b1, b2, b);
    y = a2 * b2 - (((x - a1 * b1) - a2 * b1) - a1 * b2);
  }

  HighsCDouble(double hi_, double lo_) : hi(hi_), lo(lo_) {}

 public:
  HighsCDouble() = default;

  HighsCDouble(double val) : hi(val), lo(0.0) {}

  explicit operator double() const { return hi + lo; }

  HighsCDouble& operator+=(double v) {
    double c;
    two_sum(hi, c, v, hi);
    lo += c;
    return *this;
  }

  HighsCDouble& operator+=(const HighsCDouble& v) {
    (*this) += v.hi;
    lo += v.lo;
    return *this;
  }

  HighsCDouble& operator-=(double v) {
    (*this) += -v;
    return *this;
  }

  HighsCDouble& operator-=(const HighsCDouble& v) {
    (*this) -= v.hi;
    lo -= v.lo;
    return *this;
  }

  HighsCDouble& operator*=(double v) {
    double c = lo * v;
    two_product(hi, lo, hi, v);
    *this += c;
    return *this;
  }

  HighsCDouble& operator*=(const HighsCDouble& v) {
    double c1 = hi * v.lo;
    double c2 = lo * v.hi;
    two_product(hi, lo, hi, v.hi);
    *this += c1;
    *this += c2;
    return *this;
  }

  HighsCDouble& operator/=(double v) {
    HighsCDouble d(hi / v, lo / v);
    HighsCDouble c = d * v - (*this);
    c.hi /= v;
    c.lo /= v;
    *this = d - c;
    return *this;
  }

  HighsCDouble& operator/=(const HighsCDouble& v) {
    double vdbl = v.hi + v.lo;
    HighsCDouble d(hi / vdbl, lo / vdbl);
    HighsCDouble c = d * v - (*this);
    c.hi /= vdbl;
    c.lo /= vdbl;
    *this = d - c;
    return *this;
  }

  HighsCDouble operator-() const { return HighsCDouble(-hi, -lo); }

  HighsCDouble operator+(double v) const {
    HighsCDouble res;

    two_sum(res.hi, res.lo, hi, v);
    res.lo += lo;

    return res;
  }

  HighsCDouble operator+(const HighsCDouble& v) const {
    HighsCDouble res = (*this) + v.hi;
    res.lo += v.lo;

    return res;
  }

  friend HighsCDouble operator+(double a, const HighsCDouble& b) {
    return b + a;
  }

  HighsCDouble operator-(double v) const {
    HighsCDouble res;

    two_sum(res.hi, res.lo, hi, -v);
    res.lo += lo;

    return res;
  }

  HighsCDouble operator-(const HighsCDouble& v) const {
    HighsCDouble res = (*this) - v.hi;
    res.lo -= v.lo;

    return res;
  }

  friend HighsCDouble operator-(double a, const HighsCDouble& b) {
    return -b + a;
  }

  HighsCDouble operator*(double v) const {
    HighsCDouble res;

    two_product(res.hi, res.lo, hi, v);
    res += lo * v;

    return res;
  }

  HighsCDouble operator*(const HighsCDouble& v) const {
    HighsCDouble res = (*this) * v.hi;
    res += hi * v.lo;

    return res;
  }

  friend HighsCDouble operator*(double a, const HighsCDouble& b) {
    return b * a;
  }

  HighsCDouble operator/(double v) const {
    HighsCDouble res = *this;

    res /= v;

    return res;
  }

  HighsCDouble operator/(const HighsCDouble& v) const {
    HighsCDouble res = (*this);

    res /= v;

    return res;
  }

  friend HighsCDouble operator/(double a, const HighsCDouble& b) {
    return HighsCDouble(a) / b;
  }

  bool operator>(const HighsCDouble& other) const {
    return double(*this) > double(other);
  }

  bool operator>(double other) const { return double(*this) > other; }

  friend bool operator>(double a, const HighsCDouble& b) {
    return a > double(b);
  }

  bool operator<(const HighsCDouble& other) const {
    return double(*this) < double(other);
  }

  bool operator<(double other) const { return double(*this) < other; }

  friend bool operator<(double a, const HighsCDouble& b) {
    return a < double(b);
  }

  bool operator>=(const HighsCDouble& other) const {
    return double(*this) >= double(other);
  }

  bool operator>=(double other) const { return double(*this) >= other; }

  friend bool operator>=(double a, const HighsCDouble& b) {
    return a >= double(b);
  }

  bool operator<=(const HighsCDouble& other) const {
    return double(*this) <= double(other);
  }

  bool operator<=(double other) const { return double(*this) <= other; }

  friend bool operator<=(double a, const HighsCDouble& b) {
    return a <= double(b);
  }

  bool operator==(const HighsCDouble& other) const {
    return double(*this) == double(other);
  }

  bool operator==(double other) const { return double(*this) == other; }

  friend bool operator==(double a, const HighsCDouble& b) {
    return a == double(b);
  }

  bool operator!=(const HighsCDouble& other) const {
    return double(*this) != double(other);
  }

  bool operator!=(double other) const { return double(*this) != other; }

  friend bool operator!=(double a, const HighsCDouble& b) {
    return a != double(b);
  }

  void renormalize() { two_sum(hi, lo, hi, lo); }

  friend HighsCDouble abs(const HighsCDouble& v) { return v < 0 ? -v : v; }

  friend HighsCDouble sqrt(const HighsCDouble& v) {
    double c = std::sqrt(v.hi + v.lo);

    // guard against division by zero
    if (c == 0.0) return 0.0;

    // calculate precise square root by newton step
    HighsCDouble res = v / c;
    res += c;
    // multiplication by 0.5 is exact
    res.hi *= 0.5;
    res.lo *= 0.5;
    return res;
  }

  friend HighsCDouble floor(const HighsCDouble& x) {
    // Treat |x| < 1 as special case, as per (for example)
    // https://github.com/shibatch/tlfloat: see #2041
    if (abs(x) < 1) {
      if (x == 0 || x > 0) return HighsCDouble(0.0);
      return HighsCDouble(-1.0);
    }
    double floor_x = std::floor(double(x));
    HighsCDouble res;

    two_sum(res.hi, res.lo, floor_x, std::floor(double(x - floor_x)));
    return res;
  }

  friend HighsCDouble ceil(const HighsCDouble& x) {
    // Treat |x| < 1 as special case, as per (for example)
    // https://github.com/shibatch/tlfloat: see #2041
    if (abs(x) < 1) {
      if (x == 0 || x < 0) return HighsCDouble(0.0);
      return HighsCDouble(1.0);
    }
    double ceil_x = std::ceil(double(x));
    HighsCDouble res;

    two_sum(res.hi, res.lo, ceil_x, std::ceil(double(x - ceil_x)));
    return res;
  }

  friend HighsCDouble round(const HighsCDouble& x) { return floor(x + 0.5); }

  friend HighsCDouble ldexp(const HighsCDouble& v, int exp) {
    return HighsCDouble(std::ldexp(v.hi, exp), std::ldexp(v.lo, exp));
  }

  friend const HighsCDouble& min(const HighsCDouble& a, const HighsCDouble& b) {
    return a < b ? a : b;
  }

  friend const HighsCDouble& max(const HighsCDouble& a, const HighsCDouble& b) {
    return a > b ? a : b;
  }
};

#endif
