// Numerically stable computation of iv(v+1, x) / iv(v, x)

#pragma once

#include "special/tools.h"

class IvRatioCFGenerator {
  /* This is the continued fraction for (1 + sqrt(1 + 4x))/2. Replace the
   * operator with something which generates the coefficients for Perron's
   * continued fraction.
   */

public:
  IvRatioCFGenerator(double v, double x) : v_(v), x_(x) {}

  std::pair<double, double> operator()() { return {1.0, v_ * x_}; }

private:
  double v_, x_;
};

inline double iv_ratio(double v, double x) {
  auto generator = IvRatioCFGenerator(v, x);
  return special::detail::continued_fraction_eval(
      generator, 1.0, std::numeric_limits<double>::epsilon(), 500, "iv_ratio");
}
