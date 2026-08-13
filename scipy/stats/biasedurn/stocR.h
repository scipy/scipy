#ifndef _STOCR_H_
#define _STOCR_H_

#include <cstddef>  // for NULL

struct StocRBase {

  double(*next_double)();
  double(*next_normal)(const double m, const double s);

  StocRBase() : next_double(NULL), next_normal(NULL) {}
  StocRBase(int seed) : next_double(NULL), next_normal(NULL) {}

  double Random() {
    return next_double();
  }

  double Normal(double m, double s) {
    // Also see impls.cpp for the StochasticLib1 implementation
    // (should be identical to this)
    return next_normal(m, s);
  }
};

#endif // _STOCR_H_
