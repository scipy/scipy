#ifndef _STOCR_H_
#define _STOCR_H_

#include <cstddef>  // for NULL

// Random number source for the stocc library. The generator state is held per
// instance and handed back to the callbacks, so concurrent instances cannot
// interfere with each other.
struct StocRBase {

  void *rng_state;
  double(*next_double)(void *rng_state);
  double(*next_normal)(void *rng_state, const double m, const double s);

  StocRBase() : rng_state(NULL), next_double(NULL), next_normal(NULL) {}
  StocRBase(int seed) : rng_state(NULL), next_double(NULL), next_normal(NULL) {}

  double Random() {
    return next_double(rng_state);
  }

  double Normal(double m, double s) {
    // Also see impls.cpp for the StochasticLib1 implementation
    // (should be identical to this)
    return next_normal(rng_state, m, s);
  }
};

#endif // _STOCR_H_
