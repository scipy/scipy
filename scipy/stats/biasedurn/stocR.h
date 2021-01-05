#ifndef _STOCR_H_
#define _STOCR_H_

#include <cstddef>  // for NULL

struct StocRBase {

  double(*next_double)(void*);

  StocRBase() : next_double(NULL) {}
  StocRBase(int seed) : next_double(NULL) {}

  double Random() {
    return next_double(NULL);
  }
};

#endif // _STOCR_H_
