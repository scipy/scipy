#ifndef _STOCR_H_
#define _STOCR_H_

struct StocRBase {

  double(*next_double)();

  StocRBase() : next_double(NULL) {}
  StocRBase(int seed) : next_double(NULL) {}

  double Random() {
    return next_double(NULL);
  }
};

#endif // _STOCR_H_
