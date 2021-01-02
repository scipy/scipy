#ifndef _STOCR_H_
#define _STOCR_H_

#include <cstdint>
#include <stdexcept>
#include "numpy/random/bitgen.h"

struct StocRBase {

  // bitgen_state is owned by the caller and its
  // lifetime will be assumed to be the same or
  // greater than this StocRBase object
  bitgen_t* bitgen_state;

  StocRBase() : bitgen_state(NULL) {}
  StocRBase(std::int32_t seed) : bitgen_state(NULL) {}

  void SetBitGen(bitgen_t* that_bitgen_state) {
    bitgen_state = that_bitgen_state;
  }

  void RandomInit() {
    // we can only generate random numbers if we have bitgen_t object
    if (bitgen_state == NULL) {
      throw std::runtime_error("SetBitGen(bitgen_state) has not been called!");
    }
  }

  double Random() {
    return bitgen_state->next_double(bitgen_state->state);
  }
};

#endif // _STOCR_H_