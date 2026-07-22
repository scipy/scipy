/*
Copyright (c) 2005-2022, NumPy Developers.
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
    * Redistributions of source code must retain the above copyright
       notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above
       copyright notice, this list of conditions and the following
       disclaimer in the documentation and/or other materials provided
       with the distribution.
    * Neither the name of the NumPy Developers nor the names of any
       contributors may be used to endorse or promote products derived
       from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <math.h>
#include <stdint.h>

#include "distributions.h"
#include "ziggurat_constants.h"

/* Inline generators for internal use */
static inline uint32_t next_uint32(bitgen_t *bitgen_state) {
  return bitgen_state->next_uint32(bitgen_state->state);
}

static inline uint64_t next_uint64(bitgen_t *bitgen_state) {
  return bitgen_state->next_uint64(bitgen_state->state);
}

static inline double next_double(bitgen_t *bitgen_state) {
    return bitgen_state->next_double(bitgen_state->state);
}


double random_standard_normal(bitgen_t *bitgen_state) {
  uint64_t r;
  int sign;
  uint64_t rabs;
  int idx;
  double x, xx, yy;
  for (;;) {
    /* r = e3n52sb8 */
    r = next_uint64(bitgen_state);
    idx = r & 0xff;
    r >>= 8;
    sign = r & 0x1;
    rabs = (r >> 1) & 0x000fffffffffffff;
    x = rabs * wi_double[idx];
    if (sign & 0x1)
      x = -x;
    if (rabs < ki_double[idx])
      return x; /* 99.3% of the time return here */
    if (idx == 0) {
      for (;;) {
        /* Switch to 1.0 - U to avoid log(0.0), see GH 13361 */
        xx = -ziggurat_nor_inv_r * log1p(-next_double(bitgen_state));
        yy = -log1p(-next_double(bitgen_state));
        if (yy + yy > xx * xx)
          return ((rabs >> 8) & 0x1) ? -(ziggurat_nor_r + xx)
                                     : ziggurat_nor_r + xx;
      }
    } else {
      if (((fi_double[idx - 1] - fi_double[idx]) * next_double(bitgen_state) +
           fi_double[idx]) < exp(-0.5 * x * x))
        return x;
    }
  }
}


double random_normal(bitgen_t *bitgen_state, double loc, double scale) {
  return loc + scale * random_standard_normal(bitgen_state);
}


uint64_t random_interval(bitgen_t *bitgen_state, uint64_t max) {
  uint64_t mask, value;
  if (max == 0) {
    return 0;
  }

  mask = max;

  /* Smallest bit mask >= max */
  mask |= mask >> 1;
  mask |= mask >> 2;
  mask |= mask >> 4;
  mask |= mask >> 8;
  mask |= mask >> 16;
  mask |= mask >> 32;

  /* Search a random value in [0..mask] <= max */
  if (max <= 0xffffffffUL) {
    while ((value = (next_uint32(bitgen_state) & mask)) > max)
      ;
  } else {
    while ((value = (next_uint64(bitgen_state) & mask)) > max)
      ;
  }
  return value;
}


double random_standard_uniform(bitgen_t *bitgen_state) {
    return next_double(bitgen_state);
}
