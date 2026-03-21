#include "ducc0/fft/fftnd_impl.h"

namespace ducc0{
namespace detail_fft {
#define T double
#include "ducc0/fft/fft_inst_inc.h"
#undef T

// we don't need long double when nanobind is active
#ifndef DUCC0_USE_NANOBIND
#define T long double
#include "ducc0/fft/fft_inst_inc.h"
#undef T
#endif
}
}
