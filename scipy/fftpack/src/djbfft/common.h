#ifndef _SCIPY_DJBFFT_COMMON_H
#define _SCIPY_DJBFFT_COMMON_H

#include <cycliccache.h>

#include "api.h"

namespace fft {

class DJBFFTCacheId : public CacheId {
        public:
                DJBFFTCacheId(int n) : CacheId(n) {};
};

};

#endif
