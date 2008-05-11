#ifndef _SCIPY_DJBFFT_COMMON_H
#define _SCIPY_DJBFFT_COMMON_H

class DJBFFTCacheId : public CacheId {
        public:
                DJBFFTCacheId(int n) : CacheId(n) {};
};

#endif
