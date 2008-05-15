/*
 * Last Change: Tue May 13 02:00 PM 2008 J
 *
 * FFTPACK implementation
 *
 * Original code by Pearu Peterson.
 */
#include "api.h"

#include "common.h"

using namespace fft;

static CacheManager<RFFTPackCacheId, RFFTPackCache> rfftpack_cmgr(10);

void drfft_fftpack(double *inout, int n, int direction, int howmany,
			  int normalize)
{
        int i;
        double *ptr = inout;
        RFFTPackCache* cache;

        cache = rfftpack_cmgr.get_cache(RFFTPackCacheId(n));

        switch (direction) {
                case 1:
                        for (i = 0; i < howmany; ++i, ptr += n) {
                                cache->compute_forward(ptr);
                        }
                        break;

                case -1:
                        for (i = 0; i < howmany; ++i, ptr += n) {
                                cache->compute_backward(ptr);
                        }
                        break;

                default:
                        fprintf(stderr, "drfft: invalid direction=%d\n", direction);
                        return;
        }

        if (normalize) {
                double d = 1.0 / n;
                ptr = inout;
                for (i = n * howmany - 1; i >= 0; --i) {
                        (*(ptr++)) *= d;
                }
        }
}
