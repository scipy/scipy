#ifndef _SCIPYFFT_FFTW_COMMON_H_
#define _SCIPYFFT_FFTW_COMMON_H_

#include "cycliccache.h"

#include <stddef.h>

#define COPYRFFTW2STD(SRC,DEST,N) { \
  int j,n2=(N)/2; \
  *(DEST) = *(SRC); \
  for (j=1;j<n2;++j) { \
    *(DEST+2*j-1) = *(SRC+j); \
    *(DEST+2*j) = *(SRC+(N)-j); \
  } \
  if (N>1) { \
    *(DEST+2*n2-1) = *(SRC+n2); \
    if ((N)%2) \
      *(DEST+2*n2) = *(SRC+(N)-n2); \
  } \
}
#define COPYINVRFFTW2STD(SRC,DEST,N) { \
  int j,n2=(N)/2; \
  *(DEST) = *(SRC); \
  for (j=1;j<n2;++j) { \
    *(DEST+j) = *(SRC+2*j-1); \
    *(DEST+(N)-j) = *(SRC+2*j); \
  } \
  if (N>1) {\
    *(DEST+n2) = *(SRC+2*n2-1); \
    if ((N)%2) \
      *(DEST+(N)-n2) = *(SRC+2*n2); \
  } \
}

namespace fft {

inline bool is_simd_aligned(const void * p)
{
        return (((reinterpret_cast<ptrdiff_t> (p)) & 0xF) == 0);
}

class FFTW3CacheId : public CacheId {
	public:
		FFTW3CacheId(int n, int dir, bool isalign) : 
                        CacheId(n), 
                        m_dir(dir),
                        m_isalign(isalign)
                {
                };

		virtual bool operator==(const FFTW3CacheId& other) const
		{
			return is_equal(other);
		}

		virtual bool is_equal(const FFTW3CacheId& other) const
		{
			const CacheId *ot = &other;
			const CacheId *th = this;

			return m_isalign == other.m_isalign && 
                               m_dir == other.m_dir &&  th->is_equal(*ot);
		}

	public:
		int m_dir;
                bool m_isalign;
};

}; // namespace fft

#endif
