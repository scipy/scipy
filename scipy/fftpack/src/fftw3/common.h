#ifndef _SCIPYFFT_FFTW_COMMON_H_
#define _SCIPYFFT_FFTW_COMMON_H_

#include "cycliccache.h"

#include <stddef.h>

namespace fft {

inline bool is_simd_aligned(const void * p)
{
        return ((reinterpret_cast<ptrdiff_t> (p)) & 0xF == 0);
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
