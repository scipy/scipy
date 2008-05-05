#ifndef _SCIPYFFT_FFTW_COMMON_H_
#define _SCIPYFFT_FFTW_COMMON_H_

#include "cycliccache.h"

namespace fft {

class FFTW3CacheId : public CacheId {
	public:
		FFTW3CacheId(int n, int dir) : 
                        CacheId(n), 
                        m_dir(dir)
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

			return m_dir == other.m_dir &&  th->is_equal(*ot);
		}

	public:
		int m_dir;
};

}; // namespace fft

#endif
