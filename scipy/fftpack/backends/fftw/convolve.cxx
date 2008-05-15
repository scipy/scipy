#include <new>

#include <rfftw.h>
#include <fftw.h>

#include "api.h"

#include "cycliccache.h"

using namespace fft;

class DRFFTWCacheId : public CacheId {
	public:
		DRFFTWCacheId(int n);

};

DRFFTWCacheId::DRFFTWCacheId(int n):
        CacheId(n)
{
}

class DRFFTWCache : public Cache<DRFFTWCacheId> {
	public:	
		DRFFTWCache(const DRFFTWCacheId& id);
		virtual ~DRFFTWCache();

                int convolve(double* inout, double* omega, 
                             int swap_real_imag) const;
                int convolve_z(double* inout, double* omega_real, 
                               double* omega_imag) const;

	protected:
		rfftw_plan m_plan1;
		rfftw_plan m_plan2;
};

DRFFTWCache::DRFFTWCache(const DRFFTWCacheId& id)
:	Cache<DRFFTWCacheId>(id)
{
        int flags = FFTW_ESTIMATE | FFTW_IN_PLACE;

	m_plan1 = rfftw_create_plan(id.m_n, FFTW_REAL_TO_COMPLEX, flags);
	if (m_plan1 == NULL) {
		goto fail;
	}

	m_plan2 = rfftw_create_plan(id.m_n, FFTW_COMPLEX_TO_REAL, flags);
	if (m_plan2 == NULL) {
		goto clean_plan1;
	}

        return;

clean_plan1:
	rfftw_destroy_plan(m_plan1);
fail:
	throw std::bad_alloc();
}

DRFFTWCache::~DRFFTWCache()
{
	rfftw_destroy_plan(m_plan2);
	rfftw_destroy_plan(m_plan1);
}

int DRFFTWCache::convolve(double* inout, double* omega, int swap_real_imag)
    const
{
        int n = m_id.m_n;
	int l = (n-1)/2+1;
        int i;
        double c;

	rfftw_one(m_plan1, (fftw_real *)inout, NULL);
	if (swap_real_imag) {
		inout[0] *= omega[0];
		if (!(n%2)) {
			inout[n/2] *= omega[n/2];
		}
		for(i=1;i<l;++i) {
			c = inout[i] * omega[i];
			inout[i] = omega[n-i] * inout[n-i];
			inout[n-i] = c;
		}
	} else {
		for(i=0;i<n;++i) {
			inout[i] *= omega[i];
		}
	}
	rfftw_one(m_plan2, (fftw_real *)inout, NULL);

        return 0;
}

int DRFFTWCache::convolve_z(double* inout, double* omega_real, 
                            double* omega_imag) const
{
        int n = m_id.m_n;
	int i;
	int l = (n-1)/2+1;
	double c;

	rfftw_one(m_plan1, (fftw_real *)inout, NULL);
	inout[0] *= (omega_real[0]+omega_imag[0]);

	if (!(n%2)) {
		inout[n/2] *= (omega_real[n/2]+omega_imag[n/2]);
	}
	for(i=1;i<l;++i) {
		c = inout[i] * omega_imag[i];
		inout[i] *= omega_real[i];
		inout[i] += omega_imag[n-i] * inout[n-i];
		inout[n-i] *= omega_real[n-i];
		inout[n-i] += c;
	}
	rfftw_one(m_plan2, (fftw_real *)inout, NULL);

        return 0;
}

CacheManager<DRFFTWCacheId, DRFFTWCache> drfftw_cmgr(20);

/**************** convolve **********************/
void convolve_fftw(int n,double* inout,double* omega,int swap_real_imag) 
{
        DRFFTWCache *cache;

        cache = drfftw_cmgr.get_cache(DRFFTWCacheId(n));
        cache->convolve(inout, omega, swap_real_imag);
}

/**************** convolve **********************/
void convolve_z_fftw(int n,double* inout,double* omega_real,double* omega_imag) 
{
        DRFFTWCache *cache;

        cache = drfftw_cmgr.get_cache(DRFFTWCacheId(n));
        cache->convolve_z(inout, omega_real, omega_imag);
}

void init_convolution_kernel_fftw(int n,double* omega, int d,
			     double (*kernel_func)(int),
			     int zero_nyquist) 
{
	/*
	 *  omega[k] = pow(sqrt(-1),d) * kernel_func(k)
	 *  omega[0] = kernel_func(0)
	 *  conjugate(omega[-k]) == omega[k]
	 */
	int k,l=(n-1)/2+1;
	omega[0] = (*kernel_func)(0)/n;;
	switch (d%4) {
		case 0:
			for (k=1;k<l;++k)
				omega[k] = omega[n-k] = (*kernel_func)(k)/n;
			if (!(n%2)) 
				omega[n/2] = (zero_nyquist?0.0:(*kernel_func)(n/2)/n);
			break;
		case 1:;case -3:
		       for (k=1;k<l;++k) {
			       omega[k] = (*kernel_func)(k)/n;
			       omega[n-k] = -omega[k];
		       }
		       if (!(n%2))
			       omega[n/2] = (zero_nyquist?0.0:(*kernel_func)(n/2)/n);
		       break;
		case 2:;case -2:
		       for (k=1;k<l;++k)
			       omega[k] = omega[n-k] = -(*kernel_func)(k)/n;
		       if (!(n%2))
			       omega[n/2] = (zero_nyquist?0.0:-(*kernel_func)(n/2)/n);
		       break;
		case 3:;case -1:
		       for (k=1;k<l;++k) {
			       omega[k] = -(*kernel_func)(k)/n;
			       omega[n-k] = -omega[k];
		       }
		       if (!(n%2))
			       omega[n/2] = (zero_nyquist?0.0:-(*kernel_func)(n/2)/n);
		       break;
	}
}
