#include "common.h"

extern "C" void destroy_convolve_cache_fftpack(void)
{
}

using namespace fft;

class DFFTPackCache : public RFFTPackCache {
        public:
                DFFTPackCache(const RFFTPackCacheId& id) : RFFTPackCache(id) {};
                virtual ~DFFTPackCache() {};

        public:
                int convolve(double* inout, double* omega, int swap_real_imag) const;
                int convolve_z(double* inout, double* omega_real, double* omega_imag) const;
};

int DFFTPackCache::convolve(double* inout, double* omega, int swap_real_imag) 
        const
{
	int i;
        int n = m_id.m_n;
        double* wsave = m_wsave;

	F_FUNC(dfftf,DFFTF)(&n,inout,wsave);
	if (swap_real_imag) {
		double c;
		int n1 = n-1;
		inout[0] *= omega[0];
		if (!(n%2))
			inout[n-1] *= omega[n-1];
		for(i=1;i<n1;i+=2) {
			c = inout[i] * omega[i];
			inout[i] = inout[i+1] * omega[i+1];
			inout[i+1] = c;
		}
	}
	else
		for(i=0;i<n;++i)
			inout[i] *= omega[i];
	F_FUNC(dfftb,DFFTB)(&n,inout,wsave);

        return 0;
}

int DFFTPackCache::convolve_z(double* inout, double* omega_real, 
                              double* omega_imag) const
{
	int i;
	double* wsave = m_wsave;
	double c;
        int n = m_id.m_n;
	int n1 = n-1;

	F_FUNC(dfftf,DFFTF)(&n,inout,wsave);
	inout[0] *= (omega_real[0]+omega_imag[0]);
	if (!(n%2))
		inout[n-1] *= (omega_real[n-1]+omega_imag[n-1]);
	for(i=1;i<n1;i+=2) {
		c = inout[i] * omega_imag[i];
		inout[i] *= omega_real[i];
		inout[i] += inout[i+1] * omega_imag[i+1];
		inout[i+1] *= omega_real[i+1];
		inout[i+1] += c;
	}
	F_FUNC(dfftb,DFFTB)(&n,inout,wsave);

        return 0;
}

static CacheManager<RFFTPackCacheId, DFFTPackCache> dfftpack_cmgr(20);

void convolve_fftpack(int n,double* inout,double* omega,int swap_real_imag) 
{
        DFFTPackCache* cache;

        cache = dfftpack_cmgr.get_cache(RFFTPackCacheId(n));
        cache->convolve(inout, omega, swap_real_imag);
}

void convolve_z_fftpack(int n,double* inout,double* omega_real,double* omega_imag) 
{
        DFFTPackCache* cache;

        cache = dfftpack_cmgr.get_cache(RFFTPackCacheId(n));
        cache->convolve_z(inout, omega_real, omega_imag);
}

void init_convolution_kernel_fftpack(int n,double* omega, int d,
			     double (*kernel_func)(int),
			     int zero_nyquist) 
{
	/*
	 * omega[k] = pow(sqrt(-1),d) * kernel_func(k)
	 * omega[0] = kernel_func(0)
	 * conjugate(omega[-k]) == omega[k]
	 */
	int j,k,l=(n%2?n:n-1);
	omega[0] = (*kernel_func)(0)/n;
	switch (d%4) {
		case 0:
			for (k=j=1;j<l;j+=2,++k)
				omega[j] = omega[j+1] = (*kernel_func)(k)/n;
			if (!(n%2))
				omega[n-1] = (zero_nyquist?0.0:(*kernel_func)(k)/n);
			break;
		case 1:;case -3:
		       for (k=j=1;j<l;j+=2,++k) {
			       omega[j] = (*kernel_func)(k)/n;
			       omega[j+1] = -omega[j];
		       }
		       if (!(n%2))
			       omega[n-1] = (zero_nyquist?0.0:(*kernel_func)(k)/n);
		       break;
		case 2:;case -2:
		       for (k=j=1;j<l;j+=2,++k)
			       omega[j] = omega[j+1] = -(*kernel_func)(k)/n;
		       if (!(n%2))
			       omega[n-1] = (zero_nyquist?0.0:-(*kernel_func)(k)/n);
		       break;
		case 3:;case -1:
		       for (k=j=1;j<l;j+=2,++k) {
			       omega[j] = -(*kernel_func)(k)/n;
			       omega[j+1] = -omega[j];
		       }
		       if (!(n%2))
			       omega[n-1] = (zero_nyquist?0.0:-(*kernel_func)(k)/n);
		       break;
	}
}
