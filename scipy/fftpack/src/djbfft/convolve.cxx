#ifdef WITH_FFTW
#define destroy_convolve_cache_def destroy_convolve_cache_fftw
#define convolve_def convolve_fftw
#define convolve_z_def convolve_z_fftw
#define init_convolution_kernel_def init_convolution_kernel_fftw
#else
#define destroy_convolve_cache_def destroy_convolve_cache_fftpack
#define convolve_def convolve_fftpack
#define convolve_z_def convolve_z_fftpack
#define init_convolution_kernel_def init_convolution_kernel_fftpack
#endif

GEN_CACHE(ddjbfft, (int n)
	  , double *ptr;, (caches_ddjbfft[i].n == n)
	  , caches_ddjbfft[id].ptr =
	  (double *) malloc(sizeof(double) * n);,
	  free(caches_ddjbfft[id].ptr);, 20)

extern "C" void destroy_convolve_cache(void)
{
	destroy_ddjbfft_caches();
	destroy_convolve_cache_def();
}

/**************** convolve **********************/
extern "C"
    void convolve(int n, double *inout, double *omega, int swap_real_imag)
{
	int i;
	double *ptr = NULL;
	switch (n) {
	case 2:;
	case 4:;
	case 8:;
	case 16:;
	case 32:;
	case 64:;
	case 128:;
	case 256:;
	case 512:;
	case 1024:;
	case 2048:;
	case 4096:;
	case 8192:
		i = get_cache_id_ddjbfft(n);
		ptr = caches_ddjbfft[i].ptr;
		COPYSTD2DJB(inout, ptr, n);
		switch (n) {
#define TMPCASE(N) case N: fftr8_##N(ptr); break
			TMPCASE(2);
			TMPCASE(4);
			TMPCASE(8);
			TMPCASE(16);
			TMPCASE(32);
			TMPCASE(64);
			TMPCASE(128);
			TMPCASE(256);
			TMPCASE(512);
			TMPCASE(1024);
			TMPCASE(2048);
			TMPCASE(4096);
			TMPCASE(8192);
#undef TMPCASE
		}
		if (swap_real_imag) {
			int n1 = n - 1;
			double c;
			ptr[0] *= omega[0];
			ptr[1] *= omega[1];
			for (i = 2; i < n1; i += 2) {
				c = ptr[i] * omega[i];
				ptr[i] = ptr[i + 1] * omega[i + 1];
				ptr[i + 1] = c;
			}
		} else
			for (i = 0; i < n; ++i)
				ptr[i] *= omega[i];
		switch (n) {
#define TMPCASE(N)case N:fftr8_un##N(ptr);break
			TMPCASE(2);
			TMPCASE(4);
			TMPCASE(8);
			TMPCASE(16);
			TMPCASE(32);
			TMPCASE(64);
			TMPCASE(128);
			TMPCASE(256);
			TMPCASE(512);
			TMPCASE(1024);
			TMPCASE(2048);
			TMPCASE(4096);
			TMPCASE(8192);
#undef TMPCASE
		}
		COPYINVDJB2STD2(ptr, inout, n);
		return;
	}
	{
		convolve_def(n, inout, omega, swap_real_imag);
	}
}

/**************** convolve **********************/
extern "C"
    void convolve_z(int n, double *inout, double *omega_real,
		    double *omega_imag)
{
	int i;
	double *ptr = NULL;
	switch (n) {
	case 2:;
	case 4:;
	case 8:;
	case 16:;
	case 32:;
	case 64:;
	case 128:;
	case 256:;
	case 512:;
	case 1024:;
	case 2048:;
	case 4096:;
	case 8192:
		i = get_cache_id_ddjbfft(n);
		ptr = caches_ddjbfft[i].ptr;
		COPYSTD2DJB(inout, ptr, n);
		switch (n) {
#define TMPCASE(N) case N: fftr8_##N(ptr); break
			TMPCASE(2);
			TMPCASE(4);
			TMPCASE(8);
			TMPCASE(16);
			TMPCASE(32);
			TMPCASE(64);
			TMPCASE(128);
			TMPCASE(256);
			TMPCASE(512);
			TMPCASE(1024);
			TMPCASE(2048);
			TMPCASE(4096);
			TMPCASE(8192);
#undef TMPCASE
		}
		{
			int n1 = n - 1;
			double c;
			ptr[0] *= (omega_real[0] + omega_imag[0]);
			ptr[1] *= (omega_real[1] + omega_imag[1]);
			for (i = 2; i < n1; i += 2) {
				c = ptr[i] * omega_imag[i];
				ptr[i] *= omega_real[i];
				ptr[i] += ptr[i + 1] * omega_imag[i + 1];
				ptr[i + 1] *= omega_real[i + 1];
				ptr[i + 1] += c;
			}
		}
		switch (n) {
#define TMPCASE(N)case N:fftr8_un##N(ptr);break
			TMPCASE(2);
			TMPCASE(4);
			TMPCASE(8);
			TMPCASE(16);
			TMPCASE(32);
			TMPCASE(64);
			TMPCASE(128);
			TMPCASE(256);
			TMPCASE(512);
			TMPCASE(1024);
			TMPCASE(2048);
			TMPCASE(4096);
			TMPCASE(8192);
#undef TMPCASE
		}
		COPYINVDJB2STD2(ptr, inout, n);
		return;
	}
	{
		convolve_z_def(n, inout, omega_real, omega_imag);
	}
}

extern "C"
    void init_convolution_kernel(int n, double *omega, int d,
				 double (*kernel_func) (int),
				 int zero_nyquist)
{
	/*
	   omega[k] = pow(sqrt(-1),d) * kernel_func(k)
	   omega[0] = kernel_func(0)
	   conjugate(omega[-k]) == omega[k]
	 */
	switch (n) {
	case 2:;
	case 4:;
	case 8:;
	case 16:;
	case 32:;
	case 64:;
	case 128:;
	case 256:;
	case 512:;
	case 1024:;
	case 2048:;
	case 4096:;
	case 8192:
		{
			int k, n2 = n / 2;
			unsigned int *f =
			    (unsigned int *) malloc(sizeof(int) * (n));
			fftfreq_rtable(f, n);
			for (k = 1; k < n; ++k)
				if (f[k] > n2)
					f[k] -= n;
			omega[0] = (*kernel_func) (0) / n;
			switch (d % 4) {
			case 0:
				for (k = 2; k < n - 1; k += 2) {
					omega[k] =
					    (*kernel_func) (f[k]) / n2;
					omega[k + 1] = -omega[k];
				}
				omega[1] =
				    (zero_nyquist ? 0.0
				     : (*kernel_func) (n2) / n);
				break;
			case 1:;
			case -3:
				for (k = 2; k < n - 1; k += 2)
					omega[k] = omega[k + 1] =
					    -(*kernel_func) (f[k]) / n2;
				omega[1] =
				    (zero_nyquist ? 0.0
				     : (*kernel_func) (n2) / n);
				break;
			case 2:;
			case -2:
				for (k = 2; k < n - 1; k += 2) {
					omega[k] =
					    -(*kernel_func) (f[k]) / n2;
					omega[k + 1] = -omega[k];
				}
				omega[1] =
				    (zero_nyquist ? 0.0 :
				     -(*kernel_func) (n2) / n);
				break;
			case 3:;
			case -1:
				for (k = 2; k < n - 1; k += 2)
					omega[k] = omega[k + 1] =
					    (*kernel_func) (f[k]) / n2;
				omega[1] =
				    (zero_nyquist ? 0.0 :
				     -(*kernel_func) (n2) / n);
				break;
			}
			free(f);
		}
		return;
	}
	init_convolution_kernel_def(n, omega, d, kernel_func,
				    zero_nyquist);
}
