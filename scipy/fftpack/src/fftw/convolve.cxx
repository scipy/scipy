GEN_CACHE(drfftw,(int n)
	  ,rfftw_plan plan1;
	   rfftw_plan plan2;
	  ,(caches_drfftw[i].n==n)
	  ,caches_drfftw[id].plan1 = rfftw_create_plan(n,
		FFTW_REAL_TO_COMPLEX,
		FFTW_IN_PLACE|FFTW_ESTIMATE);
	   caches_drfftw[id].plan2 = rfftw_create_plan(n,
		FFTW_COMPLEX_TO_REAL,
		FFTW_IN_PLACE|FFTW_ESTIMATE);
	  ,rfftw_destroy_plan(caches_drfftw[id].plan1);
  	   rfftw_destroy_plan(caches_drfftw[id].plan2);
	  ,20)

static void destroy_convolve_cache_fftw(void) 
{
  destroy_drfftw_caches();
}

/**************** convolve **********************/
static
void convolve_fftw(int n,double* inout,double* omega,int swap_real_imag) 
{
	int i;
	rfftw_plan plan1 = NULL;
	rfftw_plan plan2 = NULL;
	int l = (n-1)/2+1;

	i = get_cache_id_drfftw(n);
	plan1 = caches_drfftw[i].plan1;
	plan2 = caches_drfftw[i].plan2;
	rfftw_one(plan1, (fftw_real *)inout, NULL);
	if (swap_real_imag) {
		double c;
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
	rfftw_one(plan2, (fftw_real *)inout, NULL);
}

/**************** convolve **********************/
static
void convolve_z_fftw(int n,double* inout,double* omega_real,double* omega_imag) {
	int i;
	rfftw_plan plan1 = NULL;
	rfftw_plan plan2 = NULL;
	int l = (n-1)/2+1;
	double c;

	i = get_cache_id_drfftw(n);
	plan1 = caches_drfftw[i].plan1;
	plan2 = caches_drfftw[i].plan2;
	rfftw_one(plan1, (fftw_real *)inout, NULL);
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
	rfftw_one(plan2, (fftw_real *)inout, NULL);
}

static
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
