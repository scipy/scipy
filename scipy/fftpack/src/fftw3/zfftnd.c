/*
 * fftw3 backend for multi dimensional fft
 *
 * Original code by Pearu Peaterson
 *
 * Last Change: Wed Aug 08 02:00 PM 2007 J
 */

/* stub because fftw3 has no cache mechanism (yet) */
static void destroy_zfftnd_fftw3_caches(void) {}

extern void zfftnd_fftw3(complex_double * inout, int rank,
			   int *dims, int direction, int howmany,
			   int normalize)
{
    int i, sz;
    complex_double *ptr = inout;

    fftw_plan plan = NULL;
    sz = 1;
    for (i = 0; i < rank; ++i) {
        sz *= dims[i];
    }
    plan = fftw_plan_many_dft(rank, dims, howmany,
			      (fftw_complex *) ptr, NULL, 1, sz,
			      (fftw_complex *) ptr, NULL, 1, sz,
			      (direction >
			       0 ? FFTW_FORWARD : FFTW_BACKWARD),
			      FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    if (normalize) {
        ptr = inout;
        for (i = sz * howmany - 1; i >= 0; --i) {
            *((double *) (ptr)) /= sz;
            *((double *) (ptr++) + 1) /= sz;
        }
    }
}
