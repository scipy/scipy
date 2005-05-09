
#include <Python.h>
#include <stdio.h>
#include <math.h>
#include "libnumarray.h"

#ifdef MS_WIN32
#pragma warning(once : 4244)
#endif

#define logical_and(arg1, arg2) (arg1 != 0) & (arg2 != 0)
#define logical_or(arg1, arg2)  (arg1 != 0) | (arg2 != 0)
#define logical_xor(arg1, arg2) ((arg1 != 0) ^ (arg2 != 0)) & 1
#define ufmaximum(arg1, arg2) (((temp1=arg1) > (temp2=arg2)) ? temp1 : temp2)
#define ufminimum(arg1, arg2) (((temp1=arg1) < (temp2=arg2)) ? temp1 : temp2)

#define distance3d(x,y,z) sqrt(x*x + y*y + z*z)

#include "cephes.h"

double besselpoly(double a, double lambda, double nu);

Complex64 cbesi_wrap( double v, Complex64 z);
Complex64 cbesi_wrap_e( double v, Complex64 z);
Complex64 cbesj_wrap( double v, Complex64 z);
Complex64 cbesj_wrap_e( double v, Complex64 z);
Complex64 cbesy_wrap( double v, Complex64 z);
Complex64 cbesy_wrap_e( double v, Complex64 z);
Complex64 cbesk_wrap( double v, Complex64 z);
Complex64 cbesk_wrap_e( double v, Complex64 z);  
Complex64 cbesh_wrap1( double v, Complex64 z);
Complex64 cbesh_wrap1_e( double v, Complex64 z);  
Complex64 cbesh_wrap2( double v, Complex64 z);
Complex64 cbesh_wrap2_e( double v, Complex64 z);

extern double cdfbet3_wrap(double p, double x, double b);
extern double cdfbet4_wrap(double p, double x, double a);

extern double cdfbin2_wrap(double p, double xn, double pr);
extern double cdfbin3_wrap(double p, double s, double pr);

extern double cdfchi3_wrap(double p, double x);

extern double cdfchn1_wrap(double x, double df, double nc);
extern double cdfchn2_wrap(double p, double df, double nc);
extern double cdfchn3_wrap(double p, double x, double nc);
extern double cdfchn4_wrap(double p, double x, double df);

extern double cdff3_wrap(double p, double f, double dfd);
extern double cdff4_wrap(double p, double f, double dfn);

extern double cdffnc1_wrap(double f, double dfn, double dfd, double nc);
extern double cdffnc2_wrap(double p, double dfn, double dfd, double nc);
extern double cdffnc3_wrap(double p, double f, double dfd, double nc);
extern double cdffnc4_wrap(double p, double f, double dfn, double nc);
extern double cdffnc5_wrap(double p, double f, double dfn, double dfd);

extern double cdfgam1_wrap(double p, double x, double scl);
extern double cdfgam2_wrap(double p, double x, double shp);
extern double cdfgam3_wrap(double p, double x, double scl);
extern double cdfgam4_wrap(double p, double x, double shp);

extern double cdfnbn2_wrap(double p, double xn, double pr);
extern double cdfnbn3_wrap(double p, double s, double pr);

extern double cdfnor3_wrap(double p, double x, double std);
extern double cdfnor4_wrap(double p, double x, double mn);

extern double cdfpoi2_wrap(double p, double xlam);

extern double cdft1_wrap(double p, double t);
extern double cdft2_wrap(double p, double t);
extern double cdft3_wrap(double p, double t);

extern double cdftnc1_wrap(double df, double nc, double t);
extern double cdftnc2_wrap(double df, double nc, double p);
extern double cdftnc3_wrap(double p, double nc, double t);
extern double cdftnc4_wrap(double df, double p, double t);

extern double tukeylambdacdf(double x, double lambda);

Complex64 cgamma_wrap( Complex64 z);
Complex64 clngamma_wrap( Complex64 z);
Complex64 cpsi_wrap( Complex64 z);
Complex64 crgamma_wrap( Complex64 z);
Complex64 chyp2f1_wrap( double a, double b, double c, Complex64 z);
Complex64 chyp1f1_wrap( double a, double b, Complex64 z);
double hypU_wrap(double a, double b, double x);
double exp1_wrap(double x);
double expi_wrap(double x);
Complex64 cexp1_wrap( Complex64 z);
Complex64 cerf_wrap( Complex64 z);
int itairy_wrap(double x, double *apt, double *bpt, double *ant, double *bnt);

double struve_wrap(double v, double x);
double itstruve0_wrap(double x);
double it2struve0_wrap(double x);

double modstruve_wrap(double v, double x);
double itmodstruve0_wrap(double x);

double ber_wrap(double x);
double bei_wrap(double x);
double ker_wrap(double x);
double kei_wrap(double x);
double berp_wrap(double x);
double beip_wrap(double x);
double kerp_wrap(double x);
double keip_wrap(double x);

int kelvin_wrap(double x, Complex64 *Be, Complex64 *Ke, Complex64 *Bep, Complex64 *Kep);

int it1j0y0_wrap(double x, double *, double *);
int it2j0y0_wrap(double x, double *, double *);
int it1i0k0_wrap(double x, double *, double *);
int it2i0k0_wrap(double x, double *, double *);

int cfresnl_wrap(Complex64 x, Complex64 *sf, Complex64 *cf);
double cem_cva_wrap(double m, double q);
double sem_cva_wrap(double m, double q);
int cem_wrap(double m, double q, double x, double *csf, double *csd);
int sem_wrap(double m, double q, double x, double *csf, double *csd);
int mcm1_wrap(double m, double q, double x, double *f1r, double *d1r);
int msm1_wrap(double m, double q, double x, double *f1r, double *d1r);
int mcm2_wrap(double m, double q, double x, double *f2r, double *d2r);
int msm2_wrap(double m, double q, double x, double *f2r, double *d2r);
double pmv_wrap(double, double, double);
int pbwa_wrap(double, double, double *, double *);
int pbdv_wrap(double, double, double *, double *);
int pbvv_wrap(double, double, double *, double *);

int prolate_aswfa_wrap(double, double, double, double, double, double *, double *);
int prolate_radial1_wrap(double, double, double, double, double, double *, double *);
int prolate_radial2_wrap(double, double, double, double, double, double *, double *);

/*
int oblate_aswfa_wrap(double, double, double, double, double, double *, double *);
int oblate_radial1_wrap(double, double, double, double, double, double *, double *);
int oblate_radial2_wrap(double, double, double, double, double, double *, double *);
double prolate_aswfa_nocv_wrap(double, double, double, double, double *);
double prolate_radial1_nocv_wrap(double, double, double, double, double *);
double prolate_radial2_nocv_wrap(double, double, double, double, double *);
double oblate_aswfa_nocv_wrap(double, double, double, double, double *);
double oblate_radial1_nocv_wrap(double, double, double, double, double *);
double oblate_radial2_nocv_wrap(double, double, double, double, double *);
*/

double prolate_segv_wrap(double, double, double);
double oblate_segv_wrap(double, double, double);


int modified_fresnel_plus_wrap(double x, Complex64 *F, Complex64 *K);
int modified_fresnel_minus_wrap(double x, Complex64 *F, Complex64 *K);

extern Complex64 cwofz_wrap(Complex64 z);


static int airy_fxffff_vxvvvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float32    *tout1 =  (Float32 *) buffers[2];
    Float32    *tout2 =  (Float32 *) buffers[3];
    Float32    *tout3 =  (Float32 *) buffers[4];
    Float64     result0;
    Float64     result1;
    Float64     result2;
    Float64     result3;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        airy(*tin0, &result0, &result1, &result2, &result3);
                *tout0 = result0;
        *tout1 = result1;
        *tout2 = result2;
        *tout3 = result3;

	++tin0; ++tout0; ++tout1; ++tout2; ++tout3; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor airy_fxffff_vxvvvv_descr =
{ "airy_fxffff_vxvvvv", (void *) airy_fxffff_vxvvvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 4,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0 } };

static int airy_dxdddd_vxvvvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];
    Float64    *tout1 =  (Float64 *) buffers[2];
    Float64    *tout2 =  (Float64 *) buffers[3];
    Float64    *tout3 =  (Float64 *) buffers[4];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        airy(*tin0, tout0, tout1, tout2, tout3);
        
	++tin0; ++tout0; ++tout1; ++tout2; ++tout3; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor airy_dxdddd_vxvvvv_descr =
{ "airy_dxdddd_vxvvvv", (void *) airy_dxdddd_vxvvvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 4,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0 } };

static int airy_FxFFFF_vxvvvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex32  *tin0 =  (Complex32 *) buffers[0];
    Complex32  *tout0 =  (Complex32 *) buffers[1];
    Complex32  *tout1 =  (Complex32 *) buffers[2];
    Complex32  *tout2 =  (Complex32 *) buffers[3];
    Complex32  *tout3 =  (Complex32 *) buffers[4];
    Complex64     input0;
    Complex64     result0;
    Complex64     result1;
    Complex64     result2;
    Complex64     result3;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input0.r = tin0->r;
        input0.i = tin0->i;

        cairy_wrap(input0, &result0, &result1, &result2, &result3);
                tout0->r = result0.r;
        tout0->i = result0.i;
        tout1->r = result1.r;
        tout1->i = result1.i;
        tout2->r = result2.r;
        tout2->i = result2.i;
        tout3->r = result3.r;
        tout3->i = result3.i;

	++tin0; ++tout0; ++tout1; ++tout2; ++tout3; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor airy_FxFFFF_vxvvvv_descr =
{ "airy_FxFFFF_vxvvvv", (void *) airy_FxFFFF_vxvvvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 4,
  { sizeof(Complex32), sizeof(Complex32), sizeof(Complex32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0, 0, 0 } };

static int airy_DxDDDD_vxvvvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex64  *tin0 =  (Complex64 *) buffers[0];
    Complex64  *tout0 =  (Complex64 *) buffers[1];
    Complex64  *tout1 =  (Complex64 *) buffers[2];
    Complex64  *tout2 =  (Complex64 *) buffers[3];
    Complex64  *tout3 =  (Complex64 *) buffers[4];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        cairy_wrap(*tin0, tout0, tout1, tout2, tout3);
        
	++tin0; ++tout0; ++tout1; ++tout2; ++tout3; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor airy_DxDDDD_vxvvvv_descr =
{ "airy_DxDDDD_vxvvvv", (void *) airy_DxDDDD_vxvvvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 4,
  { sizeof(Complex64), sizeof(Complex64), sizeof(Complex64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0, 0, 0 } };

static int airye_fxffff_vxvvvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float32    *tout1 =  (Float32 *) buffers[2];
    Float32    *tout2 =  (Float32 *) buffers[3];
    Float32    *tout3 =  (Float32 *) buffers[4];
    Complex64     input0;
    Complex64     result0;
    Complex64     result1;
    Complex64     result2;
    Complex64     result3;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input0.r = *tin0;
        input0.i = 0;

        cairy_wrap_e(input0, &result0, &result1, &result2, &result3);
                *tout0 = result0.r;
        *tout1 = result1.r;
        *tout2 = result2.r;
        *tout3 = result3.r;

	++tin0; ++tout0; ++tout1; ++tout2; ++tout3; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor airye_fxffff_vxvvvv_descr =
{ "airye_fxffff_vxvvvv", (void *) airye_fxffff_vxvvvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 4,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0 } };

static int airye_dxdddd_vxvvvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];
    Float64    *tout1 =  (Float64 *) buffers[2];
    Float64    *tout2 =  (Float64 *) buffers[3];
    Float64    *tout3 =  (Float64 *) buffers[4];
    Complex64     input0;
    Complex64     result0;
    Complex64     result1;
    Complex64     result2;
    Complex64     result3;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input0.r = *tin0;
        input0.i = 0;

        cairy_wrap_e(input0, &result0, &result1, &result2, &result3);
                *tout0 = result0.r;
        *tout1 = result1.r;
        *tout2 = result2.r;
        *tout3 = result3.r;

	++tin0; ++tout0; ++tout1; ++tout2; ++tout3; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor airye_dxdddd_vxvvvv_descr =
{ "airye_dxdddd_vxvvvv", (void *) airye_dxdddd_vxvvvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 4,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0 } };

static int airye_FxFFFF_vxvvvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex32  *tin0 =  (Complex32 *) buffers[0];
    Complex32  *tout0 =  (Complex32 *) buffers[1];
    Complex32  *tout1 =  (Complex32 *) buffers[2];
    Complex32  *tout2 =  (Complex32 *) buffers[3];
    Complex32  *tout3 =  (Complex32 *) buffers[4];
    Complex64     input0;
    Complex64     result0;
    Complex64     result1;
    Complex64     result2;
    Complex64     result3;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input0.r = tin0->r;
        input0.i = tin0->i;

        cairy_wrap_e(input0, &result0, &result1, &result2, &result3);
                tout0->r = result0.r;
        tout0->i = result0.i;
        tout1->r = result1.r;
        tout1->i = result1.i;
        tout2->r = result2.r;
        tout2->i = result2.i;
        tout3->r = result3.r;
        tout3->i = result3.i;

	++tin0; ++tout0; ++tout1; ++tout2; ++tout3; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor airye_FxFFFF_vxvvvv_descr =
{ "airye_FxFFFF_vxvvvv", (void *) airye_FxFFFF_vxvvvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 4,
  { sizeof(Complex32), sizeof(Complex32), sizeof(Complex32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0, 0, 0 } };

static int airye_DxDDDD_vxvvvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex64  *tin0 =  (Complex64 *) buffers[0];
    Complex64  *tout0 =  (Complex64 *) buffers[1];
    Complex64  *tout1 =  (Complex64 *) buffers[2];
    Complex64  *tout2 =  (Complex64 *) buffers[3];
    Complex64  *tout3 =  (Complex64 *) buffers[4];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        cairy_wrap_e(*tin0, tout0, tout1, tout2, tout3);
        
	++tin0; ++tout0; ++tout1; ++tout2; ++tout3; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor airye_DxDDDD_vxvvvv_descr =
{ "airye_DxDDDD_vxvvvv", (void *) airye_DxDDDD_vxvvvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 4,
  { sizeof(Complex64), sizeof(Complex64), sizeof(Complex64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0, 0, 0 } };

static int bdtr_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = bdtr(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor bdtr_fffxf_vvvxf_descr =
{ "bdtr_fffxf_vvvxf", (void *) bdtr_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int bdtr_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = bdtr(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor bdtr_dddxd_vvvxf_descr =
{ "bdtr_dddxd_vvvxf", (void *) bdtr_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int bdtrc_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = bdtrc(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor bdtrc_fffxf_vvvxf_descr =
{ "bdtrc_fffxf_vvvxf", (void *) bdtrc_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int bdtrc_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = bdtrc(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor bdtrc_dddxd_vvvxf_descr =
{ "bdtrc_dddxd_vvvxf", (void *) bdtrc_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int bdtri_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = bdtri(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor bdtri_fffxf_vvvxf_descr =
{ "bdtri_fffxf_vvvxf", (void *) bdtri_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int bdtri_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = bdtri(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor bdtri_dddxd_vvvxf_descr =
{ "bdtri_dddxd_vvvxf", (void *) bdtri_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int bdtrik_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdfbin2_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor bdtrik_fffxf_vvvxf_descr =
{ "bdtrik_fffxf_vvvxf", (void *) bdtrik_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int bdtrik_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdfbin2_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor bdtrik_dddxd_vvvxf_descr =
{ "bdtrik_dddxd_vvvxf", (void *) bdtrik_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int bdtrin_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdfbin3_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor bdtrin_fffxf_vvvxf_descr =
{ "bdtrin_fffxf_vvvxf", (void *) bdtrin_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int bdtrin_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdfbin3_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor bdtrin_dddxd_vvvxf_descr =
{ "bdtrin_dddxd_vvvxf", (void *) bdtrin_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int bei_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = bei_wrap(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor bei_fxf_vxf_descr =
{ "bei_fxf_vxf", (void *) bei_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int bei_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = bei_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor bei_dxd_vxf_descr =
{ "bei_dxd_vxf", (void *) bei_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int beip_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = beip_wrap(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor beip_fxf_vxf_descr =
{ "beip_fxf_vxf", (void *) beip_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int beip_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = beip_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor beip_dxd_vxf_descr =
{ "beip_dxd_vxf", (void *) beip_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int ber_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = ber_wrap(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ber_fxf_vxf_descr =
{ "ber_fxf_vxf", (void *) ber_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int ber_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = ber_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ber_dxd_vxf_descr =
{ "ber_dxd_vxf", (void *) ber_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int berp_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = berp_wrap(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor berp_fxf_vxf_descr =
{ "berp_fxf_vxf", (void *) berp_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int berp_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = berp_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor berp_dxd_vxf_descr =
{ "berp_dxd_vxf", (void *) berp_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int besselpoly_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = besselpoly(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor besselpoly_fffxf_vvvxf_descr =
{ "besselpoly_fffxf_vvvxf", (void *) besselpoly_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int besselpoly_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = besselpoly(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor besselpoly_dddxd_vvvxf_descr =
{ "besselpoly_dddxd_vvvxf", (void *) besselpoly_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int beta_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = beta(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor beta_ffxf_vvxf_descr =
{ "beta_ffxf_vvxf", (void *) beta_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int beta_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = beta(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor beta_ddxd_vvxf_descr =
{ "beta_ddxd_vvxf", (void *) beta_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int betainc_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = incbet(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor betainc_fffxf_vvvxf_descr =
{ "betainc_fffxf_vvvxf", (void *) betainc_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int betainc_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = incbet(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor betainc_dddxd_vvvxf_descr =
{ "betainc_dddxd_vvvxf", (void *) betainc_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int betaincinv_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = incbi(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor betaincinv_fffxf_vvvxf_descr =
{ "betaincinv_fffxf_vvvxf", (void *) betaincinv_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int betaincinv_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = incbi(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor betaincinv_dddxd_vvvxf_descr =
{ "betaincinv_dddxd_vvvxf", (void *) betaincinv_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int betaln_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = lbeta(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor betaln_ffxf_vvxf_descr =
{ "betaln_ffxf_vvxf", (void *) betaln_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int betaln_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = lbeta(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor betaln_ddxd_vvxf_descr =
{ "betaln_ddxd_vvxf", (void *) betaln_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int btdtr_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = btdtr(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor btdtr_fffxf_vvvxf_descr =
{ "btdtr_fffxf_vvvxf", (void *) btdtr_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int btdtr_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = btdtr(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor btdtr_dddxd_vvvxf_descr =
{ "btdtr_dddxd_vvvxf", (void *) btdtr_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int btdtri_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = incbi(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor btdtri_fffxf_vvvxf_descr =
{ "btdtri_fffxf_vvvxf", (void *) btdtri_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int btdtri_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = incbi(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor btdtri_dddxd_vvvxf_descr =
{ "btdtri_dddxd_vvvxf", (void *) btdtri_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int btdtria_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdfbet3_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor btdtria_fffxf_vvvxf_descr =
{ "btdtria_fffxf_vvvxf", (void *) btdtria_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int btdtria_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdfbet3_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor btdtria_dddxd_vvvxf_descr =
{ "btdtria_dddxd_vvvxf", (void *) btdtria_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int btdtrib_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdfbet4_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor btdtrib_fffxf_vvvxf_descr =
{ "btdtrib_fffxf_vvvxf", (void *) btdtrib_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int btdtrib_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdfbet4_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor btdtrib_dddxd_vvvxf_descr =
{ "btdtrib_dddxd_vvvxf", (void *) btdtrib_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int cbrt_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cbrt(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor cbrt_fxf_vxf_descr =
{ "cbrt_fxf_vxf", (void *) cbrt_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int cbrt_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cbrt(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor cbrt_dxd_vxf_descr =
{ "cbrt_dxd_vxf", (void *) cbrt_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int chdtr_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = chdtr(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor chdtr_ffxf_vvxf_descr =
{ "chdtr_ffxf_vvxf", (void *) chdtr_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int chdtr_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = chdtr(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor chdtr_ddxd_vvxf_descr =
{ "chdtr_ddxd_vvxf", (void *) chdtr_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int chdtrc_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = chdtrc(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor chdtrc_ffxf_vvxf_descr =
{ "chdtrc_ffxf_vvxf", (void *) chdtrc_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int chdtrc_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = chdtrc(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor chdtrc_ddxd_vvxf_descr =
{ "chdtrc_ddxd_vvxf", (void *) chdtrc_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int chdtri_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = chdtri(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor chdtri_ffxf_vvxf_descr =
{ "chdtri_ffxf_vvxf", (void *) chdtri_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int chdtri_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = chdtri(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor chdtri_ddxd_vvxf_descr =
{ "chdtri_ddxd_vvxf", (void *) chdtri_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int chdtriv_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdfchi3_wrap(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor chdtriv_ffxf_vvxf_descr =
{ "chdtriv_ffxf_vvxf", (void *) chdtriv_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int chdtriv_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdfchi3_wrap(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor chdtriv_ddxd_vvxf_descr =
{ "chdtriv_ddxd_vvxf", (void *) chdtriv_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int chndtr_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdfchn1_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor chndtr_fffxf_vvvxf_descr =
{ "chndtr_fffxf_vvvxf", (void *) chndtr_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int chndtr_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdfchn1_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor chndtr_dddxd_vvvxf_descr =
{ "chndtr_dddxd_vvvxf", (void *) chndtr_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int chndtridf_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdfchn3_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor chndtridf_fffxf_vvvxf_descr =
{ "chndtridf_fffxf_vvvxf", (void *) chndtridf_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int chndtridf_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdfchn3_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor chndtridf_dddxd_vvvxf_descr =
{ "chndtridf_dddxd_vvvxf", (void *) chndtridf_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int chndtrinc_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdfchn4_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor chndtrinc_fffxf_vvvxf_descr =
{ "chndtrinc_fffxf_vvvxf", (void *) chndtrinc_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int chndtrinc_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdfchn4_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor chndtrinc_dddxd_vvvxf_descr =
{ "chndtrinc_dddxd_vvvxf", (void *) chndtrinc_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int chndtrix_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdfchn2_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor chndtrix_fffxf_vvvxf_descr =
{ "chndtrix_fffxf_vvvxf", (void *) chndtrix_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int chndtrix_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdfchn2_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor chndtrix_dddxd_vvvxf_descr =
{ "chndtrix_dddxd_vvvxf", (void *) chndtrix_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int cosdg_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cosdg(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor cosdg_fxf_vxf_descr =
{ "cosdg_fxf_vxf", (void *) cosdg_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int cosdg_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cosdg(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor cosdg_dxd_vxf_descr =
{ "cosdg_dxd_vxf", (void *) cosdg_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int cosm1_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cosm1(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor cosm1_fxf_vxf_descr =
{ "cosm1_fxf_vxf", (void *) cosm1_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int cosm1_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cosm1(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor cosm1_dxd_vxf_descr =
{ "cosm1_dxd_vxf", (void *) cosm1_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int cotdg_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cotdg(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor cotdg_fxf_vxf_descr =
{ "cotdg_fxf_vxf", (void *) cotdg_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int cotdg_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cotdg(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor cotdg_dxd_vxf_descr =
{ "cotdg_dxd_vxf", (void *) cotdg_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int dawsn_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = dawsn(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor dawsn_fxf_vxf_descr =
{ "dawsn_fxf_vxf", (void *) dawsn_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int dawsn_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = dawsn(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor dawsn_dxd_vxf_descr =
{ "dawsn_dxd_vxf", (void *) dawsn_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int ellipe_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = ellpe(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ellipe_fxf_vxf_descr =
{ "ellipe_fxf_vxf", (void *) ellipe_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int ellipe_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = ellpe(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ellipe_dxd_vxf_descr =
{ "ellipe_dxd_vxf", (void *) ellipe_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int ellipeinc_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = ellie(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ellipeinc_ffxf_vvxf_descr =
{ "ellipeinc_ffxf_vvxf", (void *) ellipeinc_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int ellipeinc_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = ellie(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ellipeinc_ddxd_vvxf_descr =
{ "ellipeinc_ddxd_vvxf", (void *) ellipeinc_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int ellipj_ffxffff_vvxvvvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float32    *tout1 =  (Float32 *) buffers[3];
    Float32    *tout2 =  (Float32 *) buffers[4];
    Float32    *tout3 =  (Float32 *) buffers[5];
    Float64     result0;
    Float64     result1;
    Float64     result2;
    Float64     result3;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        ellpj(*tin0, *tin1, &result0, &result1, &result2, &result3);
                *tout0 = result0;
        *tout1 = result1;
        *tout2 = result2;
        *tout3 = result3;

	++tin0; ++tin1; ++tout0; ++tout1; ++tout2; ++tout3; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ellipj_ffxffff_vvxvvvv_descr =
{ "ellipj_ffxffff_vvxvvvv", (void *) ellipj_ffxffff_vvxvvvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 4,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int ellipj_ddxdddd_vvxvvvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];
    Float64    *tout1 =  (Float64 *) buffers[3];
    Float64    *tout2 =  (Float64 *) buffers[4];
    Float64    *tout3 =  (Float64 *) buffers[5];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        ellpj(*tin0, *tin1, tout0, tout1, tout2, tout3);
        
	++tin0; ++tin1; ++tout0; ++tout1; ++tout2; ++tout3; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ellipj_ddxdddd_vvxvvvv_descr =
{ "ellipj_ddxdddd_vvxvvvv", (void *) ellipj_ddxdddd_vvxvvvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 4,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int ellipk_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = ellpk(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ellipk_fxf_vxf_descr =
{ "ellipk_fxf_vxf", (void *) ellipk_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int ellipk_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = ellpk(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ellipk_dxd_vxf_descr =
{ "ellipk_dxd_vxf", (void *) ellipk_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int ellipkinc_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = ellik(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ellipkinc_ffxf_vvxf_descr =
{ "ellipkinc_ffxf_vvxf", (void *) ellipkinc_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int ellipkinc_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = ellik(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ellipkinc_ddxd_vvxf_descr =
{ "ellipkinc_ddxd_vvxf", (void *) ellipkinc_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int erf_FxF_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex32  *tin0 =  (Complex32 *) buffers[0];
    Complex32  *tout0 =  (Complex32 *) buffers[1];
    Complex64     input0;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input0.r = tin0->r;
        input0.i = tin0->i;

        result0 = cerf_wrap(input0);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor erf_FxF_vxf_descr =
{ "erf_FxF_vxf", (void *) erf_FxF_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0 } };

static int erf_DxD_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex64  *tin0 =  (Complex64 *) buffers[0];
    Complex64  *tout0 =  (Complex64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cerf_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor erf_DxD_vxf_descr =
{ "erf_DxD_vxf", (void *) erf_DxD_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0 } };

static int erf_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = erf(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor erf_fxf_vxf_descr =
{ "erf_fxf_vxf", (void *) erf_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int erf_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = erf(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor erf_dxd_vxf_descr =
{ "erf_dxd_vxf", (void *) erf_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int erfc_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = erfc(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor erfc_fxf_vxf_descr =
{ "erfc_fxf_vxf", (void *) erfc_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int erfc_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = erfc(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor erfc_dxd_vxf_descr =
{ "erfc_dxd_vxf", (void *) erfc_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int exp1_FxF_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex32  *tin0 =  (Complex32 *) buffers[0];
    Complex32  *tout0 =  (Complex32 *) buffers[1];
    Complex64     input0;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input0.r = tin0->r;
        input0.i = tin0->i;

        result0 = cexp1_wrap(input0);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor exp1_FxF_vxf_descr =
{ "exp1_FxF_vxf", (void *) exp1_FxF_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0 } };

static int exp1_DxD_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex64  *tin0 =  (Complex64 *) buffers[0];
    Complex64  *tout0 =  (Complex64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cexp1_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor exp1_DxD_vxf_descr =
{ "exp1_DxD_vxf", (void *) exp1_DxD_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0 } };

static int exp1_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = exp1_wrap(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor exp1_fxf_vxf_descr =
{ "exp1_fxf_vxf", (void *) exp1_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int exp1_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = exp1_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor exp1_dxd_vxf_descr =
{ "exp1_dxd_vxf", (void *) exp1_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int exp10_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = exp10(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor exp10_fxf_vxf_descr =
{ "exp10_fxf_vxf", (void *) exp10_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int exp10_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = exp10(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor exp10_dxd_vxf_descr =
{ "exp10_dxd_vxf", (void *) exp10_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int exp2_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = exp2(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor exp2_fxf_vxf_descr =
{ "exp2_fxf_vxf", (void *) exp2_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int exp2_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = exp2(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor exp2_dxd_vxf_descr =
{ "exp2_dxd_vxf", (void *) exp2_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int expi_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = expi_wrap(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor expi_fxf_vxf_descr =
{ "expi_fxf_vxf", (void *) expi_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int expi_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = expi_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor expi_dxd_vxf_descr =
{ "expi_dxd_vxf", (void *) expi_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int expm1_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = expm1(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor expm1_fxf_vxf_descr =
{ "expm1_fxf_vxf", (void *) expm1_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int expm1_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = expm1(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor expm1_dxd_vxf_descr =
{ "expm1_dxd_vxf", (void *) expm1_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int expn_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = expn(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor expn_ffxf_vvxf_descr =
{ "expn_ffxf_vvxf", (void *) expn_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int expn_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = expn(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor expn_ddxd_vvxf_descr =
{ "expn_ddxd_vvxf", (void *) expn_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int fdtr_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = fdtr(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor fdtr_fffxf_vvvxf_descr =
{ "fdtr_fffxf_vvvxf", (void *) fdtr_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int fdtr_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = fdtr(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor fdtr_dddxd_vvvxf_descr =
{ "fdtr_dddxd_vvvxf", (void *) fdtr_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int fdtrc_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = fdtrc(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor fdtrc_fffxf_vvvxf_descr =
{ "fdtrc_fffxf_vvvxf", (void *) fdtrc_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int fdtrc_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = fdtrc(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor fdtrc_dddxd_vvvxf_descr =
{ "fdtrc_dddxd_vvvxf", (void *) fdtrc_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int fdtri_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = fdtri(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor fdtri_fffxf_vvvxf_descr =
{ "fdtri_fffxf_vvvxf", (void *) fdtri_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int fdtri_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = fdtri(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor fdtri_dddxd_vvvxf_descr =
{ "fdtri_dddxd_vvvxf", (void *) fdtri_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int fdtridfd_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdff4_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor fdtridfd_fffxf_vvvxf_descr =
{ "fdtridfd_fffxf_vvvxf", (void *) fdtridfd_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int fdtridfd_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdff4_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor fdtridfd_dddxd_vvvxf_descr =
{ "fdtridfd_dddxd_vvvxf", (void *) fdtridfd_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int fdtridfn_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdff3_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor fdtridfn_fffxf_vvvxf_descr =
{ "fdtridfn_fffxf_vvvxf", (void *) fdtridfn_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int fdtridfn_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdff3_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor fdtridfn_dddxd_vvvxf_descr =
{ "fdtridfn_dddxd_vvvxf", (void *) fdtridfn_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int fresnel_FxFF_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex32  *tin0 =  (Complex32 *) buffers[0];
    Complex32  *tout0 =  (Complex32 *) buffers[1];
    Complex32  *tout1 =  (Complex32 *) buffers[2];
    Complex64     input0;
    Complex64     result0;
    Complex64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input0.r = tin0->r;
        input0.i = tin0->i;

        cfresnl_wrap(input0, &result0, &result1);
                tout0->r = result0.r;
        tout0->i = result0.i;
        tout1->r = result1.r;
        tout1->i = result1.i;

	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor fresnel_FxFF_vxvv_descr =
{ "fresnel_FxFF_vxvv", (void *) fresnel_FxFF_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Complex32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0 } };

static int fresnel_DxDD_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex64  *tin0 =  (Complex64 *) buffers[0];
    Complex64  *tout0 =  (Complex64 *) buffers[1];
    Complex64  *tout1 =  (Complex64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        cfresnl_wrap(*tin0, tout0, tout1);
        
	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor fresnel_DxDD_vxvv_descr =
{ "fresnel_DxDD_vxvv", (void *) fresnel_DxDD_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Complex64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0 } };

static int fresnel_fxff_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float32    *tout1 =  (Float32 *) buffers[2];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        fresnl(*tin0, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor fresnel_fxff_vxvv_descr =
{ "fresnel_fxff_vxvv", (void *) fresnel_fxff_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int fresnel_dxdd_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];
    Float64    *tout1 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        fresnl(*tin0, tout0, tout1);
        
	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor fresnel_dxdd_vxvv_descr =
{ "fresnel_dxdd_vxvv", (void *) fresnel_dxdd_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int gamma_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = Gamma(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gamma_fxf_vxf_descr =
{ "gamma_fxf_vxf", (void *) gamma_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int gamma_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = Gamma(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gamma_dxd_vxf_descr =
{ "gamma_dxd_vxf", (void *) gamma_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int gamma_FxF_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex32  *tin0 =  (Complex32 *) buffers[0];
    Complex32  *tout0 =  (Complex32 *) buffers[1];
    Complex64     input0;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input0.r = tin0->r;
        input0.i = tin0->i;

        result0 = cgamma_wrap(input0);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gamma_FxF_vxf_descr =
{ "gamma_FxF_vxf", (void *) gamma_FxF_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0 } };

static int gamma_DxD_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex64  *tin0 =  (Complex64 *) buffers[0];
    Complex64  *tout0 =  (Complex64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cgamma_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gamma_DxD_vxf_descr =
{ "gamma_DxD_vxf", (void *) gamma_DxD_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0 } };

static int gammainc_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = igam(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gammainc_ffxf_vvxf_descr =
{ "gammainc_ffxf_vvxf", (void *) gammainc_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int gammainc_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = igam(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gammainc_ddxd_vvxf_descr =
{ "gammainc_ddxd_vvxf", (void *) gammainc_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int gammaincc_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = igamc(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gammaincc_ffxf_vvxf_descr =
{ "gammaincc_ffxf_vvxf", (void *) gammaincc_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int gammaincc_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = igamc(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gammaincc_ddxd_vvxf_descr =
{ "gammaincc_ddxd_vvxf", (void *) gammaincc_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int gammainccinv_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = igami(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gammainccinv_ffxf_vvxf_descr =
{ "gammainccinv_ffxf_vvxf", (void *) gammainccinv_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int gammainccinv_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = igami(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gammainccinv_ddxd_vvxf_descr =
{ "gammainccinv_ddxd_vvxf", (void *) gammainccinv_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int gammaln_FxF_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex32  *tin0 =  (Complex32 *) buffers[0];
    Complex32  *tout0 =  (Complex32 *) buffers[1];
    Complex64     input0;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input0.r = tin0->r;
        input0.i = tin0->i;

        result0 = clngamma_wrap(input0);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gammaln_FxF_vxf_descr =
{ "gammaln_FxF_vxf", (void *) gammaln_FxF_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0 } };

static int gammaln_DxD_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex64  *tin0 =  (Complex64 *) buffers[0];
    Complex64  *tout0 =  (Complex64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = clngamma_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gammaln_DxD_vxf_descr =
{ "gammaln_DxD_vxf", (void *) gammaln_DxD_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0 } };

static int gammaln_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = lgam(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gammaln_fxf_vxf_descr =
{ "gammaln_fxf_vxf", (void *) gammaln_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int gammaln_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = lgam(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gammaln_dxd_vxf_descr =
{ "gammaln_dxd_vxf", (void *) gammaln_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int gdtr_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = gdtr(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gdtr_fffxf_vvvxf_descr =
{ "gdtr_fffxf_vvvxf", (void *) gdtr_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int gdtr_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = gdtr(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gdtr_dddxd_vvvxf_descr =
{ "gdtr_dddxd_vvvxf", (void *) gdtr_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int gdtr2_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdfgam1_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gdtr2_fffxf_vvvxf_descr =
{ "gdtr2_fffxf_vvvxf", (void *) gdtr2_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int gdtr2_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdfgam1_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gdtr2_dddxd_vvvxf_descr =
{ "gdtr2_dddxd_vvvxf", (void *) gdtr2_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int gdtrc_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = gdtrc(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gdtrc_fffxf_vvvxf_descr =
{ "gdtrc_fffxf_vvvxf", (void *) gdtrc_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int gdtrc_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = gdtrc(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gdtrc_dddxd_vvvxf_descr =
{ "gdtrc_dddxd_vvvxf", (void *) gdtrc_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int gdtri_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = gdtri(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gdtri_fffxf_vvvxf_descr =
{ "gdtri_fffxf_vvvxf", (void *) gdtri_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int gdtri_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = gdtri(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gdtri_dddxd_vvvxf_descr =
{ "gdtri_dddxd_vvvxf", (void *) gdtri_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int gdtria_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdfgam4_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gdtria_fffxf_vvvxf_descr =
{ "gdtria_fffxf_vvvxf", (void *) gdtria_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int gdtria_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdfgam4_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gdtria_dddxd_vvvxf_descr =
{ "gdtria_dddxd_vvvxf", (void *) gdtria_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int gdtrib_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdfgam3_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gdtrib_fffxf_vvvxf_descr =
{ "gdtrib_fffxf_vvvxf", (void *) gdtrib_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int gdtrib_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdfgam3_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gdtrib_dddxd_vvvxf_descr =
{ "gdtrib_dddxd_vvvxf", (void *) gdtrib_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int gdtrix_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdfgam2_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gdtrix_fffxf_vvvxf_descr =
{ "gdtrix_fffxf_vvvxf", (void *) gdtrix_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int gdtrix_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdfgam2_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor gdtrix_dddxd_vvvxf_descr =
{ "gdtrix_dddxd_vvvxf", (void *) gdtrix_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int hankel1_fFxF_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Complex32  *tin1 =  (Complex32 *) buffers[1];
    Complex32  *tout0 =  (Complex32 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = tin1->r;
        input1.i = tin1->i;

        result0 = cbesh_wrap1(*tin0, input1);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hankel1_fFxF_vvxf_descr =
{ "hankel1_fFxF_vvxf", (void *) hankel1_fFxF_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0 } };

static int hankel1_dDxD_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Complex64  *tin1 =  (Complex64 *) buffers[1];
    Complex64  *tout0 =  (Complex64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cbesh_wrap1(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hankel1_dDxD_vvxf_descr =
{ "hankel1_dDxD_vvxf", (void *) hankel1_dDxD_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0 } };

static int hankel1e_fFxF_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Complex32  *tin1 =  (Complex32 *) buffers[1];
    Complex32  *tout0 =  (Complex32 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = tin1->r;
        input1.i = tin1->i;

        result0 = cbesh_wrap1_e(*tin0, input1);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hankel1e_fFxF_vvxf_descr =
{ "hankel1e_fFxF_vvxf", (void *) hankel1e_fFxF_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0 } };

static int hankel1e_dDxD_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Complex64  *tin1 =  (Complex64 *) buffers[1];
    Complex64  *tout0 =  (Complex64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cbesh_wrap1_e(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hankel1e_dDxD_vvxf_descr =
{ "hankel1e_dDxD_vvxf", (void *) hankel1e_dDxD_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0 } };

static int hankel2_fFxF_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Complex32  *tin1 =  (Complex32 *) buffers[1];
    Complex32  *tout0 =  (Complex32 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = tin1->r;
        input1.i = tin1->i;

        result0 = cbesh_wrap2(*tin0, input1);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hankel2_fFxF_vvxf_descr =
{ "hankel2_fFxF_vvxf", (void *) hankel2_fFxF_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0 } };

static int hankel2_dDxD_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Complex64  *tin1 =  (Complex64 *) buffers[1];
    Complex64  *tout0 =  (Complex64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cbesh_wrap2(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hankel2_dDxD_vvxf_descr =
{ "hankel2_dDxD_vvxf", (void *) hankel2_dDxD_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0 } };

static int hankel2e_fFxF_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Complex32  *tin1 =  (Complex32 *) buffers[1];
    Complex32  *tout0 =  (Complex32 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = tin1->r;
        input1.i = tin1->i;

        result0 = cbesh_wrap2_e(*tin0, input1);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hankel2e_fFxF_vvxf_descr =
{ "hankel2e_fFxF_vvxf", (void *) hankel2e_fFxF_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0 } };

static int hankel2e_dDxD_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Complex64  *tin1 =  (Complex64 *) buffers[1];
    Complex64  *tout0 =  (Complex64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cbesh_wrap2_e(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hankel2e_dDxD_vvxf_descr =
{ "hankel2e_dDxD_vvxf", (void *) hankel2e_dDxD_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0 } };

static int hyp1f1_ffFxF_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Complex32  *tin2 =  (Complex32 *) buffers[2];
    Complex32  *tout0 =  (Complex32 *) buffers[3];
    Complex64     input2;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input2.r = tin2->r;
        input2.i = tin2->i;

        result0 = chyp1f1_wrap(*tin0, *tin1, input2);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hyp1f1_ffFxF_vvvxf_descr =
{ "hyp1f1_ffFxF_vvvxf", (void *) hyp1f1_ffFxF_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0, 0 } };

static int hyp1f1_ddDxD_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Complex64  *tin2 =  (Complex64 *) buffers[2];
    Complex64  *tout0 =  (Complex64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = chyp1f1_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hyp1f1_ddDxD_vvvxf_descr =
{ "hyp1f1_ddDxD_vvvxf", (void *) hyp1f1_ddDxD_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0, 0 } };

static int hyp1f1_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = hyperg(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hyp1f1_fffxf_vvvxf_descr =
{ "hyp1f1_fffxf_vvvxf", (void *) hyp1f1_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int hyp1f1_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = hyperg(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hyp1f1_dddxd_vvvxf_descr =
{ "hyp1f1_dddxd_vvvxf", (void *) hyp1f1_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int hyp1f2_ffffxff_vvvvxfv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tout0 =  (Float32 *) buffers[4];
    Float32    *tout1 =  (Float32 *) buffers[5];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = onef2(*tin0, *tin1, *tin2, *tin3, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hyp1f2_ffffxff_vvvvxfv_descr =
{ "hyp1f2_ffffxff_vvvvxfv", (void *) hyp1f2_ffffxff_vvvvxfv, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int hyp1f2_ddddxdd_vvvvxfv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tout0 =  (Float64 *) buffers[4];
    Float64    *tout1 =  (Float64 *) buffers[5];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = onef2(*tin0, *tin1, *tin2, *tin3, tout1);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hyp1f2_ddddxdd_vvvvxfv_descr =
{ "hyp1f2_ddddxdd_vvvvxfv", (void *) hyp1f2_ddddxdd_vvvvxfv, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int hyp2f0_ffffxff_vvvvxfv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tout0 =  (Float32 *) buffers[4];
    Float32    *tout1 =  (Float32 *) buffers[5];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = hyp2f0(*tin0, *tin1, *tin2, *tin3, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hyp2f0_ffffxff_vvvvxfv_descr =
{ "hyp2f0_ffffxff_vvvvxfv", (void *) hyp2f0_ffffxff_vvvvxfv, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int hyp2f0_ddddxdd_vvvvxfv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tout0 =  (Float64 *) buffers[4];
    Float64    *tout1 =  (Float64 *) buffers[5];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = hyp2f0(*tin0, *tin1, *tin2, *tin3, tout1);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hyp2f0_ddddxdd_vvvvxfv_descr =
{ "hyp2f0_ddddxdd_vvvvxfv", (void *) hyp2f0_ddddxdd_vvvvxfv, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int hyp2f1_fffFxF_vvvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Complex32  *tin3 =  (Complex32 *) buffers[3];
    Complex32  *tout0 =  (Complex32 *) buffers[4];
    Complex64     input3;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input3.r = tin3->r;
        input3.i = tin3->i;

        result0 = chyp2f1_wrap(*tin0, *tin1, *tin2, input3);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hyp2f1_fffFxF_vvvvxf_descr =
{ "hyp2f1_fffFxF_vvvvxf", (void *) hyp2f1_fffFxF_vvvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0, 0, 0 } };

static int hyp2f1_dddDxD_vvvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Complex64  *tin3 =  (Complex64 *) buffers[3];
    Complex64  *tout0 =  (Complex64 *) buffers[4];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = chyp2f1_wrap(*tin0, *tin1, *tin2, *tin3);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hyp2f1_dddDxD_vvvvxf_descr =
{ "hyp2f1_dddDxD_vvvvxf", (void *) hyp2f1_dddDxD_vvvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0, 0, 0 } };

static int hyp2f1_ffffxf_vvvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tout0 =  (Float32 *) buffers[4];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = hyp2f1(*tin0, *tin1, *tin2, *tin3);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hyp2f1_ffffxf_vvvvxf_descr =
{ "hyp2f1_ffffxf_vvvvxf", (void *) hyp2f1_ffffxf_vvvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0 } };

static int hyp2f1_ddddxd_vvvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tout0 =  (Float64 *) buffers[4];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = hyp2f1(*tin0, *tin1, *tin2, *tin3);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hyp2f1_ddddxd_vvvvxf_descr =
{ "hyp2f1_ddddxd_vvvvxf", (void *) hyp2f1_ddddxd_vvvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0 } };

static int hyp3f0_ffffxff_vvvvxfv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tout0 =  (Float32 *) buffers[4];
    Float32    *tout1 =  (Float32 *) buffers[5];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = threef0(*tin0, *tin1, *tin2, *tin3, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hyp3f0_ffffxff_vvvvxfv_descr =
{ "hyp3f0_ffffxff_vvvvxfv", (void *) hyp3f0_ffffxff_vvvvxfv, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int hyp3f0_ddddxdd_vvvvxfv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tout0 =  (Float64 *) buffers[4];
    Float64    *tout1 =  (Float64 *) buffers[5];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = threef0(*tin0, *tin1, *tin2, *tin3, tout1);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hyp3f0_ddddxdd_vvvvxfv_descr =
{ "hyp3f0_ddddxdd_vvvvxfv", (void *) hyp3f0_ddddxdd_vvvvxfv, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int hyperu_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = hypU_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hyperu_fffxf_vvvxf_descr =
{ "hyperu_fffxf_vvvxf", (void *) hyperu_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int hyperu_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = hypU_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor hyperu_dddxd_vvvxf_descr =
{ "hyperu_dddxd_vvvxf", (void *) hyperu_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int i0_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = i0(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor i0_fxf_vxf_descr =
{ "i0_fxf_vxf", (void *) i0_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int i0_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = i0(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor i0_dxd_vxf_descr =
{ "i0_dxd_vxf", (void *) i0_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int i0e_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = i0e(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor i0e_fxf_vxf_descr =
{ "i0e_fxf_vxf", (void *) i0e_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int i0e_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = i0e(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor i0e_dxd_vxf_descr =
{ "i0e_dxd_vxf", (void *) i0e_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int i1_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = i1(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor i1_fxf_vxf_descr =
{ "i1_fxf_vxf", (void *) i1_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int i1_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = i1(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor i1_dxd_vxf_descr =
{ "i1_dxd_vxf", (void *) i1_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int i1e_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = i1e(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor i1e_fxf_vxf_descr =
{ "i1e_fxf_vxf", (void *) i1e_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int i1e_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = i1e(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor i1e_dxd_vxf_descr =
{ "i1e_dxd_vxf", (void *) i1e_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int it2i0k0_fxff_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float32    *tout1 =  (Float32 *) buffers[2];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        it2i0k0_wrap(*tin0, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor it2i0k0_fxff_vxvv_descr =
{ "it2i0k0_fxff_vxvv", (void *) it2i0k0_fxff_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int it2i0k0_dxdd_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];
    Float64    *tout1 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        it2i0k0_wrap(*tin0, tout0, tout1);
        
	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor it2i0k0_dxdd_vxvv_descr =
{ "it2i0k0_dxdd_vxvv", (void *) it2i0k0_dxdd_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int it2j0y0_fxff_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float32    *tout1 =  (Float32 *) buffers[2];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        it2j0y0_wrap(*tin0, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor it2j0y0_fxff_vxvv_descr =
{ "it2j0y0_fxff_vxvv", (void *) it2j0y0_fxff_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int it2j0y0_dxdd_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];
    Float64    *tout1 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        it2j0y0_wrap(*tin0, tout0, tout1);
        
	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor it2j0y0_dxdd_vxvv_descr =
{ "it2j0y0_dxdd_vxvv", (void *) it2j0y0_dxdd_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int it2struve0_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = it2struve0_wrap(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor it2struve0_fxf_vxf_descr =
{ "it2struve0_fxf_vxf", (void *) it2struve0_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int it2struve0_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = it2struve0_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor it2struve0_dxd_vxf_descr =
{ "it2struve0_dxd_vxf", (void *) it2struve0_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int itairy_fxffff_vxvvvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float32    *tout1 =  (Float32 *) buffers[2];
    Float32    *tout2 =  (Float32 *) buffers[3];
    Float32    *tout3 =  (Float32 *) buffers[4];
    Float64     result0;
    Float64     result1;
    Float64     result2;
    Float64     result3;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        itairy_wrap(*tin0, &result0, &result1, &result2, &result3);
                *tout0 = result0;
        *tout1 = result1;
        *tout2 = result2;
        *tout3 = result3;

	++tin0; ++tout0; ++tout1; ++tout2; ++tout3; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor itairy_fxffff_vxvvvv_descr =
{ "itairy_fxffff_vxvvvv", (void *) itairy_fxffff_vxvvvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 4,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0 } };

static int itairy_dxdddd_vxvvvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];
    Float64    *tout1 =  (Float64 *) buffers[2];
    Float64    *tout2 =  (Float64 *) buffers[3];
    Float64    *tout3 =  (Float64 *) buffers[4];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        itairy_wrap(*tin0, tout0, tout1, tout2, tout3);
        
	++tin0; ++tout0; ++tout1; ++tout2; ++tout3; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor itairy_dxdddd_vxvvvv_descr =
{ "itairy_dxdddd_vxvvvv", (void *) itairy_dxdddd_vxvvvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 4,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0 } };

static int iti0k0_fxff_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float32    *tout1 =  (Float32 *) buffers[2];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        it1i0k0_wrap(*tin0, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor iti0k0_fxff_vxvv_descr =
{ "iti0k0_fxff_vxvv", (void *) iti0k0_fxff_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int iti0k0_dxdd_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];
    Float64    *tout1 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        it1i0k0_wrap(*tin0, tout0, tout1);
        
	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor iti0k0_dxdd_vxvv_descr =
{ "iti0k0_dxdd_vxvv", (void *) iti0k0_dxdd_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int itj0y0_fxff_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float32    *tout1 =  (Float32 *) buffers[2];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        it1j0y0_wrap(*tin0, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor itj0y0_fxff_vxvv_descr =
{ "itj0y0_fxff_vxvv", (void *) itj0y0_fxff_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int itj0y0_dxdd_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];
    Float64    *tout1 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        it1j0y0_wrap(*tin0, tout0, tout1);
        
	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor itj0y0_dxdd_vxvv_descr =
{ "itj0y0_dxdd_vxvv", (void *) itj0y0_dxdd_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int itmodstruve0_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = itmodstruve0_wrap(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor itmodstruve0_fxf_vxf_descr =
{ "itmodstruve0_fxf_vxf", (void *) itmodstruve0_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int itmodstruve0_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = itmodstruve0_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor itmodstruve0_dxd_vxf_descr =
{ "itmodstruve0_dxd_vxf", (void *) itmodstruve0_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int itstruve0_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = itstruve0_wrap(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor itstruve0_fxf_vxf_descr =
{ "itstruve0_fxf_vxf", (void *) itstruve0_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int itstruve0_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = itstruve0_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor itstruve0_dxd_vxf_descr =
{ "itstruve0_dxd_vxf", (void *) itstruve0_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int iv_fFxF_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Complex32  *tin1 =  (Complex32 *) buffers[1];
    Complex32  *tout0 =  (Complex32 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = tin1->r;
        input1.i = tin1->i;

        result0 = cbesi_wrap(*tin0, input1);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor iv_fFxF_vvxf_descr =
{ "iv_fFxF_vvxf", (void *) iv_fFxF_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0 } };

static int iv_dDxD_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Complex64  *tin1 =  (Complex64 *) buffers[1];
    Complex64  *tout0 =  (Complex64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cbesi_wrap(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor iv_dDxD_vvxf_descr =
{ "iv_dDxD_vvxf", (void *) iv_dDxD_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0 } };

static int iv_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = iv(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor iv_ffxf_vvxf_descr =
{ "iv_ffxf_vvxf", (void *) iv_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int iv_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = iv(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor iv_ddxd_vvxf_descr =
{ "iv_ddxd_vvxf", (void *) iv_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int ive_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = *tin1;
        input1.i = 0;

        result0 = cbesi_wrap_e(*tin0, input1);
                *tout0 = result0.r;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ive_ffxf_vvxf_descr =
{ "ive_ffxf_vvxf", (void *) ive_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int ive_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = *tin1;
        input1.i = 0;

        result0 = cbesi_wrap_e(*tin0, input1);
                *tout0 = result0.r;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ive_ddxd_vvxf_descr =
{ "ive_ddxd_vvxf", (void *) ive_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int ive_fFxF_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Complex32  *tin1 =  (Complex32 *) buffers[1];
    Complex32  *tout0 =  (Complex32 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = tin1->r;
        input1.i = tin1->i;

        result0 = cbesi_wrap_e(*tin0, input1);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ive_fFxF_vvxf_descr =
{ "ive_fFxF_vvxf", (void *) ive_fFxF_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0 } };

static int ive_dDxD_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Complex64  *tin1 =  (Complex64 *) buffers[1];
    Complex64  *tout0 =  (Complex64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cbesi_wrap_e(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ive_dDxD_vvxf_descr =
{ "ive_dDxD_vvxf", (void *) ive_dDxD_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0 } };

static int j0_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = j0(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor j0_fxf_vxf_descr =
{ "j0_fxf_vxf", (void *) j0_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int j0_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = j0(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor j0_dxd_vxf_descr =
{ "j0_dxd_vxf", (void *) j0_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int j1_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = j1(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor j1_fxf_vxf_descr =
{ "j1_fxf_vxf", (void *) j1_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int j1_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = j1(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor j1_dxd_vxf_descr =
{ "j1_dxd_vxf", (void *) j1_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int jn_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = jn(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor jn_ffxf_vvxf_descr =
{ "jn_ffxf_vvxf", (void *) jn_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int jn_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = jn(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor jn_ddxd_vvxf_descr =
{ "jn_ddxd_vvxf", (void *) jn_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int jv_fFxF_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Complex32  *tin1 =  (Complex32 *) buffers[1];
    Complex32  *tout0 =  (Complex32 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = tin1->r;
        input1.i = tin1->i;

        result0 = cbesj_wrap(*tin0, input1);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor jv_fFxF_vvxf_descr =
{ "jv_fFxF_vvxf", (void *) jv_fFxF_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0 } };

static int jv_dDxD_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Complex64  *tin1 =  (Complex64 *) buffers[1];
    Complex64  *tout0 =  (Complex64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cbesj_wrap(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor jv_dDxD_vvxf_descr =
{ "jv_dDxD_vvxf", (void *) jv_dDxD_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0 } };

static int jv_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = jv(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor jv_ffxf_vvxf_descr =
{ "jv_ffxf_vvxf", (void *) jv_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int jv_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = jv(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor jv_ddxd_vvxf_descr =
{ "jv_ddxd_vvxf", (void *) jv_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int jve_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = *tin1;
        input1.i = 0;

        result0 = cbesj_wrap_e(*tin0, input1);
                *tout0 = result0.r;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor jve_ffxf_vvxf_descr =
{ "jve_ffxf_vvxf", (void *) jve_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int jve_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = *tin1;
        input1.i = 0;

        result0 = cbesj_wrap_e(*tin0, input1);
                *tout0 = result0.r;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor jve_ddxd_vvxf_descr =
{ "jve_ddxd_vvxf", (void *) jve_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int jve_fFxF_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Complex32  *tin1 =  (Complex32 *) buffers[1];
    Complex32  *tout0 =  (Complex32 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = tin1->r;
        input1.i = tin1->i;

        result0 = cbesj_wrap_e(*tin0, input1);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor jve_fFxF_vvxf_descr =
{ "jve_fFxF_vvxf", (void *) jve_fFxF_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0 } };

static int jve_dDxD_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Complex64  *tin1 =  (Complex64 *) buffers[1];
    Complex64  *tout0 =  (Complex64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cbesj_wrap_e(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor jve_dDxD_vvxf_descr =
{ "jve_dDxD_vvxf", (void *) jve_dDxD_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0 } };

static int k0_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = k0(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor k0_fxf_vxf_descr =
{ "k0_fxf_vxf", (void *) k0_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int k0_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = k0(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor k0_dxd_vxf_descr =
{ "k0_dxd_vxf", (void *) k0_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int k0e_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = k0e(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor k0e_fxf_vxf_descr =
{ "k0e_fxf_vxf", (void *) k0e_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int k0e_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = k0e(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor k0e_dxd_vxf_descr =
{ "k0e_dxd_vxf", (void *) k0e_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int k1_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = k1(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor k1_fxf_vxf_descr =
{ "k1_fxf_vxf", (void *) k1_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int k1_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = k1(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor k1_dxd_vxf_descr =
{ "k1_dxd_vxf", (void *) k1_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int k1e_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = k1e(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor k1e_fxf_vxf_descr =
{ "k1e_fxf_vxf", (void *) k1e_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int k1e_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = k1e(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor k1e_dxd_vxf_descr =
{ "k1e_dxd_vxf", (void *) k1e_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int kei_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = kei_wrap(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kei_fxf_vxf_descr =
{ "kei_fxf_vxf", (void *) kei_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int kei_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = kei_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kei_dxd_vxf_descr =
{ "kei_dxd_vxf", (void *) kei_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int keip_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = keip_wrap(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor keip_fxf_vxf_descr =
{ "keip_fxf_vxf", (void *) keip_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int keip_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = keip_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor keip_dxd_vxf_descr =
{ "keip_dxd_vxf", (void *) keip_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int kelvin_fxFFFF_vxvvvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Complex32  *tout0 =  (Complex32 *) buffers[1];
    Complex32  *tout1 =  (Complex32 *) buffers[2];
    Complex32  *tout2 =  (Complex32 *) buffers[3];
    Complex32  *tout3 =  (Complex32 *) buffers[4];
    Complex64     result0;
    Complex64     result1;
    Complex64     result2;
    Complex64     result3;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        kelvin_wrap(*tin0, &result0, &result1, &result2, &result3);
                tout0->r = result0.r;
        tout0->i = result0.i;
        tout1->r = result1.r;
        tout1->i = result1.i;
        tout2->r = result2.r;
        tout2->i = result2.i;
        tout3->r = result3.r;
        tout3->i = result3.i;

	++tin0; ++tout0; ++tout1; ++tout2; ++tout3; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kelvin_fxFFFF_vxvvvv_descr =
{ "kelvin_fxFFFF_vxvvvv", (void *) kelvin_fxFFFF_vxvvvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 4,
  { sizeof(Float32), sizeof(Complex32), sizeof(Complex32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0, 0, 0 } };

static int kelvin_dxDDDD_vxvvvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Complex64  *tout0 =  (Complex64 *) buffers[1];
    Complex64  *tout1 =  (Complex64 *) buffers[2];
    Complex64  *tout2 =  (Complex64 *) buffers[3];
    Complex64  *tout3 =  (Complex64 *) buffers[4];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        kelvin_wrap(*tin0, tout0, tout1, tout2, tout3);
        
	++tin0; ++tout0; ++tout1; ++tout2; ++tout3; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kelvin_dxDDDD_vxvvvv_descr =
{ "kelvin_dxDDDD_vxvvvv", (void *) kelvin_dxDDDD_vxvvvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 4,
  { sizeof(Float64), sizeof(Complex64), sizeof(Complex64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0, 0, 0 } };

static int ker_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = ker_wrap(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ker_fxf_vxf_descr =
{ "ker_fxf_vxf", (void *) ker_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int ker_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = ker_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ker_dxd_vxf_descr =
{ "ker_dxd_vxf", (void *) ker_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int kerp_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = kerp_wrap(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kerp_fxf_vxf_descr =
{ "kerp_fxf_vxf", (void *) kerp_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int kerp_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = kerp_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kerp_dxd_vxf_descr =
{ "kerp_dxd_vxf", (void *) kerp_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int kn_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = kn(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kn_ffxf_vvxf_descr =
{ "kn_ffxf_vvxf", (void *) kn_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int kn_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = kn(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kn_ddxd_vvxf_descr =
{ "kn_ddxd_vvxf", (void *) kn_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int kolmogi_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = kolmogi(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kolmogi_fxf_vxf_descr =
{ "kolmogi_fxf_vxf", (void *) kolmogi_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int kolmogi_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = kolmogi(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kolmogi_dxd_vxf_descr =
{ "kolmogi_dxd_vxf", (void *) kolmogi_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int kolmogorov_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = kolmogorov(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kolmogorov_fxf_vxf_descr =
{ "kolmogorov_fxf_vxf", (void *) kolmogorov_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int kolmogorov_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = kolmogorov(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kolmogorov_dxd_vxf_descr =
{ "kolmogorov_dxd_vxf", (void *) kolmogorov_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int kv_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = *tin1;
        input1.i = 0;

        result0 = cbesk_wrap(*tin0, input1);
                *tout0 = result0.r;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kv_ffxf_vvxf_descr =
{ "kv_ffxf_vvxf", (void *) kv_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int kv_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = *tin1;
        input1.i = 0;

        result0 = cbesk_wrap(*tin0, input1);
                *tout0 = result0.r;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kv_ddxd_vvxf_descr =
{ "kv_ddxd_vvxf", (void *) kv_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int kv_fFxF_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Complex32  *tin1 =  (Complex32 *) buffers[1];
    Complex32  *tout0 =  (Complex32 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = tin1->r;
        input1.i = tin1->i;

        result0 = cbesk_wrap(*tin0, input1);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kv_fFxF_vvxf_descr =
{ "kv_fFxF_vvxf", (void *) kv_fFxF_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0 } };

static int kv_dDxD_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Complex64  *tin1 =  (Complex64 *) buffers[1];
    Complex64  *tout0 =  (Complex64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cbesk_wrap(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kv_dDxD_vvxf_descr =
{ "kv_dDxD_vvxf", (void *) kv_dDxD_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0 } };

static int kve_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = *tin1;
        input1.i = 0;

        result0 = cbesk_wrap_e(*tin0, input1);
                *tout0 = result0.r;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kve_ffxf_vvxf_descr =
{ "kve_ffxf_vvxf", (void *) kve_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int kve_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = *tin1;
        input1.i = 0;

        result0 = cbesk_wrap_e(*tin0, input1);
                *tout0 = result0.r;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kve_ddxd_vvxf_descr =
{ "kve_ddxd_vvxf", (void *) kve_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int kve_fFxF_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Complex32  *tin1 =  (Complex32 *) buffers[1];
    Complex32  *tout0 =  (Complex32 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = tin1->r;
        input1.i = tin1->i;

        result0 = cbesk_wrap_e(*tin0, input1);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kve_fFxF_vvxf_descr =
{ "kve_fFxF_vvxf", (void *) kve_fFxF_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0 } };

static int kve_dDxD_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Complex64  *tin1 =  (Complex64 *) buffers[1];
    Complex64  *tout0 =  (Complex64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cbesk_wrap_e(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor kve_dDxD_vvxf_descr =
{ "kve_dDxD_vvxf", (void *) kve_dDxD_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0 } };

static int log1p_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = log1p(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor log1p_fxf_vxf_descr =
{ "log1p_fxf_vxf", (void *) log1p_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int log1p_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = log1p(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor log1p_dxd_vxf_descr =
{ "log1p_dxd_vxf", (void *) log1p_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int lpmv_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = pmv_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor lpmv_fffxf_vvvxf_descr =
{ "lpmv_fffxf_vvvxf", (void *) lpmv_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int lpmv_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = pmv_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor lpmv_dddxd_vvvxf_descr =
{ "lpmv_dddxd_vvvxf", (void *) lpmv_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int mathieu_a_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cem_cva_wrap(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor mathieu_a_ffxf_vvxf_descr =
{ "mathieu_a_ffxf_vvxf", (void *) mathieu_a_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int mathieu_a_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cem_cva_wrap(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor mathieu_a_ddxd_vvxf_descr =
{ "mathieu_a_ddxd_vvxf", (void *) mathieu_a_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int mathieu_b_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = sem_cva_wrap(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor mathieu_b_ffxf_vvxf_descr =
{ "mathieu_b_ffxf_vvxf", (void *) mathieu_b_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int mathieu_b_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = sem_cva_wrap(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor mathieu_b_ddxd_vvxf_descr =
{ "mathieu_b_ddxd_vvxf", (void *) mathieu_b_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int mathieu_cem_fffxff_vvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float32    *tout1 =  (Float32 *) buffers[4];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        cem_wrap(*tin0, *tin1, *tin2, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor mathieu_cem_fffxff_vvvxvv_descr =
{ "mathieu_cem_fffxff_vvvxvv", (void *) mathieu_cem_fffxff_vvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0 } };

static int mathieu_cem_dddxdd_vvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];
    Float64    *tout1 =  (Float64 *) buffers[4];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        cem_wrap(*tin0, *tin1, *tin2, tout0, tout1);
        
	++tin0; ++tin1; ++tin2; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor mathieu_cem_dddxdd_vvvxvv_descr =
{ "mathieu_cem_dddxdd_vvvxvv", (void *) mathieu_cem_dddxdd_vvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0 } };

static int mathieu_modcem1_fffxff_vvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float32    *tout1 =  (Float32 *) buffers[4];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        mcm1_wrap(*tin0, *tin1, *tin2, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor mathieu_modcem1_fffxff_vvvxvv_descr =
{ "mathieu_modcem1_fffxff_vvvxvv", (void *) mathieu_modcem1_fffxff_vvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0 } };

static int mathieu_modcem1_dddxdd_vvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];
    Float64    *tout1 =  (Float64 *) buffers[4];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        mcm1_wrap(*tin0, *tin1, *tin2, tout0, tout1);
        
	++tin0; ++tin1; ++tin2; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor mathieu_modcem1_dddxdd_vvvxvv_descr =
{ "mathieu_modcem1_dddxdd_vvvxvv", (void *) mathieu_modcem1_dddxdd_vvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0 } };

static int mathieu_modcem2_fffxff_vvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float32    *tout1 =  (Float32 *) buffers[4];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        mcm2_wrap(*tin0, *tin1, *tin2, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor mathieu_modcem2_fffxff_vvvxvv_descr =
{ "mathieu_modcem2_fffxff_vvvxvv", (void *) mathieu_modcem2_fffxff_vvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0 } };

static int mathieu_modcem2_dddxdd_vvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];
    Float64    *tout1 =  (Float64 *) buffers[4];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        mcm2_wrap(*tin0, *tin1, *tin2, tout0, tout1);
        
	++tin0; ++tin1; ++tin2; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor mathieu_modcem2_dddxdd_vvvxvv_descr =
{ "mathieu_modcem2_dddxdd_vvvxvv", (void *) mathieu_modcem2_dddxdd_vvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0 } };

static int mathieu_modsem1_fffxff_vvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float32    *tout1 =  (Float32 *) buffers[4];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        msm1_wrap(*tin0, *tin1, *tin2, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor mathieu_modsem1_fffxff_vvvxvv_descr =
{ "mathieu_modsem1_fffxff_vvvxvv", (void *) mathieu_modsem1_fffxff_vvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0 } };

static int mathieu_modsem1_dddxdd_vvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];
    Float64    *tout1 =  (Float64 *) buffers[4];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        msm1_wrap(*tin0, *tin1, *tin2, tout0, tout1);
        
	++tin0; ++tin1; ++tin2; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor mathieu_modsem1_dddxdd_vvvxvv_descr =
{ "mathieu_modsem1_dddxdd_vvvxvv", (void *) mathieu_modsem1_dddxdd_vvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0 } };

static int mathieu_modsem2_fffxff_vvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float32    *tout1 =  (Float32 *) buffers[4];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        msm2_wrap(*tin0, *tin1, *tin2, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor mathieu_modsem2_fffxff_vvvxvv_descr =
{ "mathieu_modsem2_fffxff_vvvxvv", (void *) mathieu_modsem2_fffxff_vvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0 } };

static int mathieu_modsem2_dddxdd_vvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];
    Float64    *tout1 =  (Float64 *) buffers[4];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        msm2_wrap(*tin0, *tin1, *tin2, tout0, tout1);
        
	++tin0; ++tin1; ++tin2; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor mathieu_modsem2_dddxdd_vvvxvv_descr =
{ "mathieu_modsem2_dddxdd_vvvxvv", (void *) mathieu_modsem2_dddxdd_vvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0 } };

static int mathieu_sem_fffxff_vvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float32    *tout1 =  (Float32 *) buffers[4];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        sem_wrap(*tin0, *tin1, *tin2, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor mathieu_sem_fffxff_vvvxvv_descr =
{ "mathieu_sem_fffxff_vvvxvv", (void *) mathieu_sem_fffxff_vvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0 } };

static int mathieu_sem_dddxdd_vvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];
    Float64    *tout1 =  (Float64 *) buffers[4];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        sem_wrap(*tin0, *tin1, *tin2, tout0, tout1);
        
	++tin0; ++tin1; ++tin2; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor mathieu_sem_dddxdd_vvvxvv_descr =
{ "mathieu_sem_dddxdd_vvvxvv", (void *) mathieu_sem_dddxdd_vvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0 } };

static int modfresnelm_fxFF_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Complex32  *tout0 =  (Complex32 *) buffers[1];
    Complex32  *tout1 =  (Complex32 *) buffers[2];
    Complex64     result0;
    Complex64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        modified_fresnel_minus_wrap(*tin0, &result0, &result1);
                tout0->r = result0.r;
        tout0->i = result0.i;
        tout1->r = result1.r;
        tout1->i = result1.i;

	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor modfresnelm_fxFF_vxvv_descr =
{ "modfresnelm_fxFF_vxvv", (void *) modfresnelm_fxFF_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Float32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0 } };

static int modfresnelm_dxDD_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Complex64  *tout0 =  (Complex64 *) buffers[1];
    Complex64  *tout1 =  (Complex64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        modified_fresnel_minus_wrap(*tin0, tout0, tout1);
        
	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor modfresnelm_dxDD_vxvv_descr =
{ "modfresnelm_dxDD_vxvv", (void *) modfresnelm_dxDD_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Float64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0 } };

static int modfresnelp_fxFF_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Complex32  *tout0 =  (Complex32 *) buffers[1];
    Complex32  *tout1 =  (Complex32 *) buffers[2];
    Complex64     result0;
    Complex64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        modified_fresnel_plus_wrap(*tin0, &result0, &result1);
                tout0->r = result0.r;
        tout0->i = result0.i;
        tout1->r = result1.r;
        tout1->i = result1.i;

	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor modfresnelp_fxFF_vxvv_descr =
{ "modfresnelp_fxFF_vxvv", (void *) modfresnelp_fxFF_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Float32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0 } };

static int modfresnelp_dxDD_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Complex64  *tout0 =  (Complex64 *) buffers[1];
    Complex64  *tout1 =  (Complex64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        modified_fresnel_plus_wrap(*tin0, tout0, tout1);
        
	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor modfresnelp_dxDD_vxvv_descr =
{ "modfresnelp_dxDD_vxvv", (void *) modfresnelp_dxDD_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Float64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0 } };

static int modstruve_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = modstruve_wrap(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor modstruve_ffxf_vvxf_descr =
{ "modstruve_ffxf_vvxf", (void *) modstruve_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int modstruve_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = modstruve_wrap(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor modstruve_ddxd_vvxf_descr =
{ "modstruve_ddxd_vvxf", (void *) modstruve_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int nbdtr_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = nbdtr(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nbdtr_fffxf_vvvxf_descr =
{ "nbdtr_fffxf_vvvxf", (void *) nbdtr_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int nbdtr_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = nbdtr(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nbdtr_dddxd_vvvxf_descr =
{ "nbdtr_dddxd_vvvxf", (void *) nbdtr_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int nbdtrc_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = nbdtrc(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nbdtrc_fffxf_vvvxf_descr =
{ "nbdtrc_fffxf_vvvxf", (void *) nbdtrc_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int nbdtrc_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = nbdtrc(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nbdtrc_dddxd_vvvxf_descr =
{ "nbdtrc_dddxd_vvvxf", (void *) nbdtrc_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int nbdtri_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = nbdtri(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nbdtri_fffxf_vvvxf_descr =
{ "nbdtri_fffxf_vvvxf", (void *) nbdtri_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int nbdtri_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = nbdtri(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nbdtri_dddxd_vvvxf_descr =
{ "nbdtri_dddxd_vvvxf", (void *) nbdtri_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int nbdtrik_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdfnbn2_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nbdtrik_fffxf_vvvxf_descr =
{ "nbdtrik_fffxf_vvvxf", (void *) nbdtrik_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int nbdtrik_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdfnbn2_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nbdtrik_dddxd_vvvxf_descr =
{ "nbdtrik_dddxd_vvvxf", (void *) nbdtrik_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int nbdtrin_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdfnbn3_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nbdtrin_fffxf_vvvxf_descr =
{ "nbdtrin_fffxf_vvvxf", (void *) nbdtrin_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int nbdtrin_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdfnbn3_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nbdtrin_dddxd_vvvxf_descr =
{ "nbdtrin_dddxd_vvvxf", (void *) nbdtrin_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int ncfdtr_ffffxf_vvvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tout0 =  (Float32 *) buffers[4];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdffnc1_wrap(*tin0, *tin1, *tin2, *tin3);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ncfdtr_ffffxf_vvvvxf_descr =
{ "ncfdtr_ffffxf_vvvvxf", (void *) ncfdtr_ffffxf_vvvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0 } };

static int ncfdtr_ddddxd_vvvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tout0 =  (Float64 *) buffers[4];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdffnc1_wrap(*tin0, *tin1, *tin2, *tin3);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ncfdtr_ddddxd_vvvvxf_descr =
{ "ncfdtr_ddddxd_vvvvxf", (void *) ncfdtr_ddddxd_vvvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0 } };

static int ncfdtri_ffffxf_vvvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tout0 =  (Float32 *) buffers[4];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdffnc2_wrap(*tin0, *tin1, *tin2, *tin3);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ncfdtri_ffffxf_vvvvxf_descr =
{ "ncfdtri_ffffxf_vvvvxf", (void *) ncfdtri_ffffxf_vvvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0 } };

static int ncfdtri_ddddxd_vvvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tout0 =  (Float64 *) buffers[4];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdffnc2_wrap(*tin0, *tin1, *tin2, *tin3);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ncfdtri_ddddxd_vvvvxf_descr =
{ "ncfdtri_ddddxd_vvvvxf", (void *) ncfdtri_ddddxd_vvvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0 } };

static int ncfdtridfd_ffffxf_vvvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tout0 =  (Float32 *) buffers[4];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdffnc4_wrap(*tin0, *tin1, *tin2, *tin3);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ncfdtridfd_ffffxf_vvvvxf_descr =
{ "ncfdtridfd_ffffxf_vvvvxf", (void *) ncfdtridfd_ffffxf_vvvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0 } };

static int ncfdtridfd_ddddxd_vvvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tout0 =  (Float64 *) buffers[4];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdffnc4_wrap(*tin0, *tin1, *tin2, *tin3);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ncfdtridfd_ddddxd_vvvvxf_descr =
{ "ncfdtridfd_ddddxd_vvvvxf", (void *) ncfdtridfd_ddddxd_vvvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0 } };

static int ncfdtridfn_ffffxf_vvvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tout0 =  (Float32 *) buffers[4];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdffnc3_wrap(*tin0, *tin1, *tin2, *tin3);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ncfdtridfn_ffffxf_vvvvxf_descr =
{ "ncfdtridfn_ffffxf_vvvvxf", (void *) ncfdtridfn_ffffxf_vvvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0 } };

static int ncfdtridfn_ddddxd_vvvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tout0 =  (Float64 *) buffers[4];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdffnc3_wrap(*tin0, *tin1, *tin2, *tin3);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ncfdtridfn_ddddxd_vvvvxf_descr =
{ "ncfdtridfn_ddddxd_vvvvxf", (void *) ncfdtridfn_ddddxd_vvvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0 } };

static int ncfdtrinc_ffffxf_vvvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tout0 =  (Float32 *) buffers[4];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdffnc5_wrap(*tin0, *tin1, *tin2, *tin3);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ncfdtrinc_ffffxf_vvvvxf_descr =
{ "ncfdtrinc_ffffxf_vvvvxf", (void *) ncfdtrinc_ffffxf_vvvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0 } };

static int ncfdtrinc_ddddxd_vvvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tout0 =  (Float64 *) buffers[4];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdffnc5_wrap(*tin0, *tin1, *tin2, *tin3);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ncfdtrinc_ddddxd_vvvvxf_descr =
{ "ncfdtrinc_ddddxd_vvvvxf", (void *) ncfdtrinc_ddddxd_vvvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0 } };

static int nctdtr_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdftnc1_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nctdtr_fffxf_vvvxf_descr =
{ "nctdtr_fffxf_vvvxf", (void *) nctdtr_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int nctdtr_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdftnc1_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nctdtr_dddxd_vvvxf_descr =
{ "nctdtr_dddxd_vvvxf", (void *) nctdtr_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int nctdtridf_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdftnc3_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nctdtridf_fffxf_vvvxf_descr =
{ "nctdtridf_fffxf_vvvxf", (void *) nctdtridf_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int nctdtridf_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdftnc3_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nctdtridf_dddxd_vvvxf_descr =
{ "nctdtridf_dddxd_vvvxf", (void *) nctdtridf_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int nctdtrinc_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdftnc4_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nctdtrinc_fffxf_vvvxf_descr =
{ "nctdtrinc_fffxf_vvvxf", (void *) nctdtrinc_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int nctdtrinc_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdftnc4_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nctdtrinc_dddxd_vvvxf_descr =
{ "nctdtrinc_dddxd_vvvxf", (void *) nctdtrinc_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int nctdtrit_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdftnc2_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nctdtrit_fffxf_vvvxf_descr =
{ "nctdtrit_fffxf_vvvxf", (void *) nctdtrit_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int nctdtrit_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdftnc2_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nctdtrit_dddxd_vvvxf_descr =
{ "nctdtrit_dddxd_vvvxf", (void *) nctdtrit_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int ndtr_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = ndtr(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ndtr_fxf_vxf_descr =
{ "ndtr_fxf_vxf", (void *) ndtr_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int ndtr_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = ndtr(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ndtr_dxd_vxf_descr =
{ "ndtr_dxd_vxf", (void *) ndtr_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int ndtri_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = ndtri(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ndtri_fxf_vxf_descr =
{ "ndtri_fxf_vxf", (void *) ndtri_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int ndtri_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = ndtri(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor ndtri_dxd_vxf_descr =
{ "ndtri_dxd_vxf", (void *) ndtri_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int nrdtrimn_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdfnor3_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nrdtrimn_fffxf_vvvxf_descr =
{ "nrdtrimn_fffxf_vvvxf", (void *) nrdtrimn_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int nrdtrimn_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdfnor3_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nrdtrimn_dddxd_vvvxf_descr =
{ "nrdtrimn_dddxd_vvvxf", (void *) nrdtrimn_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int nrdtrisd_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdfnor4_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nrdtrisd_fffxf_vvvxf_descr =
{ "nrdtrisd_fffxf_vvvxf", (void *) nrdtrisd_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int nrdtrisd_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdfnor4_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor nrdtrisd_dddxd_vvvxf_descr =
{ "nrdtrisd_dddxd_vvvxf", (void *) nrdtrisd_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int obl_ang1_ffffxff_vvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tout0 =  (Float32 *) buffers[4];
    Float32    *tout1 =  (Float32 *) buffers[5];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        oblate_aswfa_nocv_wrap(*tin0, *tin1, *tin2, *tin3, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor obl_ang1_ffffxff_vvvvxvv_descr =
{ "obl_ang1_ffffxff_vvvvxvv", (void *) obl_ang1_ffffxff_vvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int obl_ang1_ddddxdd_vvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tout0 =  (Float64 *) buffers[4];
    Float64    *tout1 =  (Float64 *) buffers[5];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        oblate_aswfa_nocv_wrap(*tin0, *tin1, *tin2, *tin3, tout0, tout1);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor obl_ang1_ddddxdd_vvvvxvv_descr =
{ "obl_ang1_ddddxdd_vvvvxvv", (void *) obl_ang1_ddddxdd_vvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int obl_ang1_cv_fffffxff_vvvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tin4 =  (Float32 *) buffers[4];
    Float32    *tout0 =  (Float32 *) buffers[5];
    Float32    *tout1 =  (Float32 *) buffers[6];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        oblate_aswfa_wrap(*tin0, *tin1, *tin2, *tin3, *tin4, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tin3; ++tin4; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor obl_ang1_cv_fffffxff_vvvvvxvv_descr =
{ "obl_ang1_cv_fffffxff_vvvvvxvv", (void *) obl_ang1_cv_fffffxff_vvvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 5, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0, 0, 0 } };

static int obl_ang1_cv_dddddxdd_vvvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tin4 =  (Float64 *) buffers[4];
    Float64    *tout0 =  (Float64 *) buffers[5];
    Float64    *tout1 =  (Float64 *) buffers[6];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        oblate_aswfa_wrap(*tin0, *tin1, *tin2, *tin3, *tin4, tout0, tout1);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tin4; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor obl_ang1_cv_dddddxdd_vvvvvxvv_descr =
{ "obl_ang1_cv_dddddxdd_vvvvvxvv", (void *) obl_ang1_cv_dddddxdd_vvvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 5, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0, 0, 0 } };

static int obl_cv_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = oblate_segv_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor obl_cv_fffxf_vvvxf_descr =
{ "obl_cv_fffxf_vvvxf", (void *) obl_cv_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int obl_cv_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = oblate_segv_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor obl_cv_dddxd_vvvxf_descr =
{ "obl_cv_dddxd_vvvxf", (void *) obl_cv_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int obl_rad1_ffffxff_vvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tout0 =  (Float32 *) buffers[4];
    Float32    *tout1 =  (Float32 *) buffers[5];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        oblate_radial1_nocv_wrap(*tin0, *tin1, *tin2, *tin3, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor obl_rad1_ffffxff_vvvvxvv_descr =
{ "obl_rad1_ffffxff_vvvvxvv", (void *) obl_rad1_ffffxff_vvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int obl_rad1_ddddxdd_vvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tout0 =  (Float64 *) buffers[4];
    Float64    *tout1 =  (Float64 *) buffers[5];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        oblate_radial1_nocv_wrap(*tin0, *tin1, *tin2, *tin3, tout0, tout1);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor obl_rad1_ddddxdd_vvvvxvv_descr =
{ "obl_rad1_ddddxdd_vvvvxvv", (void *) obl_rad1_ddddxdd_vvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int obl_rad1_cv_fffffxff_vvvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tin4 =  (Float32 *) buffers[4];
    Float32    *tout0 =  (Float32 *) buffers[5];
    Float32    *tout1 =  (Float32 *) buffers[6];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        oblate_radial1_wrap(*tin0, *tin1, *tin2, *tin3, *tin4, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tin3; ++tin4; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor obl_rad1_cv_fffffxff_vvvvvxvv_descr =
{ "obl_rad1_cv_fffffxff_vvvvvxvv", (void *) obl_rad1_cv_fffffxff_vvvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 5, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0, 0, 0 } };

static int obl_rad1_cv_dddddxdd_vvvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tin4 =  (Float64 *) buffers[4];
    Float64    *tout0 =  (Float64 *) buffers[5];
    Float64    *tout1 =  (Float64 *) buffers[6];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        oblate_radial1_wrap(*tin0, *tin1, *tin2, *tin3, *tin4, tout0, tout1);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tin4; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor obl_rad1_cv_dddddxdd_vvvvvxvv_descr =
{ "obl_rad1_cv_dddddxdd_vvvvvxvv", (void *) obl_rad1_cv_dddddxdd_vvvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 5, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0, 0, 0 } };

static int obl_rad2_ffffxff_vvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tout0 =  (Float32 *) buffers[4];
    Float32    *tout1 =  (Float32 *) buffers[5];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        oblate_radial2_nocv_wrap(*tin0, *tin1, *tin2, *tin3, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor obl_rad2_ffffxff_vvvvxvv_descr =
{ "obl_rad2_ffffxff_vvvvxvv", (void *) obl_rad2_ffffxff_vvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int obl_rad2_ddddxdd_vvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tout0 =  (Float64 *) buffers[4];
    Float64    *tout1 =  (Float64 *) buffers[5];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        oblate_radial2_nocv_wrap(*tin0, *tin1, *tin2, *tin3, tout0, tout1);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor obl_rad2_ddddxdd_vvvvxvv_descr =
{ "obl_rad2_ddddxdd_vvvvxvv", (void *) obl_rad2_ddddxdd_vvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int obl_rad2_cv_fffffxff_vvvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tin4 =  (Float32 *) buffers[4];
    Float32    *tout0 =  (Float32 *) buffers[5];
    Float32    *tout1 =  (Float32 *) buffers[6];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        oblate_radial2_wrap(*tin0, *tin1, *tin2, *tin3, *tin4, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tin3; ++tin4; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor obl_rad2_cv_fffffxff_vvvvvxvv_descr =
{ "obl_rad2_cv_fffffxff_vvvvvxvv", (void *) obl_rad2_cv_fffffxff_vvvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 5, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0, 0, 0 } };

static int obl_rad2_cv_dddddxdd_vvvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tin4 =  (Float64 *) buffers[4];
    Float64    *tout0 =  (Float64 *) buffers[5];
    Float64    *tout1 =  (Float64 *) buffers[6];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        oblate_radial2_wrap(*tin0, *tin1, *tin2, *tin3, *tin4, tout0, tout1);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tin4; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor obl_rad2_cv_dddddxdd_vvvvvxvv_descr =
{ "obl_rad2_cv_dddddxdd_vvvvvxvv", (void *) obl_rad2_cv_dddddxdd_vvvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 5, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0, 0, 0 } };

static int pbdv_ffxff_vvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float32    *tout1 =  (Float32 *) buffers[3];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        pbdv_wrap(*tin0, *tin1, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pbdv_ffxff_vvxvv_descr =
{ "pbdv_ffxff_vvxvv", (void *) pbdv_ffxff_vvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int pbdv_ddxdd_vvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];
    Float64    *tout1 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        pbdv_wrap(*tin0, *tin1, tout0, tout1);
        
	++tin0; ++tin1; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pbdv_ddxdd_vvxvv_descr =
{ "pbdv_ddxdd_vvxvv", (void *) pbdv_ddxdd_vvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int pbvv_ffxff_vvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float32    *tout1 =  (Float32 *) buffers[3];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        pbvv_wrap(*tin0, *tin1, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pbvv_ffxff_vvxvv_descr =
{ "pbvv_ffxff_vvxvv", (void *) pbvv_ffxff_vvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int pbvv_ddxdd_vvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];
    Float64    *tout1 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        pbvv_wrap(*tin0, *tin1, tout0, tout1);
        
	++tin0; ++tin1; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pbvv_ddxdd_vvxvv_descr =
{ "pbvv_ddxdd_vvxvv", (void *) pbvv_ddxdd_vvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int pbwa_ffxff_vvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float32    *tout1 =  (Float32 *) buffers[3];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        pbwa_wrap(*tin0, *tin1, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pbwa_ffxff_vvxvv_descr =
{ "pbwa_ffxff_vvxvv", (void *) pbwa_ffxff_vvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int pbwa_ddxdd_vvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];
    Float64    *tout1 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        pbwa_wrap(*tin0, *tin1, tout0, tout1);
        
	++tin0; ++tin1; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pbwa_ddxdd_vvxvv_descr =
{ "pbwa_ddxdd_vvxvv", (void *) pbwa_ddxdd_vvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int pdtr_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = pdtr(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pdtr_ffxf_vvxf_descr =
{ "pdtr_ffxf_vvxf", (void *) pdtr_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int pdtr_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = pdtr(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pdtr_ddxd_vvxf_descr =
{ "pdtr_ddxd_vvxf", (void *) pdtr_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int pdtrc_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = pdtrc(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pdtrc_ffxf_vvxf_descr =
{ "pdtrc_ffxf_vvxf", (void *) pdtrc_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int pdtrc_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = pdtrc(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pdtrc_ddxd_vvxf_descr =
{ "pdtrc_ddxd_vvxf", (void *) pdtrc_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int pdtri_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = pdtri(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pdtri_ffxf_vvxf_descr =
{ "pdtri_ffxf_vvxf", (void *) pdtri_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int pdtri_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = pdtri(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pdtri_ddxd_vvxf_descr =
{ "pdtri_ddxd_vvxf", (void *) pdtri_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int pdtrik_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdfpoi2_wrap(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pdtrik_ffxf_vvxf_descr =
{ "pdtrik_ffxf_vvxf", (void *) pdtrik_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int pdtrik_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdfpoi2_wrap(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pdtrik_ddxd_vvxf_descr =
{ "pdtrik_ddxd_vvxf", (void *) pdtrik_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int pro_ang1_ffffxff_vvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tout0 =  (Float32 *) buffers[4];
    Float32    *tout1 =  (Float32 *) buffers[5];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        prolate_aswfa_nocv_wrap(*tin0, *tin1, *tin2, *tin3, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pro_ang1_ffffxff_vvvvxvv_descr =
{ "pro_ang1_ffffxff_vvvvxvv", (void *) pro_ang1_ffffxff_vvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int pro_ang1_ddddxdd_vvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tout0 =  (Float64 *) buffers[4];
    Float64    *tout1 =  (Float64 *) buffers[5];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        prolate_aswfa_nocv_wrap(*tin0, *tin1, *tin2, *tin3, tout0, tout1);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pro_ang1_ddddxdd_vvvvxvv_descr =
{ "pro_ang1_ddddxdd_vvvvxvv", (void *) pro_ang1_ddddxdd_vvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int pro_ang1_cv_fffffxff_vvvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tin4 =  (Float32 *) buffers[4];
    Float32    *tout0 =  (Float32 *) buffers[5];
    Float32    *tout1 =  (Float32 *) buffers[6];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        prolate_aswfa_wrap(*tin0, *tin1, *tin2, *tin3, *tin4, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tin3; ++tin4; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pro_ang1_cv_fffffxff_vvvvvxvv_descr =
{ "pro_ang1_cv_fffffxff_vvvvvxvv", (void *) pro_ang1_cv_fffffxff_vvvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 5, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0, 0, 0 } };

static int pro_ang1_cv_dddddxdd_vvvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tin4 =  (Float64 *) buffers[4];
    Float64    *tout0 =  (Float64 *) buffers[5];
    Float64    *tout1 =  (Float64 *) buffers[6];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        prolate_aswfa_wrap(*tin0, *tin1, *tin2, *tin3, *tin4, tout0, tout1);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tin4; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pro_ang1_cv_dddddxdd_vvvvvxvv_descr =
{ "pro_ang1_cv_dddddxdd_vvvvvxvv", (void *) pro_ang1_cv_dddddxdd_vvvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 5, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0, 0, 0 } };

static int pro_cv_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = prolate_segv_wrap(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pro_cv_fffxf_vvvxf_descr =
{ "pro_cv_fffxf_vvvxf", (void *) pro_cv_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int pro_cv_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = prolate_segv_wrap(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pro_cv_dddxd_vvvxf_descr =
{ "pro_cv_dddxd_vvvxf", (void *) pro_cv_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int pro_rad1_ffffxff_vvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tout0 =  (Float32 *) buffers[4];
    Float32    *tout1 =  (Float32 *) buffers[5];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        prolate_radial1_nocv_wrap(*tin0, *tin1, *tin2, *tin3, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pro_rad1_ffffxff_vvvvxvv_descr =
{ "pro_rad1_ffffxff_vvvvxvv", (void *) pro_rad1_ffffxff_vvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int pro_rad1_ddddxdd_vvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tout0 =  (Float64 *) buffers[4];
    Float64    *tout1 =  (Float64 *) buffers[5];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        prolate_radial1_nocv_wrap(*tin0, *tin1, *tin2, *tin3, tout0, tout1);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pro_rad1_ddddxdd_vvvvxvv_descr =
{ "pro_rad1_ddddxdd_vvvvxvv", (void *) pro_rad1_ddddxdd_vvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int pro_rad1_cv_fffffxff_vvvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tin4 =  (Float32 *) buffers[4];
    Float32    *tout0 =  (Float32 *) buffers[5];
    Float32    *tout1 =  (Float32 *) buffers[6];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        prolate_radial1_wrap(*tin0, *tin1, *tin2, *tin3, *tin4, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tin3; ++tin4; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pro_rad1_cv_fffffxff_vvvvvxvv_descr =
{ "pro_rad1_cv_fffffxff_vvvvvxvv", (void *) pro_rad1_cv_fffffxff_vvvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 5, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0, 0, 0 } };

static int pro_rad1_cv_dddddxdd_vvvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tin4 =  (Float64 *) buffers[4];
    Float64    *tout0 =  (Float64 *) buffers[5];
    Float64    *tout1 =  (Float64 *) buffers[6];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        prolate_radial1_wrap(*tin0, *tin1, *tin2, *tin3, *tin4, tout0, tout1);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tin4; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pro_rad1_cv_dddddxdd_vvvvvxvv_descr =
{ "pro_rad1_cv_dddddxdd_vvvvvxvv", (void *) pro_rad1_cv_dddddxdd_vvvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 5, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0, 0, 0 } };

static int pro_rad2_ffffxff_vvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tout0 =  (Float32 *) buffers[4];
    Float32    *tout1 =  (Float32 *) buffers[5];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        prolate_radial2_nocv_wrap(*tin0, *tin1, *tin2, *tin3, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pro_rad2_ffffxff_vvvvxvv_descr =
{ "pro_rad2_ffffxff_vvvvxvv", (void *) pro_rad2_ffffxff_vvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int pro_rad2_ddddxdd_vvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tout0 =  (Float64 *) buffers[4];
    Float64    *tout1 =  (Float64 *) buffers[5];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        prolate_radial2_nocv_wrap(*tin0, *tin1, *tin2, *tin3, tout0, tout1);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pro_rad2_ddddxdd_vvvvxvv_descr =
{ "pro_rad2_ddddxdd_vvvvxvv", (void *) pro_rad2_ddddxdd_vvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 4, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0, 0 } };

static int pro_rad2_cv_fffffxff_vvvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tin3 =  (Float32 *) buffers[3];
    Float32    *tin4 =  (Float32 *) buffers[4];
    Float32    *tout0 =  (Float32 *) buffers[5];
    Float32    *tout1 =  (Float32 *) buffers[6];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        prolate_radial2_wrap(*tin0, *tin1, *tin2, *tin3, *tin4, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tin1; ++tin2; ++tin3; ++tin4; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pro_rad2_cv_fffffxff_vvvvvxvv_descr =
{ "pro_rad2_cv_fffffxff_vvvvvxvv", (void *) pro_rad2_cv_fffffxff_vvvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 5, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0, 0, 0, 0 } };

static int pro_rad2_cv_dddddxdd_vvvvvxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tin3 =  (Float64 *) buffers[3];
    Float64    *tin4 =  (Float64 *) buffers[4];
    Float64    *tout0 =  (Float64 *) buffers[5];
    Float64    *tout1 =  (Float64 *) buffers[6];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        prolate_radial2_wrap(*tin0, *tin1, *tin2, *tin3, *tin4, tout0, tout1);
        
	++tin0; ++tin1; ++tin2; ++tin3; ++tin4; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor pro_rad2_cv_dddddxdd_vvvvvxvv_descr =
{ "pro_rad2_cv_dddddxdd_vvvvvxvv", (void *) pro_rad2_cv_dddddxdd_vvvvvxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 5, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0, 0, 0, 0 } };

static int psi_FxF_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex32  *tin0 =  (Complex32 *) buffers[0];
    Complex32  *tout0 =  (Complex32 *) buffers[1];
    Complex64     input0;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input0.r = tin0->r;
        input0.i = tin0->i;

        result0 = cpsi_wrap(input0);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor psi_FxF_vxf_descr =
{ "psi_FxF_vxf", (void *) psi_FxF_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0 } };

static int psi_DxD_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex64  *tin0 =  (Complex64 *) buffers[0];
    Complex64  *tout0 =  (Complex64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cpsi_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor psi_DxD_vxf_descr =
{ "psi_DxD_vxf", (void *) psi_DxD_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0 } };

static int psi_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = psi(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor psi_fxf_vxf_descr =
{ "psi_fxf_vxf", (void *) psi_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int psi_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = psi(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor psi_dxd_vxf_descr =
{ "psi_dxd_vxf", (void *) psi_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int radian_fffxf_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tin2 =  (Float32 *) buffers[2];
    Float32    *tout0 =  (Float32 *) buffers[3];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = radian(*tin0, *tin1, *tin2);
                *tout0 = result0;

	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor radian_fffxf_vvvxf_descr =
{ "radian_fffxf_vvvxf", (void *) radian_fffxf_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0, 0 } };

static int radian_dddxd_vvvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tin2 =  (Float64 *) buffers[2];
    Float64    *tout0 =  (Float64 *) buffers[3];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = radian(*tin0, *tin1, *tin2);
        
	++tin0; ++tin1; ++tin2; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor radian_dddxd_vvvxf_descr =
{ "radian_dddxd_vvvxf", (void *) radian_dddxd_vvvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 3, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0, 0 } };

static int rgamma_FxF_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex32  *tin0 =  (Complex32 *) buffers[0];
    Complex32  *tout0 =  (Complex32 *) buffers[1];
    Complex64     input0;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input0.r = tin0->r;
        input0.i = tin0->i;

        result0 = crgamma_wrap(input0);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor rgamma_FxF_vxf_descr =
{ "rgamma_FxF_vxf", (void *) rgamma_FxF_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0 } };

static int rgamma_DxD_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex64  *tin0 =  (Complex64 *) buffers[0];
    Complex64  *tout0 =  (Complex64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = crgamma_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor rgamma_DxD_vxf_descr =
{ "rgamma_DxD_vxf", (void *) rgamma_DxD_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0 } };

static int rgamma_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = rgamma(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor rgamma_fxf_vxf_descr =
{ "rgamma_fxf_vxf", (void *) rgamma_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int rgamma_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = rgamma(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor rgamma_dxd_vxf_descr =
{ "rgamma_dxd_vxf", (void *) rgamma_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int round_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = round(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor round_fxf_vxf_descr =
{ "round_fxf_vxf", (void *) round_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int round_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = round(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor round_dxd_vxf_descr =
{ "round_dxd_vxf", (void *) round_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int shichi_fxff_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float32    *tout1 =  (Float32 *) buffers[2];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        shichi(*tin0, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor shichi_fxff_vxvv_descr =
{ "shichi_fxff_vxvv", (void *) shichi_fxff_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int shichi_dxdd_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];
    Float64    *tout1 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        shichi(*tin0, tout0, tout1);
        
	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor shichi_dxdd_vxvv_descr =
{ "shichi_dxdd_vxvv", (void *) shichi_dxdd_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int sici_fxff_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float32    *tout1 =  (Float32 *) buffers[2];
    Float64     result0;
    Float64     result1;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        sici(*tin0, &result0, &result1);
                *tout0 = result0;
        *tout1 = result1;

	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor sici_fxff_vxvv_descr =
{ "sici_fxff_vxvv", (void *) sici_fxff_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int sici_dxdd_vxvv(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];
    Float64    *tout1 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        sici(*tin0, tout0, tout1);
        
	++tin0; ++tout0; ++tout1; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor sici_dxdd_vxvv_descr =
{ "sici_dxdd_vxvv", (void *) sici_dxdd_vxvv, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 2,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int sindg_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = sindg(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor sindg_fxf_vxf_descr =
{ "sindg_fxf_vxf", (void *) sindg_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int sindg_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = sindg(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor sindg_dxd_vxf_descr =
{ "sindg_dxd_vxf", (void *) sindg_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int smirnov_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = smirnov(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor smirnov_ffxf_vvxf_descr =
{ "smirnov_ffxf_vvxf", (void *) smirnov_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int smirnov_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = smirnov(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor smirnov_ddxd_vvxf_descr =
{ "smirnov_ddxd_vvxf", (void *) smirnov_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int smirnovi_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = smirnovi(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor smirnovi_ffxf_vvxf_descr =
{ "smirnovi_ffxf_vvxf", (void *) smirnovi_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int smirnovi_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = smirnovi(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor smirnovi_ddxd_vvxf_descr =
{ "smirnovi_ddxd_vvxf", (void *) smirnovi_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int spence_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = spence(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor spence_fxf_vxf_descr =
{ "spence_fxf_vxf", (void *) spence_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int spence_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = spence(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor spence_dxd_vxf_descr =
{ "spence_dxd_vxf", (void *) spence_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int stdtr_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdft1_wrap(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor stdtr_ffxf_vvxf_descr =
{ "stdtr_ffxf_vvxf", (void *) stdtr_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int stdtr_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdft1_wrap(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor stdtr_ddxd_vvxf_descr =
{ "stdtr_ddxd_vvxf", (void *) stdtr_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int stdtri_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = stdtri(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor stdtri_ffxf_vvxf_descr =
{ "stdtri_ffxf_vvxf", (void *) stdtri_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int stdtri_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = stdtri(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor stdtri_ddxd_vvxf_descr =
{ "stdtri_ddxd_vvxf", (void *) stdtri_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int stdtridf_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdft3_wrap(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor stdtridf_ffxf_vvxf_descr =
{ "stdtridf_ffxf_vvxf", (void *) stdtridf_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int stdtridf_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdft3_wrap(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor stdtridf_ddxd_vvxf_descr =
{ "stdtridf_ddxd_vvxf", (void *) stdtridf_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int stdtrit_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = cdft2_wrap(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor stdtrit_ffxf_vvxf_descr =
{ "stdtrit_ffxf_vvxf", (void *) stdtrit_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int stdtrit_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cdft2_wrap(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor stdtrit_ddxd_vvxf_descr =
{ "stdtrit_ddxd_vvxf", (void *) stdtrit_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int struve_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = struve(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor struve_ffxf_vvxf_descr =
{ "struve_ffxf_vvxf", (void *) struve_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int struve_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = struve(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor struve_ddxd_vvxf_descr =
{ "struve_ddxd_vvxf", (void *) struve_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int tandg_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = tandg(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor tandg_fxf_vxf_descr =
{ "tandg_fxf_vxf", (void *) tandg_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int tandg_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = tandg(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor tandg_dxd_vxf_descr =
{ "tandg_dxd_vxf", (void *) tandg_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int tklmbda_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = tukeylambdacdf(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor tklmbda_ffxf_vvxf_descr =
{ "tklmbda_ffxf_vvxf", (void *) tklmbda_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int tklmbda_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = tukeylambdacdf(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor tklmbda_ddxd_vvxf_descr =
{ "tklmbda_ddxd_vvxf", (void *) tklmbda_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int wofz_FxF_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex32  *tin0 =  (Complex32 *) buffers[0];
    Complex32  *tout0 =  (Complex32 *) buffers[1];
    Complex64     input0;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input0.r = tin0->r;
        input0.i = tin0->i;

        result0 = cwofz_wrap(input0);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor wofz_FxF_vxf_descr =
{ "wofz_FxF_vxf", (void *) wofz_FxF_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0 } };

static int wofz_DxD_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Complex64  *tin0 =  (Complex64 *) buffers[0];
    Complex64  *tout0 =  (Complex64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cwofz_wrap(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor wofz_DxD_vxf_descr =
{ "wofz_DxD_vxf", (void *) wofz_DxD_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0 } };

static int y0_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = y0(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor y0_fxf_vxf_descr =
{ "y0_fxf_vxf", (void *) y0_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int y0_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = y0(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor y0_dxd_vxf_descr =
{ "y0_dxd_vxf", (void *) y0_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int y1_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = y1(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor y1_fxf_vxf_descr =
{ "y1_fxf_vxf", (void *) y1_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int y1_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = y1(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor y1_dxd_vxf_descr =
{ "y1_dxd_vxf", (void *) y1_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static int yn_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = yn(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor yn_ffxf_vvxf_descr =
{ "yn_ffxf_vvxf", (void *) yn_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int yn_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = yn(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor yn_ddxd_vvxf_descr =
{ "yn_ddxd_vvxf", (void *) yn_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int yv_fFxF_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Complex32  *tin1 =  (Complex32 *) buffers[1];
    Complex32  *tout0 =  (Complex32 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = tin1->r;
        input1.i = tin1->i;

        result0 = cbesy_wrap(*tin0, input1);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor yv_fFxF_vvxf_descr =
{ "yv_fFxF_vvxf", (void *) yv_fFxF_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0 } };

static int yv_dDxD_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Complex64  *tin1 =  (Complex64 *) buffers[1];
    Complex64  *tout0 =  (Complex64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cbesy_wrap(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor yv_dDxD_vvxf_descr =
{ "yv_dDxD_vvxf", (void *) yv_dDxD_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0 } };

static int yv_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = yv(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor yv_ffxf_vvxf_descr =
{ "yv_ffxf_vvxf", (void *) yv_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int yv_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = yv(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor yv_ddxd_vvxf_descr =
{ "yv_ddxd_vvxf", (void *) yv_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int yve_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = *tin1;
        input1.i = 0;

        result0 = cbesy_wrap_e(*tin0, input1);
                *tout0 = result0.r;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor yve_ffxf_vvxf_descr =
{ "yve_ffxf_vvxf", (void *) yve_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int yve_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = *tin1;
        input1.i = 0;

        result0 = cbesy_wrap_e(*tin0, input1);
                *tout0 = result0.r;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor yve_ddxd_vvxf_descr =
{ "yve_ddxd_vvxf", (void *) yve_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int yve_fFxF_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Complex32  *tin1 =  (Complex32 *) buffers[1];
    Complex32  *tout0 =  (Complex32 *) buffers[2];
    Complex64     input1;
    Complex64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
                input1.r = tin1->r;
        input1.i = tin1->i;

        result0 = cbesy_wrap_e(*tin0, input1);
                tout0->r = result0.r;
        tout0->i = result0.i;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor yve_fFxF_vvxf_descr =
{ "yve_fFxF_vvxf", (void *) yve_fFxF_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Complex32), sizeof(Complex32) }, { 0, 0, 0, 0 } };

static int yve_dDxD_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Complex64  *tin1 =  (Complex64 *) buffers[1];
    Complex64  *tout0 =  (Complex64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = cbesy_wrap_e(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor yve_dDxD_vvxf_descr =
{ "yve_dDxD_vvxf", (void *) yve_dDxD_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Complex64), sizeof(Complex64) }, { 0, 0, 0, 0 } };

static int zeta_ffxf_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tin1 =  (Float32 *) buffers[1];
    Float32    *tout0 =  (Float32 *) buffers[2];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = zeta(*tin0, *tin1);
                *tout0 = result0;

	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor zeta_ffxf_vvxf_descr =
{ "zeta_ffxf_vvxf", (void *) zeta_ffxf_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float32), sizeof(Float32), sizeof(Float32) }, { 0, 0, 0, 0 } };

static int zeta_ddxd_vvxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tin1 =  (Float64 *) buffers[1];
    Float64    *tout0 =  (Float64 *) buffers[2];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = zeta(*tin0, *tin1);
        
	++tin0; ++tin1; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor zeta_ddxd_vvxf_descr =
{ "zeta_ddxd_vvxf", (void *) zeta_ddxd_vvxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 2, 1,
  { sizeof(Float64), sizeof(Float64), sizeof(Float64) }, { 0, 0, 0, 0 } };

static int zetac_fxf_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float32    *tin0 =  (Float32 *) buffers[0];
    Float32    *tout0 =  (Float32 *) buffers[1];
    Float64     result0;

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        result0 = zetac(*tin0);
                *tout0 = result0;

	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor zetac_fxf_vxf_descr =
{ "zetac_fxf_vxf", (void *) zetac_fxf_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float32), sizeof(Float32) }, { 0, 0, 0 } };

static int zetac_dxd_vxf(long niter, long ninargs, long noutargs, void **buffers, long *bsizes) {
    long i;
    Float64    *tin0 =  (Float64 *) buffers[0];
    Float64    *tout0 =  (Float64 *) buffers[1];

    BEGIN_THREADS
    for (i=0; i<niter; i++) {
        
        *tout0 = zetac(*tin0);
        
	++tin0; ++tout0; 
    }
    END_THREADS
    return 0;
}

static CfuncDescriptor zetac_dxd_vxf_descr =
{ "zetac_dxd_vxf", (void *) zetac_dxd_vxf, CFUNC_UFUNC, 0, CHECK_ALIGN, 1, 1,
  { sizeof(Float64), sizeof(Float64) }, { 0, 0, 0 } };

static PyMethodDef _na_cephesMethods[] = {

	{NULL,      NULL}        /* Sentinel */
};

static PyObject *init_funcDict(void) {
    PyObject *dict, *keytuple;
    dict = PyDict_New();
    /* airy_fxffff_vxvvvv */
    keytuple=Py_BuildValue("ss((s)(ssss))","airy","v","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&airy_fxffff_vxvvvv_descr));

    /* airy_dxdddd_vxvvvv */
    keytuple=Py_BuildValue("ss((s)(ssss))","airy","v","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&airy_dxdddd_vxvvvv_descr));

    /* airy_FxFFFF_vxvvvv */
    keytuple=Py_BuildValue("ss((s)(ssss))","airy","v","Complex32","Complex32","Complex32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&airy_FxFFFF_vxvvvv_descr));

    /* airy_DxDDDD_vxvvvv */
    keytuple=Py_BuildValue("ss((s)(ssss))","airy","v","Complex64","Complex64","Complex64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&airy_DxDDDD_vxvvvv_descr));

    /* airye_fxffff_vxvvvv */
    keytuple=Py_BuildValue("ss((s)(ssss))","airye","v","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&airye_fxffff_vxvvvv_descr));

    /* airye_dxdddd_vxvvvv */
    keytuple=Py_BuildValue("ss((s)(ssss))","airye","v","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&airye_dxdddd_vxvvvv_descr));

    /* airye_FxFFFF_vxvvvv */
    keytuple=Py_BuildValue("ss((s)(ssss))","airye","v","Complex32","Complex32","Complex32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&airye_FxFFFF_vxvvvv_descr));

    /* airye_DxDDDD_vxvvvv */
    keytuple=Py_BuildValue("ss((s)(ssss))","airye","v","Complex64","Complex64","Complex64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&airye_DxDDDD_vxvvvv_descr));

    /* bdtr_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","bdtr","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&bdtr_fffxf_vvvxf_descr));

    /* bdtr_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","bdtr","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&bdtr_dddxd_vvvxf_descr));

    /* bdtrc_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","bdtrc","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&bdtrc_fffxf_vvvxf_descr));

    /* bdtrc_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","bdtrc","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&bdtrc_dddxd_vvvxf_descr));

    /* bdtri_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","bdtri","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&bdtri_fffxf_vvvxf_descr));

    /* bdtri_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","bdtri","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&bdtri_dddxd_vvvxf_descr));

    /* bdtrik_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","bdtrik","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&bdtrik_fffxf_vvvxf_descr));

    /* bdtrik_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","bdtrik","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&bdtrik_dddxd_vvvxf_descr));

    /* bdtrin_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","bdtrin","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&bdtrin_fffxf_vvvxf_descr));

    /* bdtrin_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","bdtrin","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&bdtrin_dddxd_vvvxf_descr));

    /* bei_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","bei","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&bei_fxf_vxf_descr));

    /* bei_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","bei","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&bei_dxd_vxf_descr));

    /* beip_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","beip","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&beip_fxf_vxf_descr));

    /* beip_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","beip","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&beip_dxd_vxf_descr));

    /* ber_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","ber","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ber_fxf_vxf_descr));

    /* ber_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","ber","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ber_dxd_vxf_descr));

    /* berp_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","berp","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&berp_fxf_vxf_descr));

    /* berp_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","berp","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&berp_dxd_vxf_descr));

    /* besselpoly_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","besselpoly","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&besselpoly_fffxf_vvvxf_descr));

    /* besselpoly_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","besselpoly","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&besselpoly_dddxd_vvvxf_descr));

    /* beta_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","beta","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&beta_ffxf_vvxf_descr));

    /* beta_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","beta","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&beta_ddxd_vvxf_descr));

    /* betainc_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","betainc","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&betainc_fffxf_vvvxf_descr));

    /* betainc_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","betainc","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&betainc_dddxd_vvvxf_descr));

    /* betaincinv_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","betaincinv","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&betaincinv_fffxf_vvvxf_descr));

    /* betaincinv_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","betaincinv","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&betaincinv_dddxd_vvvxf_descr));

    /* betaln_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","betaln","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&betaln_ffxf_vvxf_descr));

    /* betaln_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","betaln","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&betaln_ddxd_vvxf_descr));

    /* btdtr_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","btdtr","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&btdtr_fffxf_vvvxf_descr));

    /* btdtr_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","btdtr","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&btdtr_dddxd_vvvxf_descr));

    /* btdtri_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","btdtri","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&btdtri_fffxf_vvvxf_descr));

    /* btdtri_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","btdtri","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&btdtri_dddxd_vvvxf_descr));

    /* btdtria_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","btdtria","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&btdtria_fffxf_vvvxf_descr));

    /* btdtria_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","btdtria","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&btdtria_dddxd_vvvxf_descr));

    /* btdtrib_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","btdtrib","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&btdtrib_fffxf_vvvxf_descr));

    /* btdtrib_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","btdtrib","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&btdtrib_dddxd_vvvxf_descr));

    /* cbrt_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","cbrt","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&cbrt_fxf_vxf_descr));

    /* cbrt_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","cbrt","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&cbrt_dxd_vxf_descr));

    /* chdtr_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","chdtr","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&chdtr_ffxf_vvxf_descr));

    /* chdtr_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","chdtr","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&chdtr_ddxd_vvxf_descr));

    /* chdtrc_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","chdtrc","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&chdtrc_ffxf_vvxf_descr));

    /* chdtrc_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","chdtrc","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&chdtrc_ddxd_vvxf_descr));

    /* chdtri_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","chdtri","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&chdtri_ffxf_vvxf_descr));

    /* chdtri_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","chdtri","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&chdtri_ddxd_vvxf_descr));

    /* chdtriv_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","chdtriv","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&chdtriv_ffxf_vvxf_descr));

    /* chdtriv_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","chdtriv","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&chdtriv_ddxd_vvxf_descr));

    /* chndtr_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","chndtr","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&chndtr_fffxf_vvvxf_descr));

    /* chndtr_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","chndtr","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&chndtr_dddxd_vvvxf_descr));

    /* chndtridf_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","chndtridf","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&chndtridf_fffxf_vvvxf_descr));

    /* chndtridf_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","chndtridf","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&chndtridf_dddxd_vvvxf_descr));

    /* chndtrinc_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","chndtrinc","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&chndtrinc_fffxf_vvvxf_descr));

    /* chndtrinc_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","chndtrinc","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&chndtrinc_dddxd_vvvxf_descr));

    /* chndtrix_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","chndtrix","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&chndtrix_fffxf_vvvxf_descr));

    /* chndtrix_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","chndtrix","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&chndtrix_dddxd_vvvxf_descr));

    /* cosdg_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","cosdg","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&cosdg_fxf_vxf_descr));

    /* cosdg_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","cosdg","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&cosdg_dxd_vxf_descr));

    /* cosm1_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","cosm1","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&cosm1_fxf_vxf_descr));

    /* cosm1_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","cosm1","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&cosm1_dxd_vxf_descr));

    /* cotdg_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","cotdg","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&cotdg_fxf_vxf_descr));

    /* cotdg_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","cotdg","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&cotdg_dxd_vxf_descr));

    /* dawsn_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","dawsn","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&dawsn_fxf_vxf_descr));

    /* dawsn_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","dawsn","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&dawsn_dxd_vxf_descr));

    /* ellipe_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","ellipe","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ellipe_fxf_vxf_descr));

    /* ellipe_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","ellipe","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ellipe_dxd_vxf_descr));

    /* ellipeinc_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","ellipeinc","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ellipeinc_ffxf_vvxf_descr));

    /* ellipeinc_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","ellipeinc","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ellipeinc_ddxd_vvxf_descr));

    /* ellipj_ffxffff_vvxvvvv */
    keytuple=Py_BuildValue("ss((ss)(ssss))","ellipj","vv","Float32","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ellipj_ffxffff_vvxvvvv_descr));

    /* ellipj_ddxdddd_vvxvvvv */
    keytuple=Py_BuildValue("ss((ss)(ssss))","ellipj","vv","Float64","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ellipj_ddxdddd_vvxvvvv_descr));

    /* ellipk_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","ellipk","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ellipk_fxf_vxf_descr));

    /* ellipk_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","ellipk","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ellipk_dxd_vxf_descr));

    /* ellipkinc_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","ellipkinc","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ellipkinc_ffxf_vvxf_descr));

    /* ellipkinc_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","ellipkinc","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ellipkinc_ddxd_vvxf_descr));

    /* erf_FxF_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","erf","v","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&erf_FxF_vxf_descr));

    /* erf_DxD_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","erf","v","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&erf_DxD_vxf_descr));

    /* erf_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","erf","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&erf_fxf_vxf_descr));

    /* erf_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","erf","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&erf_dxd_vxf_descr));

    /* erfc_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","erfc","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&erfc_fxf_vxf_descr));

    /* erfc_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","erfc","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&erfc_dxd_vxf_descr));

    /* exp1_FxF_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","exp1","v","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&exp1_FxF_vxf_descr));

    /* exp1_DxD_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","exp1","v","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&exp1_DxD_vxf_descr));

    /* exp1_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","exp1","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&exp1_fxf_vxf_descr));

    /* exp1_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","exp1","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&exp1_dxd_vxf_descr));

    /* exp10_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","exp10","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&exp10_fxf_vxf_descr));

    /* exp10_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","exp10","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&exp10_dxd_vxf_descr));

    /* exp2_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","exp2","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&exp2_fxf_vxf_descr));

    /* exp2_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","exp2","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&exp2_dxd_vxf_descr));

    /* expi_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","expi","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&expi_fxf_vxf_descr));

    /* expi_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","expi","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&expi_dxd_vxf_descr));

    /* expm1_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","expm1","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&expm1_fxf_vxf_descr));

    /* expm1_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","expm1","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&expm1_dxd_vxf_descr));

    /* expn_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","expn","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&expn_ffxf_vvxf_descr));

    /* expn_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","expn","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&expn_ddxd_vvxf_descr));

    /* fdtr_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","fdtr","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&fdtr_fffxf_vvvxf_descr));

    /* fdtr_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","fdtr","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&fdtr_dddxd_vvvxf_descr));

    /* fdtrc_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","fdtrc","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&fdtrc_fffxf_vvvxf_descr));

    /* fdtrc_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","fdtrc","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&fdtrc_dddxd_vvvxf_descr));

    /* fdtri_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","fdtri","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&fdtri_fffxf_vvvxf_descr));

    /* fdtri_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","fdtri","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&fdtri_dddxd_vvvxf_descr));

    /* fdtridfd_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","fdtridfd","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&fdtridfd_fffxf_vvvxf_descr));

    /* fdtridfd_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","fdtridfd","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&fdtridfd_dddxd_vvvxf_descr));

    /* fdtridfn_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","fdtridfn","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&fdtridfn_fffxf_vvvxf_descr));

    /* fdtridfn_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","fdtridfn","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&fdtridfn_dddxd_vvvxf_descr));

    /* fresnel_FxFF_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","fresnel","v","Complex32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&fresnel_FxFF_vxvv_descr));

    /* fresnel_DxDD_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","fresnel","v","Complex64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&fresnel_DxDD_vxvv_descr));

    /* fresnel_fxff_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","fresnel","v","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&fresnel_fxff_vxvv_descr));

    /* fresnel_dxdd_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","fresnel","v","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&fresnel_dxdd_vxvv_descr));

    /* gamma_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","gamma","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gamma_fxf_vxf_descr));

    /* gamma_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","gamma","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gamma_dxd_vxf_descr));

    /* gamma_FxF_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","gamma","v","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gamma_FxF_vxf_descr));

    /* gamma_DxD_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","gamma","v","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gamma_DxD_vxf_descr));

    /* gammainc_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","gammainc","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gammainc_ffxf_vvxf_descr));

    /* gammainc_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","gammainc","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gammainc_ddxd_vvxf_descr));

    /* gammaincc_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","gammaincc","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gammaincc_ffxf_vvxf_descr));

    /* gammaincc_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","gammaincc","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gammaincc_ddxd_vvxf_descr));

    /* gammainccinv_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","gammainccinv","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gammainccinv_ffxf_vvxf_descr));

    /* gammainccinv_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","gammainccinv","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gammainccinv_ddxd_vvxf_descr));

    /* gammaln_FxF_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","gammaln","v","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gammaln_FxF_vxf_descr));

    /* gammaln_DxD_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","gammaln","v","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gammaln_DxD_vxf_descr));

    /* gammaln_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","gammaln","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gammaln_fxf_vxf_descr));

    /* gammaln_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","gammaln","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gammaln_dxd_vxf_descr));

    /* gdtr_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","gdtr","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gdtr_fffxf_vvvxf_descr));

    /* gdtr_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","gdtr","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gdtr_dddxd_vvvxf_descr));

    /* gdtr2_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","gdtr2","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gdtr2_fffxf_vvvxf_descr));

    /* gdtr2_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","gdtr2","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gdtr2_dddxd_vvvxf_descr));

    /* gdtrc_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","gdtrc","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gdtrc_fffxf_vvvxf_descr));

    /* gdtrc_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","gdtrc","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gdtrc_dddxd_vvvxf_descr));

    /* gdtri_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","gdtri","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gdtri_fffxf_vvvxf_descr));

    /* gdtri_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","gdtri","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gdtri_dddxd_vvvxf_descr));

    /* gdtria_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","gdtria","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gdtria_fffxf_vvvxf_descr));

    /* gdtria_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","gdtria","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gdtria_dddxd_vvvxf_descr));

    /* gdtrib_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","gdtrib","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gdtrib_fffxf_vvvxf_descr));

    /* gdtrib_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","gdtrib","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gdtrib_dddxd_vvvxf_descr));

    /* gdtrix_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","gdtrix","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gdtrix_fffxf_vvvxf_descr));

    /* gdtrix_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","gdtrix","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&gdtrix_dddxd_vvvxf_descr));

    /* hankel1_fFxF_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","hankel1","vv","Float32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hankel1_fFxF_vvxf_descr));

    /* hankel1_dDxD_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","hankel1","vv","Float64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hankel1_dDxD_vvxf_descr));

    /* hankel1e_fFxF_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","hankel1e","vv","Float32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hankel1e_fFxF_vvxf_descr));

    /* hankel1e_dDxD_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","hankel1e","vv","Float64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hankel1e_dDxD_vvxf_descr));

    /* hankel2_fFxF_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","hankel2","vv","Float32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hankel2_fFxF_vvxf_descr));

    /* hankel2_dDxD_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","hankel2","vv","Float64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hankel2_dDxD_vvxf_descr));

    /* hankel2e_fFxF_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","hankel2e","vv","Float32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hankel2e_fFxF_vvxf_descr));

    /* hankel2e_dDxD_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","hankel2e","vv","Float64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hankel2e_dDxD_vvxf_descr));

    /* hyp1f1_ffFxF_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","hyp1f1","vvv","Float32","Float32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hyp1f1_ffFxF_vvvxf_descr));

    /* hyp1f1_ddDxD_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","hyp1f1","vvv","Float64","Float64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hyp1f1_ddDxD_vvvxf_descr));

    /* hyp1f1_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","hyp1f1","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hyp1f1_fffxf_vvvxf_descr));

    /* hyp1f1_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","hyp1f1","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hyp1f1_dddxd_vvvxf_descr));

    /* hyp1f2_ffffxff_vvvvxfv */
    keytuple=Py_BuildValue("ss((ssss)(ss))","hyp1f2","vvvv","Float32","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hyp1f2_ffffxff_vvvvxfv_descr));

    /* hyp1f2_ddddxdd_vvvvxfv */
    keytuple=Py_BuildValue("ss((ssss)(ss))","hyp1f2","vvvv","Float64","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hyp1f2_ddddxdd_vvvvxfv_descr));

    /* hyp2f0_ffffxff_vvvvxfv */
    keytuple=Py_BuildValue("ss((ssss)(ss))","hyp2f0","vvvv","Float32","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hyp2f0_ffffxff_vvvvxfv_descr));

    /* hyp2f0_ddddxdd_vvvvxfv */
    keytuple=Py_BuildValue("ss((ssss)(ss))","hyp2f0","vvvv","Float64","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hyp2f0_ddddxdd_vvvvxfv_descr));

    /* hyp2f1_fffFxF_vvvvxf */
    keytuple=Py_BuildValue("ss((ssss)(s))","hyp2f1","vvvv","Float32","Float32","Float32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hyp2f1_fffFxF_vvvvxf_descr));

    /* hyp2f1_dddDxD_vvvvxf */
    keytuple=Py_BuildValue("ss((ssss)(s))","hyp2f1","vvvv","Float64","Float64","Float64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hyp2f1_dddDxD_vvvvxf_descr));

    /* hyp2f1_ffffxf_vvvvxf */
    keytuple=Py_BuildValue("ss((ssss)(s))","hyp2f1","vvvv","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hyp2f1_ffffxf_vvvvxf_descr));

    /* hyp2f1_ddddxd_vvvvxf */
    keytuple=Py_BuildValue("ss((ssss)(s))","hyp2f1","vvvv","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hyp2f1_ddddxd_vvvvxf_descr));

    /* hyp3f0_ffffxff_vvvvxfv */
    keytuple=Py_BuildValue("ss((ssss)(ss))","hyp3f0","vvvv","Float32","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hyp3f0_ffffxff_vvvvxfv_descr));

    /* hyp3f0_ddddxdd_vvvvxfv */
    keytuple=Py_BuildValue("ss((ssss)(ss))","hyp3f0","vvvv","Float64","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hyp3f0_ddddxdd_vvvvxfv_descr));

    /* hyperu_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","hyperu","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hyperu_fffxf_vvvxf_descr));

    /* hyperu_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","hyperu","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&hyperu_dddxd_vvvxf_descr));

    /* i0_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","i0","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&i0_fxf_vxf_descr));

    /* i0_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","i0","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&i0_dxd_vxf_descr));

    /* i0e_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","i0e","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&i0e_fxf_vxf_descr));

    /* i0e_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","i0e","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&i0e_dxd_vxf_descr));

    /* i1_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","i1","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&i1_fxf_vxf_descr));

    /* i1_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","i1","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&i1_dxd_vxf_descr));

    /* i1e_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","i1e","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&i1e_fxf_vxf_descr));

    /* i1e_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","i1e","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&i1e_dxd_vxf_descr));

    /* it2i0k0_fxff_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","it2i0k0","v","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&it2i0k0_fxff_vxvv_descr));

    /* it2i0k0_dxdd_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","it2i0k0","v","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&it2i0k0_dxdd_vxvv_descr));

    /* it2j0y0_fxff_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","it2j0y0","v","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&it2j0y0_fxff_vxvv_descr));

    /* it2j0y0_dxdd_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","it2j0y0","v","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&it2j0y0_dxdd_vxvv_descr));

    /* it2struve0_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","it2struve0","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&it2struve0_fxf_vxf_descr));

    /* it2struve0_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","it2struve0","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&it2struve0_dxd_vxf_descr));

    /* itairy_fxffff_vxvvvv */
    keytuple=Py_BuildValue("ss((s)(ssss))","itairy","v","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&itairy_fxffff_vxvvvv_descr));

    /* itairy_dxdddd_vxvvvv */
    keytuple=Py_BuildValue("ss((s)(ssss))","itairy","v","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&itairy_dxdddd_vxvvvv_descr));

    /* iti0k0_fxff_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","iti0k0","v","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&iti0k0_fxff_vxvv_descr));

    /* iti0k0_dxdd_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","iti0k0","v","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&iti0k0_dxdd_vxvv_descr));

    /* itj0y0_fxff_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","itj0y0","v","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&itj0y0_fxff_vxvv_descr));

    /* itj0y0_dxdd_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","itj0y0","v","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&itj0y0_dxdd_vxvv_descr));

    /* itmodstruve0_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","itmodstruve0","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&itmodstruve0_fxf_vxf_descr));

    /* itmodstruve0_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","itmodstruve0","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&itmodstruve0_dxd_vxf_descr));

    /* itstruve0_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","itstruve0","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&itstruve0_fxf_vxf_descr));

    /* itstruve0_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","itstruve0","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&itstruve0_dxd_vxf_descr));

    /* iv_fFxF_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","iv","vv","Float32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&iv_fFxF_vvxf_descr));

    /* iv_dDxD_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","iv","vv","Float64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&iv_dDxD_vvxf_descr));

    /* iv_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","iv","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&iv_ffxf_vvxf_descr));

    /* iv_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","iv","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&iv_ddxd_vvxf_descr));

    /* ive_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","ive","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ive_ffxf_vvxf_descr));

    /* ive_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","ive","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ive_ddxd_vvxf_descr));

    /* ive_fFxF_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","ive","vv","Float32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ive_fFxF_vvxf_descr));

    /* ive_dDxD_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","ive","vv","Float64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ive_dDxD_vvxf_descr));

    /* j0_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","j0","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&j0_fxf_vxf_descr));

    /* j0_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","j0","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&j0_dxd_vxf_descr));

    /* j1_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","j1","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&j1_fxf_vxf_descr));

    /* j1_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","j1","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&j1_dxd_vxf_descr));

    /* jn_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","jn","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&jn_ffxf_vvxf_descr));

    /* jn_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","jn","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&jn_ddxd_vvxf_descr));

    /* jv_fFxF_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","jv","vv","Float32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&jv_fFxF_vvxf_descr));

    /* jv_dDxD_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","jv","vv","Float64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&jv_dDxD_vvxf_descr));

    /* jv_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","jv","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&jv_ffxf_vvxf_descr));

    /* jv_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","jv","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&jv_ddxd_vvxf_descr));

    /* jve_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","jve","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&jve_ffxf_vvxf_descr));

    /* jve_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","jve","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&jve_ddxd_vvxf_descr));

    /* jve_fFxF_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","jve","vv","Float32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&jve_fFxF_vvxf_descr));

    /* jve_dDxD_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","jve","vv","Float64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&jve_dDxD_vvxf_descr));

    /* k0_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","k0","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&k0_fxf_vxf_descr));

    /* k0_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","k0","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&k0_dxd_vxf_descr));

    /* k0e_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","k0e","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&k0e_fxf_vxf_descr));

    /* k0e_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","k0e","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&k0e_dxd_vxf_descr));

    /* k1_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","k1","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&k1_fxf_vxf_descr));

    /* k1_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","k1","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&k1_dxd_vxf_descr));

    /* k1e_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","k1e","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&k1e_fxf_vxf_descr));

    /* k1e_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","k1e","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&k1e_dxd_vxf_descr));

    /* kei_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","kei","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kei_fxf_vxf_descr));

    /* kei_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","kei","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kei_dxd_vxf_descr));

    /* keip_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","keip","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&keip_fxf_vxf_descr));

    /* keip_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","keip","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&keip_dxd_vxf_descr));

    /* kelvin_fxFFFF_vxvvvv */
    keytuple=Py_BuildValue("ss((s)(ssss))","kelvin","v","Float32","Complex32","Complex32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kelvin_fxFFFF_vxvvvv_descr));

    /* kelvin_dxDDDD_vxvvvv */
    keytuple=Py_BuildValue("ss((s)(ssss))","kelvin","v","Float64","Complex64","Complex64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kelvin_dxDDDD_vxvvvv_descr));

    /* ker_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","ker","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ker_fxf_vxf_descr));

    /* ker_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","ker","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ker_dxd_vxf_descr));

    /* kerp_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","kerp","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kerp_fxf_vxf_descr));

    /* kerp_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","kerp","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kerp_dxd_vxf_descr));

    /* kn_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","kn","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kn_ffxf_vvxf_descr));

    /* kn_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","kn","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kn_ddxd_vvxf_descr));

    /* kolmogi_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","kolmogi","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kolmogi_fxf_vxf_descr));

    /* kolmogi_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","kolmogi","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kolmogi_dxd_vxf_descr));

    /* kolmogorov_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","kolmogorov","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kolmogorov_fxf_vxf_descr));

    /* kolmogorov_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","kolmogorov","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kolmogorov_dxd_vxf_descr));

    /* kv_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","kv","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kv_ffxf_vvxf_descr));

    /* kv_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","kv","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kv_ddxd_vvxf_descr));

    /* kv_fFxF_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","kv","vv","Float32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kv_fFxF_vvxf_descr));

    /* kv_dDxD_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","kv","vv","Float64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kv_dDxD_vvxf_descr));

    /* kve_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","kve","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kve_ffxf_vvxf_descr));

    /* kve_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","kve","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kve_ddxd_vvxf_descr));

    /* kve_fFxF_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","kve","vv","Float32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kve_fFxF_vvxf_descr));

    /* kve_dDxD_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","kve","vv","Float64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&kve_dDxD_vvxf_descr));

    /* log1p_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","log1p","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&log1p_fxf_vxf_descr));

    /* log1p_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","log1p","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&log1p_dxd_vxf_descr));

    /* lpmv_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","lpmv","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&lpmv_fffxf_vvvxf_descr));

    /* lpmv_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","lpmv","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&lpmv_dddxd_vvvxf_descr));

    /* mathieu_a_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","mathieu_a","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&mathieu_a_ffxf_vvxf_descr));

    /* mathieu_a_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","mathieu_a","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&mathieu_a_ddxd_vvxf_descr));

    /* mathieu_b_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","mathieu_b","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&mathieu_b_ffxf_vvxf_descr));

    /* mathieu_b_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","mathieu_b","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&mathieu_b_ddxd_vvxf_descr));

    /* mathieu_cem_fffxff_vvvxvv */
    keytuple=Py_BuildValue("ss((sss)(ss))","mathieu_cem","vvv","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&mathieu_cem_fffxff_vvvxvv_descr));

    /* mathieu_cem_dddxdd_vvvxvv */
    keytuple=Py_BuildValue("ss((sss)(ss))","mathieu_cem","vvv","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&mathieu_cem_dddxdd_vvvxvv_descr));

    /* mathieu_modcem1_fffxff_vvvxvv */
    keytuple=Py_BuildValue("ss((sss)(ss))","mathieu_modcem1","vvv","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&mathieu_modcem1_fffxff_vvvxvv_descr));

    /* mathieu_modcem1_dddxdd_vvvxvv */
    keytuple=Py_BuildValue("ss((sss)(ss))","mathieu_modcem1","vvv","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&mathieu_modcem1_dddxdd_vvvxvv_descr));

    /* mathieu_modcem2_fffxff_vvvxvv */
    keytuple=Py_BuildValue("ss((sss)(ss))","mathieu_modcem2","vvv","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&mathieu_modcem2_fffxff_vvvxvv_descr));

    /* mathieu_modcem2_dddxdd_vvvxvv */
    keytuple=Py_BuildValue("ss((sss)(ss))","mathieu_modcem2","vvv","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&mathieu_modcem2_dddxdd_vvvxvv_descr));

    /* mathieu_modsem1_fffxff_vvvxvv */
    keytuple=Py_BuildValue("ss((sss)(ss))","mathieu_modsem1","vvv","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&mathieu_modsem1_fffxff_vvvxvv_descr));

    /* mathieu_modsem1_dddxdd_vvvxvv */
    keytuple=Py_BuildValue("ss((sss)(ss))","mathieu_modsem1","vvv","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&mathieu_modsem1_dddxdd_vvvxvv_descr));

    /* mathieu_modsem2_fffxff_vvvxvv */
    keytuple=Py_BuildValue("ss((sss)(ss))","mathieu_modsem2","vvv","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&mathieu_modsem2_fffxff_vvvxvv_descr));

    /* mathieu_modsem2_dddxdd_vvvxvv */
    keytuple=Py_BuildValue("ss((sss)(ss))","mathieu_modsem2","vvv","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&mathieu_modsem2_dddxdd_vvvxvv_descr));

    /* mathieu_sem_fffxff_vvvxvv */
    keytuple=Py_BuildValue("ss((sss)(ss))","mathieu_sem","vvv","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&mathieu_sem_fffxff_vvvxvv_descr));

    /* mathieu_sem_dddxdd_vvvxvv */
    keytuple=Py_BuildValue("ss((sss)(ss))","mathieu_sem","vvv","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&mathieu_sem_dddxdd_vvvxvv_descr));

    /* modfresnelm_fxFF_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","modfresnelm","v","Float32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&modfresnelm_fxFF_vxvv_descr));

    /* modfresnelm_dxDD_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","modfresnelm","v","Float64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&modfresnelm_dxDD_vxvv_descr));

    /* modfresnelp_fxFF_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","modfresnelp","v","Float32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&modfresnelp_fxFF_vxvv_descr));

    /* modfresnelp_dxDD_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","modfresnelp","v","Float64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&modfresnelp_dxDD_vxvv_descr));

    /* modstruve_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","modstruve","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&modstruve_ffxf_vvxf_descr));

    /* modstruve_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","modstruve","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&modstruve_ddxd_vvxf_descr));

    /* nbdtr_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nbdtr","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nbdtr_fffxf_vvvxf_descr));

    /* nbdtr_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nbdtr","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nbdtr_dddxd_vvvxf_descr));

    /* nbdtrc_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nbdtrc","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nbdtrc_fffxf_vvvxf_descr));

    /* nbdtrc_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nbdtrc","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nbdtrc_dddxd_vvvxf_descr));

    /* nbdtri_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nbdtri","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nbdtri_fffxf_vvvxf_descr));

    /* nbdtri_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nbdtri","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nbdtri_dddxd_vvvxf_descr));

    /* nbdtrik_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nbdtrik","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nbdtrik_fffxf_vvvxf_descr));

    /* nbdtrik_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nbdtrik","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nbdtrik_dddxd_vvvxf_descr));

    /* nbdtrin_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nbdtrin","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nbdtrin_fffxf_vvvxf_descr));

    /* nbdtrin_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nbdtrin","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nbdtrin_dddxd_vvvxf_descr));

    /* ncfdtr_ffffxf_vvvvxf */
    keytuple=Py_BuildValue("ss((ssss)(s))","ncfdtr","vvvv","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ncfdtr_ffffxf_vvvvxf_descr));

    /* ncfdtr_ddddxd_vvvvxf */
    keytuple=Py_BuildValue("ss((ssss)(s))","ncfdtr","vvvv","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ncfdtr_ddddxd_vvvvxf_descr));

    /* ncfdtri_ffffxf_vvvvxf */
    keytuple=Py_BuildValue("ss((ssss)(s))","ncfdtri","vvvv","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ncfdtri_ffffxf_vvvvxf_descr));

    /* ncfdtri_ddddxd_vvvvxf */
    keytuple=Py_BuildValue("ss((ssss)(s))","ncfdtri","vvvv","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ncfdtri_ddddxd_vvvvxf_descr));

    /* ncfdtridfd_ffffxf_vvvvxf */
    keytuple=Py_BuildValue("ss((ssss)(s))","ncfdtridfd","vvvv","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ncfdtridfd_ffffxf_vvvvxf_descr));

    /* ncfdtridfd_ddddxd_vvvvxf */
    keytuple=Py_BuildValue("ss((ssss)(s))","ncfdtridfd","vvvv","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ncfdtridfd_ddddxd_vvvvxf_descr));

    /* ncfdtridfn_ffffxf_vvvvxf */
    keytuple=Py_BuildValue("ss((ssss)(s))","ncfdtridfn","vvvv","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ncfdtridfn_ffffxf_vvvvxf_descr));

    /* ncfdtridfn_ddddxd_vvvvxf */
    keytuple=Py_BuildValue("ss((ssss)(s))","ncfdtridfn","vvvv","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ncfdtridfn_ddddxd_vvvvxf_descr));

    /* ncfdtrinc_ffffxf_vvvvxf */
    keytuple=Py_BuildValue("ss((ssss)(s))","ncfdtrinc","vvvv","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ncfdtrinc_ffffxf_vvvvxf_descr));

    /* ncfdtrinc_ddddxd_vvvvxf */
    keytuple=Py_BuildValue("ss((ssss)(s))","ncfdtrinc","vvvv","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ncfdtrinc_ddddxd_vvvvxf_descr));

    /* nctdtr_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nctdtr","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nctdtr_fffxf_vvvxf_descr));

    /* nctdtr_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nctdtr","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nctdtr_dddxd_vvvxf_descr));

    /* nctdtridf_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nctdtridf","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nctdtridf_fffxf_vvvxf_descr));

    /* nctdtridf_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nctdtridf","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nctdtridf_dddxd_vvvxf_descr));

    /* nctdtrinc_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nctdtrinc","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nctdtrinc_fffxf_vvvxf_descr));

    /* nctdtrinc_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nctdtrinc","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nctdtrinc_dddxd_vvvxf_descr));

    /* nctdtrit_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nctdtrit","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nctdtrit_fffxf_vvvxf_descr));

    /* nctdtrit_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nctdtrit","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nctdtrit_dddxd_vvvxf_descr));

    /* ndtr_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","ndtr","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ndtr_fxf_vxf_descr));

    /* ndtr_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","ndtr","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ndtr_dxd_vxf_descr));

    /* ndtri_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","ndtri","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ndtri_fxf_vxf_descr));

    /* ndtri_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","ndtri","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&ndtri_dxd_vxf_descr));

    /* nrdtrimn_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nrdtrimn","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nrdtrimn_fffxf_vvvxf_descr));

    /* nrdtrimn_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nrdtrimn","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nrdtrimn_dddxd_vvvxf_descr));

    /* nrdtrisd_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nrdtrisd","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nrdtrisd_fffxf_vvvxf_descr));

    /* nrdtrisd_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","nrdtrisd","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&nrdtrisd_dddxd_vvvxf_descr));

    /* obl_ang1_ffffxff_vvvvxvv */
    keytuple=Py_BuildValue("ss((ssss)(ss))","obl_ang1","vvvv","Float32","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&obl_ang1_ffffxff_vvvvxvv_descr));

    /* obl_ang1_ddddxdd_vvvvxvv */
    keytuple=Py_BuildValue("ss((ssss)(ss))","obl_ang1","vvvv","Float64","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&obl_ang1_ddddxdd_vvvvxvv_descr));

    /* obl_ang1_cv_fffffxff_vvvvvxvv */
    keytuple=Py_BuildValue("ss((sssss)(ss))","obl_ang1_cv","vvvvv","Float32","Float32","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&obl_ang1_cv_fffffxff_vvvvvxvv_descr));

    /* obl_ang1_cv_dddddxdd_vvvvvxvv */
    keytuple=Py_BuildValue("ss((sssss)(ss))","obl_ang1_cv","vvvvv","Float64","Float64","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&obl_ang1_cv_dddddxdd_vvvvvxvv_descr));

    /* obl_cv_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","obl_cv","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&obl_cv_fffxf_vvvxf_descr));

    /* obl_cv_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","obl_cv","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&obl_cv_dddxd_vvvxf_descr));

    /* obl_rad1_ffffxff_vvvvxvv */
    keytuple=Py_BuildValue("ss((ssss)(ss))","obl_rad1","vvvv","Float32","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&obl_rad1_ffffxff_vvvvxvv_descr));

    /* obl_rad1_ddddxdd_vvvvxvv */
    keytuple=Py_BuildValue("ss((ssss)(ss))","obl_rad1","vvvv","Float64","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&obl_rad1_ddddxdd_vvvvxvv_descr));

    /* obl_rad1_cv_fffffxff_vvvvvxvv */
    keytuple=Py_BuildValue("ss((sssss)(ss))","obl_rad1_cv","vvvvv","Float32","Float32","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&obl_rad1_cv_fffffxff_vvvvvxvv_descr));

    /* obl_rad1_cv_dddddxdd_vvvvvxvv */
    keytuple=Py_BuildValue("ss((sssss)(ss))","obl_rad1_cv","vvvvv","Float64","Float64","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&obl_rad1_cv_dddddxdd_vvvvvxvv_descr));

    /* obl_rad2_ffffxff_vvvvxvv */
    keytuple=Py_BuildValue("ss((ssss)(ss))","obl_rad2","vvvv","Float32","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&obl_rad2_ffffxff_vvvvxvv_descr));

    /* obl_rad2_ddddxdd_vvvvxvv */
    keytuple=Py_BuildValue("ss((ssss)(ss))","obl_rad2","vvvv","Float64","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&obl_rad2_ddddxdd_vvvvxvv_descr));

    /* obl_rad2_cv_fffffxff_vvvvvxvv */
    keytuple=Py_BuildValue("ss((sssss)(ss))","obl_rad2_cv","vvvvv","Float32","Float32","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&obl_rad2_cv_fffffxff_vvvvvxvv_descr));

    /* obl_rad2_cv_dddddxdd_vvvvvxvv */
    keytuple=Py_BuildValue("ss((sssss)(ss))","obl_rad2_cv","vvvvv","Float64","Float64","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&obl_rad2_cv_dddddxdd_vvvvvxvv_descr));

    /* pbdv_ffxff_vvxvv */
    keytuple=Py_BuildValue("ss((ss)(ss))","pbdv","vv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pbdv_ffxff_vvxvv_descr));

    /* pbdv_ddxdd_vvxvv */
    keytuple=Py_BuildValue("ss((ss)(ss))","pbdv","vv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pbdv_ddxdd_vvxvv_descr));

    /* pbvv_ffxff_vvxvv */
    keytuple=Py_BuildValue("ss((ss)(ss))","pbvv","vv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pbvv_ffxff_vvxvv_descr));

    /* pbvv_ddxdd_vvxvv */
    keytuple=Py_BuildValue("ss((ss)(ss))","pbvv","vv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pbvv_ddxdd_vvxvv_descr));

    /* pbwa_ffxff_vvxvv */
    keytuple=Py_BuildValue("ss((ss)(ss))","pbwa","vv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pbwa_ffxff_vvxvv_descr));

    /* pbwa_ddxdd_vvxvv */
    keytuple=Py_BuildValue("ss((ss)(ss))","pbwa","vv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pbwa_ddxdd_vvxvv_descr));

    /* pdtr_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","pdtr","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pdtr_ffxf_vvxf_descr));

    /* pdtr_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","pdtr","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pdtr_ddxd_vvxf_descr));

    /* pdtrc_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","pdtrc","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pdtrc_ffxf_vvxf_descr));

    /* pdtrc_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","pdtrc","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pdtrc_ddxd_vvxf_descr));

    /* pdtri_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","pdtri","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pdtri_ffxf_vvxf_descr));

    /* pdtri_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","pdtri","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pdtri_ddxd_vvxf_descr));

    /* pdtrik_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","pdtrik","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pdtrik_ffxf_vvxf_descr));

    /* pdtrik_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","pdtrik","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pdtrik_ddxd_vvxf_descr));

    /* pro_ang1_ffffxff_vvvvxvv */
    keytuple=Py_BuildValue("ss((ssss)(ss))","pro_ang1","vvvv","Float32","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pro_ang1_ffffxff_vvvvxvv_descr));

    /* pro_ang1_ddddxdd_vvvvxvv */
    keytuple=Py_BuildValue("ss((ssss)(ss))","pro_ang1","vvvv","Float64","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pro_ang1_ddddxdd_vvvvxvv_descr));

    /* pro_ang1_cv_fffffxff_vvvvvxvv */
    keytuple=Py_BuildValue("ss((sssss)(ss))","pro_ang1_cv","vvvvv","Float32","Float32","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pro_ang1_cv_fffffxff_vvvvvxvv_descr));

    /* pro_ang1_cv_dddddxdd_vvvvvxvv */
    keytuple=Py_BuildValue("ss((sssss)(ss))","pro_ang1_cv","vvvvv","Float64","Float64","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pro_ang1_cv_dddddxdd_vvvvvxvv_descr));

    /* pro_cv_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","pro_cv","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pro_cv_fffxf_vvvxf_descr));

    /* pro_cv_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","pro_cv","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pro_cv_dddxd_vvvxf_descr));

    /* pro_rad1_ffffxff_vvvvxvv */
    keytuple=Py_BuildValue("ss((ssss)(ss))","pro_rad1","vvvv","Float32","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pro_rad1_ffffxff_vvvvxvv_descr));

    /* pro_rad1_ddddxdd_vvvvxvv */
    keytuple=Py_BuildValue("ss((ssss)(ss))","pro_rad1","vvvv","Float64","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pro_rad1_ddddxdd_vvvvxvv_descr));

    /* pro_rad1_cv_fffffxff_vvvvvxvv */
    keytuple=Py_BuildValue("ss((sssss)(ss))","pro_rad1_cv","vvvvv","Float32","Float32","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pro_rad1_cv_fffffxff_vvvvvxvv_descr));

    /* pro_rad1_cv_dddddxdd_vvvvvxvv */
    keytuple=Py_BuildValue("ss((sssss)(ss))","pro_rad1_cv","vvvvv","Float64","Float64","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pro_rad1_cv_dddddxdd_vvvvvxvv_descr));

    /* pro_rad2_ffffxff_vvvvxvv */
    keytuple=Py_BuildValue("ss((ssss)(ss))","pro_rad2","vvvv","Float32","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pro_rad2_ffffxff_vvvvxvv_descr));

    /* pro_rad2_ddddxdd_vvvvxvv */
    keytuple=Py_BuildValue("ss((ssss)(ss))","pro_rad2","vvvv","Float64","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pro_rad2_ddddxdd_vvvvxvv_descr));

    /* pro_rad2_cv_fffffxff_vvvvvxvv */
    keytuple=Py_BuildValue("ss((sssss)(ss))","pro_rad2_cv","vvvvv","Float32","Float32","Float32","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pro_rad2_cv_fffffxff_vvvvvxvv_descr));

    /* pro_rad2_cv_dddddxdd_vvvvvxvv */
    keytuple=Py_BuildValue("ss((sssss)(ss))","pro_rad2_cv","vvvvv","Float64","Float64","Float64","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&pro_rad2_cv_dddddxdd_vvvvvxvv_descr));

    /* psi_FxF_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","psi","v","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&psi_FxF_vxf_descr));

    /* psi_DxD_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","psi","v","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&psi_DxD_vxf_descr));

    /* psi_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","psi","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&psi_fxf_vxf_descr));

    /* psi_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","psi","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&psi_dxd_vxf_descr));

    /* radian_fffxf_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","radian","vvv","Float32","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&radian_fffxf_vvvxf_descr));

    /* radian_dddxd_vvvxf */
    keytuple=Py_BuildValue("ss((sss)(s))","radian","vvv","Float64","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&radian_dddxd_vvvxf_descr));

    /* rgamma_FxF_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","rgamma","v","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&rgamma_FxF_vxf_descr));

    /* rgamma_DxD_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","rgamma","v","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&rgamma_DxD_vxf_descr));

    /* rgamma_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","rgamma","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&rgamma_fxf_vxf_descr));

    /* rgamma_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","rgamma","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&rgamma_dxd_vxf_descr));

    /* round_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","round","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&round_fxf_vxf_descr));

    /* round_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","round","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&round_dxd_vxf_descr));

    /* shichi_fxff_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","shichi","v","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&shichi_fxff_vxvv_descr));

    /* shichi_dxdd_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","shichi","v","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&shichi_dxdd_vxvv_descr));

    /* sici_fxff_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","sici","v","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&sici_fxff_vxvv_descr));

    /* sici_dxdd_vxvv */
    keytuple=Py_BuildValue("ss((s)(ss))","sici","v","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&sici_dxdd_vxvv_descr));

    /* sindg_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","sindg","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&sindg_fxf_vxf_descr));

    /* sindg_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","sindg","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&sindg_dxd_vxf_descr));

    /* smirnov_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","smirnov","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&smirnov_ffxf_vvxf_descr));

    /* smirnov_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","smirnov","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&smirnov_ddxd_vvxf_descr));

    /* smirnovi_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","smirnovi","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&smirnovi_ffxf_vvxf_descr));

    /* smirnovi_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","smirnovi","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&smirnovi_ddxd_vvxf_descr));

    /* spence_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","spence","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&spence_fxf_vxf_descr));

    /* spence_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","spence","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&spence_dxd_vxf_descr));

    /* stdtr_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","stdtr","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&stdtr_ffxf_vvxf_descr));

    /* stdtr_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","stdtr","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&stdtr_ddxd_vvxf_descr));

    /* stdtri_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","stdtri","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&stdtri_ffxf_vvxf_descr));

    /* stdtri_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","stdtri","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&stdtri_ddxd_vvxf_descr));

    /* stdtridf_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","stdtridf","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&stdtridf_ffxf_vvxf_descr));

    /* stdtridf_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","stdtridf","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&stdtridf_ddxd_vvxf_descr));

    /* stdtrit_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","stdtrit","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&stdtrit_ffxf_vvxf_descr));

    /* stdtrit_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","stdtrit","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&stdtrit_ddxd_vvxf_descr));

    /* struve_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","struve","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&struve_ffxf_vvxf_descr));

    /* struve_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","struve","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&struve_ddxd_vvxf_descr));

    /* tandg_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","tandg","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&tandg_fxf_vxf_descr));

    /* tandg_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","tandg","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&tandg_dxd_vxf_descr));

    /* tklmbda_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","tklmbda","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&tklmbda_ffxf_vvxf_descr));

    /* tklmbda_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","tklmbda","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&tklmbda_ddxd_vvxf_descr));

    /* wofz_FxF_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","wofz","v","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&wofz_FxF_vxf_descr));

    /* wofz_DxD_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","wofz","v","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&wofz_DxD_vxf_descr));

    /* y0_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","y0","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&y0_fxf_vxf_descr));

    /* y0_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","y0","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&y0_dxd_vxf_descr));

    /* y1_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","y1","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&y1_fxf_vxf_descr));

    /* y1_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","y1","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&y1_dxd_vxf_descr));

    /* yn_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","yn","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&yn_ffxf_vvxf_descr));

    /* yn_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","yn","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&yn_ddxd_vvxf_descr));

    /* yv_fFxF_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","yv","vv","Float32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&yv_fFxF_vvxf_descr));

    /* yv_dDxD_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","yv","vv","Float64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&yv_dDxD_vvxf_descr));

    /* yv_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","yv","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&yv_ffxf_vvxf_descr));

    /* yv_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","yv","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&yv_ddxd_vvxf_descr));

    /* yve_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","yve","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&yve_ffxf_vvxf_descr));

    /* yve_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","yve","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&yve_ddxd_vvxf_descr));

    /* yve_fFxF_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","yve","vv","Float32","Complex32","Complex32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&yve_fFxF_vvxf_descr));

    /* yve_dDxD_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","yve","vv","Float64","Complex64","Complex64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&yve_dDxD_vvxf_descr));

    /* zeta_ffxf_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","zeta","vv","Float32","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&zeta_ffxf_vvxf_descr));

    /* zeta_ddxd_vvxf */
    keytuple=Py_BuildValue("ss((ss)(s))","zeta","vv","Float64","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&zeta_ddxd_vvxf_descr));

    /* zetac_fxf_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","zetac","v","Float32","Float32");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&zetac_fxf_vxf_descr));

    /* zetac_dxd_vxf */
    keytuple=Py_BuildValue("ss((s)(s))","zetac","v","Float64","Float64");
    PyDict_SetItem(dict,keytuple,
                   (PyObject*)NA_new_cfunc((void*)&zetac_dxd_vxf_descr));

    return dict;
}

/* platform independent*/
#ifdef MS_WIN32
__declspec(dllexport)
#endif
void init_na_cephes(void) {
    PyObject *m, *d;
    m = Py_InitModule("_na_cephes", _na_cephesMethods);
    d = PyModule_GetDict(m);
    import_libnumarray();
    PyDict_SetItemString(d, "functionDict", init_funcDict());
    ADD_VERSION(m);
}
