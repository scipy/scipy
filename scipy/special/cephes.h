#ifndef CEPHES_H
#define CEPHES_H

#include "cephes/cephes_names.h"

#ifdef __cplusplus
extern "C" {
#endif

extern int airy(double x, double *ai, double *aip, double *bi, double *bip);

extern double bdtrc(double k, int n, double p);
extern double bdtr(double k, int n, double p);
extern double bdtri(double k, int n, double y);

extern double besselpoly(double a, double lambda, double nu);

extern double beta(double a, double b);
extern double lbeta(double a, double b);

extern double btdtr(double a, double b, double x);

extern double cbrt(double x);
extern double chbevl(double x, double array[], int n);
extern double chdtrc(double df, double x);
extern double chdtr(double df, double x);
extern double chdtri(double df, double y);
extern double dawsn(double xx);

extern double ellie(double phi, double m);
extern double ellik(double phi, double m);
extern double ellpe(double x);

extern int ellpj(double u, double m, double *sn, double *cn, double *dn, double *ph);
extern double ellpk(double x);
extern double exp10(double x);
extern double exp2(double x);

extern double expn(int n, double x);

extern double fdtrc(double a, double b, double x);
extern double fdtr(double a, double b, double x);
extern double fdtri(double a, double b, double y);

extern int fresnl(double xxa, double *ssa, double *cca);
extern double Gamma(double x);
extern double lgam(double x);
extern double lgam_sgn(double x, int *sign);
extern double gammasgn(double x);

extern double gdtr(double a, double b, double x);
extern double gdtrc(double a, double b, double x);
extern double gdtri(double a, double b, double y);

extern double hyp2f1(double a, double b, double c, double x);
extern double hyperg(double a, double b, double x);
extern double threef0(double a, double b, double c, double x, double *err);

extern double i0(double x);
extern double i0e(double x);
extern double i1(double x);
extern double i1e(double x);
extern double igamc(double a, double x);
extern double igam(double a, double x);
extern double igam_fac(double a, double x);
extern double igamci(double a, double q);
extern double igami(double a, double p);

extern double incbet(double aa, double bb, double xx);
extern double incbi(double aa, double bb, double yy0);

extern double iv(double v, double x);
extern double j0(double x);
extern double y0(double x);
extern double j1(double x);
extern double y1(double x);

extern double jn(int n, double x);
extern double jv(double n, double x);
extern double k0(double x);
extern double k0e(double x);
extern double k1(double x);
extern double k1e(double x);
extern double kn(int nn, double x);

extern double nbdtrc(int k, int n, double p);
extern double nbdtr(int k, int n, double p);
extern double nbdtri(int k, int n, double p);

extern double ndtr(double a);
extern double log_ndtr(double a);
extern double erfc(double a);
extern double erf(double x);
extern double erfinv(double y);
extern double erfcinv(double y);
extern double ndtri(double y0);

extern double pdtrc(double k, double m);
extern double pdtr(double k, double m);
extern double pdtri(int k, double y);

extern double poch(double x, double m);

extern double psi(double x);

extern double rgamma(double x);
extern double round(double x);

extern int shichi(double x, double *si, double *ci);
extern int sici(double x, double *si, double *ci);

extern double radian(double d, double m, double s);
extern double sindg(double x);
extern double sinpi(double x);
extern double cosdg(double x);
extern double cospi(double x);

extern double spence(double x);

extern double stdtr(int k, double t);
extern double stdtri(int k, double p);

extern double struve_h(double v, double x);
extern double struve_l(double v, double x);
extern double struve_power_series(double v, double x, int is_h, double *err);
extern double struve_asymp_large_z(double v, double z, int is_h, double *err);
extern double struve_bessel_series(double v, double z, int is_h, double *err);

extern double yv(double v, double x);

extern double tandg(double x);
extern double cotdg(double x);

extern double log1p(double x);
extern double log1pmx(double x);
extern double expm1(double x);
extern double cosm1(double x);
extern double lgam1p(double x);

extern double yn(int n, double x);
extern double zeta(double x, double q);
extern double zetac(double x);

extern double smirnov(int n, double d);
extern double smirnovi(int n, double p);
extern double smirnovp(int n, double d);
extern double smirnovc(int n, double d);
extern double smirnovci(int n, double p);
extern double kolmogorov(double x);
extern double kolmogi(double p);
extern double kolmogp(double x);
extern double kolmogc(double x);
extern double kolmogci(double p);

extern double lanczos_sum_expg_scaled(double x);

extern double owens_t(double h, double a);

extern double tukeylambdacdf(double x, double lambda);

#ifdef __cplusplus
}
#endif

#endif /* CEPHES_H */
