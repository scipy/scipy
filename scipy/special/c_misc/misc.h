#ifndef C_MISC_MISC_H
#define C_MISC_MISC_H

#if __GNUC__ >= 4
#pragma GCC visibility push(hidden)
#endif

double besselpoly(double a, double lambda, double nu);
double gammasgn(double x);
double poch(double x, double m);

double struve_h(double v, double x);
double struve_l(double v, double x);
double struve_power_series(double v, double x, int is_h, double *err);
double struve_asymp_large_z(double v, double z, int is_h, double *err);
double struve_bessel_series(double v, double z, int is_h, double *err);

#if __GNUC__ >= 4
#pragma GCC visibility pop
#endif

#endif /* C_MISC_MISC_H */
