#ifndef _DD_REAL_H
#error "dd_real.h needs to be included"
#endif

DD_EXTERN_INLINE double quick_two_sum(double a, double b, double *err);
DD_EXTERN_INLINE double quick_two_diff(double a, double b, double *err);
DD_EXTERN_INLINE double two_sum(double a, double b, double *err);
DD_EXTERN_INLINE double two_diff(double a, double b, double *err);
DD_EXTERN_INLINE void two_split(double a, double *hi, double *lo);
DD_EXTERN_INLINE double two_prod(double a, double b, double *err);
DD_EXTERN_INLINE double two_sqr(double a, double *err);
DD_EXTERN_INLINE double two_div(double a, double b, double *err);
DD_EXTERN_INLINE double two_nint(double d);
DD_EXTERN_INLINE double two_aint(double d);
DD_EXTERN_INLINE int two_comp(const double a, const double b);

