#ifndef _DD_REAL_H
#error "dd_real.h needs to be included"
#endif

/*********** Inline basic functions ************/

DD_EXTERN_INLINE double dd_hi(const double2 a);
DD_EXTERN_INLINE double dd_lo(const double2 a);
DD_EXTERN_INLINE int dd_isnan(const double2 a);
DD_EXTERN_INLINE int dd_isfinite(const double2 a);
DD_EXTERN_INLINE int dd_isinf(const double2 a);
DD_EXTERN_INLINE int dd_is_zero(const double2 a);
DD_EXTERN_INLINE int dd_is_one(const double2 a);
DD_EXTERN_INLINE int dd_is_positive(const double2 a);
DD_EXTERN_INLINE int dd_is_negative(const double2 a);
DD_EXTERN_INLINE double dd_to_double(const double2 a);
DD_EXTERN_INLINE int dd_to_int(const double2 a);

/*********** Equality Comparisons ************/

DD_EXTERN_INLINE int dd_comp(const double2 a, const double2 b);
DD_EXTERN_INLINE int dd_comp_dd_d(const double2 a, double b);
DD_EXTERN_INLINE int dd_comp_d_dd(double a, const double2 b);

/*********** Creation ************/
DD_EXTERN_INLINE double2 dd_create(double hi, double lo);
DD_EXTERN_INLINE double2 dd_zero(void);
DD_EXTERN_INLINE double2 dd_create_d(double h);
DD_EXTERN_INLINE double2 dd_create_i(int h);
DD_EXTERN_INLINE double2 dd_create_dp(const double *d);

/*********** Unary Minus ***********/
DD_EXTERN_INLINE double2 dd_neg(const double2 a);

/*********** Rounding to integer ***********/

DD_EXTERN_INLINE double2 dd_nint(const double2 a);
DD_EXTERN_INLINE double2 dd_floor(const double2 a);
DD_EXTERN_INLINE double2 dd_ceil(const double2 a);
DD_EXTERN_INLINE double2 dd_aint(const double2 a);

DD_EXTERN_INLINE double2 dd_abs(const double2 a);
DD_EXTERN_INLINE double2 dd_fabs(const double2 a);

/*********** Normalizing ***********/

DD_EXTERN_INLINE double2 dd_ldexp(const double2 a, int expt);
DD_EXTERN_INLINE double2 dd_frexp(const double2 a, int *expt);

/*********** Additions/Subtractions ************/
DD_EXTERN_INLINE double2 dd_add_d_d(double a, double b);
DD_EXTERN_INLINE double2 dd_ieee_add(const double2 a, const double2 b);
DD_EXTERN_INLINE double2 dd_sloppy_add(const double2 a, const double2 b);
DD_EXTERN_INLINE double2 dd_add(const double2 a, const double2 b);
DD_EXTERN_INLINE double2 dd_add_dd_d(const double2 a, double b);
DD_EXTERN_INLINE double2 dd_add_d_dd(double a, const double2 b);

DD_EXTERN_INLINE double2 dd_sub_d_d(double a, double b);
DD_EXTERN_INLINE double2 dd_sub(const double2 a, const double2 b);
DD_EXTERN_INLINE double2 dd_sub_dd_d(const double2 a, double b);
DD_EXTERN_INLINE double2 dd_sub_d_dd(double a, const double2 b);

/*********** Multiplications ************/
DD_EXTERN_INLINE double2 dd_mul_d_d(double a, double b);
DD_EXTERN_INLINE double2 dd_mul_pwr2(const double2 a, double b);
DD_EXTERN_INLINE double2 dd_mul_dd_d(const double2 a, double b);
DD_EXTERN_INLINE double2 dd_mul_d_dd(double a, const double2 b);
DD_EXTERN_INLINE double2 dd_mul(const double2 a, const double2 b);

/*********** Divisions ************/

DD_EXTERN_INLINE double2 dd_sloppy_div(const double2 a, const double2 b);
DD_EXTERN_INLINE double2 dd_accurate_div(const double2 a, const double2 b);
DD_EXTERN_INLINE double2 dd_div(const double2 a, const double2 b);
DD_EXTERN_INLINE double2 dd_inv(const double2 a);

DD_EXTERN_INLINE double2 dd_div_d_dd(double a, const double2 b);
DD_EXTERN_INLINE double2 dd_div_dd_d(const double2 a, double b);
DD_EXTERN_INLINE double2 dd_div_d_d(double a, double b);

/********** Remainder **********/
DD_EXTERN_INLINE double2 dd_drem(const double2 a, const double2 b);
DD_EXTERN_INLINE double2 dd_divrem(const double2 a, const double2 b, double2 *r);
DD_EXTERN_INLINE double2 dd_fmod(const double2 a, const double2 b);

/*********** Squaring **********/
DD_EXTERN_INLINE double2 dd_sqr(const double2 a);
DD_EXTERN_INLINE double2 dd_sqr_d(double a);

