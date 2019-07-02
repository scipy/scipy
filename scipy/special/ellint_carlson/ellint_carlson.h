#ifndef _ELLINT_CARLSON_H_INCLUDED
#define _ELLINT_CARLSON_H_INCLUDED


#ifdef __cplusplus
extern "C" {
#endif


typedef enum ELLINT_RETURN_VALUE
{
    ELLINT_STATUS_SUCCESS = 0,
    ELLINT_STATUS_SINGULAR,
    ELLINT_STATUS_UNDERFLOW,
    ELLINT_STATUS_OVERFLOW,
    ELLINT_STATUS_NITER,
    ELLINT_STATUS_PRECLOSS,
    ELLINT_STATUS_NORESULT,
    ELLINT_STATUS_BAD_ARGS,
    ELLINT_STATUS_BAD_RERR,
    ELLINT_STATUS_OTHER,
    ELLINT_STATUS_UNUSED
} EllInt_Status_t;


#include "ellint_compat.h"


extern EllInt_Status_t fellint_RF(double x, double y, double z, double rerr,
                                  double * restrict res);

extern EllInt_Status_t cellint_RF(double_complex x, double_complex y,
				  double_complex z,
				  double rerr, double_complex * restrict res);

extern EllInt_Status_t fellint_RD(double x, double y, double z, double rerr,
                                  double * restrict res);

extern EllInt_Status_t cellint_RD(double_complex x, double_complex y,
                                  double_complex z,
				  double rerr, double_complex * restrict res);

extern EllInt_Status_t fellint_RJ(double x, double y, double z, double p,
                                  double rerr, double * restrict res);

extern EllInt_Status_t cellint_RJ(double_complex x, double_complex y,
                                  double_complex z, double_complex p,
                                  double rerr, double_complex * restrict res);

extern EllInt_Status_t fellint_RC(double x, double y,
                                  double rerr, double * restrict res);

extern EllInt_Status_t cellint_RC(double_complex x, double_complex y,
                                  double rerr, double_complex * restrict res);

extern EllInt_Status_t fellint_RG(double x, double y, double z, double rerr,
                                  double * restrict res);

extern EllInt_Status_t cellint_RG(double_complex x, double_complex y,
                                  double_complex z,
				  double rerr, double_complex * restrict res);


#ifdef __cplusplus
}
#endif


#endif  /* _ELLINT_CARLSON_H_INCLUDED */
