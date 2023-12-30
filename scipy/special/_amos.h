#ifndef _AMOS_H
#define _AMOS_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#include <math.h>
#include <complex.h>

double complex amos_airy(double complex, int, int, int *, int *);
int amos_besh(double complex, double, int, int, int, double complex *, int *);
int amos_besi(double complex, double, int, int, double complex *, int *);
int amos_besj(double complex, double, int, int, double complex *, int *);
int amos_besk(double complex, double, int, int, double complex *, int *);
int amos_besy(double complex, double, int, int, double complex *, int *);
double complex amos_biry(double complex,int, int, int *);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* ifndef */