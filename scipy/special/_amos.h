#ifndef _AMOS_H
#define _AMOS_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#include <math.h>
#include <complex.h>

#ifndef CMPLX
#define CMPLX(x, y) ((double complex)((double)(x) + I * (double)(y)))
#endif /* CMPLX */

int amos_acai(double complex, double, int, int, int, double complex *, double, double, double, double);
int amos_acon(double complex, double, int, int, int, double complex *, double, double, double, double, double);
double complex amos_airy(double complex, int, int, int *, int *);
int amos_asyi(double complex, double, int, int, double complex *, double, double, double, double);
int amos_besh(double complex, double, int, int, int, double complex *, int *);
int amos_besi(double complex, double, int, int, double complex *, int *);
int amos_besj(double complex, double, int, int, double complex *, int *);
int amos_besk(double complex, double, int, int, double complex *, int *);
int amos_besy(double complex, double, int, int, double complex *, int *);
int amos_binu(double complex, double fnu, int, int, double complex *, double, double, double, double, double);
int amos_bknu(double complex, double, int, int, double complex *, double, double, double);
double complex amos_biry(double complex,int, int, int *);
int amos_buni(double complex, double, int, int, double complex *, int, int *, double, double, double, double);
int amos_bunk(double complex, double, int, int, int, double complex *, double, double, double);
double amos_gamln(double);
int amos_kscl(double complex, double, int, double complex *, double complex, double *, double, double);
int amos_mlri(double complex, double, int, int, double complex *, double);
void amos_rati(double complex, double, int, double complex *, double);
int amos_seri(double complex, double, int, int, double complex *, double, double, double);
int amos_s1s2(double complex, double complex *, double complex *, double, double, int *);
int amos_uchk(double complex, double, double);
void amos_unhj(double complex, double, int, double, double complex *, double complex *, double complex *, double complex *, double complex *, double complex *);
void amos_uni1(double complex, double, int, int, double complex *, int *, int *, double, double, double, double);
void amos_uni2(double complex, double, int, int, double complex *, int *, int *, double, double, double, double);
void amos_unik(double complex, double, int, int, double, int *, double complex *, double complex *, double complex *, double complex *, double complex *);
int amos_unk1(double complex, double, int, int, int, double complex *, double, double, double);
int amos_unk2(double complex, double, int, int, int, double complex *, double, double, double);
int amos_uoik(double complex, double, int, int, int, double complex *, double, double, double);
int amos_wrsk(double complex, double, int, int, double complex *, double complex *, double, double, double);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* ifndef */