#ifndef WRIGHT_H
#define WRIGHT_H 
#include<stdlib.h>
#include<complex.h>
#include<math.h>
#include<fenv.h>
#include<float.h>
#define TWOITERTOL DBL_EPSILON
double complex wrightomega(double complex z);
int wrightomega_ext(double complex z,double complex *w,\
                    double complex *e,double complex *r,\
                    double complex *condest);
#endif /* wright.h */
