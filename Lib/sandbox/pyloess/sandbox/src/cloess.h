#include <stdio.h>
#include <stdlib.h>


// from loess.c
void loess_setup(double *x, double *y, int n, int p, loess *lo);
void loess_fit(loess *lo);
void loess_free_mem(loess *lo);
void loess_summary(loess *lo);
               
// from misc.c
void anova(loess *one, loess *two, anova_struct *out);
void pointwise(prediction *pre, int m, double coverage, conf_inv *ci);     
void pw_free_mem(conf_inv *ci);
////
// from predict.c
void predict(double *eval, int m, loess *lo, prediction *pre, int se);
void pred_free_mem(prediction *pre);
