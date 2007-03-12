#include <stdio.h>
#include <stdlib.h>
#include "loess.h"
/*
loess.c
*/
void loess_setup(double *x, double *y, int n, int p, loess *lo);
void loess_fit(loess *lo);
void
loess_(double *y, double *x_, int *size_info, double *weights, double *span,
       int *degree, int *parametric, int *drop_square, int *normalize,
       char **statistics, char **surface, double *cell, char **trace_hat_in,
       int *iterations, double *fitted_values, double *fitted_residuals,
       double *enp, double *s, double *one_delta, double *two_delta,
       double *pseudovalues, double *trace_hat_out, double *diagonal,
       double *robust, double *divisor, int *parameter, int *a, double *xi,
       double *vert, double *vval);
void loess_free_mem(loess *lo);
void loess_summary(loess *lo);
void condition(char **surface, char *new_stat, char **trace_hat_in);
int comp(double *d1, double *d2);

/*
loessc.c
*/
void
loess_raw(double *y, double *x, double *weights, double *robust, int *d,
          int*n, double *span, int *degree, int *nonparametric,
          int *drop_square, int *sum_drop_sqr, double *cell, char **surf_stat,
          double *surface, int *parameter, int *a, double *xi, double *vert,
          double *vval, double *diagonal, double *trL, double *one_delta,
          double *two_delta, int *setLf);
void
loess_dfit(double *y, double *x, double *x_evaluate, double *weights,
           double *span, int *degree, int *nonparametric,
           int *drop_square, int *sum_drop_sqr, int *d, int *n, int *m,
           double *fit);
void
loess_dfitse(double *y, double *x, double *x_evaluate, double *weights,
             double *robust, int *family, double *span, int *degree,
             int *nonparametric, int *drop_square, int *sum_drop_sqr,
             int *d, int *n, int *m, double *fit, double *L);
void
loess_ifit(int *parameter, int *a, double *xi, double *vert, double *vval,
           int *m, double *x_evaluate, double *fit);
void
loess_ise(double *y, double *x, double *x_evaluate, double *weights,
          double *span, int *degree, int *nonparametric, int *drop_square,
          int *sum_drop_sqr, double *cell, int *d, int *n, int *m,
          double *fit, double *L);
void
loess_workspace(int *d, int *n, double *span, int *degree,
                int *nonparametric, int *drop_square, int *sum_drop_sqr,
                int *setLf);
void
loess_prune(int *parameter, int *a, double *xi, double *vert, double *vval);
void
loess_grow(int *parameter, int *a, double *xi, double *vert, double *vval);
void
F77_SUB(lowesw)(double *res, int *n, double *rw, int *pi);
void
F77_SUB(lowesp)(int *n, double *y, double *yhat, double *pwgts,
                double *rwgts, int *pi, double *ytilde);
                
//                
///*
//misc.c
//*/
//void
//anova(struct loess_struct *one, struct loess_struct *two,
//      struct anova_struct *out);
//
///*
//predict.c
//*/
//void
//predict(double *eval, int m, struct loess_struct *lo, struct pred_struct *pre,
//        int se);
//void
//pred_(double *y, double *x_, double *new_x, int *size_info, double *s,
//      double *weights, double *robust, double *span, int *degree,
//      int *normalize, int *parametric, int *drop_square, char **surface,
//      double *cell, char **family, int *parameter, int *a, double *xi,
//      double *vert, double *vval, double *divisor, int *se, double *fit,
//      double *se_fit);
//void
//pred_free_mem(struct pred_struct *pre);
