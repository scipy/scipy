# -*- Mode: Python -*-  

cdef extern from "loess.h":
    ctypedef struct c_loess_inputs "loess_inputs":
        long    n
        long    p
        double  *y
        double  *x
        double  *weights
    ctypedef struct c_loess_model "loess_model":
        double  span
        long    degree
        long    normalize
        long    parametric[8]
        long    drop_square[8]
        char    *family
    ctypedef struct c_loess_control "loess_control":
        char    *surface
        char    *statistics
        double  cell
        char    *trace_hat
        long    iterations
    ctypedef struct c_loess_kd_tree "loess_kd_tree":
        long    *parameter
        long    *a
        double  *xi
        double  *vert
        double  *vval
    ctypedef struct c_loess_outputs "loess_outputs":
        double  *fitted_values
        double  *fitted_residuals
        double  enp
        double  s
        double  one_delta
        double  two_delta
        double  *pseudovalues
        double  trace_hat
        double  *diagonal
        double  *robust
        double  *divisor
    ctypedef struct c_loess "loess":
        c_loess_inputs inputs
        c_loess_model model
        c_loess_control control
        c_loess_kd_tree kd_tree
        c_loess_outputs outputs
    #typedef struct {
    #    double  *fit;
    #    double  *se_fit;
    #    double  residual_scale;
    #    double  df;
    #} predicted;
    #
    #struct anova_struct {
    #    double  dfn;
    #    double  dfd;
    #    double  F_value;
    #    double  Pr_F;
    #};
    #
    #struct ci_struct {
    #    double    *fit;
    #    double    *upper;
    #    double  *lower;
    #};
cdef extern from "loess.h":    
    void loess_setup(double *x, double *y, long n, long p, c_loess *lo)
    void loess_fit(c_loess *lo)
    void loess_(double *y, double *x_, int *size_info, double *weights, 
        double *span, int *degree, int *parametric, int *drop_square, 
        int *normalize, char **statistics, char **surface, double *cell, 
        char **trace_hat_in, int *iterations, double *fitted_values, 
        double *fitted_residuals, double *enp, double *s, double *one_delta, 
        double *two_delta, double *pseudovalues, double *trace_hat_out, 
        double *diagonal, double *robust, double *divisor, long *parameter, 
        long *a, double *xi, double *vert, double *vval)
    void loess_free_mem(c_loess *lo)
    void loess_summary(c_loess *lo)
    void condition(char **surface, char *new_stat, char **trace_hat_in)
    int comp(double *d1, double *d2)

    void loess_raw(double *y, double *x, double *weights, double *robust, int *d,
         int*n, double *span, int *degree, int *nonparametric,
         int *drop_square, int *sum_drop_sqr, double *cell, char **surf_stat,
         double *surface, long *parameter, long *a, double *xi, double *vert,
         double *vval, double *diagonal, double *trL, double *one_delta,
         double *two_delta, int *setLf)
