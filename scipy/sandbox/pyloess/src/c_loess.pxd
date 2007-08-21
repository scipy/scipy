# -*- Mode: Python -*-  

cdef extern from "loess.h":
    ctypedef struct c_loess_errstatus "loess_errstatus":
        int err_status
        char *err_msg
    ctypedef struct c_loess_inputs "loess_inputs":
        long   n
        long   p
        double *y
        double *x
        double *weights
    ctypedef struct c_loess_model "loess_model":
        double span
        int    degree
        int    normalize
        int    parametric[8]
        int    drop_square[8]
        char   *family
    ctypedef struct c_loess_control "loess_control":
        char   *surface
        char   *statistics
        double cell
        char   *trace_hat
        int    iterations
    ctypedef struct c_loess_kd_tree "loess_kd_tree":
        pass
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
        c_loess_errstatus status
    ctypedef struct c_prediction "prediction": 
        double  *fit
        double  *se_fit
        double  residual_scale
        double  df
#    ctypedef struct c_anova "anova_struct":
#        double  dfn
#        double  dfd
#        double  F_value
#        double  Pr_F
    ctypedef struct c_conf_inv "conf_inv":
        double  *fit
        double  *upper
        double  *lower
    
cdef extern from "cloess.h":    
    void loess_setup(double *x, double *y, long n, long p, c_loess *lo)
    void loess_fit(c_loess *lo)
    void loess_free_mem(c_loess *lo)
    void loess_summary(c_loess *lo)
    #
    void c_predict "predict" (double *eval, int m, c_loess *lo, c_prediction *pre, int se) 
    void pred_free_mem(c_prediction *pre)
    #
    void c_pointwise "pointwise" (c_prediction *pre, int m, double coverage, c_conf_inv *ci)
    double pf (double q, double df1, double df2)
    double ibeta (double x, double a, double b)
    void pw_free_mem (c_conf_inv *ci)
    
