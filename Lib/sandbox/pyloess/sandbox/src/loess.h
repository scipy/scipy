/* for the meaning of these fields, see struct.m */
/* longs are used here so that the codes can be called from S */

#define TRUE  1
#define FALSE 0

typedef struct {
    int err_status;
    char *err_msg;
    } loess_errstatus;

typedef struct {
    long    n;
    long    p;
    double  *y;
    double  *x;
    double  *weights;
    } loess_inputs;

typedef struct {
    double span;
    int    degree;
    int    normalize;
    int    parametric[8];
    int    drop_square[8];
    char   *family;
    } loess_model;

typedef struct {
    char   *surface;
    char   *statistics;
    double cell;
    char   *trace_hat;
    int    iterations;
    } loess_control;

typedef struct {
    long   *parameter;
    long   *a;
    double *xi;
    double *vert;
    double *vval;
    } loess_kd_tree;

typedef struct {
    double *fitted_values;
    double *fitted_residuals;
    double enp;
    double s;
    double one_delta;
    double two_delta;
    double *pseudovalues;
    double trace_hat;
    double *diagonal;
    double *robust;
    double *divisor;
    } loess_outputs;

typedef struct {
    loess_inputs inputs;
    loess_model model;
    loess_control control;
    loess_kd_tree kd_tree;
    loess_outputs outputs;
    loess_errstatus status;
} loess;

typedef struct {
    double *fit;
    double *se_fit;
    double residual_scale;
    double df;
    } prediction;

typedef struct {
    double dfn;
    double dfd;
    double F_value;
    double Pr_F;
    } anova_struct;

typedef struct {
    double *fit;
    double *upper;
    double *lower;
    } conf_inv;

