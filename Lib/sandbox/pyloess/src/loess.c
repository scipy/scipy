#include "S.h"
#include "loess.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


static char *surf_stat;

int error_status = 0;
char *error_message = NULL;

void
loess_setup(double *x, double *y, long n, long p, loess *lo)
{
    int i, max_kd;

    max_kd = n > 200 ? n : 200;

    lo->inputs.y = (double *) malloc(n * sizeof(double));
    lo->inputs.x = (double *) malloc(n * p * sizeof(double));
    lo->inputs.weights = (double *) malloc(n * sizeof(double));
    for(i = 0; i < (n * p); i++)
        lo->inputs.x[i] = x[i];
    for(i = 0; i < n; i++) {
        lo->inputs.y[i] = y[i];
        lo->inputs.weights[i] = 1;
    }
    lo->inputs.n = n;
    lo->inputs.p = p;
    lo->model.span = 0.75;
    lo->model.degree = 2;
    lo->model.normalize = TRUE;
    for(i = 0; i < 8; i++)
        lo->model.parametric[i] = lo->model.drop_square[i] = FALSE;
    lo->model.family = "gaussian";

    lo->control.surface = "interpolate";
    lo->control.statistics = "approximate";
    lo->control.cell = 0.2;
    lo->control.trace_hat = "wait.to.decide";
    lo->control.iterations = 4;

    lo->outputs.fitted_values = (double *) malloc(n * sizeof(double));
    lo->outputs.fitted_residuals = (double *) malloc(n * sizeof(double));
    lo->outputs.pseudovalues = (double *) malloc(n * sizeof(double));
    lo->outputs.diagonal = (double *) malloc(n * sizeof(double));
    lo->outputs.robust = (double *) malloc(n * sizeof(double));
    lo->outputs.divisor = (double *) malloc(p * sizeof(double));

    lo->kd_tree.parameter = (long *) malloc(7 * sizeof(long));
    lo->kd_tree.a = (long *) malloc(max_kd * sizeof(long));
    lo->kd_tree.xi = (double *) malloc(max_kd * sizeof(double));
    lo->kd_tree.vert = (double *) malloc(p * 2 * sizeof(double));
    lo->kd_tree.vval = (double *) malloc((p + 1) * max_kd * sizeof(double));
}

void
loess_fit(loess *lo)
{
    int    size_info[2], iterations;
    void    loess_();

    size_info[0] = lo->inputs.p;
    size_info[1] = lo->inputs.n;
    
    //Reset the default error status...
    error_status = 0;
    lo->status.err_status = 0;
    lo->status.err_msg = NULL;

    iterations = (!strcmp(lo->model.family, "gaussian")) ? 0 :
    lo->control.iterations;
    if(!strcmp(lo->control.trace_hat, "wait.to.decide")) {
        if(!strcmp(lo->control.surface, "interpolate"))
            lo->control.trace_hat = (lo->inputs.n < 500) ? "exact" : "approximate";
        else
            lo->control.trace_hat = "exact";
        }
    loess_(lo->inputs.y, lo->inputs.x, size_info, lo->inputs.weights,
           &lo->model.span,
           &lo->model.degree,
           lo->model.parametric,
           lo->model.drop_square,
           &lo->model.normalize,
           &lo->control.statistics,
           &lo->control.surface,
           &lo->control.cell,
           &lo->control.trace_hat,
           &iterations,
           lo->outputs.fitted_values,
           lo->outputs.fitted_residuals,
           &lo->outputs.enp,
           &lo->outputs.s,
           &lo->outputs.one_delta,
           &lo->outputs.two_delta,
           lo->outputs.pseudovalues,
           &lo->outputs.trace_hat,
           lo->outputs.diagonal,
           lo->outputs.robust,
           lo->outputs.divisor,
           lo->kd_tree.parameter,
           lo->kd_tree.a,
           lo->kd_tree.xi,
           lo->kd_tree.vert,
           lo->kd_tree.vval);

    if(error_status){
        lo->status.err_status = error_status;
        lo->status.err_msg = error_message;
    }
}

void
loess_(double *y, double *x_, int *size_info, double *weights, double *span,
       int *degree, int *parametric, int *drop_square, int *normalize,
       char **statistics, char **surface, double *cell, char **trace_hat_in,
       int *iterations, double *fitted_values, double *fitted_residuals,
       double *enp, double *s, double *one_delta, double *two_delta,
       double *pseudovalues, double *trace_hat_out, double *diagonal,
       double *robust, double *divisor, int *parameter, int *a, double *xi,
       double *vert, double *vval)
{
    double  *x, *x_tmp, new_cell, trL, delta1, delta2, sum_squares = 0,
            *pseudo_resid, *temp, *xi_tmp, *vert_tmp, *vval_tmp,
            *diag_tmp, trL_tmp = 0, d1_tmp = 0, d2_tmp = 0, sum, mean;
    int    i, j, k, p, N, D, sum_drop_sqr = 0, sum_parametric = 0,
            setLf, nonparametric = 0, *order_parametric,
            *order_drop_sqr, zero = 0, max_kd, *a_tmp, *param_tmp;
    int     cut, comp();
    char    *new_stat, *mess;
    void    condition();

    D = size_info[0];
    N = size_info[1];
    max_kd = (N > 200 ? N : 200);
    *one_delta = *two_delta = *trace_hat_out = 0;

    x = (double *) malloc(D * N * sizeof(double));
    x_tmp = (double *) malloc(D * N * sizeof(double));
    temp = (double *) malloc(N * sizeof(double));
    a_tmp = (int *) malloc(max_kd * sizeof(int));
    xi_tmp = (double *) malloc(max_kd * sizeof(double));
    vert_tmp = (double *) malloc(D * 2 * sizeof(double));
    vval_tmp = (double *) malloc((D + 1) * max_kd * sizeof(double));
    diag_tmp = (double *) malloc(N * sizeof(double));
    param_tmp = (int *) malloc(N * sizeof(int));
    order_parametric = (int *) malloc(D * sizeof(int));
    order_drop_sqr = (int *) malloc(D * sizeof(int));
    if((*iterations) > 0)
        pseudo_resid = (double *) malloc(N * sizeof(double));

    new_cell = (*span) * (*cell);
    for(i = 0; i < N; i++)
        robust[i] = 1;
        for(i = 0; i < (N * D); i++)
            x_tmp[i] = x_[i];
    if((*normalize) && (D > 1)) {
        cut = ceil(0.100000000000000000001 * N);
        for(i = 0; i < D; i++) {
            k = i * N;
            for(j = 0; j < N; j++)
                temp[j] = x_[k + j];
            qsort(temp, N, sizeof(double), comp);
            sum = 0;
            for(j = cut; j <= (N - cut - 1); j++)
                sum = sum + temp[j];
            mean = sum / (N - 2 * cut);
            sum = 0;
            for(j = cut; j <= (N - cut - 1); j++) {
                temp[j] = temp[j] - mean;
                sum = sum + temp[j] * temp[j];
            }
            divisor[i] = sqrt(sum / (N - 2 * cut - 1));
            for(j = 0; j < N; j++) {
                p = k + j;
                x_tmp[p] = x_[p] / divisor[i];
            }
        }
    }
    else
        for(i = 0; i < D; i++) divisor[i] = 1;
    
    j = D - 1;
    for(i = 0; i < D; i++) {
        sum_drop_sqr = sum_drop_sqr + drop_square[i];
        sum_parametric = sum_parametric + parametric[i];
        if(parametric[i])
            order_parametric[j--] = i;
        else
            order_parametric[nonparametric++] = i;
    }
    //Reorder the predictors w/ the non-parametric first
    for(i = 0; i < D; i++) {
        order_drop_sqr[i] = 2 - drop_square[order_parametric[i]];
        k = i * N;
        p = order_parametric[i] * N;
        for(j = 0; j < N; j++)
            x[k + j] = x_tmp[p + j];
    }
    
    // Misc. checks .............................
    if((*degree) == 1 && sum_drop_sqr) {
    	error_status = 1;
    	error_message = "Specified the square of a factor predictor to be "\
               			"dropped when degree = 1";
        return;
    }
    
    if(D == 1 && sum_drop_sqr) {
    	error_status = 1;
        error_message = "Specified the square of a predictor to be dropped "\
                        "with only one numeric predictor";
        return;
    }
    
    if(sum_parametric == D) {
    	error_status = 1;
        error_message = "Specified parametric for all predictors";
        return;
        }
        
    // Start the iterations .....................
    for(j = 0; j <= (*iterations); j++) {
        new_stat = j ? "none" : *statistics;
        for(i = 0; i < N; i++)
            robust[i] = weights[i] * robust[i];
        condition(surface, new_stat, trace_hat_in);
        setLf = !strcmp(surf_stat, "interpolate/exact");
        loess_raw(y, x, weights, robust, &D, &N, span, degree,
                  &nonparametric, order_drop_sqr, &sum_drop_sqr,
                  &new_cell, &surf_stat, fitted_values, parameter, a,
                  xi, vert, vval, diagonal, &trL, &delta1, &delta2,
                  &setLf);
        if(j == 0) {
            *trace_hat_out = trL;
            *one_delta = delta1;
            *two_delta = delta2;
        }
        for(i = 0; i < N; i++){
            fitted_residuals[i] = y[i] - fitted_values[i];
        };
        if(j < (*iterations))
            F77_SUB(lowesw)(fitted_residuals, &N, robust, temp);
    }
    if((*iterations) > 0) {
        F77_SUB(lowesp)(&N, y, fitted_values, weights, robust, temp,
						pseudovalues);
        loess_raw(pseudovalues, x, weights, weights, &D, &N, span,
                  degree, &nonparametric, order_drop_sqr, &sum_drop_sqr,
                  &new_cell, &surf_stat, temp, param_tmp, a_tmp, xi_tmp,
                  vert_tmp, vval_tmp, diag_tmp, &trL_tmp, &d1_tmp, &d2_tmp,
                  &zero);
        for(i = 0; i < N; i++)
            pseudo_resid[i] = pseudovalues[i] - temp[i];
    }
    if((*iterations) == 0)
        for(i = 0; i < N; i++)
            sum_squares = sum_squares + weights[i] *
                    fitted_residuals[i] * fitted_residuals[i];
    else
        for(i = 0; i < N; i++)
            sum_squares = sum_squares + weights[i] *
                    pseudo_resid[i] * pseudo_resid[i];
    *enp = (*one_delta) + 2 * (*trace_hat_out) - N;
    *s = sqrt(sum_squares / (*one_delta));
    
    //Clean the mess and leave ..................
    free(x);
    free(x_tmp);
    free(temp);
    free(xi_tmp);
    free(vert_tmp);
    free(vval_tmp);
    free(diag_tmp);
    free(a_tmp);
    free(param_tmp);
    free(order_parametric);
    free(order_drop_sqr);
        if((*iterations) > 0)
                free(pseudo_resid);
}

void
loess_free_mem(loess *lo)
{
    free(lo->inputs.x);
    free(lo->inputs.y);
    free(lo->inputs.weights);
    free(lo->outputs.fitted_values);
    free(lo->outputs.fitted_residuals);
    free(lo->outputs.pseudovalues);
    free(lo->outputs.diagonal);
    free(lo->outputs.robust);
    free(lo->outputs.divisor);
    free(lo->kd_tree.parameter);
    free(lo->kd_tree.a);
    free(lo->kd_tree.xi);
    free(lo->kd_tree.vert);
    free(lo->kd_tree.vval);
}

void
loess_summary(loess *lo)
{
    printf("Number of Observations         : %d\n", lo->inputs.n);
    printf("Equivalent Number of Parameters: %.1f\n", lo->outputs.enp);
    if(!strcmp(lo->model.family, "gaussian"))
        printf("Residual Standard Error        : ");
    else
        printf("Residual Scale Estimate: ");
    printf("%.4f\n", lo->outputs.s);
}

void
condition(char **surface, char *new_stat, char **trace_hat_in)
{
    if(!strcmp(*surface, "interpolate")) {
        if(!strcmp(new_stat, "none"))
            surf_stat = "interpolate/none";
        else if(!strcmp(new_stat, "exact"))
            surf_stat = "interpolate/exact";
        else if(!strcmp(new_stat, "approximate"))
        {
            if(!strcmp(*trace_hat_in, "approximate"))
                surf_stat = "interpolate/2.approx";
            else if(!strcmp(*trace_hat_in, "exact"))
                surf_stat = "interpolate/1.approx";
        }
    }
    else if(!strcmp(*surface, "direct")) {
        if(!strcmp(new_stat, "none"))
            surf_stat = "direct/none";
        else if(!strcmp(new_stat, "exact"))
            surf_stat = "direct/exact";
        else if(!strcmp(new_stat, "approximate"))
            surf_stat = "direct/approximate";
    }
}

int
comp(double *d1, double *d2)
{
    if(*d1 < *d2)
        return(-1);
    else if(*d1 == *d2)
        return(0);
    else
        return(1);
}
