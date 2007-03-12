import os
from math import ceil

import numpy
from numpy import bool_, int_, float_

narray = numpy.array
nempty = numpy.empty
nzeros = numpy.zeros

import logging
logging.basicConfig(level=logging.DEBUG)
loesslogging = logging.getLogger('loesss')

import cloess
reload(cloess)

class paramdict(dict):
    def __init__(self, **kwargs):
        dict.__init__(self, **kwargs)
    
    def __getattr__(self, attr):
        _got = getattr(super(paramdict,self), attr, None)
        if _got is None:
            _got = self.get(attr, None)
        return _got

    def __setattr__(self, attr, value):
        if attr in self.keys():
            self[attr] = value
        else:
            dict.__setattr__(self, attr, value)

#-------------------------------------------------------------------------------    
class _loess_control(dict):

    _optiondict = dict(surface=('interpolate', 'direct'),
                       statistics=('approximate', 'exact'),
                       trace_hat=('wait.to.decide', 'exact', 'approximate'),)
    
    def __init__(self, surface="interpolate", statistics="approximate",
                 trace_hat="wait.to.decide", iterations=4, cell=0.2):
        dict.__init__(self, surface=None, statistics=None, 
                      trace_hat=None, iterations=iterations, 
                      cell=cell)
        self.surface = surface.lower()
        self.statistics = statistics.lower()
        self.trace_hat = trace_hat.lower()
    #
    def __getattr__(self, attr):
        _got = getattr(super(dict,self), attr, None)
        if _got is None:
            _got = self.get(attr, None)
        return _got
    #
    def __setattr__(self, attr, value):
        if attr in self.keys():
            self.__setitem__(attr, value)
        else:
            dict.__setattr__(self, attr, value)
    #
    def __setitem__(self, attr, value):
        if attr in _loess_control._optiondict.keys():
            self.__setoption__(attr, value)
        else:
            dict.__setitem__(self, attr, value)
    #
    def __setoption__(self, attr, value):
        msg = "Invalid argument: %s must be in %s (got %s)"
        attrlist = _loess_control._optiondict[attr]
        if (not isinstance(value, str)) or \
            (value.lower() not in attrlist):
            raise ValueError, msg % (attr, attrlist, value)
        dict.__setitem__(self, attr, value.lower())
    #
#    def set_surface_status(self, new_stat):
#        if new_stat.lower() not in ('none','exact','approximate'):
#            raise ValueError,"Invalid value for the new_stat parameter: "\
#                  "should be in ('none','exact','approximate'), got %s" % new_stat
#        if self.surface == "interpolate" and new_stat == "approximate":
#            if self.trace_hat == "approximate":
#                new_stat = "2.approx"
#            else:
#                trace_hat='exact' in that case
#                new_stat = "1.approx"
#        return "/".join([self.surface, new_stat])
    def set_surface_status(self, new_stat):
        if new_stat.lower() not in ('none','exact','approximate'):
            raise ValueError,"Invalid value for the new_stat parameter: "\
                  "should be in ('none','exact','approximate'), got %s" % new_stat
        statusdict = {"interpolate":{"none": 10,
                                     "exact":11,
                                     "1.approx":12,
                                     "2.approx":13 },
                      "direct":{"none" : 20,
                                "exact": 21,
                                "approximate": 22}}
        if self.surface == "interpolate" and new_stat == "approximate":
            if self.trace_hat == "approximate":
                status = 13
            else:
                # trace_hat='exact' in that case
                status = 12
        else:
            status = statusdict[self.surface][new_stat]
        return status
        
#      
#            
#
#class loess(object):
#    """
#    
#:Ivariables:
#    x : ndarray
#        Independent variables
#    y : ndarray
#        Dependent variables
#    weights : ndarray
#    """
#    def __init__(self, x, y, weights=None,
#                 span=0.75, degree=2, normalize=True, family="gaussian",
#                 surface="interpolate", statistics="approximate", 
#                 cell=0.2, trace_hat="wait.to.decide",iterations=4):
#        x = narray(x, copy=False, subok=True, order='F')
#        if x.ndim == 2:
#            (n,d) = x.shape
#        elif x.ndim == 1:
#            (n,d) = (len(x),1)
#        else:
#            raise ValueError,"Argument 'x' should be 2D at most!"
#        max_kd = max(n, 200)
#        #    
#        y = narray(y, copy=False, subok=True, order='F')
#        if weights is None:
#            weights = numpy.ones_like(y)
#        self.inputs = paramdict(x=x, y=y, weights=weights, nobs=n, nvars=d)
#        self.model = paramdict(span=span, degree=degree, 
#                               normalize=normalize, family=family.lower(),
#                               parametric=nzeros((d,), bool_, 'F'),
#                               drop_square=nzeros((d,), bool_, 'F')
#                               )        
#        self.control = _loess_control(surface=surface.lower(), 
#                                      statistics=statistics.lower(),
#                                      trace_hat=trace_hat.lower(),
#                                      iterations=iterations, cell=cell, 
#                                      )
#        self.outputs = paramdict(fitted_values=nempty((n,), float_, 'F'),
#                                 fitted_residuals=nempty((n,), float_, 'F'),
#                                 pseudovalues=nempty((n,), float_, 'F'),
#                                 diagonal=nempty((n,), float_, 'F'),
#                                 robust=numpy.ones((n,), float_, 'F'),
#                                 divisor=nempty((d,), float_, 'F'),
#                                 enp=0, s=0, one_delta=0, two_delta=0,
#                                 trace_hat=0
#                                 )
#        self.kd_tree = paramdict(parameter=nempty((7,), int_, 'F'),
#                                 a=nempty((max_kd,), int_, 'F'),
#                                 xi=nempty((max_kd,), float_, 'F'),
#                                 vert=nempty((d*2,), float_, 'F'),
#                                 vval=nempty(((d+1)*max_kd,), float_, 'F')
#                                 )
#        #
#        if self.model.family == "gaussian":
#            self.control['iterations'] = 0
#        if self.control.trace_hat == "wait.to.decide":
#            if (self.control.surface == "interpolate") and n >= 500:
#                self.control.trace_hat = "approximate"
#            else:
#                self.control.trace_hat = "exact"
#    #......................................................            
#    def fit(self, span = 0.75, degree = 2, parametric = False,
#            drop_square = False, normalize = True,
#            statistics = "approximate", surface = "interpolate",
#            cell = 0.2, iterations = 1, trace_hat = "exact"):
#        
#        # Get input....................
#        inputs = self.inputs
#        (x, y, n, d, weights) = [inputs[k] 
#                                 for k in ('x','y','nobs','nvars', 'weights')]
#
#        max_kd = max(n, 200) 
##        a = max_kd
#        one_delta = two_delta = trace_hat_out = 0
#        
#        # Set temporary ...............
#        kd_tree = self.kd_tree
#        (a, xi, vert, vval, parameter) = [kd_tree[k] 
#                                          for k in ('a', 'xi', 'vert', 'vval', 
#                                                    'parameter')]
#        a_tmp = nempty((max_kd,), int_)
#        xi_tmp = nempty((max_kd,), float_)
#        vert_tmp = nempty((2*d,), float_)
#        vval_tmp = nempty(((d+1)*max_kd,), float_)
#        
#        # Check control ................
#        control = self.control
#        surface = control.surface
#        statistics = control.statistics
#        iterations = control.iterations
#        trace_hat = control.trace_hat
#        # Get model ....................
#        model = self.model
#        family = model.family
#        parametric = model.parametric
#        drop_square = model.drop_square
#        (span, degree, normalize) = [model[k] 
#                                     for k in ('span', 'degree', 'normalize')]                         
#        #
#        outputs = self.outputs 
#        fitted_values = outputs.fitted_values
#        fitted_residuals = outputs.fitted_residuals
#        pseudovalues = outputs.pseudovalues 
#        diagonal = outputs.diagonal
#        robust = outputs.robust
#        (enp, s, one_delta, two_delta) = [outputs[k]
#                                          for k in ('enp','s','one_delta','two_delta')]
#        trace_hat = outputs.trace_hat
##        parameter = 7
#        (d1_tmp, d2_tmp, trL_tmp, zero,) = (0., 0., 0., 0.)
#        (delta1, delta2, trL, trace_hat_in) = (0, 0, 0, 0)
#        temp = nempty((n,), float_)
#        diag_tmp = nempty((n,), float_)
#        param_tmp = nempty((n,), int_)
#        if iterations > 0:
#            pseudo_resid = nempty((n,), float_)
#        #
#        new_cell = span * cell
#        #
#        loesslogging.debug("initial divisor: %s" % self.outputs.divisor)
#        if normalize and d > 1:
#            cut = int(ceil(0.1*n))
#            x_trimmed = numpy.sort(x, axis=0)[cut:-cut]
#            outputs.divisor = x_trimmed.std(axis=0)
#            outputs.divisor *= numpy.sqrt(n/float(n-1))
#            x = x / outputs.divisor
#        else:
#            outputs.divisor = numpy.ones(d, float_)
#        loesslogging.debug("final divisor: %s" % self.outputs.divisor)
#        #
#        sum_drop_sqr = sum(drop_square)
#        #
#        parametric = narray(parametric, copy=True)
#        parametric.resize((d,))
#        sum_parametric = parametric.sum()
#        nonparametric = numpy.logical_not(parametric).sum()
#        order_parametric = numpy.argsort(parametric)
#        #
#        order_drop_sqr = 2 - drop_square[order_parametric]
#        x = x[:,order_parametric].ravel()
#        #
#        if degree == 1 and sum_drop_sqr:
#            msg = "Specified the square of a factor predictor to be dropped"\
#                  " when degree = 1"
#            raise ValueError, msg
#        if d == 1 and sum_drop_sqr:
#            msg = "Specified the square of a predictor to be dropped "\
#                  "with only one numeric predictor"
#            raise ValueError, msg  
#        if sum_parametric == d:
#            raise ValueError, "Specified parametric for all predictors"
#        #
#        new_stat = statistics.lower()
#        loesslogging.debug("start iteration: %s" % new_stat)
#        for j in range(iterations+1):
#            if j > 0:
#                new_stat = "none"
#            robust = weights * robust
#            surf_stat = control.set_surface_status(new_stat)
#            #setLf = (surf_stat.lower() == "interpolate/exact")
#            setLf = int(surf_stat == 11)
#            loesslogging.debug("iteration: %i: %s - %s" % (j, surf_stat, setLf))
#            #
#            (surf_stat, fitted_values, parameter, a, 
#             xi, vert, vval, diagonal, trL, delta1, delta2, 
#             ) = loess_raw(y, x, weights, robust, d, n, span, degree, 
#                      nonparametric, order_drop_sqr, sum_drop_sqr, 
#                      new_cell, surf_stat, fitted_values, parameter, a, 
#                      xi, vert, vval, diagonal, trL, delta1, delta2, 
#                      setLf)
#            #
#            if j == 0:
#                trace_hat_out = trL
#                one_delta = delta1
#                two_delta = delta2
#            fitted_residuals = y - fitted_values
#            if j < iterations:
#                (fitted_residuals, n, robust, temp) = lowesw(fitted_residuals, n, robust, temp)
#        #
#        if (iterations > 0):
#            lowesp(n, y, fitted_values, weights, robust, temp, pseudovalues)
#            (temp, param_tmp, a_tmp, xi_tmp,
#             vert_tmp, vval_tmp, diag_tmp, trL_tmp, d1_tmp, d2_tmp, 
#             ) = loess_raw(pseudovalues, x, weights, weights, d, n, span, 
#                               degree, nonparametric, order_drop_sqr, sum_drop_sqr,
#                               new_cell, surf_stat, temp, param_tmp, a_tmp, xi_tmp,
#                               vert_tmp, vval_tmp, diag_tmp, trL_tmp, d1_tmp, d2_tmp, 
#                               zero)
#            pseudo_resid = pseudovalues - temp
#        #
#        if (iterations == 0):
#            sum_squares = numpy.sum(weights * fitted_residuals**2)
#        else:
#            sum_squares = numpy.sum(weights * pseudo_resid**2)
#        #
#        loesslogging.debug("setLf:%s" % setLf)
#        loesslogging.debug("SSR:%s" % sum_squares)
#        outputs.enp = (one_delta) + 2 * (trace_hat_out) - n;
#        loesslogging.debug("one_delta:%s-trace_out:%s" % (one_delta, trace_hat_out))
#        outputs.s = numpy.sqrt(sum_squares / (one_delta))
#        for attr in ('one_delta','two_delta','diagonal','trace_hat',
#                     'fitted_values','fitted_residuals','pseudovalues'):
#            setattr(outputs,attr,eval(attr))
##        (outputs.one_delta, outputs.two_delta) = (one_delta, two_delta)
##        outputs.diagonal = diagonal
##        outputs.
#        #
#    #......................................................
#    def summary(self):
#        print "Number of Observations         : %d" % self.inputs.nobs
#        print "Equivalent Number of Parameters: %.1f" % self.outputs.enp
#        if self.model.family == "gaussian":
#            print "Residual Standard Error        : %.4f" % self.outputs.s
#        else:
#            print "Residual Scale Estimate        : %.4f" % self.outputs.s
#    #.......................................................
#    def predict(self):
#        outputs = self.outputs
#        self.prediction = paramdict(fit=None,
#                                    se_fit=None,
#                                    residual_scale=outputs.s,
#                                    df=outputs.one_delta**2 / outputs.two_delta,
#                                    )
#        raise NotImplementedError
#    
#    
#
#
#    
#
##    size_info[0] = lo->in.p;
##    size_info[1] = lo->in.n;
##    size_info[2] = m;
##    
##    pred_(lo->in.y, lo->in.x, eval, size_info, &lo->out.s,
##        lo->in.weights,
##        lo->out.robust,
##        &lo->model.span,
##        &lo->model.degree,
##        &lo->model.normalize,
##        lo->model.parametric,
##        lo->model.drop_square,
##        &lo->control.surface,
##        &lo->control.cell,
##        &lo->model.family,
##        lo->kd_tree.parameter,
##        lo->kd_tree.a,
##        lo->kd_tree.xi,
##        lo->kd_tree.vert,
##        lo->kd_tree.vval,
##        lo->out.divisor,
##        &se,
##        pre->fit,
##        pre->se_fit);
##}
#
##void
##pred_(y, x_, new_x, size_info, s, weights, robust, span, degree, 
##    normalize, parametric, drop_square, surface, cell, family,
##    parameter, a, xi, vert, vval, divisor, se, fit, se_fit)
##double  *y, *x_, *new_x, *weights, *robust, *span, *cell, *fit, *s,
##        *xi, *vert, *vval, *divisor, *se_fit;
##long    *size_info, *degree, *normalize, *parametric, *drop_square, 
##        *parameter, *a, *se;
##char    **surface, **family;
##{     
##    double  *x, *x_tmp, *x_evaluate, *L, new_cell, z, tmp, *fit_tmp,
##            *temp, sum, mean;
##    long    N, D, M, sum_drop_sqr = 0, sum_parametric = 0,
##            nonparametric = 0, *order_parametric, *order_drop_sqr;
##    int     i, j, k, p, cut, comp();
##
##    D = size_info[0];
##    N = size_info[1];
##    M = size_info[2];
##
##    x = (double *) malloc(N * D * sizeof(double));
##    x_tmp = (double *) malloc(N * D * sizeof(double));
##    x_evaluate = (double *) malloc(M * D * sizeof(double));
##    L = (double *) malloc(N * M * sizeof(double));
##    order_parametric = (long *) malloc(D * sizeof(long));
##    order_drop_sqr = (long *) malloc(D * sizeof(long));
##    temp = (double *) malloc(N * D * sizeof(double));
##
##    for(i = 0; i < (N * D); i++)
##        x_tmp[i] = x_[i];
##    for(i = 0; i < D; i++) {
##        k = i * M;
##        for(j = 0; j < M; j++) {
##            p = k + j;
##            new_x[p] = new_x[p] / divisor[i];
##        }
##    }
##    if(!strcmp(*surface, "direct") || se) {
##        for(i = 0; i < D; i++) {
##            k = i * N;
##            for(j = 0; j < N; j++) {
##                p = k + j;
##                x_tmp[p] = x_[p] / divisor[i];
##                }
##        }
##    }
##    j = D - 1;
##    for(i = 0; i < D; i++) {
##            sum_drop_sqr = sum_drop_sqr + drop_square[i];
##            sum_parametric = sum_parametric + parametric[i];
##            if(parametric[i])
##                order_parametric[j--] = i;
##        else
##                order_parametric[nonparametric++] = i;
##    }
##    for(i = 0; i < D; i++) {
##        order_drop_sqr[i] = 2 - drop_square[order_parametric[i]];
##        k = i * M;
##        p = order_parametric[i] * M;
##        for(j = 0; j < M; j++)
##            x_evaluate[k + j] = new_x[p + j];
##        k = i * N;
##        p = order_parametric[i] * N;
##        for(j = 0; j < N; j++)
##            x[k + j] = x_tmp[p + j];
##        }
##    for(i = 0; i < N; i++)
##        robust[i] = weights[i] * robust[i];
##
##    if(!strcmp(*surface, "direct")) {
##        if(*se) {
##            loess_dfitse(y, x, x_evaluate, weights, robust,
##            !strcmp(*family, "gaussian"), span, degree,
##                    &nonparametric, order_drop_sqr, &sum_drop_sqr,
##                    &D, &N, &M, fit, L);
##            }
##        else {
##            loess_dfit(y, x, x_evaluate, robust, span, degree,
##                       &nonparametric, order_drop_sqr, &sum_drop_sqr,
##                       &D, &N, &M, fit);
##            }
##        }
##    else {
##        loess_ifit(parameter, a, xi, vert, vval, &M, x_evaluate, fit);
##        if(*se) {
##                new_cell = (*span) * (*cell);
##                fit_tmp = (double *) malloc(M * sizeof(double));
##                loess_ise(y, x, x_evaluate, weights, span, degree,
##                &nonparametric, order_drop_sqr, &sum_drop_sqr,
##                &new_cell, &D, &N, &M, fit_tmp, L);
##                free(fit_tmp);
##                }
##        }
##    if(*se) {
##        for(i = 0; i < N; i++) {
##            k = i * M;
##            for(j = 0; j < M; j++) {
##                p = k + j;
##                L[p] = L[p] / weights[i];
##                L[p] = L[p] * L[p];
##            }
##        }
##        for(i = 0; i < M; i++) {
##            tmp = 0;
##            for(j = 0; j < N; j++)
##                tmp = tmp + L[i + j * M];
##            se_fit[i] = (*s) * sqrt(tmp);
##        }
##    }
##    free(x);
##    free(x_tmp);
##    free(x_evaluate);
##    free(L);
##    free(order_parametric);
##    free(order_drop_sqr);
##    free(temp);
##}
##
##void
##pred_free_mem(pre)
##struct    pred_struct    *pre;
##{
##    free(pre->fit);
##    free(pre->se_fit);
##}



################################################################################
if __name__ == '__main__':
    import numpy as N
    _data = open(os.path.join('examples','sin_data'), 'r')
    _result = open(os.path.join('examples','sin_result'), 'r')
    x = N.arange(1.,121.)
    y = N.concatenate([N.fromiter((float(v) for v in L.rstrip().split()), float_) 
                       for L in _data.readlines()])
    z = N.concatenate([N.fromiter((float(v) for v in L.rstrip().split()), float_) 
                       for L in _result.readlines()])
#    x = N.concatenate([N.fromiter((float(v) for v in L.rstrip().split()), float_) 
#                       for L in open('_data','r').readlines()])
#    x.shape = (-1,2)
#    y = N.concatenate([N.fromiter((float(v) for v in L.rstrip().split()), float_) 
#                       for L in open('_response','r').readlines()])
    tester = cloess.loess(x,y)
    enp_theo = 4.34
    rse_theo = 0.579
    trc_smoother = 4.73
    print "OK"
    tester.fit()
    tester.summary()
    print "Fit OK"
