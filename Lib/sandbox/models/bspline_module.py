import numpy as N
from scipy.weave import ext_tools
import scipy.special.orthogonal

def build_bspline_module():
    """ Builds an extension module with Bspline basis calculators.
    """

    mod = ext_tools.ext_module('_bspline', compiler='gcc')
    knots = N.linspace(0,1,11).astype(N.Float64)
    nknots = knots.shape[0]
    x = N.array([0.4,0.5], N.Float64)
    nx = x.shape[0]
    m = 4
    d = 0
    lower = 0
    upper = 13
    
    # Bspline code in C
    eval_code = '''
    double *bspline(double *x, int nx, 
                    double *knots, int nknots,
                    int m, int d, int lower, int upper)
    {                   
       int nbasis;
       int index, i, j, k;
       double *result, *b, *b0, *b1;
       double *f0, *f1;
       double denom;

       nbasis = upper - lower;

       result = (double *) malloc(sizeof(*result) * nbasis * nx);
       f0 = (double *) malloc(sizeof(*f0) * nx);
       f1 = (double *) malloc(sizeof(*f1) * nx);

       if (m == 1) {
           for(i=0; i<nbasis; i++) {
               index = i + lower;

               if(index < nknots - 1) {
                   if ((knots[index] != knots[index+1]) && (d <= 0)) {
                       for (k=0; k<nx; k++) {

                           *result = (double) (x[k] >= knots[index]) * (x[k] < knots[index+1]);
                           result++;
                       }
                   } 
                   else {
                       for (k=0; k<nx; k++) {
                           *result = 0.;
                           result++;
                       }
                   }
                }
                else {
                   for (k=0; k<nx; k++) {
                       *result = 0.;
                       result++;
                   }
               }
            }
        }    
        else {
            b = bspline(x, nx, knots, nknots, m-1, d-1, lower, upper+1);

            for(i=0; i<nbasis; i++) {
                b0 = b + nx*i;
                b1 = b + nx*(i+1);

                index = i+lower;

                if ((knots[index] != knots[index+m-1]) && (index+m-1 < nknots)) {
                    denom = knots[index+m-1] - knots[index];   
                    if (d <= 0) {
                        for (k=0; k<nx; k++) {
                            f0[k] = (x[k] - knots[index]) / denom;
                        }
                    }
                    else {
                        for (k=0; k<nx; k++) {
                            f0[k] = (m-1) / (knots[index+m-1] - knots[index]);
                        }
                    }
                }
                else {
                    for (k=0; k<nx; k++) {
                        f0[k] = 0.;
                    }
                }

                index = i+lower+1;
                if ((knots[index] != knots[index+m-1]) && (index+m-1 < nknots)) {
                    denom = knots[index+m-1] - knots[index];
                    if (d <= 0) {
                        for (k=0; k<nx; k++) { 
                            f1[k] = (knots[index+m-1] - x[k]) / denom;
                        }
                    }
                    else {
                        for (k=0; k<nx; k++) { 
                            f1[k] = -(m-1) / (knots[index+m-1] - knots[index]);
                        }
                    }
                }  
                else {
                    for (k=0; k<nx; k++) {
                        f1[k] = 0.;
                    }
                }

                for (k=0; k<nx; k++) {
                    *result = f0[k]*(*b0) + f1[k]*(*b1);
                    b0++; b1++; result++;
                }
            }
        free(b);
        }    
        free(f0); free(f1); 
        result = result - nx * nbasis;

        return(result);
    }
    '''

    eval_ext_code = '''

    int dim[2] = {upper-lower, Nx[0]};

    return_val = PyArray_SimpleNewFromData(2, dim, PyArray_DOUBLE,
        (void *) bspline(x, Nx[0], knots, Nknots[0], m, d, lower, upper));

    '''    

    bspline_eval = ext_tools.ext_function('evaluate',
                                          eval_ext_code,
                                          ['x', 'knots',
                                           'm', 'd', 'lower', 'upper'])
    mod.add_function(bspline_eval)
    bspline_eval.customize.add_support_code(eval_code)

    nq = 5
    qx, qw = scipy.special.orthogonal.p_roots(nq)   
    dl = dr = 2

    gram_code = '''

    double *bspline_prod(double *x, int nx, double *knots, int nknots,
                        int m, int l, int r, int dl, int dr) 
    {
        double *result, *bl, *br;
        int k;

        if (fabs(r - l) <= m) {
            result = (double *) malloc(sizeof(*result) * nx);

            bl = bspline(x, nx, knots, nknots, m, dl, l, l+1);
            br = bspline(x, nx, knots, nknots, m, dr, r, r+1);

            for (k=0; k<nx; k++) {
                result[k] = bl[k] * br[k];
            }
        }
        else {
            for (k=0; k<nx; k++) {
                result[k] = 0.;
            }
        }    

        return(result);
    }

    
    double bspline_quad(double *knots, int nknots,
                        int m, int l, int r, int dl, int dr) 
    {
        double *y;
        double qx[%(nq)d]={%(qx)s};
        double qw[%(nq)d]={%(qw)s};
        double x[%(nq)d];
        int nq=%(nq)d;
        int k, kk;
        int lower, upper;
        double result, a, b, partial;

        result = 0;

        if (l <= r) {
            lower = l - m / 2 - 1;
            upper = r + m / 2 + 1;
        }
        else {
            lower = r - m / 2 - 1;
            upper = l + m / 2 + 1;
        }        

        if (lower < m) {
            lower = m - 1;
        }
        if (upper > nknots-m-1) {
            upper = nknots-m;
        }

        for (k=lower; k<=upper; k++) {
            partial = 0.;
            a = knots[k]; b=knots[k+1];
            for (kk=0; kk<nq; kk++) {
               x[kk] = (b - a) * (qx[kk] + 1) / 2. + a;
            }

            y = bspline_prod(x, nq, knots, nknots, m, l, r, dl, dr);
            for (kk=0; kk<nq; kk++) {
                partial += y[kk] * qw[kk];
            }
            result += (b - a) * partial / 2.;
        }
        return(result);
    }    

    double *bspline_gram(double *knots, int nknots,
                        int m, int dl, int dr) 

    /* Presumes that the first m and last m knots are to be ignored, i.e.
    the interior knots are knots[(m+1):-(m+1)] and the boundary knots are
    knots[m] and knots[-m]. In this setting the first basis element of interest
    is the 1st not the 0th. Should maybe be fixed? */

    {
        double *result;
        int l, r, i, j;
        int nbasis, nband;

        nbasis = nknots - m;
        nband = 2*m - 1;

        result = (double *) malloc(sizeof(*result) * nbasis * nband);
        for (i=0; i<nbasis; i++) {
            l = i;
            for (j=0; j<nband; j++) {
                r = l-m+1+j;
                *result = bspline_quad(knots, nknots, m, l, r, dl, dr);
                result++;
            }
        } 
        result -= nbasis * nband;
        return(result);
    }    

    ''' % {'qx':`[q for q in N.real(qx)]`[1:-1], 'qw':`[q for q in qw]`[1:-1], 'nq':nq}

    gram_ext_code = '''

    int dim[2] = {Nknots[0]-m, 2*m-1};

    return_val = PyArray_SimpleNewFromData(2, dim, PyArray_DOUBLE,
        (void *) bspline_gram(knots, Nknots[0], m, dl, dr));

    '''    

    bspline_gram = ext_tools.ext_function('gram',
                                          gram_ext_code,
                                          ['knots',
                                           'm', 'dl', 'dr'])

    bspline_gram.customize.add_support_code(gram_code)
    mod.add_function(bspline_gram)

    mod.compile()

try:
    import _bspline
except ImportError:
    build_bspline_module()
    import _bspline

## if __name__ == '__main__':
##     knots = N.hstack([[0]*3, N.linspace(0,1,11).astype(N.Float64), [1]*3])
##     x = N.array([0.4,0.5])
##     print bspline_ext.bspline_eval(x, knots, 4, 2, 0, 13)

##     knots = N.hstack([[0]*3, N.linspace(0,1,501).astype(N.Float64), [1]*3])
##     nknots = knots.shape[0]
##     x = N.linspace(0,1,1000)
##     m = 4
##     d = 0


##     import time, gc
##     t = 0
##     for i in range(100):
##         lower = i
##         toc = time.time()
##         gc.collect()
##         y = bspline_ext.bspline_eval(x, knots, m, 2, 0, 503)
##         z = bspline_ext.bspline_prod(x, knots, m, 2, 1, 0, 0)
##         tic = time.time()
##         t += tic-toc
##         del(y); del(z)
    
##     print t / 100
