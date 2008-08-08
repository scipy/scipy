
#include <stdlib.h>

/*  function prototypes */

double *bspline(double **, double *, int, double *, int, int, int, int, int); 
double bspline_quad(double *, int, int, int, int, int, int);
double *bspline_prod(double *, int, double *, int, int, int, int, int, int);
void bspline_gram(double **, double *, int, int, int, int);
void invband_compute(double **, double *, int, int);

                
double *bspline(double **output, double *x, int nx,
                double *knots, int nknots,
                int m, int d, int lower, int upper){
    
       int nbasis;
       int index, i, j, k;
       double *result, *b, *b0, *b1;
       double *f0, *f1;
       double denom;

       nbasis = upper - lower;

       result = *((double **) output);
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
            b = (double *) malloc(sizeof(*b) * (nbasis+1) * nx);
            bspline(&b, x, nx, knots, nknots, m-1, d-1, lower, upper+1);

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


double bspline_quad(double *knots, int nknots,
                    int m, int l, int r, int dl, int dr)

        /* This is based on scipy.integrate.fixed_quad */

    {
        double *y;
        double qx[18]={-0.9915651684209309, -0.95582394957139838,
                       -0.89260246649755604, -0.80370495897252303, -0.69168704306035333,
                       -0.55977083107394743, -0.41175116146284346, -0.25188622569150576,
                       -0.084775013041735417, 0.084775013041735306, 0.25188622569150554,
                       0.41175116146284246, 0.55977083107394743, 0.69168704306035189,
                       0.80370495897252314, 0.89260246649755637, 0.95582394957139616,
                       0.9915651684209319};
        double qw[18]={0.021616013526480963, 0.049714548894972385,
                       0.076425730254889301, 0.10094204410628659, 0.12255520671147889,
                       0.14064291467065104, 0.15468467512626605, 0.16427648374583206,
                       0.16914238296314324, 0.16914238296314299, 0.16427648374583295,
                       0.1546846751262658, 0.14064291467065093, 0.12255520671147752,
                       0.10094204410628753, 0.076425730254888483, 0.049714548894967854,
                       0.021616013526484387};
        double x[18];
        int nq=18;
        int k, kk;
        int lower, upper;
        double result, a, b, partial;

        result = 0;

        /* TO DO: figure out knot span more efficiently */

        lower = l - m - 1;
        if (lower < 0) { lower = 0;}
        upper = lower + 2 * m + 4;
        if (upper > nknots - 1) { upper = nknots-1; }

        for (k=lower; k<upper; k++) {
            partial = 0.;
            a = knots[k]; b=knots[k+1];
            for (kk=0; kk<nq; kk++) {
               x[kk] = (b - a) * (qx[kk] + 1) / 2. + a;
            }

            y = bspline_prod(x, nq, knots, nknots, m, l, r, dl, dr);

            for (kk=0; kk<nq; kk++) {
                partial += y[kk] * qw[kk];
            }
            free(y); /* bspline_prod malloc's memory, but does not free it */

            result += (b - a) * partial / 2.;

        }

        return(result);

}



double *bspline_prod(double *x, int nx, double *knots, int nknots,
                     int m, int l, int r, int dl, int dr){
    
        double *result, *bl, *br;
        int k;

        if (abs(r - l) <= m) {
            result = (double *) malloc(sizeof(*result) * nx);
            bl = (double *) malloc(sizeof(*bl) * nx);
            br = (double *) malloc(sizeof(*br) * nx);

            bl = bspline(&bl, x, nx, knots, nknots, m, dl, l, l+1);
            br = bspline(&br, x, nx, knots, nknots, m, dr, r, r+1);

            for (k=0; k<nx; k++) {
                result[k] = bl[k] * br[k];
            }
            free(bl); free(br);
        }
        else {
            for (k=0; k<nx; k++) {
                result[k] = 0.;
            }
        }

        return(result);
}

void bspline_gram(double **output, double *knots, int nknots,
                  int m, int dl, int dr){

    /* Presumes that the first m and last m knots are to be ignored, i.e.
    the interior knots are knots[(m+1):-(m+1)] and the boundary knots are
    knots[m] and knots[-m]. In this setting the first basis element of interest
    is the 1st not the 0th. Should maybe be fixed? */

    
        double *result;
        int l, r, i, j;
        int nbasis;

        nbasis = nknots - m;

        result = *((double **) output);
        for (i=0; i<nbasis; i++) {
            for (j=0; j<m; j++) {
                l = i;
                r = l+j;
                *result = bspline_quad(knots, nknots, m, l, r, dl, dr);
                result++;
            }
        }

}


void invband_compute(double **dataptr, double *L, int n, int m) {

        /* Note: m is number of bands not including the diagonal so L is of size (m+1)xn */

        int i,j,k;
        int idx, idy;
        double *data, *odata;
        double diag;

        data = *((double **) dataptr);

        for (i=0; i<n; i++) {
             diag = L[i];
             data[i] = 1.0 / (diag*diag) ;

             for (j=0; j<=m; j++) {
                 L[j*n+i] /= diag;
                 if (j > 0) { data[j*n+i] = 0;}
             }
         }

        for (i=n-1; i>=0; i--) {
             for (j=1; j <= (m<n-1-i ? m:n-1-i); j++) {
                  for (k=1; k<=(n-1-i<m ? n-1-i:m); k++) {
                      idx = (j<k ? k-j:j-k); idy = (j<k ? i+j:i+k);
                      data[j*n+i] -= L[k*n+i] * data[idx*n+idy];
                  }
             }

             for (k=1; k<=(n-1-i<m ? n-1-i:m); k++) {
                  data[i] -= L[k*n+i] * data[k*n+i];
             }
        }

    return;

}



