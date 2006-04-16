#include <math.h>
#include <stdio.h>
#include "Python.h"
#include "pysparse/fortran.h"
#include "pysparse/blas.h"
#define SPMATRIX_UNIQUE_SYMBOL itsolvers_spmatrix
#include "pysparse/spmatrix.h"
#include "pysparse/gmres.h"

#define SpMatrix_PRECON(prec_obj, n, x, y) \
        {if (SpMatrix_Precon((prec_obj),(n),(x),(y))) return -1;}
#define SpMatrix_MATVEC(mat_obj, n, x, m, y) \
        {if (SpMatrix_Matvec((mat_obj), (n), (x), (m), (y))) return -1;}
	

static double InnerProd(int n, double *x, double *y)
{
    double result;

    int one = 1;
    result = F77(ddot)(&n, x, &one, y, &one);

    return result;
}

static void CopyVector(int n, double *x, double *y)
{
    int one = 1;
    F77(dcopy)(&n, x, &one, y, &one);
}

static void ScaleVector(int n, double alpha, double *x)
{
    int one = 1;
    F77(dscal)(&n, &alpha, x, &one);
}

static void Axpy(int n, double alpha, double *x, double *y)
{
    int one = 1;
    F77(daxpy)(&n, &alpha, x, &one, y, &one);
}

/* simulate 2-D arrays at the cost of some arithmetic */
#define V(i) (&V[(i)*n])
#define W(i) (&W[(i)*n])
#define H(i,j) (H[(j)*m1+(i)])

static void
GeneratePlaneRotation(double dx, double dy, double *cs, double *sn)
{
  if (dy == 0.0) {
    *cs = 1.0;
    *sn = 0.0;
  } else if (fabs(dy) > fabs(dx)) {
    double temp = dx / dy;
    *sn = 1.0 / sqrt( 1.0 + temp*temp );
    *cs = temp * *sn;
  } else {
    double temp = dy / dx;
    *cs = 1.0 / sqrt( 1.0 + temp*temp );
    *sn = temp * *cs;
  }
}

static void ApplyPlaneRotation(double *dx, double *dy, double cs, double sn)
{
  double temp  =  cs * *dx + sn * *dy;
  *dy = -sn * *dx + cs * *dy;
  *dx = temp;
}

int Itsolvers_gmres_kernel(int n, double errtol, int it_max,
			    int *it, double *relres, int dim,
			    double *x, double *b, double *work,
			    PyObject *mat_obj,
			    PyObject *prec_obj) {
    int ONE = 1;

    int mype = 1; /* use 0 for monitoring results */
    int iter;
    double rel_resid;

    double *H  = (double *) malloc(dim*(dim+1) * sizeof(double));


    int m1 = dim+1; /* used inside H macro */
    int i, j, k;
    double beta, resid0, n2b;

    double *s  = (double *) malloc((dim+1) * sizeof(double));
    double *cs = (double *) malloc(dim * sizeof(double));
    double *sn = (double *) malloc(dim * sizeof(double));

    double *V  = (double *) malloc(n*(dim+1) * sizeof(double));
    double *W  = (double *) malloc(n*dim * sizeof(double));

   /* Check for all zero right hand side vector => all zero solution */

    n2b = F77(dnrm2)(&n, b, &ONE);/* Norm of rhs vector, b */
    if (n2b == 0.0) {		/* if rhs vector is all zeros */
    for (i = 0; i < n; i ++)	/* then  solution is all zeros */
      x[i] = 0.0;
     			        /* a valid solution has been obtained */
    *relres = 0.0;		/* the relative residual is actually 0/0 */
    *it = 0;			/* no iterations need be performed */
    return(0);
    }

    iter = 0;
    do
    {
        /* compute initial residual and its norm */
	SpMatrix_MATVEC(mat_obj, n, x, n, V(0));         /* V(0) = A*x        */
        Axpy(n, -1.0, b, V(0));                          /* V(0) = V(0) - b   */
        beta = sqrt(InnerProd(n, V(0), V(0)));     /* beta = norm(V(0)) */
        ScaleVector(n, -1.0/beta, V(0));                 /* V(0) = -V(0)/beta */

        /* save very first residual norm */
        if (iter == 0)
            resid0 = beta;

        for (i = 1; i < dim+1; i++)
            s[i] = 0.0;
        s[0] = beta;

        i = -1;
        do
        {
            i++;
            iter++;


	    if (prec_obj){
               SpMatrix_PRECON(prec_obj, n, V(i), W(i));
            } else {
               CopyVector(n, V(i), W(i));
            }
	    SpMatrix_MATVEC(mat_obj, n, W(i), n, V(i+1));

            for (k = 0; k <= i; k++)
            {
                H(k, i) = InnerProd(n, V(i+1), V(k));
                /* V(i+1) -= H(k, i) * V(k); */
                Axpy(n, -H(k,i), V(k), V(i+1));
            }

            H(i+1, i) = sqrt(InnerProd(n, V(i+1), V(i+1)));
            /* V(i+1) = V(i+1) / H(i+1, i) */
            ScaleVector(n, 1.0 / H(i+1, i), V(i+1));

            for (k = 0; k < i; k++)
                ApplyPlaneRotation(&H(k,i), &H(k+1,i), cs[k], sn[k]);

            GeneratePlaneRotation(H(i,i), H(i+1,i), &cs[i], &sn[i]);
            ApplyPlaneRotation(&H(i,i), &H(i+1,i), cs[i], sn[i]);
            ApplyPlaneRotation(&s[i], &s[i+1], cs[i], sn[i]);

            rel_resid = fabs(s[i+1]) / resid0;
            if (mype == 0)
               printf("Iter (%d): rel. resid. norm: %e\n", iter, rel_resid);

            if (rel_resid <= errtol)
                break;
        }
        while (i+1 < dim && iter+1 <= it_max);

        /* solve upper triangular system in place */
        for (j = i; j >= 0; j--)
        {
            s[j] /= H(j,j);
            for (k = j-1; k >= 0; k--)
                s[k] -= H(k,j) * s[j];
        }

        /* update the solution */
        for (j = 0; j <= i; j++)
        {
            /* x = x + s[j] * W(j) */
            Axpy(n, s[j], W(j), x);
        }
    }
    while (rel_resid > errtol && iter+1 <= it_max);

    /* compute exact residual norm reduction */

    SpMatrix_MATVEC(mat_obj, n, x, n, V(0));             /* V(0) = A*x        */
    Axpy(n, -1.0, b, V(0));                             /* V(0) = V(0) - b   */
    beta = sqrt(InnerProd(n, V(0), V(0)));        /* beta = norm(V(0)) */
    rel_resid = beta / resid0;

    if (mype == 0)
        printf("Iter (%d): computed true norm    : %e\n", iter, beta);

    *it = iter;
    *relres = rel_resid;

    free(H);
    free(s);
    free(cs);
    free(sn);
    free(V);
    free(W);
    
    return 0;
}


