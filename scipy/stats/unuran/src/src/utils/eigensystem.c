/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE: eigensystem.c                                                     *
 *                                                                           *
 *   Routines for calculating eigenvalues and eigenvectors of real           *
 *   symmetric matrices.                                                     *
 *                                                                           *
 *   Mainly adapted and modified from FORTRAN77 code POLYEN created by       *
 *   G. Derflinger                                                           *
 *                                                                           *
 *****************************************************************************
 *   FORTRAN version by Gerhard Derflinger.                                  *
 *   Translated into C and adapted for UNU.RAN by Roman Karawatzki.          *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/
					 
#include <unur_source.h>
#include "matrix_source.h"

/* ---------------------------------------------------------------------------- */

static const double EPS = DBL_EPSILON;
static const double EPS2 = (DBL_EPSILON * DBL_EPSILON);

/* ---------------------------------------------------------------------------- */

static int _unur_eigensystem_house ( int dim, double *A, double *d, double *e, double *e2 );
/* householder reduction of a real symmetric matrix A to its tridiagonal form */

static int _unur_eigensystem_newqr(int dim, double *a, double *b, double *b2, double *g);
/* computes the eigenvalues of a real symmetric tri-diagonal matrix */

static int _unur_eigensystem_trinv(int dim, double *a, double *b, double *g, double *c,
                            double *p, double *q, double *r, double *w, double *y,
  	                    int *in);
/* computes the eigenvectors of a real symmetric tri-diagonal matrix */

static int _unur_eigensystem_back(int dim, double *a, double *e, double *c);
/* backtransformation of eigenvectors of tridiagonal matrix to original matrix */

/* ---------------------------------------------------------------------------- */


int _unur_matrix_eigensystem(int dim, const double *M, double *values, double *vectors)
     /*----------------------------------------------------------------------*/
     /* Computes the eigenvalues and the corresponding eigenvectors          */
     /* of a real symmetric matrix M. 					     */
     /* The eigenvectors are normalized and (almost) orthognal.              */
     /* The eigenvectors are stored consecutively in vectors.                */
     /*	                                                                     */
     /* Method:                                                              */
     /*   Householder tri-diagonalization                                    */
     /*   Secant-QR method                                                   */
     /*   Inverse iteration                                                  */
     /*	                                                                     */
     /* Parameters:                                                          */
     /*   dim     : dimension                                                */
     /*   M       : real symmetric dim x dim matrix                          */
     /*	                                                                     */
     /* Output:                                                              */
     /*   values  : eigenvalues of M in ascending order                      */
     /*   vectors : normalized eigenvectors of M stored row-wise             */
     /*	                                                                     */
     /* Return:                                                              */
     /*   UNUR_SUCCESS : eignesystem successfully computed                   */
     /*   UNUR_FAILURE : qr algorithm did not converge                       */
     /*   UNUR_ERR_NULL: invalid NULL pointer occurred                       */
     /*----------------------------------------------------------------------*/
{
  double *A; /* local working copy of M (elements of A will be overwritten) */
  double *diag; 
  double *codiag; 
  double *wk;
  int *in;

  int i;
  int ret = 0; /* UNUR_SUCCESS */
 
  /* Check arguments */
  CHECK_NULL(M,UNUR_ERR_NULL);

  /* Special case when M is one-dimensional */
  if (dim==1) {
    values[0]=M[0];
    vectors[0]=1.;
    return ret;
  }

  /* make a local copy of the matrix M -> A */
  A = _unur_xmalloc(dim*dim*sizeof(double));
  memcpy(A, M, dim*dim*sizeof(double));

  /* alocate working arrays */
  diag = _unur_xmalloc(dim*sizeof(double));   
  codiag = _unur_xmalloc(dim*sizeof(double)); /* stored in 0..(dim-2) */
  wk = _unur_xmalloc((5*dim+2)*sizeof(double)); /*working array */
  in = _unur_xmalloc(dim*sizeof(int)); /*working array */

  /* calculate tridiagonal Householder matrix */
  _unur_eigensystem_house(dim, A, diag, codiag, &wk[0]);
  for (i=1; i<dim; i++) { 
    wk[  dim+i-1] = diag[i-1];
    wk[2*dim+i-1] = codiag[i-1];
    wk[3*dim+i-1] = wk[i-1];
  }
  wk[dim+dim-1]=diag[dim-1];

  /* obtain the eigenvalues */
  ret = _unur_eigensystem_newqr(dim, &wk[dim], &wk[2*dim], &wk[3*dim], values);
  if (ret != UNUR_SUCCESS) {
    /* could not compute eigenvalues */
    goto free_memory;
  }

  /* obtain the eigenvectors of the tri-diagonal matrix */
  _unur_eigensystem_trinv(dim, diag, codiag, values, vectors,
     &wk[0], &wk[dim], &wk[2*dim], &wk[3*dim], &wk[4*dim], in);
  
  /* obtain eigenvectors of original matrix A */
  _unur_eigensystem_back(dim, A, codiag, vectors);
  
free_memory:

  free(A);
  free(diag);
  free(codiag);
  free(wk);
  free(in);

  return ret; 
} /* end of _unur_matrix_eigensystem() */

/* ---------------------------------------------------------------------------- */
int 
_unur_eigensystem_house(int dim, double *A, double *d, double *e, double *e2)
	/*----------------------------------------------------------------------*/
	/* Householder reduction to tri-diagonal form                          	*/
	/* 									*/
	/* Parameters	 							*/
	/*    dim : dimension of the input matrix A				*/
	/*    A	  : real symmetrix dim x dim matrix				*/
	/*    d   : resultant diagonal						*/
	/*    e   : resultant co-diagonal (in e[0]..e[dim-2])			*/
	/*    e2  : squares of e						*/
	/*									*/
	/* This routine is a modified version of the ALGOL procedure TRED1	*/
	/* of Reinsch, Bauer and Wilkinson : Handbook for automatic computation */
	/* vol. 2, Springer 1971, p. 215					*/
	/*----------------------------------------------------------------------*/
{
#define idx1(a,b) (((a)-1)*dim+((b)-1))

  int i,j,k,k1;
  double s,f,g,h,ek,fk,uj; /* parameter for Householder algorithm */
  
  for (k=1; k<=dim; k++) { 
    k1=k+1;
    if (k1<dim) {
      /* obtain parameters for Householder step */
      s=0.;
      for (i=k1; i<=dim; i++) {
        s += A[idx1(i,k)]*A[idx1(i,k)];
      }
      e2[k-1]=s;
      if (_unur_iszero(s)) e[k-1]=0.;
      else {
        f = A[idx1(k1,k)];
        ek = sqrt(s) * ((f <= 0.) ? 1. : -1.);
        A[idx1(k1,k)] = f-ek;
	h = s-f*ek;

	/* perform Householder step */
	fk = 0.;
	for (j=k1; j<=dim; j++) {
          g=0.;
	  for (i=k1;  i<=j;   i++) g += A[idx1(i,k)]*A[idx1(j,i)];
	  for (i=j+1; i<=dim; i++) g += A[idx1(i,k)]*A[idx1(i,j)];
	  g /= h; /* can h be zero ? */
	  e[j-2] = g;
	  fk += g*A[idx1(j,k)];
        }
	fk = fk/(h+h);
	for (j=k1; j<=dim; j++) {
          uj=A[idx1(j,k)];
          s=e[j-2]-2.*fk*uj;
	  for (i=j; i<=dim; i++) A[idx1(i,j)] -= (uj*e[i-2] + A[idx1(i,k)]*s);
        }
	e[k-1]=ek;
      }
    }
    d[k-1]=A[idx1(k,k)];
  } /* next k */
  e[dim-2]=A[idx1(dim,dim-1)];
  e2[dim-2]=e[dim-2]*e[dim-2];

  return UNUR_SUCCESS;

#undef idx1
} /* end of _unur_eigensystem_house() */

/*---------------------------------------------------------------------------*/

int _unur_eigensystem_newqr(int dim, double *a, double *b, double *b2, double *g)
     /*----------------------------------------------------------------------*/
     /* Computes the eigenvalues of a real symmetric tridiagonal matrix      */
     /* using newton iteration and the secant-QR method                      */
     /*                                                                      */
     /* See also G. Derflinger : The NEWQL Algorithm for Obtaining           */
     /*    Eigenvalues of Symmetric Tridiagonal matrices.                    */
     /*    Statistik, Informatik und Oekonomie (hrsg. W. Janko),             */
     /*    Springer 1988, p. 22-33                                           */
     /*                                                                      */
     /* Parameters:                                                          */
     /*   dim ... dimension                                                  */
     /*   a   ... diagonal of trigonal matrix in a[0]..a[dim-1]              */
     /*   b   ... co-diagonal in b[0]..b[dim-2]                              */
     /*   b2  ... squares of b                                               */
     /*   g   ... resultant eigenvalues in ascending order                   */
     /*                                                                      */
     /* Note : the values of a[], b2[] and b[dim-1] are overwritten          */
     /*                                                                      */
     /* Return value:                                                        */
     /*      UNUR_SUCCESS : normal case                                      */
     /*      UNUR_FAILURE : qr algorithm did not converge within MAXIT steps */
     /*----------------------------------------------------------------------*/
{
#define ZERO 0.
#define ONE  1.
#define TWO  2.
#define MAXIT 15

  double shift, crit, crit2, gamm, oldgam;
  double xnull, aa, p, s, d, tanx, sec, cosec, u, w, w1, oldw;

  int i, it, k, l, ll;

  b[dim-1]=ZERO;
  b2[dim-1]=ZERO;

  /* Gershgorin */
  xnull = a[0]-fabs(b[0]);
  for (i=2; i<=dim; i++) xnull = _unur_min( xnull, a[i-1]-fabs(b[i-2])-fabs(b[i-1]) );
  for (i=1; i<=dim; i++) a[i-1] -= xnull;
  
  shift = ZERO; 
  crit = ZERO;
  sec = ZERO;

  for (ll=1; ll<dim; ll++) {
    l = dim+1-ll;
    crit = _unur_max(crit, EPS*(a[l-1]+fabs(b[l-2])) );
    crit2 = crit * crit;
    /* Newton */
    w = ZERO;
    do {
      oldw = w;
      d=a[0]-shift;
      if ( _unur_iszero(d) ) goto label_9; 
      u = ONE / d;
      w = u;
      for (i=2; i<=l; i++) {
        s = b2[i-2]/d;
	d = a[i-1]-shift-s;
        if ( _unur_iszero(d) ) goto label_9;
        u = (u*s+ONE)/d;
	w += u;
      }
      w1 = ONE / w;
      shift += w1;
    } while(w1>crit && w>oldw);
   
label_9:
    g[ll-1]=xnull+shift;
    
    if (ll<dim) { /* here this is always true ... */
      /* seqant-qr algorithm */
      for (it=1; it<=MAXIT; it++) {
        gamm = a[0]-shift;
	p = gamm*gamm;
	k=1;
	for (i=1; i<l; i++) {
	  oldgam = gamm;
	  if ( ! _unur_iszero(p) ) {
            aa = a[i] - shift;
	    tanx = b2[i-1] / p;
	    sec = tanx + ONE;
	    if (tanx > EPS2) {
	      cosec = p / b2[i-1] + ONE;
	      if (sec > TWO) gamm = (aa+oldgam)/sec - oldgam;
	      else gamm = aa - (aa+oldgam)/cosec;
	      p = gamm*gamm*sec;
	      b2[i-1] = (p + b2[i]) / cosec;
	    }
	    else {
	      gamm = aa;
	      p = gamm * gamm;
	      b2[i-1] = ZERO;
	    }
	  }
	  else {
            gamm = -oldgam;
	    if (i==1) sec = ONE;
	    if ( ! _unur_iszero(sec) ) p = b2[i-1] / sec;
	    sec = ZERO;
	    b2[i-1] = p + b2[i];
	  }
          a[i-1] = a[i] - gamm + oldgam;
	  if (a[i-1] <= a[k-1]) k=i;
	}

	a[l-1] = gamm + shift;
	if (a[l-1] <= a[k-1]) k=l;

        /* convergence */
	if (k>1) {
          if (b2[k-2]>crit2) continue;
	  b2[k-2] = ZERO;
	}
	if (k<l) {
          if (b2[k-1]>crit2) continue;
          for (i=k; i<l; i++) {
            a[i-1]=a[i];
	    b2[i-1]=b2[i];
	  }
	}
	
	goto next_eigenvalue;

      } /* next iteration */

      /* here the maximum number of iterations is exceeded */
      /* ll contain the actual index of the eigenvalue */
      return UNUR_FAILURE;
    }
    
next_eigenvalue: ;
  }
  
  g[dim-1] = xnull + a[0];

  return UNUR_SUCCESS;

#undef ZERO
#undef ONE
#undef TWO
#undef MAXIT
} /* end of _unur_eigensystem_newqr() */

/* ---------------------------------------------------------------------------- */

int
_unur_eigensystem_trinv(int dim, double *a, double *b, double *g, double *c,
                        double *p, double *q, double *r, double *w, double *y, int *in)
	/*----------------------------------------------------------------------*/
	/* Computes the eigenvectors to specified eigenvalues of a real 	*/
	/* symmetric tri-diagonal matrix by the method of inverse iteration	*/
	/* The routine is a modified version of the corresponding part of the   */
	/* subroutine GIVHO of J. Ortega in Ralston Wilf. : Mathematical 	*/
	/* methods for digital computers, vol 2, Wiley 1967			*/
	/*									*/
	/* Parameters								*/
	/*   dim : dimension							*/
	/*   a   : diagonal of tridiagonal matrix in a[0]..a[dim-1]		*/
	/*   b   : co-diagonal in b[0]..b[dim-2]				*/
	/*       : with exception of b[dim-1], which is set to 0		*/
	/*	 : a and b remain unchanged					*/
	/*   g   : eigenvalues in g[0]..g[dim-1]				*/
	/*   c   : resultant eigenvectors (stored row-wise)			*/
	/*   p, q, r, w, y, in : working storage (the length of these arrays    */
	/*                       is dim, with exception of y which is dim+2)	*/
	/*									*/
	/*----------------------------------------------------------------------*/
{
#define idx1(a,b) ((a-1)*dim+(b-1))
  int i,j,k,l=0;
  int close;
  
  double epsc = 1.e-5;
  double reps, epsa, f, t, pi, pidim;

  reps = _unur_max(epsc, 16000.*EPS);
  epsa = EPS2*EPS2;
  close = 0;
  b[dim-1]=0.;
  y[dim]=0.;
  y[dim+1]=0.;
  pi=3.141592653589793;
  pidim=pi/(dim+1);

  for (i=dim; i>=1; i--) {
    for (j=1; j<=dim; j++) {
      p[j-1]=0.;
      q[j-1]=b[j-1];
      r[j-1]=a[j-1]-g[i-1];
      y[j-1]=epsa;
    }
    if (i<dim) close = (fabs(g[i]-g[i-1]) < reps);
    if (close) {
      l=l+1;
      for (j=1; j<=dim; j++) {
        y[j-1]=epsa*sin((l+1)*j*pidim);
      }
    }
    else l=0;

    /* Gauss algorithm */

    for (j=1; j<dim; j++) {
      if(fabs(r[j-1]) >= fabs(b[j-1])) {
        if ( _unur_iszero(r[j-1]) ) r[j-1]=EPS2;
        in[j-1]=0;
        f=b[j-1]/r[j-1];
      }
      else {
        in[j-1]=1;
        f=r[j-1]/b[j-1];
        r[j-1]=b[j-1];
        t=r[j];
        r[j]=q[j-1];
        q[j-1]=t;
        p[j-1]=q[j];
        q[j]=0.;
      }
      w[j-1]=f;
      q[j] -= f*p[j-1];
      r[j] -= f*q[j-1];
      if ( _unur_iszero(r[j-1]) ) r[j-1]=EPS2;
    }
    if ( _unur_iszero(r[dim-1]) ) r[dim-1]=EPS2;

    /* two iterations */ 
    for (j=dim; j>=1; j--) {
      y[j-1]=(y[j-1]-y[j]*q[j-1]-y[j+1]*p[j-1])/r[j-1];
      c[idx1(1,j)]=y[j-1];
    }
    for (j=1; j<dim; j++) {
      if (in[j-1]==0) y[j] -= w[j-1]*y[j-1];
      else {
        t=y[j-1];
        y[j-1]=y[j];
        y[j]=t-w[j-1]*y[j];
      }
    }
    t=0.;
    for (j=dim; j>=1; j--) {
	y[j-1]=(y[j-1]-y[j]*q[j-1]-y[j+1]*p[j-1])/r[j-1];
        t += y[j-1]*y[j-1];
    }
    /* Orthogonalize in case of close eigenvalues */
    if (close) {
	  for (k=i+1; k<=i+l; k++) {
        t=0.;
        for (j=1; j<=dim; j++) t += y[j-1]*c[idx1(k,j)];
        for (j=1; j<=dim; j++) y[j-1] -= t*c[idx1(k,j)];
	  }
      t=0.;
      for (j=1; j<=dim; j++) t += y[j-1]*y[j-1];
    }

    /* normalize and store eigenvector */ 
    t=sqrt(t);
    for (j=1; j<=dim; j++) {
      c[idx1(i,j)]=y[j-1]/t;
    }

  }
  
  return UNUR_SUCCESS;

#undef idx1  
} /* end of _unur_eigensystem_trinv() */

/* ---------------------------------------------------------------------------- */

int _unur_eigensystem_back(int dim, double *a, double *e, double *c)
{
	/*----------------------------------------------------------------------*/
	/* Performas a back-transformation on the eigenvectors of the tri-	*/
	/* diagonal matrix in order to obtain the eigenvectors of the original	*/
	/* matrix.								*/
	/*									*/
	/* Parameters								*/
	/*    dim  : dimension							*/
	/*    a, e : the strictly lower triangle of a and e contain the 	*/
	/*	   : necessary information on the householder reduction		*/
	/*	   : (diagonal and co-diagonal elements)			*/
	/*         : a and e are left unaltered by this routine			*/
	/*    c	   : eigenvectors of tridiagonal matrix - overwritten by 	*/
	/*         : resultant eigenvectors					*/
	/*----------------------------------------------------------------------*/
#define idx1(a,b) ((a-1)*dim+(b-1))

  int i,j,k,k1;
  double h, s, s2;

  for (k=dim-2; k>=1; k--) {
    if( ! _unur_iszero(e[k-1]) ) {
      k1=k+1;
      h=-e[k-1]*a[idx1(k1,k)];
      for (j=1; j<=dim; j++) {
        s=0.;
        for (i=k1; i<=dim; i++) s += a[idx1(i,k)]*c[idx1(j,i)];
        s=s/h;
        for (i=k1; i<=dim; i++) c[idx1(j,i)] -= s*a[idx1(i,k)];

        /* non-negative sum of eigenvector components */ 
        s=0.;
        for (i=1; i<=dim; i++) s += c[idx1(j,i)];
        if (s < 0.) {
          for (i=1; i<=dim; i++) c[idx1(j,i)] = -c[idx1(j,i)];
        }			   

        /* eigenvector normalization */ 
        s2=0.;
        for (i=1; i<=dim; i++) s2 += c[idx1(j,i)]*c[idx1(j,i)] ;
        s2 = sqrt(s2);
        for (i=1; i<=dim; i++) c[idx1(j,i)] /= s2 ;

      }
    }
  }

  return UNUR_SUCCESS;

#undef idx1  
} /* end of _unur_eigensystem_back() */

/* ------------------------------------------------------------------------- */
