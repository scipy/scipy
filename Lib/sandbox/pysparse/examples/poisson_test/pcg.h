#ifndef PCG_H
#define PCG_H

void pcg(int n, 
	 double *x, 
	 double *b,
	 double tol, 
	 int maxit,
	 int clvl,
	 int *iter, 
	 double *relres, 
	 int *flag,
	 double *work,
	 void (*matvec)(double *, double *),
	 void (*precon)(double *, double *));

void F77(pcg_f77)(int *n, 
		  double *x, 
		  double *b,
		  double *tol, 
		  int *maxit,
		  int *clvl,
		  int *iter, 
		  double *relres, 
		  int *flag,
		  double *work,
		  void (*matvec)(double *, double *),
		  void (*precon)(double *, double *));

#endif
