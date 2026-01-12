#ifndef __MINPACK_H
#define __MINPACK_H


void CHKDER(
  const int m,
  const int n,
  double* x,
  const double* fvec,
  const double* fjac,
  const int ldfjac,
  double* xp,
  const double* fvecp,
  const int mode,
  double* err
);

void HYBRD(
  int(*fcn)(int* n, double* x, double* fvec, int* iflag),
  const int n,
  double* x,
  double* fvec,
  const double xtol,
  const int maxfev,
  const int ml,
  const int mu,
  const double epsfcn,
  double* diag,
  const int mode,
  const double factor,
  const int nprint,
  int* info,
  int* nfev,
  double* fjac,
  const int ldfjac,
  double* r,
  const int lr,
  double* qtf,
  double* wa1,
  double* wa2,
  double* wa3,
  double* wa4
);


void HYBRJ(
  int(*fcn)(int* n, double* x, double* fvec, double* fjac, int* ldfjac, int* iflag),
  const int n,
  double* x,
  double* fvec,
  double* fjac,
  const int ldfjac,
  const double xtol,
  const int maxfev,
  double* diag,
  const int mode,
  const double factor,
  const int nprint,
  int* info,
  int* nfev,
  int* njev,
  double* r,
  const int ldr,
  double* qtf,
  double* wa1,
  double* wa2,
  double* wa3,
  double* wa4
);

void LMDIF(
  int(*fcn)(int* m, int* n, double* x, double* fvec, int* iflag),
  const int m,
  const int n,
  double* x,
  double* fvec,
  const double ftol,
  const double xtol,
  const double gtol,
  const int maxfev,
  const double epsfcn,
  double* diag,
  const int mode,
  const double factor,
  const int nprint,
  int* info,
  int* nfev,
  double* fjac,
  const int ldfjac,
  int* ipvt,
  double* qtf,
  double* wa1,
  double* wa2,
  double* wa3,
  double* wa4
);

void LMDER(
  int(*fcn)(int* m, int* n, double* x, double* fvec, double* fjac, int* ldfjac, int* iflag),
  const int m,
  const int n,
  double* x,
  double* fvec,
  double* fjac,
  const int ldfjac,
  const double ftol,
  const double xtol,
  const double gtol,
  const int maxfev,
  double* diag,
  const int mode,
  const double factor,
  const int nprint,
  int* info,
  int* nfev,
  int* njev,
  int* ipvt,
  double* qtf,
  double* wa1,
  double* wa2,
  double* wa3,
  double* wa4
);

void LMSTR(
  int(*fcn)(int* m, int* n, double* x, double* fvec, double* wa3, int* iflag),
  const int m,
  const int n,
  double* x,
  double* fvec,
  double* fjac,
  const int ldfjac,
  const double ftol,
  const double xtol,
  const double gtol,
  const int maxfev,
  double* diag,
  const int mode,
  const double factor,
  const int nprint,
  int* info,
  int* nfev,
  int* njev,
  int* ipvt,
  double* qtf,
  double* wa1,
  double* wa2,
  double* wa3,
  double* wa4
);

#endif
