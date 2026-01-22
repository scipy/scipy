#ifndef DFITPACK_H
#define DFITPACK_H

#include <math.h>

// Public API functions exposed to Python layer

void   bispeu(const double *tx, int nx, const double *ty, int ny, const double *c, int kx, int ky,
              const double *x, const double *y, double *z, int m, double *wrk, int lwrk, int *ier);
void   bispev(const double *tx, int nx, const double *ty, int ny, const double *c, int kx, int ky,
              const double *x, int mx, const double *y, int my, double *z, double *wrk, int lwrk, int *iwrk, int kwrk, int *ier);
void   curfit(const int iopt, const int m, const double *x, const double *y, const double *w,
              const double xb, const double xe, const int k, const double s, const int nest,
              int *n, double *t, double *c, double *fp, double *wrk, const int lwrk,
              int *iwrk, int *ier);
double dblint(double* tx, const int nx, double* ty, const int ny, double* c, const int kx, const int ky, const double xb,
              const double xe, const double yb, const double ye, double* wrk);
void   fpchec(const double* x, const int m, const double* t, const int n, const int k, int* ier);
void   insert(const int iopt, const double* t, const int n, const double* c, const int k, const double x, double* tt, int* nn, double* cc, const int nest, int* ier);
void   parcur(const int iopt, const int ipar, const int idim, const int m, double *u, const int mx, const double *x,
              const double *w, double *ub, double *ue, const int k, const double s, const int nest, int *n, double *t,
              const int nc, double *c, double *fp, double *wrk, const int lwrk, int *iwrk, int *ier);
void   clocur(const int iopt, const int ipar, const int idim, const int m, double *u, const int mx,
              const double *x, const double *w, const int k, const double s, const int nest,
              int *n, double *t, const int nc, double *c, double *fp, double *wrk, const int lwrk,
              int *iwrk, int *ier);
void   parder(const double *tx, int nx, const double *ty, int ny, double *c, int kx, int ky, int nux, int nuy, const double *x, int mx,
              const double *y, int my, double *z, double *wrk, int lwrk, int *iwrk, int kwrk, int *ier);
void   pardeu(const double *tx, int nx, const double *ty, int ny, double *c, int kx, int ky, int nux, int nuy,
              const double *x, const double *y, double *z, int m, double *wrk, int lwrk, int *iwrk, int kwrk, int *ier);
void   pardtc(const double* tx, const int nx, const double* ty, const int ny, const double* c, const int kx, const int ky,
              const int nux, const int nuy, double* newc, int* ier);
void   percur(const int iopt, const int m, const double *x, const double *y, const double *w, const int k, const double s,
              const int nest, int *n, double *t, double *c, double *fp, double *wrk, const int lwrk, int *iwrk, int *ier);
void   regrid(const int iopt, const int mx, const double *x, const int my, const double *y, const double *z,
              const double xb, const double xe, const double yb, const double ye, const int kx, const int ky, const double s,
              const int nxest, const int nyest, const int maxit, int *nx, double *tx, int *ny, double *ty, double *c, double *fp,
              double *wrk, const int lwrk, int *iwrk, const int kwrk, int *ier);
void   spalde(const double *t, const int n, const double *c, const int nc, const int k1, const double x, double* d, int* ier);
void   spgrid(const int *iopt, const int *ider, const int mu, const double *u, const int mv, const double *v,
              const double *r, const double r0, const double r1, const double s, const int nuest, const int nvest,
              int *nu, double *tu, int *nv, double *tv, double *c, double *fp, double *wrk, const int lwrk, int *iwrk, const int kwrk, int *ier);
void   sphere(const int iopt, const int m, const double *teta, const double *phi, const double *r, const double *w,
              const double s, const int ntest, const int npest, const double eps, int *nt, double *tt, int *np, double *tp,
              double *c, double *fp, double *wrk1, const int lwrk1, double *wrk2, const int lwrk2, int *iwrk, const int kwrk, int *ier);
void   splder(const double *t, const int n, const double *c, const int nc, const int k, const int nu,
              const double *x, double *y, const int m, const int e, double *wrk, int *ier);
void   splev(const double *t, const int n, const double *c, const int nc, const int k, const double *x, double *y, const int m, const int e, int *ier);
double splint(const double *t, const int n, const double *c, const int nc, const int k, const double a, const double b, double* wrk);
void   sproot(const double *t, const int n, const double *c, const int nc, double *zero, const int mest, int *m, int *ier);
void   surfit(int iopt, int m, double* x, double* y, double* z, double* w, double xb, double xe, double yb, double ye, int kx, int ky, double s,
              int nxest, int nyest, int nmax, double eps, int* nx, double* tx, int* ny, double* ty, double* c, double* fp, double* wrk1, int lwrk1,
              double* wrk2, int lwrk2, int* iwrk, int kwrk, int* ier);

#endif // DFITPACK_H
