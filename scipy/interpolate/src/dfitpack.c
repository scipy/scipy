#include "dfitpack.h"

// The following are not yet translated from the original FITPACK Fortran77 code.
// cocosp
// concon
// concur
// cualde
// curev
// evapol
// fourco
// fpadno
// fpadpo
// fpbfout
// fpched
// fpcoco
// fpcons
// fpcosp
// fpcsin
// fpdeno
// fpfrno
// fpgrdi
// fpgrpa
// fpopdi
// fppasu
// fppocu
// fppogr
// fppola
// fprppo
// fpseno
// fpsuev
// fptrnp
// fptrpe
// parsur
// pogrid
// polar
// profil
// surev

static const double PI = 3.1415926535897932384626433832795;

/* Forward declarations for internal (static) functions */

static void   fpader(const double* t, const int n, const double* c, const int k1, const double x, const int l, double* d);
static void   fpback(double* a, double* z, const int n, const int k, double* c, const int nest);
static void   fpbacp(const double *a, const double *b, const double *z, const int n, const int k, double *c, const int k1, const int nest);
static void   fpbisp(const double *tx, const int nx, const double *ty, const int ny, const double *c, const int kx, const int ky, const double *x, const int mx,
                     const double *y, const int my, double *z, double *wx, double *wy, int *lx, int *ly);
static void   fpbspl(const double* t, const int n, const int k, const double x, const int l, double* h);
static void   fpchep(const double* x, const int m, double* t, const int n, const int k, int* ier);
static void   fpclos(const int iopt, const int idim, const int m, const double *u, const int mx,
                     const double *x, const double *w, const int k, const double s, const int nest,
                     const double tol, const int maxit, const int k1, const int k2, int *n, double *t,
                     const int nc, double *c, double *fp, double *fpint, double *z, double *a1,
                     double *a2, double *b, double *g1, double *g2, double *q, int *nrdata, int *ier);
static void   fpcurf(const int iopt, const double *x, const double *y, const double *w, const int m, const double xb, const double xe,
                     const int k, const double s, const int nest, const double tol, const int maxit, const int k1, const int k2,
                     int *n, double *t, double *c, double *fp, double *fpint, double *z, double *a, double *b, double *g,
                     double *q, int *nrdata, int *ier);
static void   fpcuro(const double a, const double b, const double c, const double d, double *x, int *n);
static void   fpcyt1(double *a, const int n, const int nn);
static void   fpcyt2(const double *a, const int n, const double *b, double *c, const int nn);
static void   fpdisc(const double* t, const int n, const int k2, double* b, const int nest);
static void   fpgivs(const double piv, double* ww, double* cos, double* sin);
static void   fpgrre(int ifsx, int ifsy, int ifbx, int ifby, const double *x, const int mx, const double *y, const int my,
                     const double *z, const int mz, const int kx, const int ky, const double *tx, const int nx, const double *ty,
                     const int ny, const double p, double *c, const int nc, double *fp, double *fpx, double *fpy, const int mm,
                     const int mynx, const int kx1, const int kx2, const int ky1, const int ky2, double *spx, double *spy,
                     double *right, double *q, double *ax, double *ay, double *bx, double *by, int *nrx, int *nry);
static void   fpgrsp(int ifsu, int ifsv, int ifbu, int ifbv, int iback, const double *u, const int mu, const double *v,
                     const int mv, const double *r, const int mr, const double *dr, int iop0, int iop1, const double *tu, const int nu,
                     const double *tv, const int nv, const double p, double *c, const int nc, double *sq, double *fp, double *fpu,
                     double *fpv, const int mm, const int mvnu, double *spu, double *spv, double *right, double *q, double *au,
                     double *av1, double *av2, double *bu, double *bv, double *a0, double *a1, double *b0, double *b1, double *c0,
                     double *c1, double *cosi, int *nru, int *nrv);
static void   fpinst(const int iopt, const double* t, const int n, const double* c, const int k, const double x, const int l, double* tt, int* nn, double* cc, const int nest);
static void   fpintb(const double *t, int n, double *bint, const int nk1, const double x, const double y);
static void   fpknot(const double* x, int m, double* t, int* n, double* fpint, int* nrdata, int* nrint, const int nest, const int istart);
static void   fpopsp(const int ifsu, const int ifsv, const int ifbu, const int ifbv, const double *u, const int mu, const double *v, const int mv,
                     const double *r, const int mr, const double r0, const double r1, double *dr, const int *iopt, const int *ider, const double *tu,
                     const int nu, const double *tv, const int nv, const int nuest, const int nvest, const double p, const double *step, double *c,
                     const int nc, double *fp, double *fpu, double *fpv, int *nru, int *nrv, double *wrk, const int lwrk);
static void   fporde(const double* x, const double* y, const int m, const int kx, const int ky, const double* tx, const int nx, const double* ty,
                     const int ny, int* nummer, int* index, const int nreg);
static void   fppara(const int iopt, const int idim, const int m, const double *u, const int mx, const double *x, const double *w,
                     const double ub, const double ue, const int k, const double s, const int nest, const double tol, const int maxit,
                     const int k1, const int k2, int *n, double *t, const int nc, double *c, double *fp, double *fpint, double *z,
                     double *a, double *b, double *g, double *q, int *nrdata, int *ier);
static void   fpperi(const int iopt, const double *x, const double *y, const double *w, const int m,
                     const int k, const double s, const int nest, const double tol, const int maxit,
                     const int k1, const int k2, int *n, double *t, double *c, double *fp, double *fpint,
                     double *z, double *a1, double *a2, double *b, double *g1, double *g2, double *q,
                     int *nrdata, int *ier);
static void   fprank(double* a, double* f, const int n, const int m, const int na, const double tol, double* c, double* sq, int* rank,
                     double* aa, double* ff, double* h);
static double fprati(double* p1, double* f1, double* p2, double* f2, double* p3, double* f3);
static void   fpregr(const int iopt, const double *x, const int mx, const double *y, const int my, const double *z, const int mz,
                     const double xb, const double xe, const double yb, const double ye, const int kx, const int ky, const double s,
                     const int nxest, const int nyest, const double tol, const int maxit, const int nc, int *nx, double *tx, int *ny,
                     double *ty, double *c, double *fp, double *fp0, double *fpold, double *reducx, double *reducy, double *fpintx,
                     double *fpinty, int *lastdi, int *nplusx, int *nplusy, int *nrx, int *nry, int *nrdatx, int *nrdaty,
                     double *wrk, const int lwrk, int *ier);
static void   fprota(const double c, const double s, double* restrict a, double* restrict b);
static void   fprpsp(const int nt, const int np, const double *co, const double *si, double *c, double *f, const int ncoff);
static void   fpsphe(const int iopt, const int m, const double *teta, const double *phi, const double *r, const double *w,
                     const double s, const int ntest, const int npest, const double eta, const double tol, const int maxit,
                     const int ib1, const int ib3, const int nc, const int ncc, const int intest, const int nrest,
                     int *nt, double *tt, int *np, double *tp, double *c, double *fp, double *sup, double *fpint, double *coord,
                     double *f, double *ff, double *row, double *coco, double *cosi, double *a, double *q, double *bt, double *bp,
                     double *spt, double *spp, double *h, int *index, int *nummer, double *wrk, const int lwrk, int *ier);
static void   fpspgr(const int* iopt, const int* ider, const double* u, const int mu, const double* v, const int mv, const double* r, const int mr,
                     const double r0, const double r1, const double s, const int nuest, const int nvest, const double tol, const int maxit,
                     const int nc, int* nu, double* tu, int* nv, double* tv, double* c, double* fp, double* fp0, double* fpold, double* reducu, double* reducv,
                     double* fpintu, double* fpintv, double* dr, double* step, int* lastdi, int* nplusu, int* nplusv, int* lastu0, int* lastu1, int* nru, int* nrv,
                     int* nrdatu, int* nrdatv, double* wrk, const int lwrk, int* ier);
static void   fpsurf(int iopt, int m, double* x, double* y, double* z, double* w, double xb, double xe, double yb, double ye, int kxx, int kyy,
                     double s, int nxest, int nyest, double eta, double tol, int maxit, int nmax, int km1, int km2, int ib1, int ib3, int nc, int intest,
                     int nrest, int* nx0, double* tx, int* ny0, double* ty, double* c, double* fp, double* fp0, double* fpint, double* coord, double* f,
                     double* ff, double* a, double* q, double* bx, double* by, double* spx, double* spy, double* h, int* index, int* nummer, double* wrk,
                     int lwrk, int* ier);
static void   fpsysy(double* restrict a, const int n, double* restrict g);


/****************************************************************************/


void
bispeu(const double *tx, int nx, const double *ty, int ny, const double *c,
       int kx, int ky, const double *x, const double *y, double *z, int m,
       double *wrk, int lwrk, int *ier)
{
    int iwrk[2];
    int lwest = kx + ky + 2;

    // before starting computations a data check is made. if the input data
    // are invalid control is immediately repassed to the calling program.
    *ier = 10;
    if (lwrk < lwest) { return; }
    if (m < 1) { return; }

    *ier = 0;
    for (int i = 0; i < m; i++) {
        fpbisp(tx, nx, ty, ny, c, kx, ky, &x[i], 1, &y[i], 1, &z[i], wrk, &wrk[kx + 1], iwrk, &iwrk[1]);
    }
}


void
bispev(const double *tx, int nx, const double *ty, int ny, const double *c,
       int kx, int ky, const double *x, int mx, const double *y, int my,
       double *z, double *wrk, int lwrk, int *iwrk, int kwrk, int *ier)
{
    int i, iw, lwest;

    // before starting computations a data check is made. if the input data
    // are invalid control is immediately repassed to the calling program.
    *ier = 10;
    lwest = (kx + 1) * mx + (ky + 1) * my;
    if ((lwrk < lwest) || (kwrk < (mx + my)) || (mx < 1) || (my < 1)) { return; }

    if (mx > 1) {
        for (i = 1; i < mx; i++) {
            if (x[i] < x[i - 1]) { return; }
        }
    }

    if (my > 1) {
        for (i = 1; i < my; i++) {
            if (y[i] < y[i - 1]) { return; }
        }
    }

    *ier = 0;
    // Partition the working space and evaluate the bivariate spline
    iw = mx * (kx + 1);
    fpbisp(tx, nx, ty, ny, c, kx, ky, x, mx, y, my, z, wrk, &wrk[iw], iwrk, &iwrk[mx]);
}


void
clocur(const int iopt, const int ipar, const int idim, const int m, double *u, const int mx,
       const double *x, const double *w, const int k, const double s, const int nest,
       int *n, double *t, const int nc, double *c, double *fp, double *wrk, const int lwrk,
       int *iwrk, int *ier)
{
    // subroutine clocur determines a smooth approximating closed spline
    // curve s(u), i.e.
    //     x1 = s1(u)
    //     x2 = s2(u)       u(1) <= u <= u(m)
    //     .........
    //     xidim = sidim(u)
    // with sj(u),j=1,2,...,idim periodic spline functions of degree k with
    // common knots t(j),j=1,2,...,n.
    double per, tol, dist;
    int i, ia1, ia2, ib, ifp, ig1, ig2, iq, iz, i1, i2, j1, j2, k1, k2, lwest;
    int maxit, m1, nmin, ncc, j;

    // we set up the parameters tol and maxit
    maxit = 20;
    tol = 0.001;

    // before starting computations a data check is made. if the input data
    // are invalid, control is immediately repassed to the calling program.
    *ier = 10;
    if ((iopt < -1) || (iopt > 1)) { return; }
    if ((ipar < 0) || (ipar > 1)) { return; }
    if ((idim <= 0) || (idim > 10)) { return; }
    if ((k <= 0) || (k > 5)) { return; }
    k1 = k + 1;
    k2 = k1 + 1;
    nmin = 2 * k1;
    if ((m < 2) || (nest < nmin)) { return; }
    ncc = nest * idim;
    if ((mx < m * idim) || (nc < ncc)) { return; }
    lwest = m * k1 + nest * (7 + idim + 5 * k);
    if (lwrk < lwest) { return; }
    i1 = idim;
    i2 = m * idim;
    for (j = 1; j <= idim; j++) {
        if (x[i1 - 1] != x[i2 - 1]) { return; }
        i1--;
        i2--;
    }

    if ((ipar == 0) && (iopt <= 0)) {
        // calculate parameter values if ipar=0
        i1 = 0;
        i2 = idim;
        u[0] = 0.0;
        for (i = 2; i <= m; i++) {
            dist = 0.0;
            for (j1 = 1; j1 <= idim; j1++) {
                i1 = i1 + 1;
                i2 = i2 + 1;
                dist = dist + (x[i2 - 1] - x[i1 - 1]) * (x[i2 - 1] - x[i1 - 1]);
            }
            u[i - 1] = u[i - 2] + sqrt(dist);
        }
        if (u[m - 1] <= 0.0) { return; }
        for (i = 2; i <= m; i++) {
            u[i - 1] = u[i - 1] / u[m - 1];
        }
        u[m - 1] = 1.0;
    }

    if (w[0] <= 0.0) { return; }
    m1 = m - 1;
    for (i = 1; i <= m1; i++) {
        if ((u[i - 1] >= u[i]) || (w[i - 1] <= 0.0)) { return; }
    }

    if (iopt >= 0) {
        if (s < 0.0) { return; }
        if ((s == 0.0) && (nest < (m + 2 * k))) { return; }
        *ier = 0;
    } else {
        if ((*n <= nmin) || (*n > nest)) { return; }
        per = u[m - 1] - u[0];
        j1 = k1;
        t[j1 - 1] = u[0];
        i1 = *n - k;
        t[i1 - 1] = u[m - 1];
        j2 = j1;
        i2 = i1;
        for (i = 1; i <= k; i++) {
            i1 = i1 + 1;
            i2 = i2 - 1;
            j1 = j1 + 1;
            j2 = j2 - 1;
            t[j2 - 1] = t[i2 - 1] - per;
            t[i1 - 1] = t[j1 - 1] + per;
        }
        fpchep(u, m, t, *n, k, ier);
        if (*ier != 0) { return; }
    }

    // we partition the working space and determine the spline approximation.
    ifp = 0;
    iz = ifp + nest;
    ia1 = iz + ncc;
    ia2 = ia1 + nest * k1;
    ib = ia2 + nest * k;
    ig1 = ib + nest * k2;
    ig2 = ig1 + nest * k2;
    iq = ig2 + nest * k1;

    fpclos(iopt, idim, m, u, mx, x, w, k, s, nest, tol, maxit, k1, k2, n, t,
           ncc, c, fp, &wrk[ifp], &wrk[iz], &wrk[ia1], &wrk[ia2], &wrk[ib],
           &wrk[ig1], &wrk[ig2], &wrk[iq], iwrk, ier);
}


void
curfit(const int iopt, const int m, const double *x, const double *y, const double *w,
       const double xb, const double xe, const int k, const double s, const int nest,
       int *n, double *t, double *c, double *fp, double *wrk, const int lwrk,
       int *iwrk, int *ier)
{
    int i, ia, ib, ifp, ig, iq, iz, j, k1, k2, lwest, maxit, nmin;

    // we set up the parameters tol and maxit
    double tol = 1e-3;
    maxit = 20;

    // before starting computations a data check is made. if the input data
    // are invalid, control is immediately repassed to the calling program.
    *ier = 10;
    if ((k <= 0) || (k > 5)) { return; }
    k1 = k + 1;
    k2 = k1 + 1;
    if ((iopt < -1) || (iopt > 1)) { return; }
    nmin = 2 * k1;
    if ((m < k1) || (nest < nmin)) { return; }
    lwest = m * k1 + nest * (7 + 3 * k);
    if (lwrk < lwest) { return; }
    if ((xb > x[0]) || (xe < x[m - 1])) { return; }
    for (i = 1; i < m; i++) {
        if (x[i - 1] > x[i]) { return; }
    }
    if (iopt >= 0) {
        if (s < 0.0) { return; }
        if ((s == 0.0) && (nest < (m + k1))) { return; }
        *ier = 0;
    } else {
        if ((*n < nmin) || (*n > nest)) { return; }
        j = *n;
        for (i = 0; i < k1; i++) {
            t[i] = xb;
            t[j - 1] = xe;
            j = j - 1;
        }
        fpchec(x, m, t, *n, k, ier);
        if (*ier != 0) { return; }
    }

    // we partition the working space and determine the spline approximation.
    ifp = 0;
    iz = ifp + nest;
    ia = iz + nest;
    ib = ia + nest * k1;
    ig = ib + nest * k2;
    iq = ig + nest * k2;

    fpcurf(iopt, x, y, w, m, xb, xe, k, s, nest, tol, maxit, k1, k2, n, t, c, fp,
           &wrk[ifp], &wrk[iz], &wrk[ia], &wrk[ib], &wrk[ig], &wrk[iq], iwrk, ier);

    return;
}


double
dblint(double* tx, const int nx, double* ty, const int ny, double* c, const int kx, const int ky,
       const double xb, const double xe, const double yb, const double ye, double* wrk)
{
    int nkx1 = nx - kx - 1;
    int nky1 = ny - ky - 1;
    double res = 0.0;
    double result = 0.0;
    // we calculate the integrals of the normalized b-splines ni,kx+1(x)
    fpintb(tx, nx, wrk, nkx1, xb, xe);
    // we calculate the integrals of the normalized b-splines nj,ky+1(y)
    fpintb(ty, ny, &wrk[nkx1], nky1, yb, ye);
    // calculate the integral of s(x,y)
    for (int i = 1; i <= nkx1; i++) {
        res = wrk[i - 1];
        if (res == 0.0) { continue; }
        int m = (i - 1) * nky1;
        int l = nkx1;
        for (int j = 1; j <= nky1; j++) {
            m++;
            l++;
            result = result + res*wrk[l - 1] * c[m - 1];
        }
    }
    return result;
}


static void
fpader(const double* t, const int n, const double* c, const int k1, const double x, const int l, double* d)
{
    // subroutine fpader calculates the derivatives
    //           (j-1)
    //   d(j) = s     (x) , j=1,2,...,k1
    // of a spline of order k1 at the point t(l)<=x<t(l+1), using the
    // stable recurrence scheme of de boor
    (void)n;  // Unused
    double h[20];
    int lk = l - k1;
    for (int i = 1; i <= k1; i++) {
        int ik = i + lk;
        h[i - 1] = c[ik - 1];
    }
    int kj = k1;
    double fac = 1.0;
    for (int j = 1; j <= k1; j++) {
        int ki = kj;
        int j1 = j + 1;

        if (j != 1) {
            int i = k1;
            for (int jj = j; jj <= k1; jj++) {
                int li = i + lk;
                int lj = li + kj;
                h[i - 1] = (h[i - 1] - h[i - 2]) / (t[lj - 1] - t[li - 1]);
                i--;
            }
        }

        for (int i = j; i <= k1; i++) {
            d[i - 1] = h[i - 1];
        }

        if (j != k1) {
            for (int jj = j1; jj <= k1; jj++) {
                ki--;
                int i = k1;
                for (int j2 = jj; j2 <= k1; j2++) {
                    int li = i + lk;
                    int lj = li + ki;
                    d[i - 1] = ((x - t[li - 1]) * d[i - 1] + (t[lj - 1] - x) * d[i - 2]) / (t[lj - 1] - t[li - 1]);
                    i--;
                }
            }
        }
        d[j - 1] = d[k1 - 1] * fac;
        fac *= (double)(k1 - j);
        kj--;
    }

}


static void
fpback(double* a, double* z, const int n, const int k, double* c, const int nest)
{
    c[n - 1] = z[n - 1] / a[n - 1];
    if (n == 1) { return; }

    int k1 = k - 1;
    int i = n - 1;
    for (int j = 2; j <= n; j++) {
        double store = z[i - 1];
        int i1 = k1;
        if (j <= k1) { i1 = j - 1; }
        int m = i;
        for (int l = 1; l <= i1; l++) {
            m++;
            store = store - c[m - 1] * a[(i - 1) + l * nest];
        }
        c[i - 1] = store / a[i - 1];
        i--;
    }
}


static void
fpbacp(const double *a, const double *b, const double *z, const int n, const int k, double *c,
       const int k1, const int nest)
{
    (void)k1;  // Unused
    int n2 = n - k;
    int l = n;
    for (int i = 1; i <= k; i++) {
        double store = z[l - 1];
        int j = k + 2 - i;
        if (i != 1) {
            int l0 = l;
            for (int l1 = j; l1 <= k; l1++) {
                l0++;
                store = store - c[l0 - 1] * b[(l - 1) + (l1 - 1) * nest];
            }
        }
        c[l - 1] = store / b[(l - 1) + (j - 2) * nest];
        l--;
        if (l == 0) return;
    }
    for (int i = 1; i <= n2; i++) {
        double store = z[i - 1];
        int l = n2;
        for (int j = 1; j <= k; j++) {
            l++;
            store = store - c[l - 1] * b[(i - 1) + (j - 1) * nest];
        }
        c[i - 1] = store;
    }
    int i = n2;
    c[i - 1] = c[i - 1] / a[(i - 1)];
    if (i == 1) return;
    for (int j = 2; j <= n2; j++) {
        i--;
        double store = c[i - 1];
        int i1 = (j <= k) ? (j - 1) : k;
        int l = i;
        for (int l0 = 1; l0 <= i1; l0++) {
            l++;
            store = store - c[l - 1] * a[(i - 1) + l0 * nest];
        }
        c[i - 1] = store / a[(i - 1)];
    }
}


static void
fpbisp(const double* tx, const int nx, const double* ty, const int ny, const double* c,
       const int kx, const int ky, const double* x, const int mx, const double* y, const int my,
       double* z, double* wx, double* wy, int* lx, int* ly)
{
    int kx1, ky1, l, l1, l2, nkx1, nky1, i1, m;
    double arg, sp, tb, te, h[6];

    kx1 = kx + 1;
    nkx1 = nx - kx1;
    tb = tx[kx1 - 1];
    te = tx[nkx1];
    l  = kx1;
    l1 = l + 1;

    for (int i = 1; i <= mx; i++) {
        arg = x[i - 1];
        if (arg < tb) { arg = tb; }
        if (arg > te) { arg = te; }

        while (!((arg < tx[l1 - 1]) || (l == nkx1))) {
            l = l1;
            l1 = l + 1;
        }

        fpbspl(tx, nx, kx, arg, l, h);
        lx[i - 1] = l - kx1;
        for (int j = 1; j <= kx1; j++) {
            wx[(i - 1) + (j - 1)*mx] = h[j - 1];
        }
    }

    ky1 = ky + 1;
    nky1 = ny - ky1;
    tb = ty[ky1 - 1];
    te = ty[nky1];
    l  = ky1;
    l1 = l + 1;

    for (int i = 1; i <= my; i++) {
        arg = y[i - 1];
        if (arg < tb) { arg = tb; }
        if (arg > te) { arg = te; }

        while (!((arg < ty[l1 - 1]) || (l == nky1))) {
            l = l1;
            l1 = l + 1;
        }

        fpbspl(ty, ny, ky, arg, l, h);
        ly[i - 1] = l - ky1;
        for (int j = 1; j <= ky1; j++) {
            wy[(i - 1) + (j - 1) * my] = h[j - 1];
        }
    }

    m = 0;
    for (int i = 1; i <= mx; i++) {
        l = lx[i - 1] * nky1;
        for (int i1 = 1; i1 <= kx1; i1++) {
            h[i1 - 1] = wx[(i - 1) + (i1 - 1) * mx];
        }
        for (int j = 1; j <= my; j++) {
            l1 = l + ly[j - 1];
            sp = 0.0;
            for (i1 = 1; i1 <= kx1; i1++) {
                l2 = l1;
                for (int j1 = 1; j1 <= ky1; j1++) {
                    l2++;
                    sp = sp + c[l2 - 1] * h[i1 - 1] * wy[(j - 1) + (j1 - 1) * my];
                }
                l1 += nky1;
            }
            m++;
            z[m - 1] = sp;
        }
    }
}


static void
fpbspl(const double* t, const int n, const int k, const double x, const int l, double* h)
{
    (void)n;  // Unused
    // subroutine fpbspl evaluates the (k+1) non-zero b-splines of degree k
    // at t(l) <= x < t(l+1) using the stable recurrence relation of de boor and cox.
    // Travis Oliphant  2007
    //    changed so that weighting of 0 is used when knots with multiplicity are present.
    //    Also, notice that l+k <= n and 1 <= l+1-k or else the routine will be
    //    accessing memory outside t.  Thus it is imperative that k <= l <= n-k
    //    but this is not checked.
    int li, lj;
    double f;
    double hh[19] = { 0.0 };

    h[0] = 1.0;
    for (int j = 1; j <= k; j++) {
        for (int i = 1; i <= j; i++)
        {
            hh[i - 1] = h[i - 1];
        }

        h[0] = 0.0;
        for (int i = 1; i <= j; i++)
        {
            li = l + i;
            lj  = li - j;
            if (t[li - 1] != t[lj - 1])
            {
                f = hh[i - 1] / (t[li - 1] - t[lj - 1]);
                h[i - 1] = h[i - 1] + f * (t[li - 1] - x);
                h[i] = f * (x - t[lj - 1]);
            }
            else {
                h[i] = 0.0;
            }
        }
    }
}


void
fpchec(const double* x, const int m, const double* t, const int n, const int k, int* ier)
{
    int k1 = k + 1;
    int k2 = k1 + 1;
    int nk1 = n - k1;
    int nk2 = nk1 + 1;
    *ier = 10;

    // check condition no 1
    if ((nk1 < k1) || (nk1 > m)) { *ier = 10; return; }
    // check condition no 2
    int j = n;
    for (int i = 1; i <= k; i++) {
        if (t[i - 1] > t[i]) { *ier = 20; return; }
        if (t[j - 1] < t[j - 2]) { *ier = 20; return; }
        j--;
    }
    // check condition no 3
    for (int i = k2; i <= nk2; i++) {
        if (t[i - 1] <= t[i - 2]) { *ier = 30; return; }
    }
    // check condition no 4
    if ((x[0] < t[k]) || (x[m - 1] > t[nk2 - 1])) { *ier = 40; return; }
    // check condition no 5
    if ((x[0] >= t[k2 - 1]) || (x[m - 1] <= t[nk1 - 1])) { *ier = 50; return; }
    int i = 1;
    int l = k2;
    int nk3 = nk1 - 1;
    if (nk3 < 2) { *ier = 0; return; }
    for (int j = 2; j <= nk3; j++) {
        double tj = t[j - 1];
        l++;
        double tl = t[l - 1];
        do {
            i++;
            if (i >= m) { *ier = 50; return; }
        } while (x[i - 1] <= tj);
        if (x[i - 1] >= tl) { *ier = 50; return; }
    }
    *ier = 0;
    return;
}


static void
fpchep(const double* x, const int m, double* t, const int n, const int k, int* ier)
{
    int k1 = k + 1;
    int k2 = k1 + 1;
    int nk1 = n - k1;
    int nk2 = nk1 + 1;
    int m1 = m - 1;
    *ier = 10;

    // check condition no 1
    if ((nk1 < k1) || (n > (m + 2*k))) { return; }
    // check condition no 2
    int j = n;
    for (int i = 1; i <= k; i++) {
        if (t[i - 1] > t[i]) { return; }
        if (t[j - 1] < t[j - 2]) { return; }
        j--;
    }
    // check condition no 3
    for (int i = k2; i <= nk2; i++) {
        if (t[i - 1] <= t[i - 2]) { return; }
    }
    // check condition no 4
    if ((x[0] < t[k1 - 1]) || (x[m - 1] > t[nk2 - 1])) { return; }
    // check condition no 5
    int l;
    int l1 = k1;
    int l2 = 1;
    for (l = 1; l <= m; l++) {
        double xi = x[l - 1];
        int break_twice = 0;
        while (!((xi < t[l1]) || (l == nk1))) {
            l1++;
            l2++;
            if (l2 > k1) { break_twice = 1; break; }
        }
        if (break_twice) { break; }
        if (l == m) { break; } // Avoid incrementing l if no break.
    }

    double per = t[nk2 - 1] - t[k];
    for (int i1 = 2; i1 <= l; i1++) {
        int i = i1 - 1;
        int mm = i + m1;
        int spin_i1_loop = 0;
        for (int j = k1; j <= nk1; j++) {
            double tj = t[j - 1];
            int j1 = j + k1;
            double tl = t[j1 - 1];
            double xi = 0.0;
            do {
                i++;
                if (i > mm) { spin_i1_loop = 1; break; }
                int i2 = i - m1;
                if (i2 <= 0) {
                    xi = x[i - 1];
                } else {
                    xi = x[i2 - 1] + per;
                }
            } while (xi <= tj);
            if (xi >= tl) { spin_i1_loop = 1; }
            if (spin_i1_loop) { break; }
        }
        if (spin_i1_loop) continue;
        *ier = 0;
        return;
    }
}


static void
fpclos(const int iopt, const int idim, const int m, const double *u, const int mx,
       const double *x, const double *w, const int k, const double s, const int nest,
       const double tol, const int maxit, const int k1, const int k2, int *n, double *t,
       const int nc, double *c, double *fp, double *fpint, double *z, double *a1,
       double *a2, double *b, double *g1, double *g2, double *q, int *nrdata, int *ier)
{
    (void)mx;  // Unused
    // Subroutine fpclos determines the least-squares closed curve.
    double acc=0.0, cos, d1, fac, fpart, fpms=0.0, fpold=0.0, fp0 = 0.0, f1, f2, f3, p, per, pinv, piv;
    double p1, p2, p3, sin, store, term, ui, wi, rn, one, con1, con4, con9, half;
    int i, ich1, ich3, ij, ik, it, iter, i1=0, i2, i3, j, jj, jk, jper, j1, j2, kk, kk1;
    int l, l0, l1, l5, mm, m1, new, nk1, nk2, nmax = 0, nmin, nplus=0, npl1;
    int nrint=0, n10=0, n11=0, n7=0, n8=0;
    double h[6], h1[7], h2[6], xi[10];

    // set constants
    one = 1.0;
    con1 = 0.1;
    con9 = 0.9;
    con4 = 0.04;
    half = 0.5;

    /////////////////////////////////////////////////////////////////////////
    // part 1: determination of the number of knots and their position.    //
    // **************************************************************      //
    // given a set of knots we compute the least-squares closed curve      //
    // sinf(u). if the sum f(p=inf) <= s we accept the choice of knots.    //
    // if iopt=-1 sinf(u) is the requested curve                           //
    // if iopt=0 or iopt=1 we check whether we can accept the knots:       //
    //   if fp <=s we will continue with the current set of knots.         //
    //   if fp > s we will increase the number of knots and compute the    //
    //      corresponding least-squares curve until finally fp<=s.         //
    // the initial choice of knots depends on the value of s and iopt.     //
    //   if s=0 we have spline interpolation; in that case the number of   //
    //   knots equals nmax = m+2*k.                                        //
    //   if s > 0 and                                                      //
    //     iopt=0 we first compute the least-squares polynomial curve of   //
    //     degree k; n = nmin = 2*k+2. since s(u) must be periodic we      //
    //     find that s(u) reduces to a fixed point.                        //
    //     iopt=1 we start with the set of knots found at the last         //
    //     call of the routine, except for the case that s > fp0; then     //
    //     we compute directly the least-squares polynomial curve.         //
    /////////////////////////////////////////////////////////////////////////

    m1 = m - 1;
    kk = k;
    kk1 = k1;
    nmin = 2 * k1;

    // determine the length of the period of the splines.
    per = u[m - 1] - u[0];

    // We have to carry over the Fortran goto flow since if/else based conversion
    // forces the outer goto to jump inside an if statement which is not allowed in C.
    if (iopt < 0) { goto restart_iteration; }
    // calculation of acc, the absolute tolerance for the root of f(p)=s.
    acc = tol * s;
    // determine nmax, the number of knots for periodic spline interpolation
    nmax = m + 2 * k;
    if ((s > 0.0) || (nmax == nmin)) { goto L30; }
    // if s = 0, s(u) is an interpolating curve.
    *n = nmax;
    // test whether the required storage space exceeds the available one.
    if (*n > nest) {
        *ier = 1;
        return;
    }
L5:
    // find the position of the interior knots in case of interpolation.
    if (k % 2 == 1) {
        // k is odd
        for (i = 2; i <= m1; i++) {
            j = i + k;
            t[j - 1] = u[i - 1];
        }
        if (s <= 0.0) {
            kk = k - 1;
            kk1 = k;
            if (kk <= 0) {
                // Special case k=1
                t[0] = t[m - 1] - per;
                t[1] = u[0];
                t[m] = u[m - 1];
                t[m + 1] = t[2] + per;
                jj = 0;
                for (i = 1; i <= m1; i++) {
                    j = i;
                    for (j1 = 1; j1 <= idim; j1++) {
                        jj++;
                        c[j - 1] = x[jj - 1];
                        j += *n;
                    }
                }
                jj = 1;
                j = m;
                for (j1 = 1; j1 <= idim; j1++) {
                    c[j - 1] = c[jj - 1];
                    j += *n;
                    jj += *n;
                }
                *fp = 0.0;
                fpint[*n - 1] = fp0;
                fpint[*n - 2] = 0.0;
                nrdata[*n - 1] = 0;
                *ier = -1;
                return;
            }
        }
    } else {
        // k is even
        for (i = 2; i <= m1; i++) {
            j = i + k;
            t[j - 1] = (u[i - 1] + u[i - 2]) * half;
        }

    }
    goto restart_iteration;
L30:
    // if s > 0 our initial choice depends on the value of iopt.
    // if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
    // periodic polynomial. (i.e. a fixed point).
    // if iopt=1 and fp0>s we start computing the least-squares periodic
    // spline according the set of knots found at the last call of the
    // routine.

    if ((iopt != 0) && (*n != nmin)) {
        fp0 = fpint[*n - 1];
        fpold = fpint[*n - 2];
        nplus = nrdata[*n - 1];
        if (fp0 > s) { goto restart_iteration; }
    }

    // the case that s(u) is a fixed point is treated separetely.
    // fp0 denotes the corresponding sum of squared residuals.
    fp0 = 0.0;
    d1 = 0.0;
    for (j = 1; j <= idim; j++) {
        z[j - 1] = 0.0;
    }
    jj = 0;
    for (it = 1; it <= m1; it++) {
        wi = w[it - 1];
        fpgivs(wi, &d1, &cos, &sin);
        for (j = 1; j <= idim; j++) {
            jj++;
            fac = wi * x[jj - 1];
            fprota(cos, sin, &fac, &z[j - 1]);
            fp0 = fp0 + (fac * fac);
        }
    }
    for (j = 1; j <= idim; j++) {
        z[j - 1] = z[j - 1] / d1;
    }
    // test whether that fixed point is a solution of our problem.
    fpms = fp0 - s;
    if ((fpms < acc) || (nmax == nmin)) {
        *ier = -2;
        // the point (z(1),z(2),...,z(idim)) is a solution of our problem.
        // a constant function is a spline of degree k with all b-spline
        // coefficients equal to that constant.
        for (i = 1; i <= k1; i++) {
            rn = k1 - i;
            t[i - 1] = u[0] - rn * per;
            j = i + k1;
            rn = i - 1;
            t[j - 1] = u[m - 1] + rn * per;
        }
        *n = nmin;
        j1 = 0;
        for (j = 1; j <= idim; j++) {
            fac = z[j - 1];
            j2 = j1;
            for (i = 1; i <= k1; i++) {
                j2++;
                c[j2 - 1] = fac;
            }
            j1 = j1 + *n;
        }
        *fp = fp0;
        fpint[*n - 1] = fp0;
        fpint[*n - 2] = 0.0;
        nrdata[*n - 1] = 0;
        return;
    }
    fpold = fp0;

    // test whether the required storage space exceeds the available one.
    if (*n >= nest) {
        *ier = 1;
        return;
    }

    // start computing the least-squares closed curve with one
    // interior knot.
    nplus = 1;
    *n = nmin + 1;
    mm = (m + 1) / 2;
    t[k2 - 1] = u[mm - 1];
    nrdata[0] = mm - 2;
    nrdata[1] = m1 - mm;

restart_iteration:
    // main loop for the different sets of knots. m is a save upper
    // bound for the number of trials.
    for (iter = 1; iter <= m; iter++) {
        // find nrint, the number of knot intervals.
        nrint = *n - nmin + 1;

        // find the position of the additional knots which are needed for
        // the b-spline representation of s(u). if we take
        //     t(k+1) = u(1), t(n-k) = u(m)
        //     t(k+1-j) = t(n-k-j) - per, j=1,2,...k
        //     t(n-k+j) = t(k+1+j) + per, j=1,2,...k
        // then s(u) will be a smooth closed curve if the b-spline
        // coefficients satisfy the following conditions
        //     c((i-1)*n+n7+j) = c((i-1)*n+j), j=1,...k,i=1,2,...,idim (**)
        // with n7=n-2*k-1.
        t[k1 - 1] = u[0];
        nk1 = *n - k1;
        nk2 = nk1 + 1;
        t[nk2 - 1] = u[m - 1];
        for (j = 1; j <= k; j++) {
            i1 = nk2 + j;
            i2 = nk2 - j;
            j1 = k1 + j;
            j2 = k1 - j;
            t[i1 - 1] = t[j1 - 1] + per;
            t[j2 - 1] = t[i2 - 1] - per;
        }

        // compute the b-spline coefficients of the least-squares closed curve
        // sinf(u). the observation matrix a is built up row by row while
        // taking into account condition (**) and is reduced to triangular
        // form by givens transformations .
        // at the same time fp=f(p=inf) is computed.
        // the n7 x n7 triangularised upper matrix a has the form
        //           ! a1 '    !
        //       a = !    ' a2 !
        //           ! 0  '    !
        // with a2 a n7 x k matrix and a1 a n10 x n10 upper triangular
        // matrix of bandwidth k+1 ( n10 = n7-k).
        // initialization.
        for (i = 1; i <= nc; i++) {
            z[i - 1] = 0.0;
        }
        for (i = 0; i < nk1; i++) {
            for (j = 0; j < kk1; j++) {
                a1[i + j * nest] = 0.0;
            }
        }
        n7 = nk1 - k;
        n10 = n7 - kk;
        jper = 0;
        *fp = 0.0;
        l = k1;
        jj = 0;

        for (it = 1; it <= m1; it++) {
            // fetch the current data point u(it),x(it)
            ui = u[it - 1];
            wi = w[it - 1];
            for (j = 1; j <= idim; j++) {
                jj++;
                xi[j - 1] = x[jj - 1] * wi;
            }

            // search for knot interval t(l) <= ui < t(l+1).
            while (ui >= t[l]) {
                l++;
            }

            // evaluate the (k+1) non-zero b-splines at ui and store them in q.
            fpbspl(t, *n, k, ui, l, h);
            for (i = 1; i <= k1; i++) {
                q[(it - 1) + (i - 1) * m] = h[i - 1];
                h[i - 1] = h[i - 1] * wi;
            }
            l5 = l - k1;

            // test whether the b-splines nj,k+1(u),j=1+n7,...nk1 are all zero at ui
            if (l5 < n10) {
                // rotation of the new row of the observation matrix into
                // triangle in case the b-splines nj,k+1(u),j=n7+1,...n-k-1 are all zero
                // at ui.
                j = l5;
                for (i = 1; i <= kk1; i++) {
                    j++;
                    piv = h[i - 1];
                    if (piv == 0.0) {
                        continue;
                    }

                    // calculate the parameters of the givens transformation.
                    fpgivs(piv, &a1[j - 1], &cos, &sin);

                    // transformations to right hand side.
                    j1 = j;
                    for (j2 = 1; j2 <= idim; j2++) {
                        fprota(cos, sin, &xi[j2 - 1], &z[j1 - 1]);
                        j1 += *n;
                    }

                    if (i == kk1) {
                        break;
                    }

                    i2 = 1;
                    i3 = i + 1;
                    // transformations to left hand side.
                    for (i1 = i3; i1 <= kk1; i1++) {
                        i2++;
                        fprota(cos, sin, &h[i1 - 1], &a1[(j - 1) + (i2 - 1) * nest]);
                    }
                }

                // add contribution of this row to the sum of squares of residual
                // right hand sides.
                for (j2 = 1; j2 <= idim; j2++) {
                    *fp = *fp + xi[j2 - 1] * xi[j2 - 1];
                }
            } else {
                // test whether the b-splines nj,k+1(u),j=1+n7,...nk1 are all zero at ui
                if (jper == 0) {
                    // initialize the matrix a2.
                    for (i = 1; i <= n7; i++) {
                        for (j = 1; j <= kk; j++) {
                            a2[(i - 1) + (j - 1) * nest] = 0.0;
                        }
                    }
                    jk = n10 + 1;
                    for (i = 1; i <= kk; i++) {
                        ik = jk;
                        for (j = 1; j <= kk1; j++) {
                            if (ik <= 0) { break; }
                            a2[(ik - 1) + (i - 1) * nest] = a1[(ik - 1) + (j - 1) * nest];
                            ik--;
                        }
                        jk++;
                    }
                    jper = 1;
                }

                // if one of the b-splines nj,k+1(u),j=n7+1,...nk1 is not zero at ui
                // we take account of condition (**) for setting up the new row
                // of the observation matrix a. this row is stored in the arrays h1
                // (the part with respect to a1) and h2 (the part with
                // respect to a2).
                for (i = 1; i <= kk; i++) {
                    h1[i - 1] = 0.0;
                    h2[i - 1] = 0.0;
                }
                h1[kk1 - 1] = 0.0;
                j = l5 - n10;
                for (i = 1; i <= kk1; i++) {
                    j++;
                    l0 = j;
                    while (1) {
                        l1 = l0 - kk;
                        if (l1 <= 0) {
                            h2[l0 - 1] += h[i - 1];
                            break;
                        }
                        if (l1 <= n10) {
                            h1[l1 - 1] = h[i - 1];
                            break;
                        }
                        l0 = l1 - n10;
                    }
                }

                // rotate the new row of the observation matrix into triangle
                // by givens transformations.
                if (n10 > 0) {
                    // rotation with the rows 1,2,...n10 of matrix a.
                    for (j = 1; j <= n10; j++) {
                        piv = h1[0];
                        if (piv == 0.0) {
                            for (i = 1; i <= kk; i++) {
                                h1[i - 1] = h1[i];
                            }
                            h1[kk1 - 1] = 0.0;
                            continue;
                        }

                        // calculate the parameters of the givens transformation.
                        fpgivs(piv, &a1[(j - 1) + 0 * nest], &cos, &sin);

                        // transformation to the right hand side.
                        j1 = j;
                        for (j2 = 1; j2 <= idim; j2++) {
                            fprota(cos, sin, &xi[j2 - 1], &z[j1 - 1]);
                            j1 = j1 + *n;
                        }

                        // transformations to the left hand side with respect to a2.
                        for (i = 1; i <= kk; i++) {
                            fprota(cos, sin, &h2[i - 1], &a2[(j - 1) + (i - 1) * nest]);
                        }

                        if (j == n10) { break; }

                        i2 = (n10 - j < kk) ? (n10 - j) : kk;
                        // transformations to the left hand side with respect to a1.
                        for (i = 1; i <= i2; i++) {
                            i1 = i + 1;
                            fprota(cos, sin, &h1[i], &a1[(j - 1) + i * nest]);
                            h1[i - 1] = h1[i];
                        }
                        h1[i1 - 1] = 0.0;
                    }
                }

                // rotation with the rows n10+1,...n7 of matrix a.
                for (j = 1; j <= kk; j++) {
                    ij = n10 + j;
                    if (ij <= 0) { continue; }

                    piv = h2[j - 1];

                    if (piv == 0.0) { continue; }

                    // calculate the parameters of the givens transformation.
                    fpgivs(piv, &a2[(ij - 1) + (j - 1) * nest], &cos, &sin);

                    // transformations to right hand side.
                    j1 = ij;
                    for (j2 = 1; j2 <= idim; j2++) {
                        fprota(cos, sin, &xi[j2 - 1], &z[j1 - 1]);
                        j1 = j1 + *n;
                    }

                    if (j == kk) { break; }

                    j1 = j + 1;
                    // transformations to left hand side.
                    for (i = j1; i <= kk; i++) {
                        fprota(cos, sin, &h2[i - 1], &a2[(ij - 1) + (i - 1) * nest]);
                    }

                }

                // add contribution of this row to the sum of squares of residual
                // right hand sides.
                for (j2 = 1; j2 <= idim; j2++) {
                    *fp += xi[j2 - 1] * xi[j2 - 1];
                }
            }
        }
        // 290

        fpint[*n - 1] = fp0;
        fpint[*n - 2] = fpold;
        nrdata[*n - 1] = nplus;

        // backward substitution to obtain the b-spline coefficients.
        j1 = 1;
        for (j2 = 1; j2 <= idim; j2++) {
            fpbacp(&a1[0], &a2[0], &z[j1 - 1], n7, kk, &c[j1 - 1], kk1, nest);
            j1 = j1 + *n;
        }

        // calculate from condition (**) the remaining coefficients.
        for (i = 1; i <= k; i++) {
            j1 = i;
            for (j = 1; j <= idim; j++) {
                j2 = j1 + n7;
                c[j2 - 1] = c[j1 - 1];
                j1 = j1 + *n;
            }
        }

        if (iopt < 0) { return; }

        // test whether the approximation sinf(u) is an acceptable solution.
        fpms = *fp - s;
        if (fabs(fpms) < acc) { return; }

        // if f(p=inf) < s accept the choice of knots.
        if (fpms < 0.0) {
            break;
        }

        // if n=nmax, sinf(u) is an interpolating curve.
        if (*n == nmax) {
            *ier = -1;
            return;
        }

        // increase the number of knots.
        // if n=nest we cannot increase the number of knots because of the
        // storage capacity limitation.
        if (*n == nest) {
            *ier = 1;
            return;
        }

        // determine the number of knots nplus we are going to add.
        npl1 = nplus * 2;
        rn = nplus;
        if ((fpold - *fp) > acc) {
            npl1 = (int)(rn * fpms / (fpold - *fp));
        }

        {
            // nplus = min0(nplus*2,max0(npl1,nplus/2,1))
            int temp_max = npl1;
            int temp = nplus / 2;
            if (temp > temp_max) {
                temp_max = temp;
            }
            if (1 > temp_max) {
                temp_max = 1;
            }
            int temp_min = nplus * 2;
            if (temp_max < temp_min) {
                nplus = temp_max;
            } else {
                nplus = temp_min;
            }
        }
        fpold = *fp;

        // compute the sum of squared residuals for each knot interval
        // t(j+k) <= ui <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = 0.0;
        i = 1;
        l = k1;
        jj = 0;
        new = 0;
        for (it = 1; it <= m1; it++) {
            if (u[it - 1] >= t[l - 1]) {
                new = 1;
                l = l + 1;
            }
            term = 0.0;
            l0 = l - k2;
            for (j2 = 1; j2 <= idim; j2++) {
                fac = 0.0;
                j1 = l0;
                for (j = 1; j <= k1; j++) {
                    j1++;
                    fac = fac + c[j1 - 1] * q[(it - 1) + (j - 1) * m];
                }
                jj = jj + 1;
                term = term + pow(w[it - 1] * (fac - x[jj - 1]), 2);
                l0 = l0 + *n;
            }
            fpart = fpart + term;
            if (new != 0) {
                if (l > k2) {
                    store = term * half;
                    fpint[i - 1] = fpart - store;
                    i++;
                    fpart = store;
                    new = 0;
                } else {
                    fpint[nrint - 1] = term;
                    new = 0;
                }
            }
        }
        fpint[nrint - 1] += fpart;

        for (l = 1; l <= nplus; l++) {
            // add a new knot
            fpknot(u, m, t, n, fpint, nrdata, &nrint, nest, 1);

            // if n=nmax we locate the knots as for interpolation
            if (*n == nmax) { goto L5; }

            // test whether we cannot further increase the number of knots.
            if (*n == nest) {
                break;
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////
    // part 2: determination of the smoothing closed curve sp(u).          //
    // **********************************************************          //
    // we have determined the number of knots and their position.          //
    // we now compute the b-spline coefficients of the smoothing curve     //
    // sp(u). the observation matrix a is extended by the rows of matrix   //
    // b expressing that the kth derivative discontinuities of sp(u) at    //
    // the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     //
    // ponding weights of these additional rows are set to 1/p.            //
    // iteratively we then have to determine the value of p such that f(p),//
    // the sum of squared residuals be = s. we already know that the least-//
    // squares polynomial curve corresponds to p=0, and that the least-    //
    // squares periodic spline curve corresponds to p=infinity. the        //
    // iteration process which is proposed here, makes use of rational     //
    // interpolation. since f(p) is a convex and strictly decreasing       //
    // function of p, it can be approximated by a rational function        //
    // r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  //
    // ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      //
    // to calculate the new value of p such that r(p)=s. convergence is    //
    // guaranteed by taking f1>0 and f3<0.                                 //
    /////////////////////////////////////////////////////////////////////////

    // evaluate the discontinuity jump of the kth derivative of the
    // b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
    fpdisc(t, *n, k2, b, nest);

    // initial value for p.
    p1 = 0.0;
    f1 = fp0 - s;
    p3 = -one;
    f3 = fpms;
    n11 = n10 - 1;
    n8 = n7 - 1;
    p = 0.0;
    l = n7;
    int skip_352 = 0;
    for (i = 1; i <= k; i++) {
        j = k + 1 - i;
        p = p + a2[(l - 1) + (j - 1) * nest];
        l--;
        if (l == 0) {
            skip_352 = 1;  // go over both loops
            break;
        }
    }
    if (!skip_352) {
        // 352
        for (i = 1; i <= n10; i++) {
            p = p + a1[i - 1];
        }
    }
    // 356
    rn = n7;
    p = rn / p;
    ich1 = 0;
    ich3 = 0;

    // iteration process to find the root of f(p) = s.
    for (iter = 1; iter <= maxit; iter++) {
        // form the matrix g  as the matrix a extended by the rows of matrix b.
        // the rows of matrix b with weight 1/p are rotated into
        // the triangularised observation matrix a.
        // after triangularisation our n7 x n7 matrix g takes the form
        //           ! g1 '    !
        //       g = !    ' g2 !
        //           ! 0  '    !
        // with g2 a n7 x (k+1) matrix and g1 a n11 x n11 upper triangular
        // matrix of bandwidth k+2. ( n11 = n7-k-1)
        pinv = one / p;

        // store matrix a into g
        for (i = 1; i <= nc; i++) {
            c[i - 1] = z[i - 1];
        }
        for (i = 1; i <= n7; i++) {
            g1[(i - 1) + (k1 - 1) * nest] = a1[(i - 1) + (k1 - 1) * nest];
            g1[(i - 1) + (k2 - 1) * nest] = 0.0;
            g2[i - 1] = 0.0;
            for (j = 1; j <= k; j++) {
                g1[(i - 1) + (j - 1) * nest] = a1[(i - 1) + (j - 1) * nest];
                g2[(i - 1) + j * nest] = a2[(i - 1) + (j - 1) * nest];
            }
        }
        l = n10;
        for (j = 1; j <= k1; j++) {
            if (l <= 0) { break; }
            g2[l - 1] = a1[(l - 1) + (j - 1) * nest];
            l--;
        }

        for (it = 1; it <= n8; it++) {
            // fetch a new row of matrix b and store it in the arrays h1 (the part
            // with respect to g1) and h2 (the part with respect to g2).
            for (j = 1; j <= idim; j++) {
                xi[j - 1] = 0.0;
            }
            for (i = 1; i <= k1; i++) {
                h1[i - 1] = 0.0;
                h2[i - 1] = 0.0;
            }
            h1[k2 - 1] = 0.0;

            if (it > n11) {
                l = 1;
                i = it - n10;
                for (j = 1; j <= k2; j++) {
                    i++;
                    l0 = i;
                    while (1) {
                        l1 = l0 - k1;
                        if (l1 <= 0) {
                            h2[l0 - 1] = h2[l0 - 1] + b[(it - 1) + (j - 1) * nest] * pinv;
                            break;
                        }
                        if (l1 <= n11) {
                            h1[l1 - 1] = b[(it - 1) + (j - 1) * nest] * pinv;
                            break;
                        }
                        l0 = l1 - n11;
                    }
                }
            } else {
                l = it;
                l0 = it;
                for (j = 1; j <= k2; j++) {
                    if (l0 == n10) {
                        l0 = 1;
                        for (l1 = j; l1 <= k2; l1++) {
                            h2[l0 - 1] = b[(it - 1) + (l1 - 1) * nest] * pinv;
                            l0++;
                        }
                        break;
                    }
                    h1[j - 1] = b[(it - 1) + (j - 1) * nest] * pinv;
                    l0 = l0 + 1;
                }
            }

            if (n11 > 0) {
                // 470
                // rotate this row into triangle by givens transformations
                // rotation with the rows l,l+1,...n11.
                for (j = l; j <= n11; j++) {
                    piv = h1[0];

                    // calculate the parameters of the givens transformation.
                    fpgivs(piv, &g1[j - 1], &cos, &sin);

                    // transformation to right hand side.
                    j1 = j;
                    for (j2 = 1; j2 <= idim; j2++) {
                        fprota(cos, sin, &xi[j2 - 1], &c[j1 - 1]);
                        j1 += *n;
                    }

                    // transformation to the left hand side with respect to g2.
                    for (i = 1; i <= k1; i++) {
                        fprota(cos, sin, &h2[i - 1], &g2[(j - 1) + (i - 1) * nest]);
                    }

                    if (j == n11) { break; }
                    i2 = (n11 - j < k1) ? (n11 - j) : k1;
                    // transformation to the left hand side with respect to g1.
                    for (i = 1; i <= i2; i++) {
                        i1 = i + 1;
                        fprota(cos, sin, &h1[i], &g1[(j - 1) + i * nest]);
                        h1[i - 1] = h1[i];
                    }
                    h1[i1 - 1] = 0.0;
                }
            }
            // 510

            // rotation with the rows n11+1,...n7
            for (j = 1; j <= k1; j++) {
                ij = n11 + j;
                if (ij <= 0) {
                    continue;
                }
                piv = h2[j - 1];

                // calculate the parameters of the givens transformation
                fpgivs(piv, &g2[(ij - 1) + (j - 1) * nest], &cos, &sin);

                // transformation to the right hand side.
                j1 = ij;
                for (j2 = 1; j2 <= idim; j2++) {
                    fprota(cos, sin, &xi[j2 - 1], &c[j1 - 1]);
                    j1 += *n;
                }
                if (j == k1) { break; }
                j1 = j + 1;

                // transformation to the left hand side.
                for (i = j1; i <= k1; i++) {
                    fprota(cos, sin, &h2[i - 1], &g2[(ij - 1) + (i - 1) * nest]);
                }
            }
        }
        // 540

        // backward substitution to obtain the b-spline coefficients
        j1 = 1;
        for (j2 = 1; j2 <= idim; j2++) {
            fpbacp(&g1[0], &g2[0], &c[j1 - 1], n7, k1, &c[j1 - 1], k2, nest);
            j1 = j1 + *n;
        }
        // 542

        // calculate from condition (**) the remaining b-spline coefficients.
        for (i = 1; i <= k; i++) {
            j1 = i;
            for (j = 1; j <= idim; j++) {
                j2 = j1 + n7;
                c[j2 - 1] = c[j1 - 1];
                j1 += *n;
            }
        }
        // 547

        // computation of f(p).
        *fp = 0.0;
        l = k1;
        jj = 0;
        for (it = 1; it <= m1; it++) {
            if (u[it - 1] >= t[l - 1]) { l++; }
            l0 = l - k2;
            term = 0.0;
            for (j2 = 1; j2 <= idim; j2++) {
                fac = 0.0;
                j1 = l0;
                for (j = 1; j <= k1; j++) {
                    j1++;
                    fac = fac + c[j1 - 1] * q[(it - 1) + (j - 1) * m];
                }
                jj++;
                term = term + (fac - x[jj - 1]) * (fac - x[jj - 1]);
                l0 += *n;
            }
            *fp = *fp + term * w[it - 1] * w[it - 1];
        }
        // 570

        // test whether the approximation sp(u) is an acceptable solution.
        fpms = *fp - s;
        if (fabs(fpms) < acc) {
            return;
        }

        // test whether the maximal number of iterations is reached.
        if (iter == maxit) {
            *ier = 3;
            return;
        }

        // carry out one more step of the iteration process.
        p2 = p;
        f2 = fpms;
        if (ich3 == 0) {
            if ((f2 - f3) > acc) {
                if (f2 < 0.0) {
                    ich3 = 1;
                }
            } else {
                // our initial choice of p is too large.
                p3 = p2;
                f3 = f2;
                p = p * con4;
                if (p <= p1) {
                    p = p1 * con9 + p2 * con1;
                }
                continue;  // spin main iter loop
            }
        }
        if (ich1 == 0) {
            if ((f1 - f2) > acc) {
                if (f2 > 0.0) {
                    ich1 = 1;
                }
            } else {
                // our initial choice of p is too small
                p1 = p2;
                f1 = f2;
                p = p / con4;
                if (p3 < 0.0) {
                    continue;  // spin main iter loop
                }
                if (p >= p3) {
                    p = p2 * con1 + p3 * con9;
                }
                continue;  // spin main iter loop
            }
        }

        // test whether the iteration process proceeds as theoretically
        // expected.
        if ((f2 >= f1) || (f2 <= f3)) {
            *ier = 2;
            return;
        }

        // find the new value for p.
        p = fprati(&p1, &f1, &p2, &f2, &p3, &f3);
    }
}

static void
fpcurf(const int iopt, const double *x, const double *y, const double *w, const int m,
       const double xb, const double xe, const int k, const double s, const int nest,
       const double tol, const int maxit, const int k1, const int k2, int *n, double *t,
       double *c, double *fp, double *fpint, double *z, double *a, double *b, double *g,
       double *q, int *nrdata, int *ier)
{
    double acc=0.0, con1, con4, con9, ccos, half, fpart, fpms, fpold=0.0, fp0=0.0, one, p, pinv, piv;
    double p1, p2, p3, rn, ssin, store, term, wi, xi, yi;
    int i, ich1, ich3, it, iter, i1, i2, i3, j, k3, l, l0;
    int mk1, new, nk1, nmax=0, nmin, nplus=0, npl1, nrint, n8;
    double h[7];

    // set constants
    one = 1.0;
    con1 = 0.1;
    con9 = 0.9;
    con4 = 0.04;
    half = 0.5;

    /////////////////////////////////////////////////////////////////////////
    //  part 1: determination of the number of knots and their position    //
    //  **************************************************************     //
    //  given a set of knots we compute the least-squares spline sinf(x),  //
    //  and the corresponding sum of squared residuals fp=f(p=inf).        //
    //  if iopt=-1 sinf(x) is the requested approximation.                 //
    //  if iopt=0 or iopt=1 we check whether we can accept the knots:      //
    //    if fp <=s we will continue with the current set of knots.        //
    //    if fp > s we will increase the number of knots and compute the   //
    //       corresponding least-squares spline until finally fp<=s.       //
    //    the initial choice of knots depends on the value of s and iopt.  //
    //    if s=0 we have spline interpolation; in that case the number of  //
    //    knots equals nmax = m+k+1.                                       //
    //    if s > 0 and                                                     //
    //      iopt=0 we first compute the least-squares polynomial of        //
    //      degree k; n = nmin = 2*k+2                                     //
    //      iopt=1 we start with the set of knots found at the last        //
    //      call of the routine, except for the case that s > fp0; then    //
    //      we compute directly the least-squares polynomial of degree k.  //
    /////////////////////////////////////////////////////////////////////////
    // determine nmin, the number of knots for polynomial approximation.
    nmin = 2 * k1;
    if (iopt < 0) { goto restart_iteration; }

    // calculation of acc, the absolute tolerance for the root of f(p)=s.
    acc = tol * s;
    // determine nmax, the number of knots for spline interpolation.
    nmax = m + k1;

    if (s > 0.0) {
        // if s>0 our initial choice of knots depends on the value of iopt.
        // if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
        // polynomial of degree k which is a spline without interior knots.
        // if iopt=1 and fp0>s we start computing the least squares spline
        // according to the set of knots found at the last call of the routine.
        if ((iopt != 0) && (*n != nmin)) {
            fp0 = fpint[*n - 1];
            fpold = fpint[*n - 2];
            nplus = nrdata[*n - 1];
            if (fp0 > s) { goto restart_iteration; }
        }
        *n = nmin;
        fpold = 0.0;
        nplus = 0;
        nrdata[0] = m - 2;
    } else {
        // if s=0, s(x) is an interpolating spline.
        // test whether the required storage space exceeds the available one.
        *n = nmax;
        if (nmax > nest) {
            *ier = 1;
            return;
        }
        // find the position of the interior knots in case of interpolation.
        mk1 = m - k1;
        if (mk1 != 0) {
            k3 = k / 2;
            i = k2;
            j = k3 + 2;
            if (k3 * 2 == k) {
                for (l = 1; l <= mk1; l++) {
                    t[i - 1] = (x[j - 1] + x[j - 2]) * half;
                    i++;
                    j++;
                }
            } else {
                for (l = 1; l <= mk1; l++) {
                    t[i - 1] = x[j - 1];
                    i++;
                    j++;
                }
            }
        }
    }

restart_iteration:
    // main loop for the different sets of knots. m is a save upper bound
    // for the number of trials.
    for (iter = 1; iter <= m; iter++) {
        if (*n == nmin) { *ier = -2; }
        // find nrint, tne number of knot intervals.
        nrint = *n - nmin + 1;
        // find the position of the additional knots which are needed for
        // the b-spline representation of s(x).
        nk1 = *n - k1;
        i = *n;
        for (j = 1; j <= k1; j++) {
            t[j - 1] = xb;
            t[i - 1] = xe;
            i--;
        }
        // compute the b-spline coefficients of the least-squares spline
        // sinf(x). the observation matrix a is built up row by row and
        // reduced to upper triangular form by givens transformations.
        // at the same time fp=f(p=inf) is computed.
        *fp = 0.0;
        // initialize the observation matrix a.
        for (i = 1; i <= nk1; i++) {
            z[i - 1] = 0.0;
            for (j = 1; j <= k1; j++) {
                a[(i - 1) + (j - 1) * nest] = 0.0;
            }
        }
        l = k1;
        for (it = 1; it <= m; it++) {
            // fetch the current data point x(it),y(it).
            xi = x[it - 1];
            wi = w[it - 1];
            yi = y[it - 1] * wi;
            // search for knot interval t(l) <= xi < t(l+1).
            while (!((xi < t[l]) || (l == nk1))) {
                l++;
            }
            // evaluate the (k+1) non-zero b-splines at xi and store them in q.
            fpbspl(t, *n, k, xi, l, h);
            for (i = 1; i <= k1; i++) {
                q[(it - 1) + (i - 1) * m] = h[i - 1];
                h[i - 1] = h[i - 1] * wi;
            }
            // rotate the new row of the observation matrix into triangle.
            j = l - k1;
            for (i = 1; i <= k1; i++) {
                j++;
                piv = h[i - 1];
                if (piv == 0.0) { continue; }
                // calculate the parameters of the givens transformation.
                fpgivs(piv, &a[j - 1], &ccos, &ssin);
                // transformations to right hand side.
                fprota(ccos, ssin, &yi, &z[j - 1]);
                if (i == k1) { break; }
                i2 = 1;
                i3 = i + 1;
                for (i1 = i3; i1 <= k1; i1++) {
                    i2++;
                    // transformations to left hand side.
                    fprota(ccos, ssin, &h[i1 - 1], &a[(j - 1) + (i2 - 1) * nest]);
                }
            }
            // add contribution of this row to the sum of squares of residual
            // right hand sides.
            *fp = *fp + yi * yi;
        }
        if (*ier == -2) { fp0 = *fp; }
        fpint[*n - 1] = fp0;
        fpint[*n - 2] = fpold;
        nrdata[*n - 1] = nplus;
        // backward substitution to obtain the b-spline coefficients.
        fpback(a, z, nk1, k1, c, nest);
        // test whether the approximation sinf(x) is an acceptable solution.
        if (iopt < 0) { return; }
        fpms = *fp - s;
        if (fabs(fpms) < acc) { return; }
        // if f(p=inf) < s accept the choice of knots.
        if (fpms < 0.0) { break; }
        // if n = nmax, sinf(x) is an interpolating spline.
        if (*n == nmax) {
            *ier = -1;
            return;
        }
        // increase the number of knots.
        // if n=nest we cannot increase the number of knots because of
        // the storage capacity limitation.
        if (*n == nest) {
            *ier = 1;
            return;
        }
        // determine the number of knots nplus we are going to add.
        if (*ier == 0) {
            npl1 = nplus * 2;
            rn = (double)nplus;
            if (fpold - *fp > acc) { npl1 = (int)(rn * fpms / (fpold - *fp)); }
            // nplus = min(nplus*2, max(npl1, nplus/2, 1))
            int temp1 = npl1;
            int temp2 = nplus / 2;
            if (temp2 > temp1) { temp1 = temp2; }
            if (1 > temp1) { temp1 = 1; }
            int temp3 = nplus * 2;
            if (temp1 < temp3) { temp3 = temp1; }
            nplus = temp3;
        } else {
            nplus = 1;
            *ier = 0;
        }
        fpold = *fp;
        // compute the sum((w(i)*(y(i)-s(x(i))))**2) for each knot interval
        // t(j+k) <= x(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = 0.0;
        i = 1;
        l = k2;
        new = 0;
        for (it = 1; it <= m; it++) {
            if (!((x[it - 1] < t[l - 1]) || (l > nk1))) {
                new = 1;
                l++;
            }
            term = 0.0;
            l0 = l - k2;
            for (j = 1; j <= k1; j++) {
                l0++;
                term = term + c[l0 - 1] * q[(it - 1) + (j - 1) * m];
            }
            term = pow(w[it - 1] * (term - y[it - 1]), 2);
            fpart = fpart + term;
            if (new == 0) { continue; }
            store = term * half;
            fpint[i - 1] = fpart - store;
            i++;
            fpart = store;
            new = 0;
        }
        fpint[nrint - 1] = fpart;
        for (int ll = 1; ll <= nplus; ll++) {
            // add a new knot.
            fpknot(x, m, t, n, fpint, nrdata, &nrint, nest, 1);
            // if n=nmax we locate the knots as for interpolation.
            if (*n == nmax) {
                // find the position of the interior knots in case of interpolation.
                mk1 = m - k1;
                if (mk1 == 0) { goto restart_iteration; }
                k3 = k / 2;
                i = k2;
                j = k3 + 2;
                if (k3 * 2 == k) {
                    for (l = 1; l <= mk1; l++) {
                        t[i - 1] = (x[j - 1] + x[j - 2]) * half;
                        i++;
                        j++;
                    }
                } else {
                    for (l = 1; l <= mk1; l++) {
                        t[i - 1] = x[j - 1];
                        i++;
                        j++;
                    }
                }
                goto restart_iteration;
            }
            // test whether we cannot further increase the number of knots.
            if (*n == nest) { break; }
        }
        // restart the computations with the new set of knots.
    }
    // test whether the least-squares kth degree polynomial is a solution
    // of our approximation problem.
    if (*ier == -2) { return; }

    /////////////////////////////////////////////////////////////////////////
    //  part 2: determination of the smoothing spline sp(x).               //
    //  ***************************************************                //
    //  we have determined the number of knots and their position.         //
    //  we now compute the b-spline coefficients of the smoothing spline   //
    //  sp(x). the observation matrix a is extended by the rows of matrix  //
    //  b expressing that the kth derivative discontinuities of sp(x) at   //
    //  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-    //
    //  ponding weights of these additional rows are set to 1/p.           //
    //  iteratively we then have to determine the value of p such that     //
    //  f(p)=sum((w(i)*(y(i)-sp(x(i))))**2) be = s. we already know that   //
    //  the least-squares kth degree polynomial corresponds to p=0, and    //
    //  that the least-squares spline corresponds to p=infinity. the       //
    //  iteration process which is proposed here, makes use of rational    //
    //  interpolation. since f(p) is a convex and strictly decreasing      //
    //  function of p, it can be approximated by a rational function       //
    //  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond- //
    //  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used     //
    //  to calculate the new value of p such that r(p)=s. convergence is   //
    //  guaranteed by taking f1>0 and f3<0.                                //
    /////////////////////////////////////////////////////////////////////////
    // evaluate the discontinuity jump of the kth derivative of the
    // b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
    fpdisc(t, *n, k2, b, nest);
    // initial value for p.
    p1 = 0.0;
    double f1 = fp0 - s;
    p3 = -one;
    double f3 = fpms;
    p = 0.0;
    for (i = 1; i <= nk1; i++) {
        p = p + a[i - 1];
    }
    rn = (double)nk1;
    p = rn / p;
    ich1 = 0;
    ich3 = 0;
    n8 = *n - nmin;
    // iteration process to find the root of f(p) = s.
    for (iter = 1; iter <= maxit; iter++) {
        // the rows of matrix b with weight 1/p are rotated into the
        // triangularised observation matrix a which is stored in g.
        pinv = one / p;
        for (i = 1; i <= nk1; i++) {
            c[i - 1] = z[i - 1];
            g[(i - 1) + (k2 - 1) * nest] = 0.0;
            for (j = 1; j <= k1; j++) {
                g[(i - 1) + (j - 1) * nest] = a[(i - 1) + (j - 1) * nest];
            }
        }
        for (it = 1; it <= n8; it++) {
            // the row of matrix b is rotated into triangle by givens transformation
            for (i = 1; i <= k2; i++) {
                h[i - 1] = b[(it - 1) + (i - 1) * nest] * pinv;
            }
            yi = 0.0;
            for (j = it; j <= nk1; j++) {
                piv = h[0];
                // calculate the parameters of the givens transformation.
                fpgivs(piv, &g[j - 1], &ccos, &ssin);
                // transformations to right hand side.
                fprota(ccos, ssin, &yi, &c[j - 1]);
                if (j == nk1) { break; }
                i2 = k1;
                if (j > n8) { i2 = nk1 - j; }
                for (i = 1; i <= i2; i++) {
                    // transformations to left hand side.
                    i1 = i + 1;
                    fprota(ccos, ssin, &h[i], &g[(j - 1) + i * nest]);
                    h[i - 1] = h[i];
                }
                h[i2] = 0.0;
            }
        }
        // backward substitution to obtain the b-spline coefficients.
        fpback(g, c, nk1, k2, c, nest);
        // computation of f(p).
        *fp = 0.0;
        l = k2;
        for (it = 1; it <= m; it++) {
            if (!((x[it - 1] < t[l - 1]) || (l > nk1))) {
                l++;
            }
            l0 = l - k2;
            term = 0.0;
            for (j = 1; j <= k1; j++) {
                l0++;
                term = term + c[l0 - 1] * q[(it - 1) + (j - 1) * m];
            }
            *fp = *fp + pow(w[it - 1] * (term - y[it - 1]), 2);
        }
        // test whether the approximation sp(x) is an acceptable solution.
        fpms = *fp - s;
        if (fabs(fpms) < acc) { return; }
        // test whether the maximal number of iterations is reached.
        if (iter == maxit) {
            *ier = 3;
            return;
        }
        // carry out one more step of the iteration process.
        p2 = p;
        double f2 = fpms;
        if (ich3 == 0) {
            if ((f2 - f3) > acc) {
                if (f2 < 0.0) { ich3 = 1; }
            } else {
                // our initial choice of p is too large.
                p3 = p2;
                f3 = f2;
                p = p * con4;
                if (p <= p1) { p = p1 * con9 + p2 * con1; }
                continue;
            }
        }
        if (ich1 == 0) {
            if ((f1 - f2) > acc) {
                if (f2 > 0.0) { ich1 = 1; }
            } else {
                // our initial choice of p is too small
                p1 = p2;
                f1 = f2;
                p = p / con4;
                if (p3 < 0.0) { continue; }
                if (p >= p3) { p = p2 * con1 + p3 * con9; }
                continue;
            }
        }
        // test whether the iteration process proceeds as theoretically
        // expected.
        if ((f2 >= f1) || (f2 <= f3)) {
            *ier = 2;
            return;
        }
        // find the new value for p.
        p = fprati(&p1, &f1, &p2, &f2, &p3, &f3);
    }
}


static void
fpcuro(const double a, const double b, const double c, const double d, double *x, int *n)
{
    // subroutine fpcuro finds the real zeros of a cubic polynomial
    // p(x) = a*x**3+b*x**2+c*x+d.
    int i;
    double a1, b1, c1, df, disc, d1, e3, f, pi3, p3, q, r;
    double step, u, u1, u2, y;

    // set constants
    e3 = 0.1 / 0.3;
    pi3 = PI / 3.0;
    a1 = fabs(a);
    b1 = fabs(b);
    c1 = fabs(c);
    d1 = fabs(d);

    // test whether p(x) is a third degree polynomial.
    if (fmax(fmax(b1, c1), d1) < a1 * 1.0e4) {
        // p(x) is a third degree polynomial.
        b1 = b / a * e3;
        c1 = c / a;
        d1 = d / a;
        q = c1 * e3 - b1 * b1;
        r = b1 * b1 * b1 + (d1 - b1 * c1) * 0.5;
        disc = q * q * q + r * r;
        if (disc > 0.0) {
            // one real root
            u = sqrt(disc);
            u1 = -r + u;
            u2 = -r - u;
            *n = 1;
            x[0] = copysign(pow(fabs(u1), e3), u1) + copysign(pow(fabs(u2), e3), u2) - b1;
        } else {
            // three real roots
            u = sqrt(fabs(q));
            if (r < 0.0) { u = -u; }
            p3 = atan2(sqrt(-disc), fabs(r)) * e3;
            u2 = u + u;
            *n = 3;
            x[0] = -u2 * cos(p3) - b1;
            x[1] = u2 * cos(pi3 - p3) - b1;
            x[2] = u2 * cos(pi3 + p3) - b1;
        }
    } else if (fmax(c1, d1) < b1 * 1.0e4) {
        // p(x) is a second degree polynomial.
        disc = c * c - 4.0 * b * d;
        *n = 0;
        if (disc < 0.0) { return; }
        *n = 2;
        u = sqrt(disc);
        b1 = b + b;
        x[0] = (-c + u) / b1;
        x[1] = (-c - u) / b1;
    } else if (d1 < c1 * 1.0e4) {
        // p(x) is a first degree polynomial.
        *n = 1;
        x[0] = -d / c;
    } else {
        // p(x) is a constant function.
        *n = 0;
        return;
    }

    // apply a newton iteration to improve the accuracy of the roots.
    for (i = 0; i < *n; i++) {
        y = x[i];
        f = ((a * y + b) * y + c) * y + d;
        df = (3.0 * a * y + 2.0 * b) * y + c;
        step = 0.0;
        if (fabs(f) < fabs(df) * 0.1) { step = f / df; }
        x[i] = y - step;
    }
    return;
}


static void fpcyt1(double *a, const int n, const int nn)
{
    // (l u)-decomposition of a cyclic tridiagonal matrix with the non-zero
    // elements stored as follows
    //
    //    | a(1,2) a(1,3)                                    a(1,1)  |
    //    | a(2,1) a(2,2) a(2,3)                                     |
    //    |        a(3,1) a(3,2) a(3,3)                              |
    //    |               ...............                            |
    //    |                               a(n-1,1) a(n-1,2) a(n-1,3) |
    //    | a(n,3)                                  a(n,1)   a(n,2)  |
    //

    double aa, beta, gamma, sum, teta, v;

    beta = 1.0 / a[nn];
    gamma = a[(n - 1) + 2 * nn];
    teta = a[0];
    teta *= beta;
    a[3 * nn] = beta;
    a[4 * nn] = gamma;
    a[5 * nn] = teta;
    sum = gamma * teta;

    for (int i = 2; i <= n - 2; i++) {
        v = a[(i - 2) + 2 * nn] * beta;
        aa = a[i - 1];
        beta = 1.0 / (a[(i - 1) + nn] - aa * v);
        gamma = -gamma * v;
        teta = -teta * aa * beta;
        a[(i - 1) + 3 * nn] = beta;
        a[(i - 1) + 4 * nn] = gamma;
        a[(i - 1) + 5 * nn] = teta;
        sum += gamma * teta;
    }

    v = a[(n - 3) + 2 * nn] * beta;
    aa = a[n - 2];
    beta = 1.0 / (a[(n - 2) + nn] - aa * v);
    gamma = a[n - 1] - gamma * v;
    teta = (a[(n - 2) + 2 * nn] - teta * aa) * beta;
    a[(n - 2) + 3 * nn] = beta;
    a[(n - 2) + 4 * nn] = gamma;
    a[(n - 2) + 5 * nn] = teta;
    a[(n - 1) + 3 * nn] = 1.0 / (a[(n - 1) + nn] - (sum + gamma * teta));
}


static void fpcyt2(const double *a, const int n, const double *b, double *c, const int nn)
{
    // subroutine fpcyt2 solves a linear n x n system
    //         a * c = b
    // where matrix a is a cyclic tridiagonal matrix, decomposed
    // using subroutine fpsyt1.
    double cc, sum;

    c[0] = b[0] * a[3 * nn];
    sum = c[0] * a[4 * nn];
    for (int i = 1; i < n - 1; i++) {
        c[i] = (b[i] - a[i] * c[i - 1]) * a[i + 3 * nn];
        sum += c[i] * a[i + 4 * nn];
    }
    cc = (b[n - 1] - sum) * a[(n - 1) + 3 * nn];
    c[n - 1] = cc;
    c[n - 2] -= cc * a[(n - 2) + 5 * nn];
    int j = n - 1;
    for (int i = 3; i <= n; i++){
        int j1 = j - 1;
        c[j1 - 1] = c[j1 - 1] - c[j - 1] * a[(j1 - 1) + 2 * nn] * a[(j1 - 1) + 3 * nn]
                    - cc * a[(j1 - 1) + 5 * nn];
        j = j1;
    }
}


static void
fpdisc(const double* t, const int n, const int k2, double* b, const int nest)
{
    // Subroutine fpdisc calculates the discontinuity jumps of the kth
    // derivative of the b-splines of degree k at the knots t(k+2)..t(n-k-1).

    int k1 = k2 - 1;
    int k = k1 - 1;
    int nk1 = n - k1;
    int nrint = nk1 - k;
    double an = (double)nrint;
    double fac = an / (t[nk1] - t[k1 - 1]);

    double h[12];

    for (int l = k2; l <= nk1; l++) {
        int lmk = l - k1;
        for (int j = 1; j <= k1; j++) {
            int ik = j + k1;
            int lj = l + j;
            int lk = lj - k2;
            h[j-1] = t[l-1] - t[lk-1];
            h[ik-1] = t[l-1] - t[lj-1];
        }

        int lp = lmk;
        for (int j = 1; j <= k2; j++) {
            int jk = j;
            double prod = h[j-1];
            for (int i = 1; i <= k; i++) {
                jk++;
                prod = prod*h[jk-1]*fac;
            }
            int lk = lp + k1;
            b[lmk-1 + (j-1)*nest] = (t[lk-1] - t[lp-1]) / prod;
            lp++;
        }
    }
}


static void
fpgivs(const double piv, double* ww, double* cos, double* sin)
{
    // Subroutine fpgivs calculates the parameters of a Givens transformation.

    double dd;
    const double one = 1.0;
    double store = fabs(piv);

    if (store >= *ww) {
        dd = store * sqrt(one + (*ww / piv) * (*ww / piv));
    } else {
        dd = *ww * sqrt(one + (piv / *ww) * (piv / *ww));
    }

    *cos = *ww / dd;
    *sin = piv / dd;
    *ww = dd;
}


static void
fpgrre(int ifsx, int ifsy, int ifbx, int ifby, const double *x, const int mx, const double *y, const int my,
       const double *z, const int mz, const int kx, const int ky, const double *tx, const int nx, const double *ty,
       const int ny, const double p, double *c, const int nc, double *fp, double *fpx, double *fpy, const int mm,
       const int mynx, const int kx1, const int kx2, const int ky1, const int ky2, double *spx, double *spy,
       double *right, double *q, double *ax, double *ay, double *bx, double *by, int *nrx, int *nry)
{
    (void)mz;   // Unused
    (void)nc;   // Unused
    (void)mm;   // Unused
    (void)mynx; // Unused
    const double one = 1.0;
    const double half = 0.5;

    int i, ibandx, ibandy, ic, iq, irot, it, iz, i1, i2, j, k, k1, k2, l, l1, l2;
    int ncof, nk1x, nk1y, nrold, nroldx, nroldy, number, numx, numx1, numy, numy1, n1;
    double arg, cos, fac, pinv, piv, sin, term;
    double h[7];

    nk1x = nx - kx1;
    nk1y = ny - ky1;
    if (p > 0.0) { pinv = one / p; }

    // calculate the non-zero elements of matrix (spx) - observation matrix in x-direction
    if (ifsx == 0) {
        l = kx1;
        l1 = kx2;
        number = 0;
        for (it = 1; it <= mx; it++) {
            arg = x[it - 1];
            while (!((arg < tx[l1 - 1]) || (l == nk1x))) {
                l = l1;
                l1 = l + 1;
                number++;
            }
            fpbspl(tx, nx, kx, arg, l, h);
            for (i = 1; i <= kx1; i++) {
                spx[(it - 1) + (i- 1) * mx] = h[i - 1];
            }
            nrx[it - 1] = number;
        }
        ifsx = 1;
    }

    // calculate the non-zero elements of matrix (spy) - observation matrix in y-direction
    if (ifsy == 0) {
        l = ky1;
        l1 = ky2;
        number = 0;
        for (it = 1; it <= my; it++) {
            arg = y[it - 1];
            while (!((arg < ty[l1 - 1]) || (l == nk1y))) {
                l = l1;
                l1 = l + 1;
                number++;
            }
            fpbspl(ty, ny, ky, arg, l, h);
            for (i = 1; i <= ky1; i++) {
                spy[(it - 1) + (i- 1) * my] = h[i - 1];
            }
            nry[it - 1] = number;
        }
        ifsy = 1;
    }

    if (p > 0.0) {
        // calculate the non-zero elements of matrix (bx)
        if ((ifbx == 0) && (nx != (2 * kx1))) {
            fpdisc(tx, nx, kx2, bx, nx);
        }
        // calculate the non-zero elements of matrix (by)
        if ((ifby == 0) && (ny != (2 * ky1))) {
            fpdisc(ty, ny, ky2, by, ny);
        }
    }

    // reduce matrix (ax) to upper triangular form (rx) using givens rotations
    // initialization
    for (i = 0; i < my * nk1x; i++) {
        q[i] = 0.0;
    }
    for (i = 0; i < nk1x; i++) {
        for (j = 0; j < kx2; j++) {
            ax[i + j * nx] = 0.0;
        }
    }

    l = 0;
    nrold = 0;
    // ibandx denotes the bandwidth of the matrices (ax) and (rx).
    ibandx = kx1;

    for (it = 1; it <= mx; it++) {
        number = nrx[it - 1];

        while (1) {

            if (nrold == number) {
                h[ibandx - 1] = 0.0;
                for (j = 1; j <= kx1; j++) {
                    h[j - 1] = spx[(it - 1) + (j - 1) * mx];
                }
                // find the appropriate column of q.
                for (j = 1; j <= my; j++) {
                    l++;
                    right[j - 1] = z[l - 1];
                }
                irot = number;
            } else {
                if (p <= 0.0) {
                    nrold++;
                    continue;
                }
                ibandx = kx2;
                // fetch a new row of matrix (bx).
                n1 = nrold + 1;
                for (j = 1; j <= kx2; j++) {
                    h[j - 1] = bx[(n1 - 1) + (j - 1) * nx] * pinv;
                }
                // find the appropriate column of q.
                for (j = 1; j <= my; j++) {
                    right[j - 1] = 0.0;
                }
                irot = nrold;
            }

            // rotate the new row of matrix (ax) into triangle.
            for (i = 1; i <= ibandx; i++) {
                irot++;
                piv = h[i - 1];
                if (piv == 0.0) { continue; }
                // calculate the parameters of the givens transformation.
                fpgivs(piv, &ax[irot - 1], &cos, &sin);
                // apply that transformation to the rows of matrix q.
                iq = (irot - 1) * my;
                for (j = 1; j <= my; j++) {
                    iq++;
                    fprota(cos, sin, &right[j - 1], &q[iq - 1]);
                }
                // apply that transformation to the columns of (ax).
                if (i == ibandx) { break; }
                i2 = 1;
                for (j = i + 1; j <= ibandx; j++) {
                    i2++;
                    fprota(cos, sin, &h[j - 1], &ax[(irot - 1) + (i2 - 1) * nx]);
                }
            }
            // 240
            if (nrold == number) { break; }
            nrold++;
        }
    }

    // reduce the matrix (ay) to upper triangular form (ry) using givens
    // rotations. apply the same transformations to the columns of matrix g
    // to obtain the (ny-ky-1) x (nx-kx-1) matrix h.
    // store matrix (ry) into (ay) and h into c.
    ncof = nk1x * nk1y;
    for (i = 0; i < ncof; i++) {
        c[i] = 0.0;
    }
    for (i = 0; i < nk1y; i++) {
        for (j = 0; j < ky2; j++) {
            ay[i + j * ny] = 0.0;
        }
    }

    nrold = 0;
    // ibandy denotes the bandwidth of the matrices (ay) and (ry).
    ibandy = ky1;

    for (it = 1; it <= my; it++) {
        number = nry[it - 1];

        while (1) {

            if (nrold == number) {
                h[ibandy - 1] = 0.0;
                for (j = 1; j <= ky1; j++) {
                    h[j - 1] = spy[(it - 1) + (j - 1) * my];
                }
                // find the appropriate row of g.
                l = it;
                for (j = 1; j <= nk1x; j++) {
                    right[j - 1] = q[l - 1];
                    l += my;
                }
                irot = number;
            } else {
                if (p <= 0.0) {
                    nrold++;
                    continue;
                }
                ibandy = ky2;
                // fetch a new row of matrix (by).
                n1 = nrold + 1;
                for (j = 1; j <= ky2; j++) {
                    h[j - 1] = by[(n1 - 1) + (j - 1) * ny] * pinv;
                }
                // find the appropriate row of g.
                for (j = 1; j <= nk1x; j++) {
                    right[j - 1] = 0.0;
                }
                irot = nrold;
            }

            // rotate the new row of matrix (ay) into triangle.
            for (i = 1; i <= ibandy; i++) {
                irot++;
                piv = h[i - 1];
                if (piv == 0.0) { continue; }
                // calculate the parameters of the givens transformation.
                fpgivs(piv, &ay[irot - 1], &cos, &sin);
                // apply that transformation to the rows of matrix g.
                ic = irot;
                for (j = 1; j <= nk1x; j++) {
                    fprota(cos, sin, &right[j - 1], &c[ic - 1]);
                    ic += nk1y;
                }
                // apply that transformation to the columns of (ay).
                if (i == ibandy) { break; }
                i2 = 1;
                for (j = i + 1; j <= ibandy; j++) {
                    i2++;
                    fprota(cos, sin, &h[j - 1], &ay[(irot - 1) + (i2 - 1) * ny]);
                }
            }
            // 390
            if (nrold == number) { break; }
            nrold++;
        }
    }

    // backward substitution to obtain the b-spline coefficients as the
    // solution of the linear system    (ry) c (rx)' = h.
    // first step: solve the system  (ry) (c1) = h.
    k = 1;
    for (i = 1; i <= nk1x; i++) {
        fpback(ay, &c[k - 1], nk1y, ibandy, &c[k - 1], ny);
        k += nk1y;
    }

    // second step: solve c(rx)' = (c1)
    k = 0;
    for (j = 1; j <= nk1y; j++) {
        k++;
        l = k;
        for (i = 1; i <= nk1x; i++) {
            right[i - 1] = c[l - 1];
            l += nk1y;
        }
        fpback(ax, right, nk1x, ibandx, right, nx);
        l = k;
        for (i = 1; i <= nk1x; i++) {
            c[l - 1] = right[i - 1];
            l += nk1y;
        }
    }

    // calculate the quantities
    //   res(i,j) = (z(i,j) - s(x(i),y(j)))**2 , i=1,2,..,mx;j=1,2,..,my
    //   fp = sumi=1,mx(sumj=1,my(res(i,j)))
    //   fpx(r) = sum''i(sumj=1,my(res(i,j))) , r=1,2,...,nx-2*kx-1
    //                 tx(r+kx) <= x(i) <= tx(r+kx+1)
    //   fpy(r) = sumi=1,mx(sum''j(res(i,j))) , r=1,2,...,ny-2*ky-1
    //                 ty(r+ky) <= y(j) <= ty(r+ky+1)
    *fp = 0.0;
    for (i = 0; i < nx; i++) {
        fpx[i] = 0.0;
    }
    for (i = 0; i < ny; i++) {
        fpy[i] = 0.0;
    }
    nk1y = ny - ky1;
    iz = 0;
    nroldx = 0;

    // main loop for the different grid points.
    for (i1 = 1; i1 <= mx; i1++) {
        numx = nrx[i1 - 1];
        numx1 = numx + 1;
        nroldy = 0;

        for (i2 = 1; i2 <= my; i2++) {
            numy = nry[i2 - 1];
            numy1 = numy + 1;
            iz++;

            // evaluate s(x,y) at the current grid point by making the sum of the
            // cross products of the non-zero b-splines at (x,y), multiplied with
            // the appropriate b-spline coefficients.
            term = 0.0;
            k1 = numx * nk1y + numy;
            for (l1 = 1; l1 <= kx1; l1++) {
                k2 = k1;
                fac = spx[(i1 - 1) + (l1 - 1) * mx];
                for (l2 = 1; l2 <= ky1; l2++) {
                    k2++;
                    term = term + fac * spy[(i2 - 1) + (l2 - 1) * my] * c[k2-1];
                }
                k1 += nk1y;
            }

            // calculate the squared residual at the current grid point.
            term = pow(z[iz] - term, 2);

            // adjust the different parameters.
            *fp = *fp + term;
            fpx[numx1 - 1] = fpx[numx1 - 1] + term;
            fpy[numy1 - 1] = fpy[numy1 - 1] + term;
            fac = term * half;

            if (numy != nroldy) {
                fpy[numy1 - 1] = fpy[numy1 - 1] - fac;
                fpy[numy - 1]  = fpy[numy - 1]  + fac;
            }
            nroldy = numy;

            if (numx != nroldx) {
                fpx[numx1 - 1] = fpx[numx1 - 1] - fac;
                fpx[numx - 1]  = fpx[numx - 1]  + fac;
            }
        }
        nroldx = numx;
    }
}


static void
fpgrsp(int ifsu, int ifsv, int ifbu, int ifbv, int iback, const double *u, const int mu, const double *v,
       const int mv, const double *r, const int mr, const double *dr, int iop0, int iop1, const double *tu, const int nu,
       const double *tv, const int nv, const double p, double *c, const int nc, double *sq, double *fp, double *fpu,
       double *fpv, const int mm, const int mvnu, double *spu, double *spv, double *right, double *q, double *au,
       double *av1, double *av2, double *bu, double *bv, double *a0, double *a1, double *b0, double *b1, double *c0,
       double *c1, double *cosi, int *nru, int *nrv)
{
    (void)mr;    // Unused
    (void)nc;    // Unused
    (void)mm;    // Unused
    (void)mvnu;  // Unused
    double arg, co, dr01, dr02, dr03, dr11, dr12, dr13, fac, fac0, fac1, pinv=0.0, piv;
    double si, term, one, three, half;
    int i, ic, ii, ij, ik, iq, irot, it, ir, i0, i1, i2, i3, j, jj, jk, jper;
    int j0, j1, k, k1, k2, l, l0, l1, l2, mvv, ncof, nrold, nroldu, nroldv, number;
    int numu, numu1, numv, numv1, nuu, nu4, nu7, nu8, nv11, nv4, nv7, nv8, n1;
    double h[5], h1[5], h2[4];

    // set constants
    one = 1.0;
    three = 3.0;
    half = 0.5;
    // initialization
    nu4 = nu - 4;
    nu7 = nu - 7;
    nu8 = nu - 8;
    nv4 = nv - 4;
    nv7 = nv - 7;
    nv8 = nv - 8;
    nv11 = nv - 11;
    nuu = nu4 - iop0 - iop1 - 2;
    if (p > 0.0) { pinv = one / p; }
    // it depends on the value of the flags ifsu,ifsv,ifbu,ifbv,iop0,iop1
    // and on the value of p whether the matrices (spu), (spv), (bu), (bv),
    // (cosi) still must be determined.
    if (ifsu == 0) {
        // calculate the non-zero elements of the matrix (spu) which is the ob-
        // servation matrix according to the least-squares spline approximation
        // problem in the u-direction.
        l = 4;
        l1 = 5;
        number = 0;
        for (it = 1; it <= mu; it++) {
            arg = u[it - 1];
            while (!((arg < tu[l1 - 1]) || (l == nu4))) {
                l = l1;
                l1 = l + 1;
                number++;
            }
            fpbspl(tu, nu, 3, arg, l, h);
            for (i = 1; i <= 4; i++) {
                spu[(it - 1) + (i - 1) * mu] = h[i - 1];
            }
            nru[it - 1] = number;
        }
        ifsu = 1;
    }
    // calculate the non-zero elements of the matrix (spv) which is the ob-
    // servation matrix according to the least-squares spline approximation
    // problem in the v-direction.
    if (ifsv == 0) {
        l = 4;
        l1 = 5;
        number = 0;
        for (it = 1; it <= mv; it++) {
            arg = v[it - 1];
            while (!((arg < tv[l1 - 1]) || (l == nv4))) {
                l = l1;
                l1 = l + 1;
                number++;
            }
            fpbspl(tv, nv, 3, arg, l, h);
            for (i = 1; i <= 4; i++) {
                spv[(it - 1) + (i - 1) * mv] = h[i - 1];
            }
            nrv[it - 1] = number;
        }
        ifsv = 1;
        if ((iop0 != 0) || (iop1 != 0)) {
            // calculate the coefficients of the interpolating splines for cos(v)
            // and sin(v).
            for (i = 1; i <= nv4; i++) {
                cosi[0 + (i - 1) * 2] = 0.0;
                cosi[1 + (i - 1) * 2] = 0.0;
            }
            if (nv7 >= 4) {
                for (i = 1; i <= nv7; i++) {
                    l = i + 3;
                    arg = tv[l - 1];
                    fpbspl(tv, nv, 3, arg, l, h);
                    for (j = 1; j <= 3; j++) {
                        av1[(i - 1) + (j - 1) * nv] = h[j - 1];
                    }
                    cosi[0 + (i - 1) * 2] = cos(arg);
                    cosi[1 + (i - 1) * 2] = sin(arg);
                }
                fpcyt1(av1, nv7, nv);
                for (j = 1; j <= 2; j++) {
                    for (i = 1; i <= nv7; i++) {
                        right[i - 1] = cosi[(j - 1) + (i - 1) * 2];
                    }
                    fpcyt2(av1, nv7, right, right, nv);
                    for (i = 1; i <= nv7; i++) {
                        cosi[(j - 1) + i * 2] = right[i - 1];
                    }
                    cosi[(j - 1) + 0 * 2] = cosi[(j - 1) + (nv7) * 2];
                    cosi[(j - 1) + (nv7 + 1) * 2] = cosi[(j - 1) + 1 * 2];
                    cosi[(j - 1) + (nv4 - 1) * 2] = cosi[(j - 1) + 2 * 2];
                }
            }
        }
    }
    if (p > 0.0) {
        // calculate the non-zero elements of the matrix (bu).
        if ((ifbu == 0) && (nu8 != 0)) {
            fpdisc(tu, nu, 5, bu, nu);
            ifbu = 1;
        }
        // calculate the non-zero elements of the matrix (bv).
        if ((ifbv == 0) && (nv8 != 0)) {
            fpdisc(tv, nv, 5, bv, nv);
            ifbv = 1;
        }
    }
    // substituting (2),(3) and (4) into (1), we obtain the overdetermined
    // system
    //      (5)  (avv) (cc) (auu)' = (qq)
    // from which the nuu*nv7 remaining coefficients
    //      c(i,j) , i=2+iop0,3+iop0,...,nu-5-iop1,j=1,2,...,nv-7.
    // the elements of (cc), are then determined in the least-squares sense.
    // simultaneously, we compute the resulting sum of squared residuals sq.
    dr01 = dr[0];
    dr11 = dr[3];
    for (i = 1; i <= mv; i++) {
        a0[(i - 1) * 2] = dr01;
        a1[(i - 1) * 2] = dr11;
    }
    if ((nv8 != 0) && (p > 0.0)) {
        for (i = 1; i <= nv8; i++) {
            b0[(0) + (i - 1) * 2] = 0.0;
            b1[(0) + (i - 1) * 2] = 0.0;
        }
    }
    mvv = mv;
    if (iop0 != 0) {
        fac = (tu[4] - tu[3]) / three;
        dr02 = dr[1] * fac;
        dr03 = dr[2] * fac;
        for (i = 1; i <= nv4; i++) {
            c0[i - 1] = dr01 + dr02 * cosi[(0) + (i - 1) * 2] + dr03 * cosi[(1) + (i - 1) * 2];
        }
        for (i = 1; i <= mv; i++) {
            number = nrv[i - 1];
            fac = 0.0;
            for (j = 1; j <= 4; j++) {
                number++;
                fac = fac + c0[number - 1] * spv[(i - 1) + (j - 1) * mv];
            }
            a0[(1) + (i - 1) * 2] = fac;
        }
        if ((nv8 != 0) && (p > 0.0)) {
            for (i = 1; i <= nv8; i++) {
                number = i;
                fac = 0.0;
                for (j = 1; j <= 5; j++) {
                    fac = fac + c0[number - 1] * bv[(i - 1) + (j - 1) * nv];
                    number++;
                }
                b0[(1) + (i - 1) * 2] = fac * pinv;
            }
        }
        mvv = mv + nv8;
    }
    if (iop1 != 0) {
        fac = (tu[nu4 - 1] - tu[nu4]) / three;
        dr12 = dr[4] * fac;
        dr13 = dr[5] * fac;
        for (i = 1; i <= nv4; i++) {
            c1[i - 1] = dr11 + dr12 * cosi[(0) + (i - 1) * 2] + dr13 * cosi[(1) + (i - 1) * 2];
        }
        for (i = 1; i <= mv; i++) {
            number = nrv[i - 1];
            fac = 0.0;
            for (j = 1; j <= 4; j++) {
                number++;
                fac = fac + c1[number - 1] * spv[(i - 1) + (j - 1) * mv];
            }
            a1[(1) + (i - 1) * 2] = fac;
        }
        if ((nv8 != 0) && (p > 0.0)) {
            for (i = 1; i <= nv8; i++) {
                number = i;
                fac = 0.0;
                for (j = 1; j <= 5; j++) {
                    fac = fac + c1[number - 1] * bv[(i - 1) + (j - 1) * nv];
                    number++;
                }
                b1[(1) + (i - 1) * 2] = fac * pinv;
            }
        }
        mvv = mv + nv8;
    }
    // we first determine the matrices (auu) and (qq). then we reduce the
    // matrix (auu) to an unit upper triangular form (ru) using givens
    // rotations without square roots. we apply the same transformations to
    // the rows of matrix qq to obtain the mv x nuu matrix g.
    // we store matrix (ru) into au and g into q.
    l = mvv * nuu;
    // initialization.
    *sq = 0.0;
    if (l != 0) {
        for (i = 1; i <= l; i++) {
            q[i - 1] = 0.0;
        }
        for (i = 1; i <= nuu; i++) {
            for (j = 1; j <= 5; j++) {
                au[(i - 1) + (j - 1) * nu] = 0.0;
            }
        }
    }
    l = 0;
    nrold = 0;
    n1 = nrold + 1;
    for (it = 1; it <= mu; it++) {
        number = nru[it - 1];
        // find the appropriate column of q.
        while (1) {
            for (j = 1; j <= mvv; j++) {
                right[j - 1] = 0.0;
            }
            if (nrold == number) {
                // fetch a new row of matrix (spu).
                for (j = 1; j <= 4; j++) {
                    h[j - 1] = spu[(it - 1) + (j - 1) * mu];
                }
                // find the appropriate column of q.
                for (j = 1; j <= mv; j++) {
                    l++;
                    right[j - 1] = r[l - 1];
                }
                i0 = 1;
                i1 = 4;
            } else {
                if (p <= 0.0) {
                    nrold = n1;
                    n1++;
                    continue;
                }
                // fetch a new row of matrix (bu).
                for (j = 1; j <= 5; j++) {
                    h[j - 1] = bu[(n1 - 1) + (j - 1) * nu] * pinv;
                }
                i0 = 1;
                i1 = 5;
            }
            j0 = n1;
            j1 = nu7 - number;
            // take into account that we eliminate the constraints (3)
            while (j0 - 1 <= iop0) {
                fac0 = h[i0 - 1];
                for (j = 1; j <= mv; j++) {
                    right[j - 1] = right[j - 1] - fac0 * a0[(j0 - 1) + (j - 1) * 2];
                }
                if (mv != mvv) {
                    j = mv;
                    for (jj = 1; jj <= nv8; jj++) {
                        j++;
                        right[j - 1] = right[j - 1] - fac0 * b0[(j0 - 1) + (jj - 1) * 2];
                    }
                }
                j0++;
                i0++;
            }
            // take into account that we eliminate the constraints (4)
            while (j1 - 1 <= iop1) {
                fac1 = h[i1 - 1];
                for (j = 1; j <= mv; j++) {
                    right[j - 1] = right[j - 1] - fac1 * a1[(j1 - 1) + (j - 1) * 2];
                }
                if (mv != mvv) {
                    j = mv;
                    for (jj = 1; jj <= nv8; jj++) {
                        j++;
                        right[j - 1] = right[j - 1] - fac1 * b1[(j1 - 1) + (jj - 1) * 2];
                    }
                }
                j1++;
                i1--;
            }
            irot = nrold - iop0 - 1;
            if (irot < 0) { irot = 0; }
            // rotate the new row of matrix (auu) into triangle.
            if (i0 <= i1) {
                for (i = i0; i <= i1; i++) {
                    irot++;
                    piv = h[i - 1];
                    if (piv == 0.0) { continue; }
                    // calculate the parameters of the givens transformation.
                    fpgivs(piv, &au[irot - 1], &co, &si);
                    // apply that transformation to the rows of matrix (qq).
                    iq = (irot - 1) * mvv;
                    for (j = 1; j <= mvv; j++) {
                        iq++;
                        fprota(co, si, &right[j - 1], &q[iq - 1]);
                    }
                    // apply that transformation to the columns of (auu).
                    if (i != i1) {
                        i2 = 1;
                        i3 = i + 1;
                        for (j = i3; j <= i1; j++) {
                            i2++;
                            fprota(co, si, &h[j - 1], &au[(irot - 1) + (i2 - 1) * nu]);
                        }
                    }
                }
            }
            // we update the sum of squared residuals.
            for (j = 1; j <= mvv; j++) {
                *sq = *sq + (right[j - 1] * right[j - 1]);
            }
            if (nrold == number) { break; }
            nrold = n1;
            n1++;
        }
    }
    if (nuu == 0) {
        // calculate from the conditions (2)-(3)-(4), the remaining b-spline
        // coefficients.
        ncof = nu4 * nv4;
        j = ncof;
        for (l = 1; l <= nv4; l++) {
            q[l - 1] = dr01;
            q[j - 1] = dr11;
            j--;
        }
        i = nv4;
        j = 0;
        if (iop0 != 0) {
            for (l = 1; l <= nv4; l++) {
                i++;
                q[i - 1] = c0[l - 1];
            }
        }
        if (iop1 != 0) {
            for (l = 1; l <= nv4; l++) {
                i++;
                q[i - 1] = c1[l - 1];
            }
        }
        for (i = 1; i <= ncof; i++) {
            c[i - 1] = q[i - 1];
        }
        // calculate the quantities
        //   res(i,j) = (r(i,j) - s(u(i),v(j)))**2 , i=1,2,..,mu;j=1,2,..,mv
        //   fp = sumi=1,mu(sumj=1,mv(res(i,j)))
        //   fpu(r) = sum''i(sumj=1,mv(res(i,j))) , r=1,2,...,nu-7
        //                 tu(r+3) <= u(i) <= tu(r+4)
        //   fpv(r) = sumi=1,mu(sum''j(res(i,j))) , r=1,2,...,nv-7
        //                 tv(r+3) <= v(j) <= tv(r+4)
        *fp = 0.0;
        for (i = 1; i <= nu; i++) {
            fpu[i - 1] = 0.0;
        }
        for (i = 1; i <= nv; i++) {
            fpv[i - 1] = 0.0;
        }
        ir = 0;
        nroldu = 0;
        // main loop for the different grid points.
        for (i1 = 1; i1 <= mu; i1++) {
            numu = nru[i1 - 1];
            numu1 = numu + 1;
            nroldv = 0;
            for (i2 = 1; i2 <= mv; i2++) {
                numv = nrv[i2 - 1];
                numv1 = numv + 1;
                ir++;
                // evaluate s(u,v) at the current grid point by making the sum of the
                // cross products of the non-zero b-splines at (u,v), multiplied with
                // the appropriate b-spline coefficients.
                term = 0.0;
                k1 = numu * nv4 + numv;
                for (l1 = 1; l1 <= 4; l1++) {
                    k2 = k1;
                    fac = spu[(i1 - 1) + (l1 - 1) * mu];
                    for (l2 = 1; l2 <= 4; l2++) {
                        k2++;
                        term = term + fac * spv[(i2 - 1) + (l2 - 1) * mv] * c[k2 - 1];
                    }
                    k1 = k1 + nv4;
                }
                // calculate the squared residual at the current grid point.
                term = (r[ir - 1] - term) * (r[ir - 1] - term);
                // adjust the different parameters.
                *fp = *fp + term;
                fpu[numu1 - 1] = fpu[numu1 - 1] + term;
                fpv[numv1 - 1] = fpv[numv1 - 1] + term;
                fac = term * half;
                if (numv != nroldv) {
                    fpv[numv1 - 1] = fpv[numv1 - 1] - fac;
                    fpv[numv - 1] = fpv[numv - 1] + fac;
                }
                nroldv = numv;
                if (numu != nroldu) {
                    fpu[numu1 - 1] = fpu[numu1 - 1] - fac;
                    fpu[numu - 1] = fpu[numu - 1] + fac;
                }
            }
            nroldu = numu;
        }
        return;
    }
    // we determine the matrix (avv) and then we reduce her to an unit
    // upper triangular form (rv) using givens rotations without square
    // roots. we apply the same transformations to the columns of matrix
    // g to obtain the (nv-7) x (nu-6-iop0-iop1) matrix h.
    // we store matrix (rv) into av1 and av2, h into c.
    // the nv7 x nv7 triangular unit upper matrix (rv) has the form
    //           | av1 '     |
    //    (rv) = |     ' av2 |
    //           |  0  '     |
    // with (av2) a nv7 x 4 matrix and (av1) a nv11 x nv11 unit upper
    // triangular matrix of bandwidth 5.
    ncof = nuu * nv7;
    // initialization.
    for (i = 1; i <= ncof; i++) {
        c[i - 1] = 0.0;
    }
    for (i = 1; i <= nv4; i++) {
        av1[(i - 1) + 4 * nv] = 0.0;
        for (j = 1; j <= 4; j++) {
            av1[(i - 1) + (j - 1) * nv] = 0.0;
            av2[(i - 1) + (j - 1) * nv] = 0.0;
        }
    }
    jper = 0;
    nrold = 0;
    for (it = 1; it <= mv; it++) {
        number = nrv[it - 1];
        while (1) {
            if (nrold == number) {
                // fetch a new row of matrix (spv)
                h[4] = 0.0;
                for (j = 1; j <= 4; j++) {
                    h[j - 1] = spv[(it - 1) + (j - 1) * mv];
                }
                // find the appropriate row of g.
                l = it;
                for (j = 1; j <= nuu; j++) {
                    right[j - 1] = q[l - 1];
                    l += mvv;
                }
            } else {
                if (p <= 0.0) {
                    nrold++;
                    continue;
                }
                // fetch a new row of matrix (bv).
                n1 = nrold + 1;
                for (j = 1; j <= 5; j++) {
                    h[j - 1] = bv[(n1 - 1) + (j - 1) * nv] * pinv;
                }
                // find the appropriate row of g.
                for (j = 1; j <= nuu; j++) {
                    right[j - 1] = 0.0;
                }
                if (mv != mvv) {
                    l = mv + n1;
                    for (j = 1; j <= nuu; j++) {
                        right[j - 1] = q[l - 1];
                        l += mvv;
                    }
                }
            }
            // test whether there are non-zero values in the new row of (avv)
            // corresponding to the b-splines n(j;v),j=nv7+1,...,nv4.
            if (nrold < nv11) {
                // rotation into triangle of the new row of (avv), in case the elements
                // corresponding to the b-splines n(j;v),j=nv7+1,...,nv4 are all zero.
                irot = nrold;
                for (i = 1; i <= 5; i++) {
                    irot++;
                    piv = h[i - 1];
                    if (piv == 0.0) continue;
                    // calculate the parameters of the givens transformation.
                    fpgivs(piv, &av1[irot - 1], &co, &si);
                    // apply that transformation to the columns of matrix g.
                    ic = irot;
                    for (j = 1; j <= nuu; j++) {
                        fprota(co, si, &right[j - 1], &c[ic - 1]);
                        ic += nv7;
                    }
                    // apply that transformation to the rows of (avv).
                    if (i == 5) { continue; }
                    i2 = 1;
                    i3 = i + 1;
                    for (j = i3; j <= 5; j++) {
                        i2++;
                        fprota(co, si, &h[j - 1], &av1[(irot - 1) + (i2 - 1) * nv]);
                    }
                }
                // we update the sum of squared residuals.
                for (i = 1; i <= nuu; i++) {
                    *sq = *sq + right[i - 1] * right[i - 1];
                }
            } else {
                if (jper == 0) {
                    // initialize the matrix (av2).
                    jk = nv11 + 1;
                    for (i = 1; i <= 4; i++) {
                        ik = jk;
                        for (j = 1; j <= 5; j++) {
                            if (ik <= 0) { break; }
                            av2[(ik - 1) + (i - 1) * nv] = av1[(ik - 1) + (j - 1) * nv];
                            ik--;
                        }
                        jk++;
                    }
                    jper = 1;
                }
                // if one of the non-zero elements of the new row corresponds to one of
                // the b-splines n(j;v),j=nv7+1,...,nv4, we take account of condition
                // (2) for setting up this row of (avv). the row is stored in h1( the
                // part with respect to av1) and h2 (the part with respect to av2).
                for (i = 1; i <= 4; i++) {
                    h1[i - 1] = 0.0;
                    h2[i - 1] = 0.0;
                }
                h1[4] = 0.0;
                j = nrold - nv11;
                for (i = 1; i <= 5; i++) {
                    j++;
                    l0 = j;
                    while (1) {
                        l1 = l0 - 4;
                        if (l1 <= 0) {
                            h2[l0 - 1] = h2[l0 - 1] + h[i - 1];
                            break;
                        }
                        if (l1 <= nv11) {
                            h1[l1 - 1] = h[i - 1];
                            break;
                        }
                        l0 = l1 - nv11;
                    }
                }
                // rotate the new row of (avv) into triangle.
                if (nv11 > 0) {
                    // rotations with the rows 1,2,...,nv11 of (avv).
                    for (j = 1; j <= nv11; j++) {
                        piv = h1[0];
                        i2 = (nv11 - j < 4) ? (nv11 - j) : 4;
                        if (piv != 0.0) {
                            // calculate the parameters of the givens transformation.
                            fpgivs(piv, &av1[j - 1], &co, &si);
                            // apply that transformation to the columns of matrix g.
                            ic = j;
                            for (i = 1; i <= nuu; i++) {
                                fprota(co, si, &right[i - 1], &c[ic - 1]);
                                ic += nv7;
                            }
                            // apply that transformation to the rows of (avv) with respect to av2.
                            for (i = 1; i <= 4; i++) {
                                fprota(co, si, &h2[i - 1], &av2[(j - 1) + (i - 1) * nv]);
                            }
                            // apply that transformation to the rows of (avv) with respect to av1.
                            if (i2 != 0) {
                                for (i = 1; i <= i2; i++) {
                                    fprota(co, si, &h1[i], &av1[(j - 1) + i * nv]);
                                }
                            }
                        }
                        for (i = 1; i <= i2; i++) {
                            h1[i - 1] = h1[i];
                        }
                        h1[i2] = 0.0;
                    }
                }
                // rotations with the rows nv11+1,...,nv7 of avv.
                for (j = 1; j <= 4; j++) {
                    ij = nv11 + j;
                    if (ij <= 0) continue;
                    piv = h2[j - 1];
                    if (piv == 0.0) continue;
                    // calculate the parameters of the givens transformation.
                    fpgivs(piv, &av2[(ij - 1) + (j - 1) * nv], &co, &si);
                    // apply that transformation to the columns of matrix g.
                    ic = ij;
                    for (i = 1; i <= nuu; i++) {
                        fprota(co, si, &right[i - 1], &c[ic - 1]);
                        ic = ic + nv7;
                    }
                    if (j != 4) {
                        // apply that transformation to the rows of (avv) with respect to av2.
                        j1 = j + 1;
                        for (i = j1; i <= 4; i++) {
                            fprota(co, si, &h2[i - 1], &av2[(ij - 1) + (i - 1) * nv]);
                        }
                    }
                }
                // we update the sum of squared residuals.
                for (i = 1; i <= nuu; i++) {
                    *sq = *sq + right[i - 1] * right[i - 1];
                }
            }
            if (nrold == number) { break; }  // Exit while loop, spin outer for loop
            nrold++;
        }
    }
    // test whether the b-spline coefficients must be determined.
    if (iback != 0) return;
    // backward substitution to obtain the b-spline coefficients as the
    // solution of the linear system    (rv) (cr) (ru)' = h.
    // first step: solve the system  (rv) (c1) = h.
    k = 1;
    for (i = 1; i <= nuu; i++) {
        fpbacp(av1, av2, &c[k - 1], nv7, 4, &c[k - 1], 5, nv);
        k += nv7;
    }
    // second step: solve the system  (cr) (ru)' = (c1).
    k = 0;
    for (j = 1; j <= nv7; j++) {
        k++;
        l = k;
        for (i = 1; i <= nuu; i++) {
            right[i - 1] = c[l - 1];
            l += nv7;
        }
        fpback(au, right, nuu, 5, right, nu);
        l = k;
        for (i = 1; i <= nuu; i++) {
            c[l - 1] = right[i - 1];
            l += nv7;
        }
    }
    // calculate from the conditions (2)-(3)-(4), the remaining b-spline
    // coefficients.
    ncof = nu4 * nv4;
    j = ncof;
    for (l = 1; l <= nv4; l++) {
        q[l - 1] = dr01;
        q[j - 1] = dr11;
        j--;
    }
    i = nv4;
    j = 0;
    if (iop0 != 0) {
        for (l = 1; l <= nv4; l++) {
            i++;
            q[i - 1] = c0[l - 1];
        }
    }
    if (nuu != 0) {
        for (l = 1; l <= nuu; l++) {
            ii = i;
            for (k = 1; k <= nv7; k++) {
                i++;
                j++;
                q[i - 1] = c[j - 1];
            }
            for (k = 1; k <= 3; k++) {
                ii++;
                i++;
                q[i - 1] = q[ii - 1];
            }
        }
    }
    if (iop1 != 0) {
        for (l = 1; l <= nv4; l++) {
            i++;
            q[i - 1] = c1[l - 1];
        }
    }
    for (i = 1; i <= ncof; i++) {
        c[i - 1] = q[i - 1];
    }
    // calculate the quantities
    //   res(i,j) = (r(i,j) - s(u(i),v(j)))**2 , i=1,2,..,mu;j=1,2,..,mv
    //   fp = sumi=1,mu(sumj=1,mv(res(i,j)))
    //   fpu(r) = sum''i(sumj=1,mv(res(i,j))) , r=1,2,...,nu-7
    //                 tu(r+3) <= u(i) <= tu(r+4)
    //   fpv(r) = sumi=1,mu(sum''j(res(i,j))) , r=1,2,...,nv-7
    //                 tv(r+3) <= v(j) <= tv(r+4)
    *fp = 0.0;
    for (i = 1; i <= nu; i++) {
        fpu[i - 1] = 0.0;
    }
    for (i = 1; i <= nv; i++) {
        fpv[i - 1] = 0.0;
    }
    ir = 0;
    nroldu = 0;
    // main loop for the different grid points.
    for (i1 = 1; i1 <= mu; i1++) {
        numu = nru[i1 - 1];
        numu1 = numu + 1;
        nroldv = 0;
        for (i2 = 1; i2 <= mv; i2++) {
            numv = nrv[i2 - 1];
            numv1 = numv + 1;
            ir++;
            // evaluate s(u,v) at the current grid point by making the sum of the
            // cross products of the non-zero b-splines at (u,v), multiplied with
            // the appropriate b-spline coefficients.
            term = 0.0;
            k1 = numu * nv4 + numv;
            for (l1 = 1; l1 <= 4; l1++) {
                k2 = k1;
                fac = spu[(i1 - 1) + (l1 - 1) * mu];
                for (l2 = 1; l2 <= 4; l2++) {
                    k2++;
                    term = term + fac * spv[(i2 - 1) + (l2 - 1) * mv] * c[k2 - 1];
                }
                k1 = k1 + nv4;
            }
            // calculate the squared residual at the current grid point.
            term = (r[ir - 1] - term) * (r[ir - 1] - term);
            // adjust the different parameters.
            *fp = *fp + term;
            fpu[numu1 - 1] = fpu[numu1 - 1] + term;
            fpv[numv1 - 1] = fpv[numv1 - 1] + term;
            fac = term * half;
            if (numv != nroldv) {
                fpv[numv1 - 1] = fpv[numv1 - 1] - fac;
                fpv[numv - 1] = fpv[numv - 1] + fac;
            }
            nroldv = numv;
            if (numu != nroldu) {
                fpu[numu1 - 1] = fpu[numu1 - 1] - fac;
                fpu[numu - 1] = fpu[numu - 1] + fac;
            }
        }
        nroldu = numu;
    }
}


static void
fpinst(const int iopt, const double* t, const int n, const double* c, const int k,
       const double x, const int l, double* tt, int* nn, double* cc, const int nest)
{
    (void)nest; // Unused
    int k1 = k + 1;
    int nk1 = n - k1;
    // the new knots.
    int ll = l + 1;
    int i = n;
    for (int j = ll; j <= n; j++) {
        tt[i] = t[i - 1];
        i--;
    }
    tt[ll - 1] = x;
    for (int j = 1; j <= l; j++) {
        tt[j - 1] = t[j - 1];
    }
    // the new b-spline coefficients
    i = nk1;
    for (int j = l; j <= nk1; j++) {
        cc[i] = c[i - 1];
        i--;
    }
    i = l;
    for (int j = 1; j <= k; j++) {
        int m = i + k1;
        double fac = (x - tt[i - 1]) / (tt[m - 1] - tt[i - 1]);
        int i1 = i - 1;
        cc[i - 1] = fac * c[i - 1] + (1.0 - fac) * c[i1 - 1];
        i--;
    }
    for (int j = 1; j <= i; j++) {
        cc[j - 1] = c[j - 1];
    }
    *nn = n + 1;
    if (iopt == 0) { return; }
    // incorporate the boundary conditions for a periodic spline.
    int nk = *nn - k;
    int nl = nk - k1;
    double per = tt[nk - 1] - tt[k1 - 1];
    i = k1;
    int j = nk;
    if (ll > nl) {
        for (int m = 1; m <= k; m++) {
            int mk = m + nl;
            cc[m - 1] = cc[mk - 1];
            i--;
            j--;
            tt[i - 1] = tt[j - 1] - per;
        }
        return;
    }
    if (ll > (k1 + k)) { return; }

    for (int m = 1; m <= k; m++) {
        int mk = m + nl;
        cc[mk - 1] = cc[m - 1];
        i++;
        j++;
        tt[j - 1] = tt[i - 1] + per;
    }

    return;
}


static void
fpintb(const double *t, int n, double *bint, const int nk1, const double x, const double y)
{
    int i, ia = 0, ib, it, j, j1, k, k1, l, li, lj, lk, l0, mmin;
    double a, ak, arg, b, f, one;
    double aint[6], h[6], h1[6];

    one = 1.0;
    k1 = n - nk1;
    ak = (double)k1;
    k = k1 - 1;

    for (i = 1; i <= nk1; ++i) {
        bint[i-1] = 0.0;
    }
    // the integration limits are arranged in increasing order.
    a = x;
    b = y;
    mmin = 0;
    if (a >= b) {
        if (a == b) { return; }
        a = y;
        b = x;
        mmin = 1;
    }
    if (a < t[k1 - 1]) { a = t[k1 - 1]; }
    if (b > t[nk1]) { b = t[nk1]; }
    if (a > b) { return; }
    // using the expression of gaffney for the indefinite integral of a
    // b-spline we find that
    // bint(j) = (t(j+k+1)-t(j))*(res(j,b)-res(j,a))/(k+1)
    //   where for t(l) <= x < t(l+1)
    //   res(j,x) = 0, j=1,2,...,l-k-1
    //            = 1, j=l+1,l+2,...,nk1
    //            = aint(j+k-l+1), j=l-k,l-k+1,...,l
    //              = sumi((x-t(j+i))*nj+i,k+1-i(x)/(t(j+k+1)-t(j+i)))
    //                i=0,1,...,k
    l = k1;
    l0 = l + 1;
    // set arg = a.
    arg = a;
    for (it = 1; it <= 2; ++it) {
        // search for the knot interval t(l) <= arg < t(l+1).
        while (!((arg < t[l0-1]) || (l == nk1))) {
            l = l0;
            l0 = l + 1;
        }
        // calculation of aint(j-1), j=1,2,...,k+1
        // initialization.
        for (j = 1; j <= k1; j++) { aint[j-1] = 0.0; }
        aint[0] = (arg - t[l - 1]) / (t[l] - t[l - 1]);
        h1[0] = one;
        for (j = 1; j <= k; ++j) {
            // evaluation of the non-zero b-splines of degree j at arg, i.e.
            //     h(i+1) = nl - j + i,j(arg), i=0,1,...,j.
            h[0] = 0.0;
            for (i = 1; i <= j; i++) {
                li = l + i;
                lj = li - j;
                f = h1[i - 1] / (t[li - 1] - t[lj - 1]);
                h[i - 1] = h[i - 1] + f * (t[li - 1] - arg);
                h[i] = f * (arg - t[lj - 1]);
            }
            // updating of the integrals aint
            j1 = j + 1;
            for (i = 1; i <= j1; ++i) {
                li = l + i;
                lj = li - j1;
                aint[i - 1] = aint[i - 1] + h[i - 1] * (arg - t[lj - 1]) / (t[li - 1] - t[lj - 1]);
                h1[i - 1] = h[i - 1];
            }
        }
        if (it == 2) { break; }
        // updating of the integrals bint
        lk = l - k;
        ia = lk;
        for (i = 1; i <= k1; ++i) {
            bint[lk - 1] = -aint[i - 1];
            lk++;
        }
        // set arg = b.
        arg = b;
    }
    // updating of the integrals bint.
    lk = l - k;
    ib = lk - 1;
    for (i = 1; i <= k1; ++i) {
        bint[lk - 1] += aint[i - 1];
        lk++;
    }
    if (ib >= ia) {
        for (i = ia; i <= ib; ++i) {
            bint[i - 1] += one;
        }
    }
    // the scaling factors are taken into account.
    f = one / ak;
    for (i = 1; i <= nk1; ++i) {
        j = i + k1;
        bint[i - 1] = bint[i - 1] * (t[j - 1] - t[i - 1]) * f;
    }
    // the order of the integration limits is taken into account.
    if (mmin != 0) {
        for (i = 1; i <= nk1; ++i) {
            bint[i - 1] = -bint[i - 1];
        }
    }
}


static void
fpknot(const double* x, int m, double* t, int* n, double* fpint, int* nrdata, int* nrint, const int nest, const int istart)
{
    (void)m;  // Unused
    (void)nest;  // Unused
    double an, am, fpmax;
    int ihalf, j, jbegin, jj, jk, jpoint, k, maxbeg = 0, maxpt = 0, next, nrx, number = 0;
    int iserr = 1;

    k = (*n - *nrint - 1) / 2;
    fpmax = 0.0;
    jbegin = istart;

    for (j = 1; j <= *nrint; j++) {
        jpoint = nrdata[j - 1];
        if ((fpmax >= fpint[j - 1]) || (jpoint == 0)) {
            ;
        } else {
            iserr = 0;
            fpmax = fpint[j - 1];
            number = j;
            maxpt = jpoint;
            maxbeg = jbegin;
        }
        jbegin += jpoint + 1;
    }
    // error condition detected, go to exit.
    if (iserr) {
        (*n)++;
        (*nrint)++;
        return;
    }

    // let coincide the new knot t(number+k+1) with a data point x(nrx)
    // inside the old knot interval t(number+k) <= x <= t(number+k+1).
    ihalf = maxpt / 2 + 1;
    nrx = maxbeg + ihalf;
    next = number + 1;
    if (next <= *nrint) {
        // adjust the different parameters.
        for (j = next; j <= *nrint; ++j) {
            jj = next + *nrint - j;
            fpint[jj] = fpint[jj - 1];
            nrdata[jj] = nrdata[jj - 1];
            jk = jj + k;
            t[jk] = t[jk - 1];
        }
    }
    nrdata[number - 1] = ihalf - 1;
    nrdata[next - 1] = maxpt - ihalf;
    am = (double)maxpt;
    an = (double)nrdata[number - 1];
    fpint[number - 1] = fpmax * an / am;
    an = (double)nrdata[next - 1];
    fpint[next - 1] = fpmax * an / am;
    jk = next + k;
    t[jk - 1] = x[nrx - 1];

    (*n)++;
    (*nrint)++;
    return;
}


static void
fpopsp(const int ifsu, const int ifsv, const int ifbu, const int ifbv, const double *u,
       const int mu, const double *v, const int mv, const double *r, const int mr,
       const double r0, const double r1, double *dr, const int *iopt, const int *ider,
       const double *tu, const int nu, const double *tv, const int nv, const int nuest,
       const int nvest, const double p, const double *step, double *c, const int nc,
       double *fp, double *fpu, double *fpv, int *nru, int *nrv, double *wrk, const int lwrk)
{
    (void)lwrk;  // Unused
    double sq, sqq, sq0, sq1, step1, step2, three;
    int i, id0, iop0, iop1, i1, j, l, lau, lav1, lav2, la0, la1, lbu, lbv, lb0;
    int lb1, lc0, lc1, lcs, lq, lri, lsu, lsv, l1, l2, mm, mvnu, number, id1;
    int nr[6];
    double delta[6], drr[6], sum[6], a[6 * 6], g[6];

    three = 3.0;
    lsu = 1;
    lsv = lsu + 4 * mu;
    lri = lsv + 4 * mv;
    mm = (nuest > (mv + nvest)) ? nuest : (mv + nvest);
    lq = lri + mm;
    mvnu = nuest * (mv + nvest - 8);
    lau = lq + mvnu;
    lav1 = lau + 5 * nuest;
    lav2 = lav1 + 6 * nvest;
    lbu = lav2 + 4 * nvest;
    lbv = lbu + 5 * nuest;
    la0 = lbv + 5 * nvest;
    la1 = la0 + 2 * mv;
    lb0 = la1 + 2 * mv;
    lb1 = lb0 + 2 * nvest;
    lc0 = lb1 + 2 * nvest;
    lc1 = lc0 + nvest;
    lcs = lc1 + nvest;

    iop0 = iopt[1];
    iop1 = iopt[2];
    id0 = ider[0];
    id1 = ider[2];
    fpgrsp(ifsu, ifsv, ifbu, ifbv, 0, u, mu, v, mv, r, mr, dr, iop0, iop1, tu, nu, tv, nv,
           p, c, nc, &sq, fp, fpu, fpv, mm, mvnu, &wrk[lsu - 1], &wrk[lsv - 1], &wrk[lri - 1],
           &wrk[lq - 1], &wrk[lau - 1], &wrk[lav1 - 1], &wrk[lav2 - 1], &wrk[lbu - 1],
           &wrk[lbv - 1], &wrk[la0 - 1], &wrk[la1 - 1], &wrk[lb0 - 1], &wrk[lb1 - 1], &wrk[lc0 - 1],
           &wrk[lc1 - 1], &wrk[lcs - 1], nru, nrv);
    sq0 = 0.0;
    sq1 = 0.0;
    if (id0 == 0) { sq0 = (r0 - dr[0]) * (r0 - dr[0]); }
    if (id1 == 0) { sq1 = (r1 - dr[3]) * (r1 - dr[3]); }
    sq = sq + sq0 + sq1;
    // in case all derivative values dr(i) are given (step<=0) or in case
    // we have spline interpolation, we accept this spline as a solution.
    if (sq <= 0.0) { return; }
    if (step[0] <= 0.0 && step[1] <= 0.0) { return; }
    for (i = 0; i < 6; i++) { drr[i] = dr[i]; }
    // number denotes the number of derivative values dr(i) that still must
    // be optimized. let us denote these parameters by g(j),j=1,...,number.
    number = 0;
    if (id0 <= 0) {
        number = 1;
        nr[0] = 1;
        delta[0] = step[0];
    }
    if (!((iop0 == 0) || (ider[1] != 0))) {
        step2 = step[0] * three / (tu[4] - tu[3]);
        nr[number] = 2;
        nr[number + 1] = 3;
        delta[number] = step2;
        delta[number + 1] = step2;
        number = number + 2;
    }
    if (id1 <= 0) {
        number = number + 1;
        nr[number - 1] = 4;
        delta[number - 1] = step[1];
    }
    if (!((iop1 == 0) || (ider[3] != 0))) {
        step2 = step[1] * three / (tu[nu - 1] - tu[nu - 5]);
        nr[number] = 5;
        nr[number + 1] = 6;
        delta[number] = step2;
        delta[number + 1] = step2;
        number = number + 2;
    }
    if (number == 0) { return; }
    // the sum of squared residulas sq is a quadratic polynomial in the
    // parameters g(j). we determine the unknown coefficients of this
    // polymomial by calculating (number+1)*(number+2)/2 different splines
    // according to specific values for g(j).
    int skip_to_110 = 0;
    for (i = 1; i <= number; i++) {
        l = nr[i - 1];
        step1 = delta[i - 1];
        drr[l - 1] = dr[l - 1] + step1;
        fpgrsp(ifsu, ifsv, ifbu, ifbv, 1, u, mu, v, mv, r, mr, drr, iop0, iop1, tu, nu, tv, nv,
               p, c, nc, &sum[i - 1], fp, fpu, fpv, mm, mvnu, &wrk[lsu - 1], &wrk[lsv - 1],
               &wrk[lri - 1], &wrk[lq - 1], &wrk[lau - 1], &wrk[lav1 - 1], &wrk[lav2 - 1],
               &wrk[lbu - 1], &wrk[lbv - 1], &wrk[la0 - 1], &wrk[la1 - 1], &wrk[lb0 - 1],
               &wrk[lb1 - 1], &wrk[lc0 - 1], &wrk[lc1 - 1], &wrk[lcs - 1], nru, nrv);
        if (id0 == 0) { sq0 = (r0 - drr[0]) * (r0 - drr[0]); }
        if (id1 == 0) { sq1 = (r1 - drr[3]) * (r1 - drr[3]); }
        sum[i - 1] = sum[i - 1] + sq0 + sq1;
        drr[l - 1] = dr[l - 1] - step1;
        fpgrsp(ifsu, ifsv, ifbu, ifbv, 1, u, mu, v, mv, r, mr, drr, iop0, iop1, tu, nu, tv, nv,
               p, c, nc, &sqq, fp, fpu, fpv, mm, mvnu, &wrk[lsu - 1], &wrk[lsv - 1], &wrk[lri - 1],
               &wrk[lq - 1], &wrk[lau - 1], &wrk[lav1 - 1], &wrk[lav2 - 1], &wrk[lbu - 1], &wrk[lbv - 1],
               &wrk[la0 - 1], &wrk[la1 - 1], &wrk[lb0 - 1], &wrk[lb1 - 1], &wrk[lc0 - 1], &wrk[lc1 - 1],
               &wrk[lcs - 1], nru, nrv);
        if (id0 == 0) { sq0 = (r0 - drr[0]) * (r0 - drr[0]); }
        if (id1 == 0) { sq1 = (r1 - drr[3]) * (r1 - drr[3]); }
        sqq = sqq + sq0 + sq1;
        drr[l - 1] = dr[l - 1];
        a[(i - 1) + (i - 1) * 6] = (sum[i - 1] + sqq - sq - sq) / (step1 * step1);
        if (a[(i - 1) + (i - 1) * 6] <= 0.0) { skip_to_110 = 1; break;}
        g[i - 1] = (sqq - sum[i - 1]) / (step1 + step1);
    }
    if (!skip_to_110) {
        if (number != 1) {
            for (i = 2; i <= number; i++) {
                l1 = nr[i - 1];
                step1 = delta[i - 1];
                drr[l1 - 1] = dr[l1 - 1] + step1;
                i1 = i - 1;
                for (j = 1; j <= i1; j++) {
                    l2 = nr[j - 1];
                    step2 = delta[j - 1];
                    drr[l2 - 1] = dr[l2 - 1] + step2;
                    fpgrsp(ifsu, ifsv, ifbu, ifbv, 1, u, mu, v, mv, r, mr, drr, iop0, iop1, tu, nu, tv, nv,
                           p, c, nc, &sqq, fp, fpu, fpv, mm, mvnu, &wrk[lsu - 1], &wrk[lsv - 1], &wrk[lri - 1],
                           &wrk[lq - 1], &wrk[lau - 1], &wrk[lav1 - 1], &wrk[lav2 - 1], &wrk[lbu - 1], &wrk[lbv - 1],
                           &wrk[la0 - 1], &wrk[la1 - 1], &wrk[lb0 - 1], &wrk[lb1 - 1], &wrk[lc0 - 1], &wrk[lc1 - 1],
                           &wrk[lcs - 1], nru, nrv);
                    if (id0 == 0) { sq0 = (r0 - drr[0]) * (r0 - drr[0]); }
                    if (id1 == 0) { sq1 = (r1 - drr[3]) * (r1 - drr[3]); }
                    sqq = sqq + sq0 + sq1;
                    a[(i - 1) + (j - 1) * 6] = (sq + sqq - sum[i - 1] - sum[j - 1]) / (step1 * step2);
                    drr[l2 - 1] = dr[l2 - 1];
                }
                drr[l1 - 1] = dr[l1 - 1];
            }
        }
        // the optimal values g(j) are found as the solution of the system
        // d (sq) / d (g(j)) = 0 , j=1,...,number.
        fpsysy(a, number, g);
        for (i = 1; i <= number; i++) {
            l = nr[i - 1];
            dr[l - 1] = dr[l - 1] + g[i - 1];
        }
    }
    fpgrsp(ifsu, ifsv, ifbu, ifbv, 0, u, mu, v, mv, r, mr, dr, iop0, iop1, tu, nu, tv, nv,
           p, c, nc, &sq, fp, fpu, fpv, mm, mvnu, &wrk[lsu - 1], &wrk[lsv - 1], &wrk[lri - 1],
           &wrk[lq - 1], &wrk[lau - 1], &wrk[lav1 - 1], &wrk[lav2 - 1], &wrk[lbu - 1],
           &wrk[lbv - 1], &wrk[la0 - 1], &wrk[la1 - 1], &wrk[lb0 - 1], &wrk[lb1 - 1],
           &wrk[lc0 - 1], &wrk[lc1 - 1], &wrk[lcs - 1], nru, nrv);
    if (id0 == 0) { sq0 = (r0 - dr[0]) * (r0 - dr[0]); }
    if (id1 == 0) { sq1 = (r1 - dr[3]) * (r1 - dr[3]); }
    sq = sq + sq0 + sq1;
}


static void
fporde(const double* x, const double* y, const int m, const int kx, const int ky,
       const double* tx, const int nx, const double* ty, const int ny,
       int* nummer, int* index, const int nreg)
{
    // subroutine fporde sorts the data points (x(i),y(i)),i=1,2,...,m
    // according to the panel tx(l)<=x<tx(l+1),ty(k)<=y<ty(k+1), they belong
    // to. for each panel a stack is constructed  containing the numbers
    // of data points lying inside; index(j),j=1,2,...,nreg points to the
    // first data point in the jth panel while nummer(i),i=1,2,...,m gives
    // the number of the next data point in the panel.

    int kx1 = kx + 1;
    int ky1 = ky + 1;
    int nk1x = nx - kx1;
    int nk1y = ny - ky1;
    int nyy = nk1y - ky;

    // Initialize the index array
    for (int i = 0; i < nreg; i++) { index[i] = 0; }

    for (int im = 1; im <= m; im++) {
        double xi = x[im-1];
        double yi = y[im-1];

        // Find the panel in the x-direction
        int l = kx1;
        int l1 = l + 1;
        while (!((xi < tx[l1-1]) || (l == nk1x))) {
            l = l1;
            l1 = l + 1;
        }

        // Find the panel in the y-direction
        int k = ky1;
        int k1 = k + 1;
        while (!((yi < ty[k1-1]) || k == nk1y)) {
            k = k1;
            k1 = k + 1;
        }

        // Compute the region number
        int num = (l - (kx1))*nyy +  k - ky;

        // Update the nummer and index arrays
        nummer[im-1] = index[num-1];
        index[num-1] = im;
    }
}


static void
fppara(const int iopt, const int idim, const int m, const double *u, const int mx, const double *x, const double *w,
       const double ub, const double ue, const int k, const double s, const int nest, const double tol, const int maxit,
       const int k1, const int k2, int *n, double *t, const int nc, double *c, double *fp, double *fpint, double *z,
       double *a, double *b, double *g, double *q, int *nrdata, int *ier)
{
    (void)mx; // unused
    double acc=0.0, con1, con4, con9, ccos, fac, fpart, fpms, fpold=0.0, fp0=0.0, f1, f2, f3;
    double half, one, p, pinv, piv, p1, p2, p3, rn, ssin, store, term, ui, wi;
    int i, ich1, ich3, it, iter, i1, i2, i3, j, jj, j1, j2, k3, l, l0;
    int mk1, new, nk1, nmax=0, nmin, nplus=0, npl1, nrint, n8;
    double h[7], xi[10];

    // set constants
    one = 1.0;
    con1 = 0.1;
    con9 = 0.9;
    con4 = 0.04;
    half = 0.5;

    /////////////////////////////////////////////////////////////////////////
    //  part 1: determination of the number of knots and their position    //
    //  **************************************************************     //
    //  given a set of knots we compute the least-squares curve sinf(u),   //
    //  and the corresponding sum of squared residuals fp=f(p=inf).        //
    //  if iopt=-1 sinf(u) is the requested curve.                         //
    //  if iopt=0 or iopt=1 we check whether we can accept the knots:      //
    //    if fp <=s we will continue with the current set of knots.        //
    //    if fp > s we will increase the number of knots and compute the   //
    //       corresponding least-squares curve until finally fp<=s.        //
    //    the initial choice of knots depends on the value of s and iopt.  //
    //    if s=0 we have spline interpolation; in that case the number of  //
    //    knots equals nmax = m+k+1.                                       //
    //    if s > 0 and                                                     //
    //      iopt=0 we first compute the least-squares polynomial curve of  //
    //      degree k; n = nmin = 2*k+2                                     //
    //      iopt=1 we start with the set of knots found at the last        //
    //      call of the routine, except for the case that s > fp0; then    //
    //      we compute directly the polynomial curve of degree k.          //
    /////////////////////////////////////////////////////////////////////////

    // determine nmin, the number of knots for polynomial approximation.
    nmin = 2 * k1;

    if (iopt >= 0) {
        // calculation of acc, the absolute tolerance for the root of f(p)=s.
        acc = tol * s;
        // determine nmax, the number of knots for spline interpolation.
        nmax = m + k1;

        if (s > 0.0) {
            // if s>0 our initial choice of knots depends on the value of iopt.
            // if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
            // polynomial curve which is a spline curve without interior knots.
            // if iopt=1 and fp0>s we start computing the least squares spline curve
            // according to the set of knots found at the last call of the routine.
            if ((iopt != 0) && (*n != nmin)) {
                fp0 = fpint[*n - 1];
                fpold = fpint[*n - 2];
                nplus = nrdata[*n - 1];
                if (fp0 <= s) {
                    *n = nmin;
                    fpold = 0.0;
                    nplus = 0;
                    nrdata[0] = m - 2;
                }
            } else {
                *n = nmin;
                fpold = 0.0;
                nplus = 0;
                nrdata[0] = m - 2;
            }
        } else {
            // if s=0, s(u) is an interpolating curve.
            // test whether the required storage space exceeds the available one.
            *n = nmax;
            if (nmax > nest) {
                *ier = 1;
                return;
            }

            // find the position of the interior knots in case of interpolation.
            mk1 = m - k1;
            if (mk1 != 0) {
                k3 = k / 2;
                i = k2;
                j = k3 + 2;
                if ((k3 * 2) == k) {
                    for (l = 1; l <= mk1; l++) {
                        t[i - 1] = (u[j - 1] + u[j - 2]) * half;
                        i++;
                        j++;
                    }
                } else {
                    for (l = 0; l < mk1; l++) {
                        t[i - 1] = u[j - 1];
                        i++;
                        j++;
                    }
                }
            }
        }
    }

restart_iteration:
    // main loop for the different sets of knots. m is a save upper bound
    // for the number of trials.
    for (iter = 1; iter <= m; iter++) {
        if (*n == nmin) {
            *ier = -2;
        }

        // find nrint, the number of knot intervals.
        nrint = *n - nmin + 1;

        // find the position of the additional knots which are needed for
        // the b-spline representation of s(u).
        nk1 = *n - k1;
        i = *n - 1;
        for (j = 0; j < k1; j++) {
            t[j] = ub;
            t[i] = ue;
            i--;
        }

        // compute the b-spline coefficients of the least-squares spline curve
        // sinf(u). the observation matrix a is built up row by row and
        // reduced to upper triangular form by givens transformations.
        // at the same time fp=f(p=inf) is computed.
        *fp = 0.0;

        // initialize the b-spline coefficients and the observation matrix a.
        for (i = 0; i < nc; i++) {
            z[i] = 0.0;
        }
        for (i = 1; i <= nk1; i++) {
            for (j = 1; j <= k1; j++) {
                a[(i - 1) + (j - 1) * nest] = 0.0;
            }
        }

        l = k1;
        jj = 0;
        for (it = 1; it <= m; it++) {
            // fetch the current data point u(it),x(it).
            ui = u[it - 1];
            wi = w[it - 1];
            for (j = 1; j <= idim; j++) {
                jj++;
                xi[j - 1] = x[jj - 1] * wi;
            }

            // search for knot interval t(l) <= ui < t(l+1).
            while ((ui >= t[l]) && (l != nk1)) {
                l++;
            }

            // evaluate the (k+1) non-zero b-splines at ui and store them in q.
            fpbspl(t, *n, k, ui, l, h);
            for (i = 1; i <= k1; i++) {
                q[(it - 1) + (i - 1) * m] = h[i - 1];
                h[i - 1] = h[i - 1] * wi;
            }

            // rotate the new row of the observation matrix into triangle.
            j = l - k1;
            for (i = 1; i <= k1; i++) {
                j++;
                piv = h[i - 1];
                if (piv == 0.0) {
                    continue;
                }

                // calculate the parameters of the givens transformation.
                fpgivs(piv, &a[j - 1], &ccos, &ssin);

                // transformations to right hand side.
                j1 = j;
                for (j2 = 1; j2 <= idim; j2++) {
                    fprota(ccos, ssin, &xi[j2 - 1], &z[j1 - 1]);
                    j1 += *n;
                }

                if (i == k1) { break; }

                i2 = 1;
                i3 = i + 1;
                for (i1 = i3; i1 <= k1; i1++) {
                    i2++;
                    // transformations to left hand side.
                    fprota(ccos, ssin, &h[i1 - 1], &a[(j - 1) + (i2 - 1) * nest]);
                }
            }

            // add contribution of this row to the sum of squares of residual
            // right hand sides.
            for (j2 = 0; j2 < idim; j2++) {
                *fp += pow(xi[j2], 2);
            }
        }

        if (*ier == -2) {
            fp0 = *fp;
        }
        fpint[*n - 1] = fp0;
        fpint[*n - 2] = fpold;
        nrdata[*n - 1] = nplus;

        // backward substitution to obtain the b-spline coefficients.
        j1 = 1;
        for (j2 = 1; j2 <= idim; j2++) {
            fpback(a, &z[j1 - 1], nk1, k1, &c[j1 - 1], nest);
            j1 += *n;
        }

        // test whether the approximation sinf(u) is an acceptable solution.
        if (iopt < 0) { return; }

        fpms = *fp - s;
        if (fabs(fpms) < acc) { return; }

        // if f(p=inf) < s accept the choice of knots.
        if (fpms < 0.0) { break; }

        // if n = nmax, sinf(u) is an interpolating spline curve.
        if (*n == nmax) {
            *ier = -1;
            return;
        }

        // increase the number of knots.
        // if n=nest we cannot increase the number of knots because of
        // the storage capacity limitation.
        if (*n == nest) {
            *ier = 1;
            return;
        }

        // determine the number of knots nplus we are going to add.
        if (*ier == 0) {
            npl1 = nplus * 2;
            rn = (double)nplus;
            if ((fpold - *fp) > acc) {
                npl1 = (int)(rn * fpms / (fpold - *fp));
            }
            // nplus = min0(nplus*2,max0(npl1,nplus/2,1))
            int temp1 = nplus * 2;
            int temp2 = nplus / 2;
            if (temp2 < 1) {
                temp2 = 1;
            }
            if (npl1 > temp2) {
                temp2 = npl1;
            }
            if (temp1 < temp2) {
                nplus = temp1;
            } else {
                nplus = temp2;
            }
        } else {
            nplus = 1;
            *ier = 0;
        }

        fpold = *fp;

        // compute the sum of squared residuals for each knot interval
        // t(j+k) <= u(i) <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = 0.0;
        i = 1;
        l = k2;
        new = 0;
        jj = 0;
        for (it = 1; it <= m; it++) {
            if ((u[it - 1] >= t[l - 1]) && (l <= nk1)) {
                new = 1;
                l++;
            }

            term = 0.0;
            l0 = l - k2;
            for (j2 = 1; j2 <= idim; j2++) {
                fac = 0.0;
                j1 = l0;
                for (j = 1; j <= k1; j++) {
                    j1++;
                    fac += c[j1 - 1] * q[(it - 1) + (j - 1) * m];
                }
                jj++;
                term += pow(w[it - 1] * (fac - x[jj - 1]), 2);
                l0 += *n;
            }

            fpart += term;
            if (new != 0) {
                store = term * half;
                fpint[i - 1] = fpart - store;
                i++;
                fpart = store;
                new = 0;
            }
        }
        fpint[nrint - 1] = fpart;

        for (l = 0; l < nplus; l++) {
            // add a new knot.
            fpknot(u, m, t, n, fpint, nrdata, &nrint, nest, 1);

            // if n=nmax we locate the knots as for interpolation
            if (*n == nmax) {
                mk1 = m - k1;
                if (mk1 == 0) { goto restart_iteration; }
                k3 = k / 2;
                i = k2;
                j = k3 + 2;
                if ((k3 * 2) == k) {
                    for (int lp = 1; lp <= mk1; lp++) {
                        t[i - 1] = (u[j - 1] + u[j - 2]) * half;
                        i++;
                        j++;
                    }
                } else {
                    for (int lp = 1; lp <= mk1; lp++) {
                        t[i - 1] = u[j - 1];
                        i++;
                        j++;
                    }
                }
                goto restart_iteration;
            }

            // test whether we cannot further increase the number of knots.
            if (*n == nest) {
                break;
            }
        }
    }

    // test whether the least-squares kth degree polynomial curve is a
    // solution of our approximation problem.
    if (*ier == -2) {
        return;
    }

    /////////////////////////////////////////////////////////////////////////
    //  part 2: determination of the smoothing spline curve sp(u).          //
    //  **********************************************************          //
    //  we have determined the number of knots and their position.          //
    //  we now compute the b-spline coefficients of the smoothing curve     //
    //  sp(u). the observation matrix a is extended by the rows of matrix   //
    //  b expressing that the kth derivative discontinuities of sp(u) at    //
    //  the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     //
    //  ponding weights of these additional rows are set to 1/p.            //
    //  iteratively we then have to determine the value of p such that f(p),//
    //  the sum of squared residuals be = s. we already know that the least //
    //  squares kth degree polynomial curve corresponds to p=0, and that    //
    //  the least-squares spline curve corresponds to p=infinity. the       //
    //  iteration process which is proposed here, makes use of rational     //
    //  interpolation. since f(p) is a convex and strictly decreasing       //
    //  function of p, it can be approximated by a rational function        //
    //  r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  //
    //  ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      //
    //  to calculate the new value of p such that r(p)=s. convergence is    //
    //  guaranteed by taking f1>0 and f3<0.                                 //
    /////////////////////////////////////////////////////////////////////////

    // evaluate the discontinuity jump of the kth derivative of the
    // b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
    fpdisc(t, *n, k2, b, nest);

    // initial value for p.
    p1 = 0.0;
    f1 = fp0 - s;
    p3 = -one;
    f3 = fpms;
    p = 0.0;
    for (i = 0; i < nk1; i++) {
        p += a[i];
    }
    rn = (double)nk1;
    p = rn / p;
    ich1 = 0;
    ich3 = 0;
    n8 = *n - nmin;

    // iteration process to find the root of f(p) = s.
    for (iter = 1; iter <= maxit; iter++) {
        // the rows of matrix b with weight 1/p are rotated into the
        // triangularised observation matrix a which is stored in g.
        pinv = one / p;
        for (i = 0; i < nc; i++) {
            c[i] = z[i];
        }
        for (i = 1; i <= nk1; i++) {
            g[(i - 1) + (k2 - 1) * nest] = 0.0;
            for (j = 1; j <= k1; j++) {
                g[(i - 1) + (j - 1) * nest] = a[(i - 1) + (j - 1) * nest];
            }
        }

        for (it = 1; it <= n8; it++) {
            // the row of matrix b is rotated into triangle by givens transformation
            for (i = 1; i <= k2; i++) {
                h[i - 1] = b[(it - 1) + (i - 1) * nest] * pinv;
            }
            for (j = 0; j < idim; j++) {
                xi[j] = 0.0;
            }

            for (j = it; j <= nk1; j++) {
                piv = h[0];

                // calculate the parameters of the givens transformation.
                fpgivs(piv, &g[j - 1], &ccos, &ssin);

                // transformations to right hand side.
                j1 = j;
                for (j2 = 1; j2 <= idim; j2++) {
                    fprota(ccos, ssin, &xi[j2 - 1], &c[j1 - 1]);
                    j1 += *n;
                }

                if (j == nk1) {
                    break;
                }

                i2 = k1;
                if (j > n8) {
                    i2 = nk1 - j;
                }
                for (i = 1; i <= i2; i++) {
                    // transformations to left hand side.
                    i1 = i + 1;
                    fprota(ccos, ssin, &h[i], &g[(j - 1) + i * nest]);
                    h[i - 1] = h[i];
                }
                h[i2] = 0.0;
            }
        }

        // backward substitution to obtain the b-spline coefficients.
        j1 = 1;
        for (j2 = 1; j2 <= idim; j2++) {
            fpback(g, &c[j1 - 1], nk1, k2, &c[j1 - 1], nest);
            j1 += *n;
        }

        // computation of f(p).
        *fp = 0.0;
        l = k2;
        jj = 0;
        for (it = 1; it <= m; it++) {
            if ((u[it - 1] >= t[l - 1]) && (l <= nk1)) {
                l++;
            }

            l0 = l - k2;
            term = 0.0;
            for (j2 = 1; j2 <= idim; j2++) {
                fac = 0.0;
                j1 = l0;
                for (j = 1; j <= k1; j++) {
                    j1++;
                    fac += c[j1 - 1] * q[(it - 1) + (j - 1) * m];
                }
                jj++;
                term += (fac - x[jj - 1]) * (fac - x[jj - 1]);
                l0 += *n;
            }
            *fp += term * w[it - 1] * w[it - 1];
        }

        // test whether the approximation sp(u) is an acceptable solution.
        fpms = *fp - s;
        if (fabs(fpms) < acc) { return; }

        // test whether the maximal number of iterations is reached.
        if (iter == maxit) {
            *ier = 3;
            return;
        }

        // carry out one more step of the iteration process.
        p2 = p;
        f2 = fpms;

        if (ich3 == 0) {
            if ((f2 - f3) > acc) {
                if (f2 < 0.0) {
                    ich3 = 1;
                }
            } else {
                // our initial choice of p is too large.
                p3 = p2;
                f3 = f2;
                p = p * con4;
                if (p <= p1) {
                    p = p1 * con9 + p2 * con1;
                }
                continue;
            }
        }

        if (ich1 == 0) {
            if ((f1 - f2) > acc) {
                if (f2 > 0.0) {
                    ich1 = 1;
                }
            } else {
                // our initial choice of p is too small
                p1 = p2;
                f1 = f2;
                p = p / con4;
                if (p3 >= 0.0) {
                    if (p >= p3) {
                        p = p2 * con1 + p3 * con9;
                    }
                }
                continue;
            }
        }

        // test whether the iteration process proceeds as theoretically
        // expected.
        if ((f2 >= f1) || (f2 <= f3)) {
            *ier = 2;
            return;
        }

        // find the new value for p.
        p = fprati(&p1, &f1, &p2, &f2, &p3, &f3);
    }
}


static void
fpperi(const int iopt, const double *x, const double *y, const double *w, const int m,
       const int k, const double s, const int nest, const double tol, const int maxit,
       const int k1, const int k2, int *n, double *t, double *c, double *fp, double *fpint,
       double *z, double *a1, double *a2, double *b, double *g1, double *g2, double *q,
       int *nrdata, int *ier)
{
    // Subroutine fpperi determines the periodic spline approximation.
    double acc=0.0, cos, c1, d1, fpart, fpms, fpold=0.0, fp0=0.0, f1, f2, f3, p, per, pinv, piv;
    double p1, p2, p3, sin, store, term, wi, xi, yi, rn, one, con1, con4, con9, half;
    int i, ich1, ich3, ij, ik, it, iter, i1, i2, i3, j, jk, jper, j1, j2, kk, kk1;
    int l, l0, l1, l5, mm, m1, new, nk1, nk2, nmax=0, nmin, nplus=0, npl1;
    int nrint=0, n10=0, n11=0, n7=0, n8=0;
    double h[6], h1[7], h2[6];

    // set constants
    one = 1.0;
    con1 = 0.1;
    con9 = 0.9;
    con4 = 0.04;
    half = 0.5;

    /////////////////////////////////////////////////////////////////////////
    // part 1: determination of the number of knots and their position.    //
    // **************************************************************      //
    // given a set of knots we compute the least-squares periodic spline   //
    // sinf(x). if the sum f(p=inf) <= s we accept the choice of knots.    //
    // the initial choice of knots depends on the value of s and iopt.     //
    //   if s=0 we have spline interpolation; in that case the number of   //
    //   knots equals nmax = m+2*k.                                        //
    //   if s > 0 and                                                      //
    //     iopt=0 we first compute the least-squares polynomial of         //
    //     degree k; n = nmin = 2*k+2. since s(x) must be periodic we      //
    //     find that s(x) is a constant function.                          //
    //     iopt=1 we start with the set of knots found at the last         //
    //     call of the routine, except for the case that s > fp0; then     //
    //     we compute directly the least-squares periodic polynomial.      //
    /////////////////////////////////////////////////////////////////////////

    m1 = m - 1;
    kk = k;
    kk1 = k1;
    nmin = 2 * k1;

    // determine the length of the period of s(x).
    per = x[m - 1] - x[0];

    // This goto can essentially be cleaned but causes a large code duplication.
    if (iopt < 0) { goto restart_iteration; }

    // calculation of acc, the absolute tolerance for the root of f(p)=s.
    acc = tol * s;

    // determine nmax, the number of knots for periodic spline interpolation
    nmax = m + 2 * k;

    if ((s > 0.0) || (nmax == nmin)) {
        // if s > 0 our initial choice depends on the value of iopt.
        // if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
        // periodic polynomial. (i.e. a constant function).
        // if iopt=1 and fp0>s we start computing the least-squares periodic
        // spline according the set of knots found at the last call of the
        // routine.
        if ((iopt != 0) && (*n != nmin)) {
            fp0 = fpint[*n - 1];
            fpold = fpint[*n - 2];
            nplus = nrdata[*n - 1];
            if (fp0 > s) { goto restart_iteration; }

        }
        // 35

        // the case that s(x) is a constant function is treated separately.
        // find the least-squares constant c1 and compute fp0 at the same time.
        fp0 = 0.0;
        d1 = 0.0;
        c1 = 0.0;
        for (it = 1; it <= m1; it++) {
            wi = w[it - 1];
            yi = y[it - 1] * wi;
            fpgivs(wi, &d1, &cos, &sin);
            fprota(cos, sin, &yi, &c1);
            fp0 += yi * yi;
        }
        c1 = c1 / d1;

        // test whether that constant function is a solution of our problem.
        fpms = fp0 - s;
        if ((fpms < acc) || (nmax == nmin)) {
            *ier = -2;
            // the least-squares constant function c1 is a solution of our problem.
            // a constant function is a spline of degree k with all b-spline
            // coefficients equal to that constant c1.
            for (i = 1; i <= k1; i++) {
                rn = k1 - i;
                t[i - 1] = x[0] - rn * per;
                c[i - 1] = c1;
                j = i + k1;
                rn = i - 1;
                t[j - 1] = x[m - 1] + rn * per;
            }
            *n = nmin;
            *fp = fp0;
            fpint[*n - 1] = fp0;
            fpint[*n - 2] = 0.0;
            nrdata[*n - 1] = 0;
            return;
        }
        fpold = fp0;

        // test whether the required storage space exceeds the available one.
        if (nmin >= nest) {
            *ier = 1;
            return;
        }

        // start computing the least-squares periodic spline with one interior knot.
        nplus = 1;
        *n = nmin + 1;
        mm = (m + 1) / 2;
        t[k2 - 1] = x[mm - 1];
        nrdata[0] = mm - 2;
        nrdata[1] = m1 - mm;
    } else {
        // if s = 0, s(x) is an interpolating spline.
        *n = nmax;

        // test whether the required storage space exceeds the available one.
        if (*n > nest) {
            *ier = 1;
            return;
        }

        // find the position of the interior knots in case of interpolation.
        if (k % 2 == 0) {
            // k is even
            for (i = 2; i <= m1; i++) {
                j = i + k;
                t[j - 1] = (x[i - 1] + x[i - 2]) * half;
            }
        } else {
            // k is odd
            for (i = 2; i <= m1; i++) {
                j = i + k;
                t[j - 1] = x[i - 1];
            }
            if (s <= 0.0) {
                kk = k - 1;
                kk1 = k;
                if (kk <= 0) {
                    t[0] = t[m - 1] - per;
                    t[1] = x[0];
                    t[m] = x[m - 1];
                    t[m + 1] = t[2] + per;
                    for (i = 1; i <= m1; i++) {
                        c[i - 1] = y[i - 1];
                    }
                    c[m - 1] = c[0];
                    *fp = 0.0;
                    fp0 = 0.0;
                    fpint[*n - 1] = fp0;
                    fpint[*n - 2] = 0.0;
                    nrdata[*n - 1] = 0;
                    *ier = -1;
                    return;
                }
            }
        }
    }

restart_iteration:
    // main loop for the different sets of knots. m is a save upper
    // bound for the number of trials.
    for (iter = 1; iter <= m; iter++) {
        // find nrint, the number of knot intervals.
        nrint = *n - nmin + 1;

        // find the position of the additional knots which are needed for
        // the b-spline representation of s(x). if we take
        //     t(k+1) = x(1), t(n-k) = x(m)
        //     t(k+1-j) = t(n-k-j) - per, j=1,2,...k
        //     t(n-k+j) = t(k+1+j) + per, j=1,2,...k
        // then s(x) is a periodic spline with period per if the b-spline
        // coefficients satisfy the following conditions
        //     c(n7+j) = c(j), j=1,...k   (**)   with n7=n-2*k-1.
        t[k1 - 1] = x[0];
        nk1 = *n - k1;
        nk2 = nk1 + 1;
        t[nk2 - 1] = x[m - 1];

        for (j = 1; j <= k; j++) {
            i1 = nk2 + j;
            i2 = nk2 - j;
            j1 = k1 + j;
            j2 = k1 - j;
            t[i1 - 1] = t[j1 - 1] + per;
            t[j2 - 1] = t[i2 - 1] - per;
        }

        // compute the b-spline coefficients c(j),j=1,...n7 of the least-squares
        // periodic spline sinf(x). the observation matrix a is built up row
        // by row while taking into account condition (**) and is reduced to
        // triangular form by givens transformations.
        // at the same time fp=f(p=inf) is computed.
        // the n7 x n7 triangularised upper matrix a has the form
        //           ! a1 '    !
        //       a = !    ' a2 !
        //           ! 0  '    !
        // with a2 a n7 x k matrix and a1 a n10 x n10 upper triangular
        // matrix of bandwidth k+1 ( n10 = n7-k).

        // initialization.
        for (i = 1; i <= nk1; i++) {
            z[i - 1] = 0.0;
            for (j = 1; j <= kk1; j++) {
                a1[(i - 1) + (j - 1) * nest] = 0.0;
            }
        }

        n7 = nk1 - k;
        n10 = n7 - kk;
        jper = 0;
        *fp = 0.0;
        l = k1;

        for (it = 1; it <= m1; it++) {
            // fetch the current data point x(it),y(it)
            xi = x[it - 1];
            wi = w[it - 1];
            yi = y[it - 1] * wi;

            // search for knot interval t(l) <= xi < t(l+1).
            while (xi >= t[l]) {
                l = l + 1;
            }

            // evaluate the (k+1) non-zero b-splines at xi and store them in q.
            fpbspl(t, *n, k, xi, l, h);
            for (i = 1; i <= k1; i++) {
                q[(it - 1) + (i - 1) * m] = h[i - 1];
                h[i - 1] = h[i - 1] * wi;
            }

            l5 = l - k1;

            // test whether the b-splines nj,k+1(x),j=1+n7,...nk1 are all zero at xi
            if (l5 < n10) {
                // rotation of the new row of the observation matrix into
                // triangle in case the b-splines nj,k+1(x),j=n7+1,...n-k-1 are all zero at xi.
                j = l5;
                for (i = 1; i <= kk1; i++) {
                    j++;
                    piv = h[i - 1];
                    if (piv == 0.0) { continue; }

                    // calculate the parameters of the givens transformation.
                    fpgivs(piv, &a1[j - 1], &cos, &sin);

                    // transformations to right hand side.
                    fprota(cos, sin, &yi, &z[j - 1]);

                    if (i == kk1) { break; }

                    i2 = 1;
                    i3 = i + 1;

                    // transformations to left hand side.
                    for (i1 = i3; i1 <= kk1; i1++) {
                        i2++;
                        fprota(cos, sin, &h[i1 - 1], &a1[(j - 1) + (i2 - 1) * nest]);
                    }
                }

                // add contribution of this row to the sum of squares of residual right hand sides.
                *fp = *fp + yi * yi;
            } else {
                if (jper == 0) {
                    // initialize the matrix a2.
                    for (i = 1; i <= n7; i++) {
                        for (j = 1; j <= kk; j++) {
                            a2[(i - 1) + (j - 1) * nest] = 0.0;
                        }
                    }

                    jk = n10 + 1;
                    for (i = 1; i <= kk; i++) {
                        ik = jk;
                        for (j = 1; j <= kk1; j++) {
                            if (ik <= 0) { break; }
                            a2[(ik - 1) + (i - 1) * nest] = a1[(ik - 1) + (j - 1) * nest];
                            ik--;
                        }
                        jk++;
                    }
                    jper = 1;
                }

                // if one of the b-splines nj,k+1(x),j=n7+1,...nk1 is not zero at xi
                // we take account of condition (**) for setting up the new row
                // of the observation matrix a. this row is stored in the arrays h1
                // (the part with respect to a1) and h2 (the part with respect to a2).
                for (i = 1; i <= kk; i++) {
                    h1[i - 1] = 0.0;
                    h2[i - 1] = 0.0;
                }
                h1[kk1 - 1] = 0.0;

                j = l5 - n10;
                for (i = 1; i <= kk1; i++) {
                    j++;
                    l0 = j;

                    while (1) {
                        l1 = l0 - kk;
                        if (l1 <= 0) {
                            h2[l0 - 1] = h2[l0 - 1] + h[i - 1];
                            break;
                        }
                        if (l1 <= n10) {
                            h1[l1 - 1] = h[i - 1];
                            break;
                        }
                        l0 = l1 - n10;
                    }
                }

                // rotate the new row of the observation matrix into triangle by givens transformations.
                if (n10 > 0) {
                    // rotation with the rows 1,2,...n10 of matrix a.
                    for (j = 1; j <= n10; j++) {
                        piv = h1[0];

                        if (piv == 0.0) {
                            for (i = 1; i <= kk; i++) {
                                h1[i - 1] = h1[i];
                            }
                            h1[kk1 - 1] = 0.0;
                            continue;
                        }

                        // calculate the parameters of the givens transformation.
                        fpgivs(piv, &a1[j - 1], &cos, &sin);

                        // transformation to the right hand side.
                        fprota(cos, sin, &yi, &z[j - 1]);

                        // transformations to the left hand side with respect to a2.
                        for (i = 1; i <= kk; i++) {
                            fprota(cos, sin, &h2[i - 1], &a2[(j - 1) + (i - 1) * nest]);
                        }

                        if (j == n10) { break; }

                        i2 = (n10 - j < kk) ? (n10 - j) : kk;

                        // transformations to the left hand side with respect to a1.
                        for (i = 1; i <= i2; i++) {
                            i1 = i + 1;
                            fprota(cos, sin, &h1[i], &a1[(j - 1) + i * nest]);
                            h1[i - 1] = h1[i];
                        }
                        h1[i1 - 1] = 0.0;
                    }
                }

                // rotation with the rows n10+1,...n7 of matrix a.
                for (j = 1; j <= kk; j++) {
                    ij = n10 + j;
                    if (ij <= 0) { continue; }

                    piv = h2[j - 1];
                    if (piv == 0.0) { continue; }

                    // calculate the parameters of the givens transformation.
                    fpgivs(piv, &a2[(ij - 1) + (j - 1) * nest], &cos, &sin);

                    // transformations to right hand side.
                    fprota(cos, sin, &yi, &z[ij - 1]);

                    if (j == kk) { break; }

                    j1 = j + 1;

                    // transformations to left hand side.
                    for (i = j1; i <= kk; i++) {
                        fprota(cos, sin, &h2[i - 1], &a2[(ij - 1) + (i - 1) * nest]);
                    }
                }

                // add contribution of this row to the sum of squares of residual right hand sides.
                *fp = *fp + yi * yi;
            }
        }

        fpint[*n - 1] = fp0;
        fpint[*n - 2] = fpold;
        nrdata[*n - 1] = nplus;

        // backward substitution to obtain the b-spline coefficients c(j),j=1,.n
        fpbacp(a1, a2, z, n7, kk, c, kk1, nest);

        // calculate from condition (**) the coefficients c(j+n7),j=1,2,...k.
        for (i = 1; i <= k; i++) {
            j = i + n7;
            c[j - 1] = c[i - 1];
        }

        if (iopt < 0) { return; }

        // test whether the approximation sinf(x) is an acceptable solution.
        fpms = *fp - s;
        if (fabs(fpms) < acc) { return; }

        // if f(p=inf) < s accept the choice of knots.
        if (fpms < 0.0) { break; }

        // if n=nmax, sinf(x) is an interpolating spline.
        if (*n == nmax) {
            *ier = -1;
            return;
        }

        // increase the number of knots.
        // if n=nest we cannot increase the number of knots because of the storage capacity limitation.
        if (*n == nest) {
            *ier = 1;
            return;
        }

        // determine the number of knots nplus we are going to add.
        npl1 = nplus * 2;
        rn = nplus;
        if ((fpold - *fp) > acc) {
            npl1 = (int)(rn * fpms / (fpold - *fp));
        }
        // min0(nplus*2,max0(npl1,nplus/2,1))
        int temp1 = npl1;
        int temp2 = nplus / 2;
        if (temp2 > temp1) { temp1 = temp2; }
        if (1 > temp1) { temp1 = 1; }
        int temp3 = nplus * 2;
        if (temp1 < temp3) { temp3 = temp1; }
        nplus = temp3;

        fpold = *fp;

        // compute the sum(wi*(yi-s(xi))**2) for each knot interval
        // t(j+k) <= xi <= t(j+k+1) and store it in fpint(j),j=1,2,...nrint.
        fpart = 0.0;
        i = 1;
        l = k1;
        new = 0;

        for (it = 1; it <= m1; it++) {
            if (x[it - 1] >= t[l - 1]) {
                new = 1;
                l++;
            }

            term = 0.0;
            l0 = l - k2;
            for (j = 1; j <= k1; j++) {
                l0++;
                term += c[l0 - 1] * q[(it - 1) + (j - 1) * m];
            }

            term = pow(w[it - 1] * (term - y[it - 1]), 2);
            fpart += term;

            if (new == 0) { continue; }

            if (l <= k2) {
                fpint[nrint - 1] = term;
                new = 0;
                continue;
            }

            store = term * half;
            fpint[i - 1] = fpart - store;
            i++;
            fpart = store;
            new = 0;
        }

        fpint[nrint - 1] = fpint[nrint - 1] + fpart;

        for (l = 1; l <= nplus; l++) {
            // add a new knot
            fpknot(x, m, t, n, fpint, nrdata, &nrint, nest, 1);

            // if n=nmax we locate the knots as for interpolation.
            if (*n == nmax) {
                // find the position of the interior knots in case of interpolation.
                if (k % 2 == 0) {
                    // k is even
                    for (i = 2; i <= m1; i++) {
                        j = i + k;
                        t[j - 1] = (x[i - 1] + x[i - 2]) * half;
                    }
                    goto restart_iteration;
                } else {
                    // k is odd
                    for (i = 2; i <= m1; i++) {
                        j = i + k;
                        t[j - 1] = x[i - 1];
                    }
                    if (s > 0.0) { goto restart_iteration; }
                    kk = k - 1;
                    kk1 = k;
                    if (kk > 0) { goto restart_iteration; }
                    t[0] = t[m - 1] - per;
                    t[1] = x[0];
                    t[m] = x[m - 1];
                    t[m + 1] = t[2] + per;
                    for (i = 1; i <= m1; i++) {
                        c[i - 1] = y[i - 1];
                    }
                    c[m - 1] = c[0];
                    *fp = 0.0;
                    fp0 = 0.0;
                    fpint[*n - 1] = fp0;
                    fpint[*n - 2] = 0.0;
                    nrdata[*n - 1] = 0;
                    *ier = -1;
                    return;
                }
            }

            // test whether we cannot further increase the number of knots.
            if (*n == nest) { break; }
        }

        // restart the computations with the new set of knots.
    }

    /////////////////////////////////////////////////////////////////////////
    // part 2: determination of the smoothing periodic spline sp(x).       //
    // *************************************************************       //
    // we have determined the number of knots and their position.          //
    // we now compute the b-spline coefficients of the smoothing spline    //
    // sp(x). the observation matrix a is extended by the rows of matrix   //
    // b expressing that the kth derivative discontinuities of sp(x) at    //
    // the interior knots t(k+2),...t(n-k-1) must be zero. the corres-     //
    // ponding weights of these additional rows are set to 1/sqrt(p).      //
    // iteratively we then have to determine the value of p such that      //
    // f(p)=sum(w(i)*(y(i)-sp(x(i)))**2) be = s. we already know that      //
    // the least-squares constant function corresponds to p=0, and that    //
    // the least-squares periodic spline corresponds to p=infinity. the    //
    // iteration process which is proposed here, makes use of rational     //
    // interpolation. since f(p) is a convex and strictly decreasing       //
    // function of p, it can be approximated by a rational function        //
    // r(p) = (u*p+v)/(p+w). three values of p(p1,p2,p3) with correspond-  //
    // ing values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used      //
    // to calculate the new value of p such that r(p)=s. convergence is    //
    // guaranteed by taking f1>0 and f3<0.                                 //
    /////////////////////////////////////////////////////////////////////////

    // evaluate the discontinuity jump of the kth derivative of the
    // b-splines at the knots t(l),l=k+2,...n-k-1 and store in b.
    fpdisc(t, *n, k2, b, nest);

    // initial value for p.
    p1 = 0.0;
    f1 = fp0 - s;
    p3 = -one;
    f3 = fpms;
    n11 = n10 - 1;
    n8 = n7 - 1;
    p = 0.0;
    l = n7;
    int skip_354 = 0;
    for (i = 1; i <= k; i++) {
        j = k + 1 - i;
        p = p + a2[(l - 1) + (j - 1) * nest];
        l--;
        if (l == 0) {
            skip_354 = 1;  // go over both loops
            break;
        }
    }

    if (!skip_354) {
        // 354
        for (i = 1; i <= n10; i++) {
            p = p + a1[i - 1];
        }
    }
    rn = n7;
    p = rn / p;
    ich1 = 0;
    ich3 = 0;

    // iteration process to find the root of f(p) = s.
    for (iter = 1; iter <= maxit; iter++) {
        // form the matrix g as the matrix a extended by the rows of matrix b.
        // the rows of matrix b with weight 1/p are rotated into
        // the triangularised observation matrix a.
        // after triangularisation our n7 x n7 matrix g takes the form
        //           ! g1 '    !
        //       g = !    ' g2 !
        //           ! 0  '    !
        // with g2 a n7 x (k+1) matrix and g1 a n11 x n11 upper triangular
        // matrix of bandwidth k+2. ( n11 = n7-k-1)
        pinv = one / p;

        // store matrix a into g
        for (i = 1; i <= n7; i++) {
            c[i - 1] = z[i - 1];
            g1[(i - 1) + (k1 - 1) * nest] = a1[(i - 1) + (k1 - 1) * nest];
            g1[(i - 1) + (k2 - 1) * nest] = 0.0;
            g2[i - 1] = 0.0;
            for (j = 1; j <= k; j++) {
                g1[(i - 1) + (j - 1) * nest] = a1[(i - 1) + (j - 1) * nest];
                g2[(i - 1) + j * nest] = a2[(i - 1) + (j - 1) * nest];
            }
        }

        l = n10;
        for (j = 1; j <= k1; j++) {
            if (l <= 0) { break; }
            g2[l - 1] = a1[(l - 1) + (j - 1) * nest];
            l--;
        }

        for (it = 1; it <= n8; it++) {
            // fetch a new row of matrix b and store it in the arrays h1 (the part
            // with respect to g1) and h2 (the part with respect to g2).
            yi = 0.0;
            for (i = 1; i <= k1; i++) {
                h1[i - 1] = 0.0;
                h2[i - 1] = 0.0;
            }
            h1[k2 - 1] = 0.0;

            if (it <= n11) {
                l = it;
                l0 = it;
                for (j = 1; j <= k2; j++) {
                    if (l0 == n10) {
                        l0 = 1;
                        for (l1 = j; l1 <= k2; l1++) {
                            h2[l0 - 1] = b[(it - 1) + (l1 - 1) * nest] * pinv;
                            l0++;
                        }
                        break;
                    }
                    h1[j - 1] = b[(it - 1) + (j - 1) * nest] * pinv;
                    l0++;
                }
            } else {
                l = 1;
                i = it - n10;
                for (j = 1; j <= k2; j++) {
                    i++;
                    l0 = i;

                    while (1) {
                        l1 = l0 - k1;
                        if (l1 <= 0) {
                            h2[l0 - 1] = h2[l0 - 1] + b[(it - 1) + (j - 1) * nest] * pinv;
                            break;
                        }
                        if (l1 <= n11) {
                            h1[l1 - 1] = b[(it - 1) + (j - 1) * nest] * pinv;
                            break;
                        }
                        l0 = l1 - n11;
                    }
                }
            }

            if (n11 > 0) {
                // rotate this row into triangle by givens transformations without square roots.
                // rotation with the rows l,l+1,...n11.
                for (j = l; j <= n11; j++) {
                    piv = h1[0];

                    // calculate the parameters of the givens transformation.
                    fpgivs(piv, &g1[(j - 1) + 0 * nest], &cos, &sin);

                    // transformation to right hand side.
                    fprota(cos, sin, &yi, &c[j - 1]);

                    // transformation to the left hand side with respect to g2.
                    for (i = 1; i <= k1; i++) {
                        fprota(cos, sin, &h2[i - 1], &g2[(j - 1) + (i - 1) * nest]);
                    }

                    if (j == n11) { break; }

                    i2 = (n11 - j < k1) ? (n11 - j) : k1;

                    // transformation to the left hand side with respect to g1.
                    for (i = 1; i <= i2; i++) {
                        i1 = i + 1;
                        fprota(cos, sin, &h1[i1 - 1], &g1[(j - 1) + (i1 - 1) * nest]);
                        h1[i - 1] = h1[i1 - 1];
                    }
                    h1[i1 - 1] = 0.0;
                }
            }

            // rotation with the rows n11+1,...n7
            for (j = 1; j <= k1; j++) {
                ij = n11 + j;
                if (ij <= 0) { continue; }

                piv = h2[j - 1];

                // calculate the parameters of the givens transformation
                fpgivs(piv, &g2[(ij - 1) + (j - 1) * nest], &cos, &sin);

                // transformation to the right hand side.
                fprota(cos, sin, &yi, &c[ij - 1]);

                if (j == k1) { break; }

                j1 = j + 1;

                // transformation to the left hand side.
                for (i = j1; i <= k1; i++) {
                    fprota(cos, sin, &h2[i - 1], &g2[(ij - 1) + (i - 1) * nest]);
                }
            }
        }

        // backward substitution to obtain the b-spline coefficients c(j),j=1,2,...n7 of sp(x).
        fpbacp(g1, g2, c, n7, k1, c, k2, nest);

        // calculate from condition (**) the b-spline coefficients c(n7+j),j=1,2,...k.
        for (i = 1; i <= k; i++) {
            j = i + n7;
            c[j - 1] = c[i - 1];
        }

        // computation of f(p).
        *fp = 0.0;
        l = k1;
        for (it = 1; it <= m1; it++) {
            if (x[it - 1] >= t[l - 1]) {
                l++;
            }

            l0 = l - k2;
            term = 0.0;
            for (j = 1; j <= k1; j++) {
                l0++;
                term = term + c[l0 - 1] * q[(it - 1) + (j - 1) * m];
            }
            *fp += pow(w[it - 1] * (term - y[it - 1]), 2);
        }

        // test whether the approximation sp(x) is an acceptable solution.
        fpms = *fp - s;
        if (fabs(fpms) < acc) {
            return;
        }

        // test whether the maximal number of iterations is reached.
        if (iter == maxit) {
            *ier = 3;
            return;
        }

        // carry out one more step of the iteration process.
        p2 = p;
        f2 = fpms;

        if (ich3 == 0) {
            if ((f2 - f3) > acc) {
                if (f2 < 0.0) { ich3 = 1; }
            } else {
                // our initial choice of p is too large.
                p3 = p2;
                f3 = f2;
                p = p * con4;
                if (p <= p1) { p = p1 * con9 + p2 * con1; }
                continue;
            }
        }

        if (ich1 == 0) {
            if ((f1 - f2) > acc) {
                if (f2 > 0.0) { ich1 = 1; }
            } else {
                // our initial choice of p is too small
                p1 = p2;
                f1 = f2;
                p = p / con4;
                if (p3 >= 0.0) {
                    if (p >= p3) { p = p2 * con1 + p3 * con9; }
                }
                continue;
            }
        }

        // test whether the iteration process proceeds as theoretically expected.
        if ((f2 >= f1) || (f2 <= f3)) {
            *ier = 2;
            return;
        }

        // find the new value for p.
        p = fprati(&p1, &f1, &p2, &f2, &p3, &f3);
    }
}


static void
fprank(double* a, double* f, const int n, const int m, const int na,
       const double tol, double* c, double* sq, int* rank, double* aa, double* ff,
       double* h)
{
    // subroutine fprank finds the minimum norm solution of a least-
    // squares problem in case of rank deficiency.
    //
    // input parameters:
    //   a : array, which contains the non-zero elements of the observation
    //       matrix after triangularization by givens transformations.
    //   f : array, which contains the transformed right hand side.
    //   n : integer,which contains the dimension of a.
    //   m : integer, which denotes the bandwidth of a.
    // tol : real value, giving a threshold to determine the rank of a.
    //
    // output parameters:
    //   c : array, which contains the minimum norm solution.
    //  sq : real value, giving the contribution of reducing the rank
    //       to the sum of squared residuals.
    // rank : integer, which contains the rank of matrix a.

    // SciPy note: This is a very naive way of solving a LSQ problem but kept
    // for numerical compatibility with the Fortran code.

    int m1 = m - 1;
    int nl = 0;
    *sq = 0.0;

    for (int i = 1; i <= n; i++) {
        if (a[i - 1] > tol) { continue; }

        // if a sufficient small diagonal element is found, we put it to
        // zero. the remainder of the row corresponding to that zero diagonal
        // element is then rotated into triangle by givens rotations .
        // the rank deficiency is increased by one.
        nl++;
        if (i == n) { continue; }
        double yi = f[i - 1];
        for (int j = 1; j <= m1; j++) {
            h[j - 1] = a[(i - 1) + j * na];
        }
        h[m - 1] = 0.0;
        int i1 = i + 1;
        for (int ii = i1; ii <= n; ii++) {
            int i2 = (n - ii < m1) ? (n - ii) : m1;
            double piv = h[0];
            if (piv == 0.0) {
                if (i2 == 0) { break; }
                for (int j = 1; j <= i2; j++) { h[j - 1] = h[j]; }
                h[i2] = 0.0;
                continue;
            }
            double cos, sin;
            fpgivs(piv, &a[ii - 1], &cos, &sin);
            fprota(cos, sin, &yi, &f[ii - 1]);
            if (i2 == 0) { break; }
            for (int j = 1; j <= i2; j++) {
                fprota(cos, sin, &h[j], &a[(ii - 1) + j*na]);
                h[j - 1] = h[j];
            }
            h[i2] = 0.0;
        }
        // add to the sum of squared residuals the contribution of deleting
        // the row with small diagonal element.
        *sq = *sq + yi * yi;
    }

    // rank denotes the rank of a.
    *rank = n - nl;

    // let b denote the (rank*n) upper trapezoidal matrix which can be
    // obtained from the (n*n) upper triangular matrix a by deleting
    // the rows and interchanging the columns corresponding to a zero
    // diagonal element. if this matrix is factorized using givens
    // transformations as  b = (r) (u)  where
    //   r is a (rank*rank) upper triangular matrix,
    //   u is a (rank*n) orthonormal matrix
    // then the minimal least-squares solution c is given by c = b' v,
    // where v is the solution of the system  (r) (r)' v = g  and
    // g denotes the vector obtained from the old right hand side f, by
    // removing the elements corresponding to a zero diagonal element of a.

    for (int i = 1; i <= *rank; i++) {
        for (int j = 1; j <= m; j++) {
            aa[(i - 1) + (j - 1) * n] = 0.0;
        }
    }
    // form in aa the upper triangular matrix obtained from a by
    // removing rows and columns with zero diagonal elements. form in ff
    // the new right hand side by removing the elements of the old right
    // hand side corresponding to a deleted row.
    int ii = 0;
    for (int i = 1; i <= n; i++) {
        if (a[i - 1] <= tol) { continue; }
        ii++;
        ff[ii - 1] = f[i - 1];
        aa[ii - 1] = a[i - 1];
        int jj = ii;
        int kk = 1;
        int j = i;
        int j1 = (j - 1 < m1) ? (j - 1) : m1;
        if (j1 == 0) { continue; }
        for (int k = 1; k <= j1; k++) {
            j--;
            if (a[j - 1] <= tol) { continue; }
            kk++;
            jj--;
            aa[(jj - 1) + (kk - 1) * n] = a[(j - 1) + k * na];
        }
    }

    // form successively in h the columns of a with a zero diagonal element.
    ii = 0;
    for (int i = 1; i <= n; i++) {
        if (a[i - 1] > tol) {
            ii++;
            continue;
        }
        if (ii == 0) { continue; }
        int jj = 1;
        int j = i;
        int j1 = (j - 1 < m1) ? (j - 1) : m1;
        for (int k = 1; k <= j1; k++) {
            j--;
            if (a[j - 1] <= tol) { continue; }
            h[jj - 1] = a[(j - 1) + k * na];
            jj++;
        }
        for (int kk = jj; kk <= m; kk++) { h[kk - 1] = 0.0; }
        // rotate this column into aa by givens transformations.
        jj = ii;
        for (int i1 = 1; i1 <= ii; i1++) {
            int j1 = (jj - 1 < m1) ? (jj - 1) : m1;
            double piv = h[0];
            if (piv != 0.0) {
                double cos, sin;
                fpgivs(piv, &aa[jj - 1], &cos, &sin);
                if (j1 == 0) { break; }
                int kk = jj;
                for (int j2 = 1; j2 <= j1; j2++) {
                    kk--;
                    fprota(cos, sin, &h[j2], &aa[(kk - 1) + j2*n]);
                    h[j2 - 1] = h[j2];
                }
            } else {
                if (j1 == 0) { break; }
                for (int j2 = 1; j2 <= j1; j2++) {
                    h[j2 - 1] = h[j2];
                }
            }
            jj--;
            h[j1] = 0.0;
        }
    }

    // solve the system (aa) (f1) = ff
    ff[*rank - 1] = ff[*rank - 1] / aa[*rank - 1];
    int i = *rank - 1;
    if (i != 0) {
        for (int j = 2; j <= *rank; j++) {
            double store = ff[i - 1];
            int i1 = (j - 1 < m1) ? (j - 1) : m1;
            int k = i;
            for (int ii = 1; ii <= i1; ii++) {
                k++;
                store -= ff[k - 1] * aa[(i - 1) + ii * n];
            }
            ff[i - 1] = store / aa[i - 1];
            i--;
        }
    }

    // solve the system  (aa)' (f2) = f1
    ff[0] = ff[0] / aa[0 + 0 * n];
    if (*rank > 1) {
        for (int j = 2; j <= *rank; j++) {
            double store = ff[j - 1];
            int i1 = (j - 1 < m1) ? (j - 1) : m1;
            int k = j;
            for (int ii = 1; ii <= i1; ii++) {
                k--;
                store -= ff[k - 1] * aa[(k - 1) + ii * n];
            }
            ff[j - 1] = store / aa[j - 1];
        }
    }

    // premultiply f2 by the transpose of a.
    int k = 0;
    for (int i = 1; i <= n; i++) {
        double store = 0.0;
        if (a[i - 1] > tol) { k++; }
        int j1 = (i < m) ? i : m;
        int kk = k;
        int ij = i + 1;
        for (int j = 1; j <= j1; j++) {
            ij--;
            if (a[ij - 1] <= tol) { continue; }
            store += a[(ij - 1) + (j - 1) * na] * ff[kk - 1];
            kk--;
        }
        c[i - 1] = store;
    }

    // add to the sum of squared residuals the contribution of putting
    // to zero the small diagonal elements of matrix (a).
    double stor3 = 0.0;
    for (int i = 1; i <= n; i++) {
        if (a[i - 1] > tol) { continue; }
        double store = f[i - 1];
        int i1 = (n - i < m1) ? (n - i) : m1;
        if (i1 != 0) {
            for (int j = 1; j <= i1; j++) {
                int ij = i + j;
                store -= c[ij - 1] * a[(i - 1) + j * na];
            }
        }
        stor3 += (a[i - 1] * c[i - 1]) * ((a[i - 1] * c[i - 1]) - store - store);
    }
    *sq += stor3;
}


static double
fprati(double* p1, double* f1, double* p2, double* f2, double* p3, double* f3)
{
    // given three points (p1,f1),(p2,f2) and (p3,f3), function fprati
    // gives the value of p such that the rational interpolating function
    // of the form r(p) = (u*p+v)/(p+w) equals zero at p.

    // SciPy note: This pure looking function modifies its input arguments! Fun...
    double p;
    if (*p3 > 0.0) {
        // value of p in case p3 != inf
        double h1 = (*f1) * ((*f2) - (*f3));
        double h2 = (*f2) * ((*f3) - (*f1));
        double h3 = (*f3) * ((*f1) - (*f2));
        p = -((*p1)*(*p2)*h3 + (*p2)*(*p3)*h1 + (*p3)*(*p1)*h2) / ((*p1)*h1 + (*p2)*h2 + (*p3)*h3);
    } else {
        // value of p in case p3 == inf.
        p = ((*p1) * ((*f1) - (*f3)) * (*f2) - (*p2) * ((*f2) - (*f3)) * (*f1)) / (((*f1) - (*f2)) * (*f3));
    }
    // adjust the value of p1,f1,p3 and f3 such that f1 > 0 and f3 < 0.
    if ((*f2) < 0) {
        *p3 = *p2;
        *f3 = *f2;
    } else {
        *p1 = *p2;
        *f1 = *f2;
    }
    return p;
}


static void
fpregr(const int iopt, const double *x, const int mx, const double *y, const int my, const double *z, const int mz,
       const double xb, const double xe, const double yb, const double ye, const int kx, const int ky, const double s,
       const int nxest, const int nyest, const double tol, const int maxit, const int nc, int *nx, double *tx, int *ny,
       double *ty, double *c, double *fp, double *fp0, double *fpold, double *reducx, double *reducy, double *fpintx,
       double *fpinty, int *lastdi, int *nplusx, int *nplusy, int *nrx, int *nry, int *nrdatx, int *nrdaty,
       double *wrk, const int lwrk, int *ier)
{
    (void)lwrk;  // Unused
    const double one = 1.0;
    const double half = 0.5;
    const double con1 = 0.1;
    const double con9 = 0.9;
    const double con4 = 0.04;
    int mk1 = 0, nmaxx=0, nmaxy=0, nxe=0, nye=0;
    double acc = 0.0;
    // we partition the working space.
    int onlyxknots = 0;
    int kx1 = kx + 1;
    int ky1 = ky + 1;
    int kx2 = kx1 + 1;
    int ky2 = ky1 + 1;
    int lsx = 1;
    int lsy = lsx + mx * kx1;
    int lri = lsy + my * ky1;
    int mm = (nxest > my) ? nxest : my;
    int lq = lri + mm;
    int mynx = nxest * my;
    int lax = lq + mynx;
    int nxk = nxest * kx2;
    int lbx = lax + nxk;
    int lay = lbx + nxk;
    int lby = lay + nyest * ky2;

    /////////////////////////////////////////////////////////////////////////
    // part 1: determination of the number of knots and their position.    //
    // ****************************************************************    //
    //  given a set of knots we compute the least-squares spline sinf(x,y),//
    //  and the corresponding sum of squared residuals fp=f(p=inf).        //
    //  if iopt=-1  sinf(x,y) is the requested approximation.              //
    //  if iopt=0 or iopt=1 we check whether we can accept the knots:      //
    //    if fp <=s we will continue with the current set of knots.        //
    //    if fp > s we will increase the number of knots and compute the   //
    //       corresponding least-squares spline until finally fp<=s.       //
    //    the initial choice of knots depends on the value of s and iopt.  //
    //    if s=0 we have spline interpolation; in that case the number of  //
    //    knots equals nmaxx = mx+kx+1  and  nmaxy = my+ky+1.              //
    //    if s>0 and                                                       //
    //     *iopt=0 we first compute the least-squares polynomial of degree //
    //      kx in x and ky in y; nx=nminx=2*kx+2 and ny=nymin=2*ky+2.      //
    //     *iopt=1 we start with the knots found at the last call of the   //
    //      routine, except for the case that s > fp0; then we can compute //
    //      the least-squares polynomial directly.                         //
    /////////////////////////////////////////////////////////////////////////
    //  determine the number of knots for polynomial approximation.
    int nminx = 2 * kx1;
    int nminy = 2 * ky1;

    if (iopt < 0) { goto restart_iteration; }
    // acc denotes the absolute tolerance for the root of f(p)=s.
    acc = tol * s;
    // find nmaxx and nmaxy which denote the number of knots in x- and y-
    // direction in case of spline interpolation.
    nmaxx = mx + kx1;
    nmaxy = my + ky1;

    // find nxe and nye which denote the maximum number of knots
    // allowed in each direction.
    nxe = (nxest < nmaxx) ? nxest : nmaxx;
    nye = (nyest < nmaxy) ? nyest : nmaxy;

    if (s > 0.0) { goto L100; }
    // if s = 0, s(x, y) is an interpolating spline.
    *nx = nmaxx;
    *ny = nmaxy;

    // test whether the required storage space exceeds the available one.
    if ((*ny > nyest) || (*nx > nxest)) {
        *ier = 1;
        return;
    }
L10:
    // find position of the interior knots in case of interpolation.
    // the knots in x-direction.
    mk1 = mx - kx1;
    if (mk1 == 0) { goto L60; }
    int k3 = kx / 2;
    int i = kx1 + 1;
    int j = k3 + 2;
    if ((k3 * 2) == kx) {
        for (int l = 1; l <= mk1; l++) {
            tx[i - 1] = (x[j - 1] + x[j - 2]) * half;
            i++;
            j++;
        }
        if (onlyxknots) {
            onlyxknots = 0;
            goto restart_iteration;
        }
    } else {
        for (int l = 1; l <= mk1; l++) {
            tx[i - 1] = x[j - 1];
            i++;
            j++;
        }
    }
L60:
    // the knots in y-direction.
    mk1 = my - ky1;
    if (mk1 == 0) { goto restart_iteration; }
    k3 = ky / 2;
    i = ky1 + 1;
    j = k3 + 2;
    if ((k3 * 2) == ky) {
        for (int l = 1; l <= mk1; l++) {
            ty[i - 1] = (y[j - 1] + y[j - 2]) * half;
            i++;
            j++;
        }
    } else {
        for (int l = 1; l <= mk1; l++) {
            ty[i - 1] = y[j - 1];
            i++;
            j++;
        }
    }
    goto restart_iteration;

L100:
    // if s > 0 our initial choice of knots depends on the value of iopt.
    if ((iopt == 0) || (*fp0 <= s)) {
        // if iopt=0 or iopt=1 and s>=fp0, we start computing the least-squares
        // polynomial of degree kx in x and ky in y (which is a spline without
        // interior knots).
        *nx = nminx;
        *ny = nminy;
        nrdatx[0] = mx - 2;
        nrdaty[0] = my - 2;
        *lastdi = 0;
        *nplusx = 0;
        *nplusy = 0;
        *fp0 = 0.0;
        *fpold = 0.0;
        *reducx = 0.0;
        *reducy = 0.0;
    } else {
        // if iopt=1 and fp0 > s we start computing the least- squares spline
        // according to the set of knots found at the last call of the routine.
        // we determine the number of grid coordinates x(i) inside each knot
        // interval (tx(l),tx(l+1)).
        int l = kx2;
        int j = 1;
        nrdatx[0] = 0;
        int mpm = mx - 1;
        for (int i = 2; i <= mpm; i++) {
            nrdatx[j - 1]++;
            if (x[i - 1] < tx[l - 1]) { continue; }
            nrdatx[j - 1]--;
            l++;
            j++;
            nrdatx[j - 1] = 0;
        }
        // we determine the number of grid coordinates y(i) inside each knot
        // interval (ty(l),ty(l+1)).
        l = ky2;
        j = 1;
        nrdaty[0] = 0;
        mpm = my - 1;
        for (int i = 2; i <= mpm; i++) {
            nrdaty[j - 1]++;
            if (y[i - 1] < ty[l - 1]) { continue; }
            nrdaty[j - 1]--;
            l++;
            j++;
            nrdaty[j - 1] = 0;
        }
    }

restart_iteration:
    ;
    int mpm = mx + my;
    int ifsx = 0;
    int ifsy = 0;
    int ifbx = 0;
    int ifby = 0;
    double p = -one;

    // main loop for the different sets of knots.mpm=mx+my is a save upper
    // bound for the number of trials.
    for (int iter = 1; iter <= mpm; iter++) {
        if ((*nx == nminx) && (*ny == nminy)) {
            *ier = -2;
        }

        // find nrintx (nrinty) which is the number of knot intervals in the
        // x-direction (y-direction).
        int nrintx = *nx - nminx + 1;
        int nrinty = *ny - nminy + 1;

        // find ncof, the number of b-spline coefficients for the current set
        // of knots.
        // int nk1x = *nx - kx1;
        // int nk1y = *ny - ky1;
        // int ncof = nk1x * nk1y;

        // find the position of the additional knots which are needed for the
        // b-spline representation of s(x,y).
        int i = *nx;
        for (int j = 1; j <= kx1; j++) {
            tx[j - 1] = xb;
            tx[i - 1] = xe;
            i--;
        }
        i = *ny;
        for (int j = 1; j <= ky1; j++) {
            ty[j - 1] = yb;
            ty[i - 1] = ye;
            i--;
        }

        // find the least-squares spline sinf(x,y) and calculate for each knot
        // interval tx(j+kx)<=x<=tx(j+kx+1) (ty(j+ky)<=y<=ty(j+ky+1)) the sum
        // of squared residuals fpintx(j),j=1,2,...,nx-2*kx-1 (fpinty(j),j=1,2,
        // ...,ny-2*ky-1) for the data points having their absciss (ordinate)-
        // value belonging to that interval.
        // fp gives the total sum of squared residuals.
        fpgrre(ifsx, ifsy, ifbx, ifby, x, mx, y, my, z, mz, kx, ky, tx, *nx, ty, *ny, p, c, nc, fp,
               fpintx, fpinty, mm, mynx, kx1, kx2, ky1, ky2, &wrk[lsx - 1], &wrk[lsy - 1], &wrk[lri - 1],
               &wrk[lq - 1], &wrk[lax - 1], &wrk[lay - 1], &wrk[lbx - 1], &wrk[lby - 1], nrx, nry);

        if (*ier == -2) {
            *fp0 = *fp;
        }

        // test whether the least-squares spline is an acceptable solution.
        if (iopt < 0) { return; }
        double fpms = *fp - s;
        if (fabs(fpms) < acc) { return; }

        // if f(p=inf) < s, we accept the choice of knots.
        if (fpms < 0.0) { break; }

        // if nx=nmaxx and ny=nmaxy, sinf(x,y) is an interpolating spline.
        if ((*nx == nmaxx) && (*ny == nmaxy)) {
            *ier = -1;
            *fp = 0.0;
            return;
        }

        // increase the number of knots.
        // if nx=nxe and ny=nye we cannot further increase the number of knots
        // because of the storage capacity limitation.
        if ((*nx == nxe) && (*ny == nye)) {
            *ier = 1;
            return;
        }

        *ier = 0;

        // adjust the parameter reducx or reducy according to the direction
        // in which the last added knots were located.
        if (*lastdi < 0) {
            *reducx = *fpold - *fp;
        } else if (*lastdi > 0) {
            *reducy = *fpold - *fp;
        }
        // store the sum of squared residuals for the current set of knots.
        *fpold = *fp;

        // find nplx, the number of knots we should add in the x-direction.
        int nplx = 1;
        if (*nx != nminx) {
            int npl1 = (*nplusx) * 2;
            double rn = (double)(*nplusx);
            if (*reducx > acc) {
                npl1 = (int)(rn * fpms / (*reducx));
            }
            // nplx = min0(nplusx*2,max0(npl1,nplusx/2,1))
            int temp1 = (*nplusx) * 2;
            int temp2 = (*nplusx) / 2;
            if (temp2 < 1) { temp2 = 1; }
            if (npl1 > temp2) {
                nplx = npl1;
            } else {
                nplx = temp2;
            }
            if (nplx > temp1) { nplx = temp1; }
        }

        // find nply, the number of knots we should add in the y-direction.
        int nply = 1;
        if (*ny != nminy) {
            int npl1 = (*nplusy) * 2;
            double rn = (double)(*nplusy);
            if (*reducy > acc) {
                npl1 = (int)(rn * fpms / (*reducy));
            }
            // nply = min0(nplusy*2,max0(npl1,nplusy/2,1))
            int temp1 = (*nplusy) * 2;
            int temp2 = (*nplusy) / 2;
            if (temp2 < 1) { temp2 = 1; }
            if (npl1 > temp2) {
                nply = npl1;
            } else {
                nply = temp2;
            }
            if (nply > temp1) { nply = temp1; }
        }

        // Super weird goto logic
        if (nplx < nply) { goto L210; }
        if (nplx == nply) {
            if (*lastdi < 0) { goto L230; }
        } else {
            goto L230;
        }
L210:
        if (*nx == nxe) { goto L230; }
        // addition in the x-direction.
        *lastdi = -1;
        *nplusx = nplx;
        ifsx = 0;
        for (int l = 1; l <= nplx; l++) {
            fpknot(x, mx, tx, nx, fpintx, nrdatx, &nrintx, nxest, 1);
            // test whether we cannot further increase the number of knots in the
            // x-direction.
            if (*nx == nmaxx) {
                onlyxknots = 1;
                goto L10;
            }
            if (*nx == nxe) { goto L250; }
        }
        goto L250;
L230:
        if (*ny == nye) { goto L210; }
        // addition in the y-direction.
        *lastdi = 1;
        *nplusy = nply;
        ifsy = 0;
        for (int l = 1; l <= nply; l++) {
            // test whether we cannot further increase the number of knots in the
            // y-direction.
            fpknot(y, my, ty, ny, fpinty, nrdaty, &nrinty, nyest, 1);
            if (*ny == nmaxy) { goto L60; }
            if (*ny == nye) { goto L250; }
        }
L250:
        ;
    }

    // test whether the least-squares polynomial is a solution of our
    // approximation problem.
    if (*ier == -2) { return; }

    //////////////////////////////////////////////////////////////////////////
    // part 2: determination of the smoothing spline sp(x,y)                //
    // *****************************************************                //
    //  we have determined the number of knots and their position. we now   //
    //  compute the b-spline coefficients of the smoothing spline sp(x,y).  //
    //  this smoothing spline varies with the parameter p in such a way that//
    //    f(p) = sumi=1,mx(sumj=1,my((z(i,j)-sp(x(i),y(j)))**2)             //
    //  is a continuous, strictly decreasing function of p. moreover the    //
    //  least-squares polynomial corresponds to p=0 and the least-squares   //
    //  spline to p=infinity. iteratively we then have to determine the     //
    //  positive value of p such that f(p)=s. the process which is proposed //
    //  here makes use of rational interpolation. f(p) is approximated by a //
    //  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  //
    //  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)//
    //  are used to calculate the new value of p such that r(p)=s.          //
    //  convergence is guaranteed by taking f1 > 0 and f3 < 0.              //
    //////////////////////////////////////////////////////////////////////////
    // initial value for p.
    double p1 = 0.0;
    double f1 = *fp0 - s;
    double p3 = -one;
    double fpms = *fp - s;
    double f3 = fpms;
    p = one;
    int ich1 = 0;
    int ich3 = 0;

    // iteration to find root of f(p) = s
    for (int iter = 1; iter <= maxit; iter++) {
        // find smoothing spline and residual
        fpgrre(ifsx, ifsy, ifbx, ifby, x, mx, y, my, z, mz, kx, ky, tx, *nx, ty, *ny, p, c, nc, fp,
               fpintx, fpinty, mm, mynx, kx1, kx2, ky1, ky2, &wrk[lsx - 1], &wrk[lsy - 1], &wrk[lri - 1],
               &wrk[lq - 1], &wrk[lax - 1], &wrk[lay - 1], &wrk[lbx - 1], &wrk[lby - 1], nrx, nry);

        // test whether the approximation sp(x,y) is an acceptable solution.
        fpms = *fp - s;
        if (fabs(fpms) < acc) { return; }
        // test whether the maximum allowable number of iterations has been
        // reached.
        if (iter == (maxit)) {
            *ier = 3;
            return;
        }

        // carry out one more step of the iteration process.
        double p2 = p;
        double f2 = fpms;

        if (ich3 == 0) {
            if ((f2 - f3) > acc) {
                if (f2 < 0.0) {
                    ich3 = 1;
                }
            } else {
                // our initial choice of p is too large.
                p3 = p2;
                f3 = f2;
                p = p * con4;
                if (p <= p1) {
                    p = p1 * con9 + p2 * con1;
                }
                continue;
            }
        }

        if (ich1 == 0) {
            if ((f1 - f2) > acc) {
                if (f2 > 0.0) {
                    ich1 = 1;
                }
            } else {
                // our initial choice of p is too small.
                p1 = p2;
                f1 = f2;
                p = p / con4;
                if (p3 >= 0.0) {
                    if (p >= p3) {
                        p = p2 * con1 + p3 * con9;
                    }
                }
                continue;
            }
        }

        // check iteration proceeds as expected
        if ((f2 >= f1) || (f2 <= f3)) {
            *ier = 2;
            return;
        }

        // find new value of p
        p = fprati(&p1, &f1, &p2, &f2, &p3, &f3);
    }
}


static void
fprota(const double c, const double s, double* restrict a, double* restrict b)
{
    // fprota applies a givens rotation to the pair (a,b).
    double temp = c * (*a) - s * (*b);
    *b          = c * (*b) + s * (*a);
    *a          = temp;
}


static void
fprpsp(const int nt, const int np, const double *co, const double *si, double *c, double *f, const int ncoff)
{
    double cn, c1, c2, c3;
    int i, ii, j, k, l, ncof, npp, np4, nt4;

    nt4 = nt - 4;
    np4 = np - 4;
    npp = np4 - 3;
    ncof = 6 + npp * (nt4 - 4);
    c1 = c[0];
    cn = c[ncof - 1];
    j = ncoff;

    for (i = 1; i <= np4; i++) {
        f[i - 1] = c1;
        f[j - 1] = cn;
        j--;
    }

    i = np4;
    j = 1;
    for (l = 3; l <= nt4; l++) {
        ii = i;
        if ((l == 3) || (l == nt4)) {
            if (l == nt4) { c1 = cn; }
            c2 = c[j];
            c3 = c[j + 1];
            j += 2;
            for (k = 1; k <= npp; k++) {
                i++;
                f[i - 1] = c1 + c2 * co[k - 1] + c3 * si[k - 1];
            }
        } else {
            for (k = 1; k <= npp; k++) {
                i++;
                j++;
                f[i - 1] = c[j - 1];
            }
        }

        for (k = 1; k <= 3; k++) {
            ii++;
            i++;
            f[i - 1] = f[ii - 1];
        }
    }

    for (i = 1; i <= ncoff; i++) {
        c[i - 1] = f[i - 1];
    }
}


static void
fpsphe(const int iopt, const int m, const double *teta, const double *phi, const double *r, const double *w,
       const double s, const int ntest, const int npest, const double eta, const double tol, const int maxit,
       const int ib1, const int ib3, const int nc, const int ncc, const int intest, const int nrest,
       int *nt, double *tt, int *np, double *tp, double *c, double *fp, double *sup, double *fpint, double *coord,
       double *f, double *ff, double *row, double *coco, double *cosi, double *a, double *q, double *bt, double *bp,
       double *spt, double *spp, double *h, int *index, int *nummer, double *wrk, const int lwrk, int *ier)
{
    (void)ib1;    // Unused
    (void)ib3;    // Unused
    (void)nc;     // Unused
    (void)intest; // Unused
    (void)nrest;  // Unused
    double aa, acc=0.0, arg, cn, co, c1, dmax, d1, d2, eps, facc=0.0, facs=0.0, fac1, fac2, fn;
    double fpmax, fpms=0.0, f1, f2, f3, hti, htj, p, pinv=0.0, piv, pi2, p1, p2, p3, ri, si;
    double sigma, sq, store, wi, rn, one, con1, con9, con4, half, ten;
    int i, iband=0, iband1=0, ii, irot, j, jlt, jrot, j1, j2, l, la, lf, lh, ll, lp, lt, lwest, l1, l2, l3, l4;
    int ncof=0, ncoff=0, npp=0, np4=0, nreg=0, nrint, nrr, ntt=0, nt4=0, num, num1, rank, in, iter, i1, i2, i3;
    double ht[4], hp[4];

    // set constants
    one = 1.0;
    con1 = 0.1;
    con9 = 0.9;
    con4 = 0.04;
    half = 0.5;
    ten = 10.0;
    pi2 = PI + PI;
    eps = sqrt(eta);

    if (iopt < 0) { goto L70; }
    acc = tol * s;
    if (s < *sup) {
        if (*np < 11) { goto L60; }
        goto L70;
    }

    // if iopt=0 we begin by computing the weighted least-squares polynomial
    // of the form
    //    s(teta,phi) = c1*f1(teta) + cn*fn(teta)
    // where f1(teta) and fn(teta) are the cubic polynomials satisfying
    //    f1(0) = 1, f1(pi) = f1'(0) = f1'(pi) = 0 ; fn(teta) = 1-f1(teta).
    // the corresponding weighted sum of squared residuals gives the upper
    // bound sup for the smoothing factor s.
    *sup = 0.0;
    d1 = 0.0;
    d2 = 0.0;
    c1 = 0.0;
    cn = 0.0;
    fac1 = PI * (one + half);
    fac2 = (one + one) / (PI * PI * PI);
    aa = 0.0;
    for (i = 1; i <= m; i++) {
        wi = w[i - 1];
        ri = r[i - 1] * wi;
        arg = teta[i - 1];
        fn = fac2 * arg * arg * (fac1 - arg);
        f1 = (one - fn) * wi;
        fn = fn * wi;
        if (fn != 0.0) {
            fpgivs(fn, &d1, &co, &si);
            fprota(co, si, &f1, &aa);
            fprota(co, si, &ri, &cn);
        }
        if (f1 != 0.0) {
            fpgivs(f1, &d2, &co, &si);
            fprota(co, si, &ri, &c1);
        }
        *sup = *sup + ri * ri;
    }
    if (d2 != 0.0) c1 = c1 / d2;
    if (d1 != 0.0) cn = (cn - aa * c1) / d1;
    // find the b-spline representation of this least-squares polynomial
    *nt = 8;
    *np = 8;
    for (i = 0; i < 4; i++) {
        c[i] = c1;
        c[i + 4] = c1;
        c[i + 8] = cn;
        c[i + 12] = cn;
        tt[i] = 0.0;
        tt[i + 4] = PI;
        tp[i] = 0.0;
        tp[i + 4] = pi2;
    }
    *fp = *sup;
    // test whether the least-squares polynomial is an acceptable solution
    fpms = *sup - s;
    if (fpms < acc){
        *ier = -2;
        return;
    }
L60:
    // test whether we cannot further increase the number of knots.
    if ((npest < 11) || (ntest < 9)) {
        *ier = 1;
        return;
    }
    *np = 11;
    tp[4] = PI * half;
    tp[5] = PI;
    tp[6] = tp[4] + PI;
    *nt = 9;
    tt[4] = tp[4];

L70:
    /////////////////////////////////////////////////////////////////////////
    // part 1 : computation of least-squares spherical splines.            //
    // ********************************************************            //
    // if iopt < 0 we compute the least-squares spherical spline according //
    // to the given set of knots.                                          //
    // if iopt >=0 we compute least-squares spherical splines with increas-//
    // ing numbers of knots until the corresponding sum f(p=inf)<=s.       //
    // the initial set of knots then depends on the value of iopt:         //
    //   if iopt=0 we start with one interior knot in the teta-direction   //
    //             (pi/2) and three in the phi-direction (pi/2,pi,3*pi/2). //
    //   if iopt>0 we start with the set of knots found at the last call   //
    //             of the routine.                                         //
    /////////////////////////////////////////////////////////////////////////
    // main loop for the different sets of knots. m is a save upper bound
    // for the number of trials.
    for (iter = 1; iter <= m; iter++) {
        // find the position of the additional knots which are needed for the
        // b-spline representation of s(teta,phi).
        l1 = 4;
        l2 = l1;
        l3 = *np - 3;
        l4 = l3;
        tp[l2 - 1] = 0.0;
        tp[l3 - 1] = pi2;
        for (i = 1; i <= 3; i++) {
            l1++;
            l2--;
            l3++;
            l4--;
            tp[l2 - 1] = tp[l4 - 1] - pi2;
            tp[l3 - 1] = tp[l1 - 1] + pi2;
        }
        l = *nt;
        for (i = 1; i <= 4; i++) {
            tt[i - 1] = 0.0;
            tt[l - 1] = PI;
            l--;
        }
        // find nrint, the total number of knot intervals and nreg, the number
        // of panels in which the approximation domain is subdivided by the
        // intersection of knots.
        ntt = *nt - 7;
        npp = *np - 7;
        nrr = npp / 2;
        nrint = ntt + npp;
        nreg = ntt * npp;

        // arrange the data points according to the panel they belong to.
        fporde(teta, phi, m, 3, 3, tt, *nt, tp, *np, nummer, index, nreg);

        // find the b-spline coefficients coco and cosi of the cubic spline
        // approximations sc(phi) and ss(phi) for cos(phi) and sin(phi).
        for (i = 1; i <= npp; i++) {
            coco[i - 1] = 0.0;
            cosi[i - 1] = 0.0;
            for (j = 1; j <= npp; j++) {
                a[(i - 1) + (j - 1)*ncc] = 0.0;
            }
        }
        // the coefficients coco and cosi are obtained from the conditions
        // sc(tp(i))=cos(tp(i)),resp. ss(tp(i))=sin(tp(i)),i=4,5,...np-4.
        for (i = 1; i <= npp; i++) {
            l2 = i + 3;
            arg = tp[l2 - 1];
            fpbspl(tp, *np, 3, arg, l2, hp);
            for (j = 1; j <= npp; j++) {
                row[j - 1] = 0.0;
            }
            ll = i;
            for (j = 1; j <= 3; j++) {
                if (ll > npp) { ll = 1; }
                row[ll - 1] += hp[j - 1];
                ll++;
            }
            facc = cos(arg);
            facs = sin(arg);
            for (j = 1; j <= npp; j++) {
                piv = row[j - 1];
                if (piv == 0.0) { continue; }
                fpgivs(piv, &a[j - 1], &co, &si);
                fprota(co, si, &facc, &coco[j - 1]);
                fprota(co, si, &facs, &cosi[j - 1]);
                if (j == npp) { break; }
                j1 = j + 1;
                i2 = 1;
                for (l = j1; l <= npp; l++) {
                    i2++;
                    fprota(co, si, &row[l - 1], &a[(j - 1) +  (i2 - 1)*ncc]);
                }
            }
        }
        fpback(a, coco, npp, npp, coco, ncc);
        fpback(a, cosi, npp, npp, cosi, ncc);
        // find ncof, the dimension of the spherical spline and ncoff, the
        // number of coefficients in the standard b-spline representation.
        nt4 = *nt - 4;
        np4 = *np - 4;
        ncoff = nt4 * np4;
        ncof = 6 + npp * (ntt - 1);
        // find the bandwidth of the observation matrix a.
        iband = 4 * npp;
        if (ntt == 4) iband = 3 * (npp + 1);
        if (ntt < 4) iband = ncof;
        iband1 = iband - 1;

        for (i = 1; i <= ncof; i++) {
            f[i - 1] = 0.0;
            for (j = 1; j <= iband; j++) {
                a[(i - 1) + (j - 1)*ncc] = 0.0;
            }
        }
        // initialize the sum of squared residuals.
        *fp = 0.0;
        // fetch the data points in the new order. main loop for the
        // different panels.
        for (num = 1; num <= nreg; num++) {
            // fix certain constants for the current panel; jrot records the column
            // number of the first non-zero element in a row of the observation
            // matrix according to a data point of the panel.
            num1 = num - 1;
            lt = num1 / npp;
            l1 = lt + 4;
            lp = num1 - (lt * npp) + 1;
            l2 = lp + 3;
            lt = lt + 1;
            jrot = 0;
            if (lt > 2) { jrot = 3 + (lt - 3) * npp; }
            // test whether there are still data points in the current panel.
            in = index[num - 1];
            while (in != 0) {
                // fetch a new data point.
                wi = w[in - 1];
                ri = r[in - 1] * wi;
                // evaluate for the teta-direction, the 4 non-zero b-splines at teta(in)
                fpbspl(tt, *nt, 3, teta[in - 1], l1, ht);
                // evaluate for the phi-direction, the 4 non-zero b-splines at phi(in)
                fpbspl(tp, *np, 3, phi[in - 1], l2, hp);
                // store the value of these b-splines in spt and spp resp.
                for (i = 1; i <= 4; i++) {
                    spp[(in - 1) + (i - 1)*m] = hp[i - 1];
                    spt[(in - 1) + (i - 1)*m] = ht[i - 1];
                }
                // initialize the new row of the observation matrix.
                for (i = 1; i <= iband; i++) { h[i - 1] = 0.0; }
                // calculate the non-zero elements of the new row by making the cross
                // products of the non-zero b-splines in teta- and phi-direction and
                // by taking into account the conditions of the spherical splines.
                for (i = 1; i <= npp; i++) { row[i - 1] = 0.0; }
                // take into account the condition (3) of the spherical splines.
                ll = lp;
                for (i = 1; i <= 4; i++) {
                    if (ll > npp) { ll = 1; }
                    row[ll - 1] += hp[i - 1];
                    ll++;
                }
                // take into account the other conditions of the spherical splines.
                if ((lt <= 2) || (lt >= (ntt - 1))) {
                    facc = 0.0;
                    facs = 0.0;
                    for (i = 1; i <= npp; i++) {
                        facc += row[i - 1] * coco[i - 1];
                        facs += row[i - 1] * cosi[i - 1];
                    }
                }
                // fill in the non-zero elements of the new row.
                j1 = 0;
                for (j = 1; j <= 4; j++) {
                    jlt = j + lt;
                    htj = ht[j - 1];
                    if (!((jlt > 2) && (jlt <= nt4))) {
                        j1++;
                        h[j1 - 1] += htj;
                        continue;
                    }
                    if ((jlt == 3) || (jlt == nt4)) {
                        if (jlt != 3) {
                            h[j1] = facc * htj;
                            h[j1 + 1] = facs * htj;
                            h[j1 + 2] = htj;
                            j1 += 2;
                        } else {
                            h[0] += htj;
                            h[1] = facc * htj;
                            h[2] = facs * htj;
                            j1 = 3;
                        }
                    } else {
                        for (i = 1; i <= npp; i++) {
                            j1++;
                            h[j1 - 1] = row[i - 1] * htj;
                        }
                    }
                }

                for (i = 1; i <= iband; i++) { h[i - 1] *= wi; }

                // rotate the row into triangle by givens transformations.
                irot = jrot;
                for (i = 1; i <= iband; i++) {
                    irot++;
                    piv = h[i - 1];
                    if (piv == 0.0) { continue; }
                    // calculate the parameters of the givens transformation.
                    fpgivs(piv, &a[irot - 1], &co, &si);
                    // apply that transformation to the right hand side.
                    fprota(co, si, &ri, &f[irot - 1]);
                    if (i == iband) { break; }
                    // apply that transformation to the left hand side.
                    i2 = 1;
                    i3 = i + 1;
                    for (j = i3; j <= iband; j++) {
                        i2++;
                        fprota(co, si, &h[j - 1], &a[(irot - 1) + (i2 - 1)*ncc]);
                    }
                }
                // add the contribution of the row to the sum of squares of residual
                // right hand sides.
                *fp += ri * ri;
                // find the number of the next data point in the panel.
                in = nummer[in - 1];
            }
        }
        // find dmax, the maximum value for the diagonal elements in the reduced triangle.
        dmax = 0.0;
        for (i = 1; i <= ncof; i++) {
            if (a[i - 1] > dmax) { dmax = a[i- 1]; }
        }
        // check whether the observation matrix is rank deficient.
        sigma = eps * dmax;
        int rank_deficient = 0;
        for (i = 1; i <= ncof; i++) {
            if (a[i - 1] <= sigma) {
                rank_deficient = 1;
                break;
            }
        }
        if (rank_deficient) {
            // in case of rank deficiency, find the minimum norm solution.
            lwest = ncof * iband + ncof + iband;
            if (lwrk < lwest) {
                *ier = lwest;
                return;
            }
            lf = 1;
            lh = lf + ncof;
            la = lh + iband;
            for (i = 0; i < ncof; i++) {
                ff[i] = f[i];
                for (j = 0; j < iband; j++) {
                    q[i + j*ncc] = a[i + j*ncc];
                }
            }
            fprank(q, ff, ncof, iband, ncc, sigma, c, &sq, &rank, &wrk[la - 1], &wrk[lf - 1], &wrk[lh - 1]);
            for (i = 0; i < ncof; i++) {
                q[i] /= dmax;
            }
            // add to the sum of squared residuals, the contribution of reducing
            // the rank.
            *fp += sq;
        } else {
            // backward substitution in case of full rank.
            fpback(a, f, ncof, iband, c, ncc);
            rank = ncof;
            for (i = 0; i < ncof; i++) {
                q[i] = a[i] / dmax;
            }
        }

        // find the coefficients in the standard b-spline representation of
        // the spherical spline.
        fprpsp(*nt, *np, coco, cosi, c, ff, ncoff);
        // test whether the least-squares spline is an acceptable solution.
        if (iopt < 0) {
            if (*fp <= 0.0) {
                *ier = -1;
                *fp = 0.0;
            }
            if (ncof != rank) {
                *ier = -rank;
            }
            return;
        }
        fpms = *fp - s;
        if (fabs(fpms) <= acc) {
            if (*fp <= 0.0) {
                *ier = -1;
                *fp = 0.0;
            }
            if (ncof != rank) {
                *ier = -rank;
            }
            return;
        }
        // if f(p=inf) < s, accept the choice of knots.
        if (fpms < 0.0) {
            break;
        }

        // test whether we cannot further increase the number of knots.
        if (ncof > m) {
            *ier = 4;
            return;
        }
        // search where to add a new knot.
        // find for each interval the sum of squared residuals fpint for the
        // data points having the coordinate belonging to that knot interval.
        // calculate also coord which is the same sum, weighted by the position
        // of the data points considered.
        for (i = 1; i <= nrint; i++) {
            fpint[i - 1] = 0.0;
            coord[i - 1] = 0.0;
        }
        for (num = 1; num <= nreg; num++) {
            num1 = num - 1;
            lt = num1 / npp;
            l1 = lt + 1;
            lp = num1 - lt * npp;
            l2 = lp + 1 + ntt;
            jrot = lt * np4 + lp;
            in = index[num - 1];
            while (in != 0) {
                store = 0.0;
                i1 = jrot;
                for (i = 1; i <= 4; i++) {
                    hti = spt[(in - 1) + (i - 1)*m];
                    j1 = i1;
                    for (j = 1; j <= 4; j++) {
                        j1++;
                        store += hti * spp[(in - 1) + (j - 1) * m] * c[j1 - 1];
                    }
                    i1 = i1 + np4;
                }
                store = pow(w[in - 1] * (r[in - 1] - store), 2);
                fpint[l1 - 1] += store;
                coord[l1 - 1] += store * teta[in - 1];
                fpint[l2 - 1] += store;
                coord[l2 - 1] += store * phi[in - 1];
                in = nummer[in - 1];
            }
        }
        // find the interval for which fpint is maximal on the condition that
        // there still can be added a knot.
        l1 = 1;
        l2 = nrint;
        if (ntest < *nt + 1) { l1 = ntt + 1; }
        if (npest < *np + 2) { l2 = ntt; }
        // test whether we cannot further increase the number of knots.
        if (l1 > l2) { *ier = 1; return; }
        while (1) {
            fpmax = 0.0;
            l = 0;
            for (i = l1; i <= l2; i++) {
                if (fpmax >= fpint[i - 1]) { continue; }
                l = i;
                fpmax = fpint[i - 1];
            }
            if (l == 0) { *ier = 5; return; }
            // calculate the position of the new knot.
            arg = coord[l - 1] / fpint[l - 1];
            // test in what direction the new knot is going to be added.
            if (l > ntt) {
                // addition in the phi-direction
                l4 = l + 4 - ntt;
                if (arg >= PI) {
                    arg = arg - PI;
                    l4 = l4 - nrr;
                }
                fpint[l - 1] = 0.0;
                fac1 = tp[l4 - 1] - arg;
                fac2 = arg - tp[l4 - 2];
                if ((fac1 > (ten * fac2)) || (fac2 > (ten * fac1))) { continue; }
                ll = nrr + 4;
                j = ll;
                for (i = l4; i <= ll; i++) {
                    tp[j] = tp[j - 1];
                    j--;
                }
                tp[l4 - 1] = arg;
                *np += 2;
                nrr++;
                for (i = 5; i <= ll; i++) {
                    j = i + nrr;
                    tp[j - 1] = tp[i - 1] + PI;
                }
                break;
            } else {
                // addition in the teta-direction
                l4 = l + 4;
                fpint[l - 1] = 0.0;
                fac1 = tt[l4 - 1] - arg;
                fac2 = arg - tt[l4 - 2];
                if ((fac1 > (ten * fac2)) || (fac2 > (ten * fac1))) { continue; }
                j = *nt;
                for (i = l4; i <= *nt; i++) {
                    tt[j] = tt[j - 1];
                    j--;
                }
                tt[l4 - 1] = arg;
                *nt += 1;
                break;
            }
        }
        // restart the computations with the new set of knots.
    }

    //////////////////////////////////////////////////////////////////////////
    // part 2: determination of the smoothing spherical spline.             //
    // ********************************************************             //
    // we have determined the number of knots and their position. we now    //
    // compute the coefficients of the smoothing spline sp(teta,phi).       //
    // the observation matrix a is extended by the rows of a matrix, expres-//
    // sing that sp(teta,phi) must be a constant function in the variable   //
    // phi and a cubic polynomial in the variable teta. the corresponding   //
    // weights of these additional rows are set to 1/(p). iteratively       //
    // we than have to determine the value of p such that f(p) = sum((w(i)* //
    // (r(i)-sp(teta(i),phi(i))))**2)  be = s.                              //
    // we already know that the least-squares polynomial corresponds to p=0,//
    // and that the least-squares spherical spline corresponds to p=infin.  //
    // the iteration process makes use of rational interpolation. since f(p)//
    // is a convex and strictly decreasing function of p, it can be approx- //
    // imated by a rational function of the form r(p) = (u*p+v)/(p+w).      //
    // three values of p (p1,p2,p3) with corresponding values of f(p) (f1=  //
    // f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the new value   //
    // of p such that r(p)=s. convergence is guaranteed by taking f1>0,f3<0.//
    //////////////////////////////////////////////////////////////////////////

    // evaluate the discontinuity jumps of the 3-th order derivative of
    // the b-splines at the knots tt(l),l=5,...,nt-4.
    fpdisc(tt, *nt, 5, bt, ntest);
    // evaluate the discontinuity jumps of the 3-th order derivative of
    // the b-splines at the knots tp(l),l=5,...,np-4.
    fpdisc(tp, *np, 5, bp, npest);

    // initial value for p.
    p1 = 0.0;
    f1 = *sup - s;
    p3 = -1.0;
    f3 = fpms;
    p = 0.0;
    for (i = 1; i <= ncof; i++) {
        p += a[i - 1];
    }
    rn = ncof;
    p = rn / p;
    // find the bandwidth of the extended observation matrix.
    int iband4 = iband + 3;
    if (ntt <= 4) { iband4 = ncof; }
    int iband3 = iband4 - 1;
    int ich1 = 0;
    int ich3 = 0;
    // iteration process to find the root of f(p)=s.
    for (iter = 1; iter <= maxit; iter++) {
        pinv = one / p;
        // store the triangularized observation matrix into q.
        for (i = 1; i <= ncof; i++) {
            ff[i - 1] = f[i - 1];
            for (j = 1; j <= iband4; j++) {
                q[(i - 1) + (j - 1)*ncc] = 0.0;
            }
            for (j = 1; j <= iband; j++) {
                q[(i - 1) + (j - 1)*ncc] = a[(i - 1) + (j - 1)*ncc];
            }
        }
        // extend the observation matrix with the rows of a matrix, expressing
        // that for teta=cst. sp(teta,phi) must be a constant function.
        int nt6 = *nt - 6;
        for (i = 5; i <= np4; i++) {
            ii = i - 4;
            for (l = 1; l <= npp; l++) {
                row[l - 1] = 0.0;
            }
            ll = ii;
            for (l = 1; l <= 5; l++) {
                if (ll > npp) { ll = 1; }
                row[ll - 1] += bp[(ii - 1) + (l - 1)*npest];
                ll++;
            }
            facc = 0.0;
            facs = 0.0;
            for (l = 1; l <= npp; l++) {
                facc += row[l - 1] * coco[l - 1];
                facs += row[l - 1] * cosi[l - 1];
            }
            for (j = 1; j <= nt6; j++) {
                // initialize the new row.
                for (l = 1; l <= iband; l++) {
                    h[l - 1] = 0.0;
                }
                // fill in the non-zero elements of the row. jrot records the column
                // number of the first non-zero element in the row.
                jrot = 4 + (j - 2) * npp;
                if ((j > 1) && (j < nt6)) {
                    for (l = 1; l <= npp; l++) {
                        h[l - 1] = row[l - 1];
                    }
                } else {
                    h[0] = facc;
                    h[1] = facs;
                    if (j == 1) { jrot = 2; }
                }
                for (l = 1; l <= iband; l++) {
                    h[l - 1] *= pinv;
                }
                ri = 0.0;
                // rotate the new row into triangle by givens transformations.
                for (irot = jrot; irot <= ncof; irot++) {
                    piv = h[0];
                    i2 = iband1;
                    if (ncof - irot < iband1) { i2 = ncof - irot; }
                    if (piv == 0.0) {
                        if (i2 <= 0) { break; }
                    } else {
                        // calculate the parameters of the givens transformation.
                        fpgivs(piv, &q[irot - 1], &co, &si);
                        // apply that givens transformation to the right hand side.
                        fprota(co, si, &ri, &ff[irot - 1]);
                        // apply that givens transformation to the left hand side.
                        for (l = 1; l <= i2; l++) {
                            fprota(co, si, &h[l], &q[(irot - 1) + l*ncc]);
                        }
                    }
                    for (l = 1; l <= i2; l++) {
                        h[l - 1] = h[l];
                    }
                    h[i2] = 0.0;
                }
            }
        }
        // extend the observation matrix with the rows of a matrix expressing
        // that for phi=cst. sp(teta,phi) must be a cubic polynomial.
        for (i = 5; i <= nt4; i++) {
            ii = i - 4;
            for (j = 1; j <= npp; j++) {
                // initialize the new row
                for (l = 1; l <= iband4; l++) {
                    h[l - 1] = 0.0;
                }
                // fill in the non-zero elements of the row. jrot records the column
                // number of the first non-zero element in the row.
                j1 = 1;
                for (l = 1; l <= 5; l++) {
                    int il = ii + l;
                    int ij = npp;
                    if ((il == 3) || (il == nt4)) {
                        j1 = j1 + 3 - j;
                        j2 = j1 - 2;
                        ij = 0;
                        if (il == 3) {
                            j1 = 1;
                            j2 = 2;
                            ij = j + 2;
                        }
                        h[j2 - 1] = bt[(ii - 1) + (l - 1)*ntest] * coco[j - 1];
                        h[j2] = bt[(ii - 1) + (l - 1)*ntest] * cosi[j - 1];
                    }
                    h[j1 - 1] += bt[(ii - 1) + (l - 1)*ntest];
                    j1 += ij;
                }
                for (l = 1; l <= iband4; l++) {
                    h[l - 1] *= pinv;
                }
                ri = 0.0;
                jrot = 1;
                if (ii > 2) { jrot = 3 + j + (ii - 3) * npp; }
                // rotate the new row into triangle by givens transformations.
                for (irot = jrot; irot <= ncof; irot++) {
                    piv = h[0];
                    i2 = iband3;
                    if (ncof - irot < iband3) { i2 = ncof - irot; }
                    if (piv == 0.0) {
                        if (i2 <= 0) { break; }
                    } else {
                        // calculate the parameters of the givens transformation.
                        fpgivs(piv, &q[irot - 1], &co, &si);
                        // apply that givens transformation to the right hand side.
                        fprota(co, si, &ri, &ff[irot - 1]);
                        // apply that givens transformation to the left hand side.
                        for (l = 1; l <= i2; l++) {
                            fprota(co, si, &h[l], &q[(irot - 1) + l*ncc]);
                        }
                    }
                    for (l = 1; l <= i2; l++) {
                        h[l - 1] = h[l];
                    }
                    h[i2] = 0.0;
                }
            }
        }
        // find dmax, the maximum value for the diagonal elements in the
        // reduced triangle.
        dmax = 0.0;
        for (i = 1; i <= ncof; i++) {
            if (q[i - 1] > dmax) { dmax = q[i - 1]; }
        }
        // check whether the matrix is rank deficient.
        sigma = eps * dmax;
        int rank_deficient = 0;
        for (i = 1; i <= ncof; i++) {
            if (q[i - 1] <= sigma) {
                rank_deficient = 1;
                break;
            }
        }
        if (rank_deficient) {
            // in case of rank deficiency, find the minimum norm solution.
            lwest = ncof * iband4 + ncof + iband4;
            if (lwrk < lwest) {
                *ier = lwest;
                return;
            }
            lf = 1;
            lh = lf + ncof;
            la = lh + iband4;
            fprank(q, ff, ncof, iband4, ncc, sigma, c, &sq, &rank, &wrk[la - 1], &wrk[lf - 1], &wrk[lh - 1]);
        } else {
            // backward substitution in case of full rank.
            fpback(q, ff, ncof, iband4, c, ncc);
            rank = ncof;
        }
        for (i = 1; i <= ncof; i++) {
            q[i - 1] /= dmax;
        }
        // find the coefficients in the standard b-spline representation of
        // the spherical spline.
        fprpsp(*nt, *np, coco, cosi, c, ff, ncoff);
        // compute f(p).
        *fp = 0.0;
        for (num = 1; num <= nreg; num++) {
            num1 = num - 1;
            lt = num1 / npp;
            lp = num1 - lt * npp;
            jrot = lt * np4 + lp;
            in = index[num - 1];
            while (in != 0) {
                store = 0.0;
                i1 = jrot;
                for (i = 1; i <= 4; i++) {
                    hti = spt[(in - 1) + (i - 1)*m];
                    j1 = i1;
                    for (j = 1; j <= 4; j++) {
                        j1++;
                        store += hti * spp[(in - 1) + (j - 1)*m] * c[j1 - 1];
                    }
                    i1 += np4;
                }
                *fp += pow(w[in - 1] * (r[in - 1] - store), 2);
                in = nummer[in - 1];
            }
        }
        // test whether the approximation sp(teta,phi) is an acceptable solution
        fpms = *fp - s;
        if (fabs(fpms) <= acc) {
            if (ncof != rank) {
                *ier = -rank;
            }
            return;
        }
        // test whether the maximum allowable number of iterations has been
        // reached.
        if (iter == maxit) {
            *ier = 3;
            return;
        }
        // carry out one more step of the iteration process.
        p2 = p;
        f2 = fpms;
        if (ich3 == 0) {
            if ((f2 - f3) > acc) {
                if (f2 < 0.0) { ich3 = 1; }
            } else {
                // our initial choice of p is too large.
                p3 = p2;
                f3 = f2;
                p = p * con4;
                if (p <= p1) { p = p1 * con9 + p2 * con1; }
                continue;
            }
        }
        if (ich1 == 0) {
            if ((f1 - f2) > acc) {
                if (f2 > 0.0) { ich1 = 1; }
            } else {
                // our initial choice of p is too small
                p1 = p2;
                f1 = f2;
                p = p / con4;
                if (p3 >= 0.0) {
                    if (p >= p3) { p = p2 * con1 + p3 * con9; }
                }
                continue;
            }
        }
        // test whether the iteration process proceeds as theoretically
        // expected.
        if ((f2 >= f1) || (f2 <= f3)) {
            *ier = 2;
            return;
        }
        // find the new value of p.
        p = fprati(&p1, &f1, &p2, &f2, &p3, &f3);
    }
}


static void
fpspgr(const int* iopt, const int* ider, const double* u, const int mu, const double* v,
       const int mv, const double* r, const int mr, const double r0, const double r1,
       const double s, const int nuest, const int nvest, const double tol, const int maxit,
       const int nc, int* nu, double* tu, int* nv, double* tv, double* c, double* fp,
       double* fp0, double* fpold, double* reducu, double* reducv, double* fpintu,
       double* fpintv, double* dr, double* step, int* lastdi, int* nplusu, int* nplusv,
       int* lastu0, int* lastu1, int* nru, int* nrv, int* nrdatu, int* nrdatv, double* wrk,
       const int lwrk,int* ier)
{
    double acc=0.0, fpms=0.0, f1, f2, f3, p, per, p1, p2, p3, vb, ve, rmax, rmin, rn, one;
    double con1, con4, con9;
    int i, ich1, ich3, ifbu, ifbv, ifsu, ifsv, istart, iter, i1, i2, j, ju;
    int ktu, l, l1, l2, l3, l4, mpm, mumin, mu0, mu1, nn, nplu=0, nplv, npl1, nrintu;
    int nrintv=0, nue=0, numax=0, nve=0, nvmax=0;
    int idd[4];
    double drr[6];

    // set constants
    one = 1.0;
    con1 = 0.1;
    con9 = 0.9;
    con4 = 0.04;
    // initialization
    ifsu = 0;
    ifsv = 0;
    ifbu = 0;
    ifbv = 0;
    p = -one;
    mumin = 4;
    if (ider[0] >= 0) { mumin = mumin - 1; }
    if ((iopt[1] == 1) && (ider[1] == 1)) { mumin = mumin - 1; }
    if (ider[2] >= 0) { mumin = mumin - 1; }
    if ((iopt[2] == 1) && (ider[3] == 1)) { mumin = mumin - 1; }
    if (mumin == 0) { mumin = 1; }
    per = PI + PI;
    vb = v[0];
    ve = vb + per;
    //////////////////////////////////////////////////////////////////////////
    // part 1: determination of the number of knots and their position.     //
    // ****************************************************************     //
    //  given a set of knots we compute the least-squares spline sinf(u,v)  //
    //  and the corresponding sum of squared residuals fp = f(p=inf).       //
    //  if iopt(1)=-1  sinf(u,v) is the requested approximation.            //
    //  if iopt(1)>=0  we check whether we can accept the knots:            //
    //    if fp <= s we will continue with the current set of knots.        //
    //    if fp >  s we will increase the number of knots and compute the   //
    //       corresponding least-squares spline until finally fp <= s.      //
    //    the initial choice of knots depends on the value of s and iopt.   //
    //    if s=0 we have spline interpolation; in that case the number of   //
    //     knots in the u-direction equals nu=numax=mu+6+iopt(2)+iopt(3)    //
    //     and in the v-direction nv=nvmax=mv+7.                            //
    //    if s>0 and                                                        //
    //      iopt(1)=0 we first compute the least-squares polynomial,i.e. a  //
    //       spline without interior knots : nu=8 ; nv=8.                   //
    //      iopt(1)=1 we start with the set of knots found at the last call //
    //       of the routine, except for the case that s > fp0; then we      //
    //       compute the least-squares polynomial directly.                 //
    //////////////////////////////////////////////////////////////////////////
    if (iopt[0] < 0) { goto L120; }

    // acc denotes the absolute tolerance for the root of f(p)=s.
    acc = tol * s;
    // numax and nvmax denote the number of knots needed for interpolation.
    numax = mu + 6 + iopt[1] + iopt[2];
    nvmax = mv + 7;
    nue = (numax < nuest) ? numax : nuest;
    nve = (nvmax < nvest) ? nvmax : nvest;

    if  (s > 0.0) { goto L100; }

    // if s = 0, s(u,v) is an interpolating spline.
    *nu = numax;
    *nv = nvmax;
    // test whether the required storage space exceeds the available one.
    if ((*nu > nuest) || (*nv > nvest)) {
        *ier = 1;
        return;
    }
    // find the position of the knots in the v-direction.
    for (l = 1; l <= mv; l++) {
        tv[l + 2] = v[l - 1];
    }
    tv[mv + 3] = ve;
    l1 = mv - 2;
    l2 = mv + 5;
    for (i = 1; i <= 3; i++) {
        tv[i - 1] = v[l1 - 1] - per;
        tv[l2 - 1] = v[i] + per;
        l1++;
        l2++;
    }

    // if not all the derivative values g(i,j) are given, we will first
    // estimate these values by computing a least-squares spline
    idd[0] = ider[0];
    if (idd[0] == 0) { idd[0] = 1; }
    if (idd[0] > 0) { dr[0] = r0; }
    idd[1] = ider[1];
    idd[2] = ider[2];
    if (idd[2] == 0) { idd[2] = 1; }
    if (idd[2] > 0) { dr[3] = r1; }
    idd[3] = ider[3];

    // Modelling some Fortran code goto olympics here to decide what falls through 30
    //  if(ider(1).lt.0 .or.  ider(3).lt.0) go to 30
    //  if(iopt(2).ne.0 .and. ider(2).eq.0) go to 30
    //  if(iopt(3).eq.0 .or.  ider(4).ne.0) go to 70
    int condition1 = ((ider[0]  < 0) || (ider[2] <  0));
    int condition2 = ((iopt[1] != 0) && (ider[1] == 0));
    int condition3 = ((iopt[2] == 0) || (ider[3] != 0));

    if (condition1 || condition2 || (!condition3)) {
        // we set up the knots in the u-direction for computing the least-squares
        // spline.
        i1 = 3;
        i2 = mu - 2;
        *nu = 4;
        for (i = 1; i <= mu; i++) {
            if (i1 > i2) { break; }
            (*nu)++;
            tu[*nu - 1] = u[i1 - 1];
            i1 += 2;
        }

        for (i = 1; i <= 4; i++) {
            tu[i - 1] = 0.0;
            (*nu)++;
            tu[*nu - 1] = PI;
        }
        // we compute the least-squares spline for estimating the derivatives.
        fpopsp(ifsu, ifsv, ifbu, ifbv, u, mu, v, mv, r, mr, r0, r1, dr, iopt, idd,
                tu, *nu, tv, *nv, nuest, nvest, p, step, c, nc, fp, fpintu, fpintv, nru, nrv,
                wrk, lwrk);
        ifsu = 0;
    }

    // Back to label 70
    // if all the derivatives at the origin are known, we compute the
    // interpolating spline.
    // we set up the knots in the u-direction, needed for interpolation.
    nn = numax - 8;
    if (nn != 0) {
        ju = 2 - iopt[1];
        for (l = 1; l <= nn; l++) {
            tu[l + 3] = u[ju - 1];
            ju++;
        }
        *nu = numax;
        l = *nu;
        for (i = 1; i <= 4; i++) {
            tu[i - 1] = 0.0;
            tu[l - 1] = PI;
            l--;
        }
    }

    // we compute the interpolating spline.
    fpopsp(ifsu, ifsv, ifbu, ifbv, u, mu, v, mv, r, mr, r0, r1, dr, iopt, idd,
            tu, *nu, tv, *nv, nuest, nvest, p, step, c, nc, fp, fpintu, fpintv, nru, nrv,
            wrk, lwrk);
    *ier = -1;
    *fp = 0.0;
    return;

L100:
    *ier = 0;
    if (iopt[0] == 0) { goto L115; }
    step[0] = -step[0];
    step[1] = -step[1];
    if (*fp0 <= 0.0) { goto L115; }
    // if iopt(1)=1 and fp0 > s we start computing the least-squares spline
    // according to the set of knots found at the last call of the routine.
    // we determine the number of grid coordinates u(i) inside each knot
    // interval (tu(l),tu(l+1)).
    l = 5;
    j = 1;
    nrdatu[0] = 0;
    mu0 = 2 - iopt[1];
    mu1 = mu - 1 + iopt[2];
    for (i = mu0; i <= mu1; i++) {
        nrdatu[j - 1] = nrdatu[j - 1] + 1;
        if (u[i - 1] >= tu[l - 1]) {
            nrdatu[j - 1] = nrdatu[j - 1] - 1;
            l = l + 1;
            j = j + 1;
            nrdatu[j - 1] = 0;
        }
    }
    // we determine the number of grid coordinates v(i) inside each knot
    // interval (tv(l),tv(l+1)).
    l = 5;
    j = 1;
    nrdatv[0] = 0;
    for (i = 2; i <= mv; i++) {
        nrdatv[j - 1] = nrdatv[j - 1] + 1;
        if (v[i - 1] >= tv[l - 1]) {
            nrdatv[j - 1] = nrdatv[j - 1] - 1;
            l = l + 1;
            j = j + 1;
            nrdatv[j - 1] = 0;
        }
    }
    idd[0] = ider[0];
    idd[1] = ider[1];
    idd[2] = ider[2];
    idd[3] = ider[3];
    goto L120;
L115:
    // if iopt(1)=0 or iopt(1)=1 and s >= fp0,we start computing the least-
    // squares polynomial (which is a spline without interior knots).
    // 115
    *ier = -2;
    idd[0] = ider[0];
    idd[1] = 1;
    idd[2] = ider[2];
    idd[3] = 1;
    *nu = 8;
    *nv = 8;
    nrdatu[0] = mu - 2 + iopt[1] + iopt[2];
    nrdatv[0] = mv - 1;
    *lastdi = 0;
    *nplusu = 0;
    *nplusv = 0;
    *fp0 = 0.0;
    *fpold = 0.0;
    *reducu = 0.0;
    *reducv = 0.0;
L120:
    // main loop for the different sets of knots.mpm=mu+mv is a save upper
    // bound for the number of trials.
    mpm = mu + mv;
    for (iter = 1; iter <= mpm; iter++) {
        // find nrintu (nrintv) which is the number of knot intervals in the
        // u-direction (v-direction).
        nrintu = *nu - 7;
        nrintv = *nv - 7;
        // find the position of the additional knots which are needed for the
        // b-spline representation of s(u,v).
        i = *nu;
        for (j = 1; j <= 4; j++) {
            tu[j - 1] = 0.0;
            tu[i - 1] = PI;
            i--;
        }
        l1 = 4;
        l2 = l1;
        l3 = *nv - 3;
        l4 = l3;
        tv[l2 - 1] = vb;
        tv[l3 - 1] = ve;
        for (j = 1; j <= 3; j++) {
            l1++;
            l2--;
            l3++;
            l4--;
            tv[l2 - 1] = tv[l4 - 1] - per;
            tv[l3 - 1] = tv[l1 - 1] + per;
        }
        // find an estimate of the range of possible values for the optimal
        // derivatives at the origin.
        ktu = nrdatu[0] + 2 - iopt[1];
        if (ktu < mumin) { ktu = mumin; }
        if (ktu != *lastu0) {
            rmin = r0;
            rmax = r0;
            l = mv * ktu;
            for (i = 1; i <= l; i++) {
                if (r[i - 1] < rmin) { rmin = r[i - 1]; }
                if (r[i - 1] > rmax) { rmax = r[i - 1]; }
            }
            step[0] = rmax - rmin;
            *lastu0 = ktu;
        }
        ktu = nrdatu[nrintu - 1] + 2 - iopt[2];
        if (ktu < mumin) { ktu = mumin; }
        if (ktu != *lastu1) {
            rmin = r1;
            rmax = r1;
            l = mv * ktu;
            j = mr;
            for (i = 1; i <= l; i++) {
                if (r[j - 1] < rmin) { rmin = r[j - 1]; }
                if (r[j - 1] > rmax) { rmax = r[j - 1]; }
                j--;
            }
            step[1] = rmax - rmin;
            *lastu1 = ktu;
        }
        // find the least-squares spline sinf(u,v).
        fpopsp(ifsu, ifsv, ifbu, ifbv, u, mu, v, mv, r, mr, r0, r1, dr, iopt, idd,
               tu, *nu, tv, *nv, nuest, nvest, p, step, c, nc, fp, fpintu, fpintv, nru, nrv,
               wrk, lwrk);
        if (step[0] < 0.0) { step[0] = -step[0]; }
        if (step[1] < 0.0) { step[1] = -step[1]; }
        if (*ier == -2) { *fp0 = *fp; }
        // test whether the least-squares spline is an acceptable solution.
        if (iopt[0] < 0) { return; }
        fpms = *fp - s;
        if (fabs(fpms) < acc) { return; }
        // if f(p=inf) < s, we accept the choice of knots.
        if (fpms < 0.0) { break; }
        // if nu=numax and nv=nvmax, sinf(u,v) is an interpolating spline
        if ((*nu == numax) && (*nv == nvmax)) {
            *ier = -1;
            *fp = 0.0;
            return;
        }
        // increase the number of knots.
        // if nu=nue and nv=nve we cannot further increase the number of knots
        // because of the storage capacity limitation.
        if ((*nu == nue) && (*nv == nve)) {
            *ier = 1;
            return;
        }
        if (ider[0] == 0) { fpintu[0] = fpintu[0] + (r0 - dr[0]) * (r0 - dr[0]); }
        if (ider[2] == 0) { fpintu[nrintu - 1] = fpintu[nrintu - 1] + (r1 - dr[3]) * (r1 - dr[3]); }
        *ier = 0;
        // adjust the parameter reducu or reducv according to the direction
        // in which the last added knots were located.

        // More goto olympics
        if (*lastdi == 0) {
            nplv = 3;
            idd[1] = ider[1];
            idd[3] = ider[3];
            *fpold = *fp;

            goto L230;

        }
        if (*lastdi < 0) {
            *reducu = *fpold - *fp;
        } else {
            *reducv = *fpold - *fp;
        }
        // store the sum of squared residuals for the current set of knots.
        *fpold = *fp;
        // find nplu, the number of knots we should add in the u-direction.
        nplu = 1;
        if (*nu != 8) {
            npl1 = *nplusu * 2;
            rn = (double)(*nplusu);
            if (*reducu > acc) { npl1 = (int)(rn * fpms / (*reducu)); }
            int max1 = (npl1 > *nplusu / 2) ? npl1 : *nplusu / 2;
            int max2 = (max1 > 1) ? max1 : 1;
            nplu = (*nplusu * 2 < max2) ? *nplusu * 2 : max2;
        }
        // find nplv, the number of knots we should add in the v-direction.
        nplv = 3;
        if (*nv != 8) {
            npl1 = *nplusv * 2;
            rn = (double)(*nplusv);
            if (*reducv > acc) { npl1 = (int)(rn * fpms / (*reducv)); }
            int max1 = (npl1 > *nplusv / 2) ? npl1 : *nplusv / 2;
            int max2 = (max1 > 1) ? max1 : 1;
            nplv = (*nplusv * 2 < max2) ? *nplusv * 2 : max2;
        }
        // test whether we are going to add knots in the u- or v-direction.
        if (nplu > nplv) { goto L230; }
        if (nplu < nplv) { goto L210; }
        if ((nplu == nplv) && (*lastdi > 0)) { goto L230; }
L210:
        if (*nu == nue) { goto L230; }
        // addition in the u-direction.
        *lastdi = -1;
        *nplusu = nplu;
        ifsu = 0;
        istart = 0;
        if (iopt[1] == 0) { istart = 1; }
        for (l = 1; l <= *nplusu; l++) {
            // add a new knot in the u-direction
            fpknot(u, mu, tu, nu, fpintu, nrdatu, &nrintu, nuest, istart);
            // test whether we cannot further increase the number of knots in the
            // u-direction.
            if (*nu == nue) { break; }
        }
        continue;
L230:
        if (*nv == nve) { goto L210; }  // What?
        // addition in the v-direction.
        *lastdi = 1;
        *nplusv = nplv;
        ifsv = 0;
        for (l = 1; l <= *nplusv; l++) {
            // add a new knot in the v-direction.
            fpknot(v, mv, tv, nv, fpintv, nrdatv, &nrintv, nvest, 1);
            // test whether we cannot further increase the number of knots in the
            // v-direction.
            if (*nv == nve) {
                break;
            }
        }
        continue;
    }
    // test whether the least-squares polynomial is a solution of our
    // approximation problem.
    if (*ier == -2) { return; }
    //////////////////////////////////////////////////////////////////////////
    // part 2: determination of the smoothing spline sp(u,v)                //
    // *****************************************************                //
    //  we have determined the number of knots and their position. we now   //
    //  compute the b-spline coefficients of the smoothing spline sp(u,v).  //
    //  this smoothing spline depends on the parameter p in such a way that //
    //    f(p) = sumi=1,mu(sumj=1,mv((z(i,j)-sp(u(i),v(j)))**2)             //
    //  is a continuous, strictly decreasing function of p. moreover the    //
    //  least-squares polynomial corresponds to p=0 and the least-squares   //
    //  spline to p=infinity. then iteratively we have to determine the     //
    //  positive value of p such that f(p)=s. the process which is proposed //
    //  here makes use of rational interpolation. f(p) is approximated by a //
    //  rational function r(p)=(u*p+v)/(p+w); three values of p (p1,p2,p3)  //
    //  with corresponding values of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s)//
    //  are used to calculate the new value of p such that r(p)=s.          //
    //  convergence is guaranteed by taking f1 > 0 and f3 < 0.              //
    //////////////////////////////////////////////////////////////////////////
    // initial value for p.
    p1 = 0.0;
    f1 = *fp0 - s;
    p3 = -one;
    f3 = fpms;
    p = one;
    for (i = 1; i <= 6; i++) {
        drr[i - 1] = dr[i - 1];
    }
    ich1 = 0;
    ich3 = 0;
    // iteration process to find the root of f(p)=s.
    for (iter = 1; iter <= maxit; iter++) {
        // find the smoothing spline sp(u,v) and the corresponding sum f(p).
        fpopsp(ifsu, ifsv, ifbu, ifbv, u, mu, v, mv, r, mr, r0, r1, drr, iopt, idd,
               tu, *nu, tv, *nv, nuest, nvest, p, step, c, nc, fp, fpintu, fpintv, nru, nrv,
               wrk, lwrk);
        // test whether the approximation sp(u,v) is an acceptable solution.
        fpms = *fp - s;
        if (fabs(fpms) < acc) {
            return;
        }
        // test whether the maximum allowable number of iterations has been
        // reached.
        if (iter == maxit) {
            *ier = 3;
            return;
        }
        // carry out one more step of the iteration process.
        p2 = p;
        f2 = fpms;
        if (ich3 == 0) {
            if ((f2 - f3) > acc) {
                if (f2 < 0.0) { ich3 = 1; }
            } else {
                // our initial choice of p is too large.
                p3 = p2;
                f3 = f2;
                p = p * con4;
                if (p <= p1) { p = p1 * con9 + p2 * con1; }
                continue;
            }
        }
        if (ich1 == 0) {
            if ((f1 - f2) > acc) {
                if (f2 > 0.0) { ich1 = 1; }
            } else {
                // our initial choice of p is too small
                p1 = p2;
                f1 = f2;
                p = p / con4;
                if (p3 >= 0.0) {
                    if (p >= p3) { p = p2 * con1 + p3 * con9; }
                }
                continue;
            }
        }
        // test whether the iteration process proceeds as theoretically
        // expected.
        if ((f2 >= f1) || (f2 <= f3)) {
            *ier = 2;
            return;
        }
        // find the new value of p.
        p = fprati(&p1, &f1, &p2, &f2, &p3, &f3);
    }
}


static void
fpsurf(int iopt, int m, double* x, double* y, double* z, double* w,
       double xb, double xe, double yb, double ye, int kxx, int kyy,
       double s, int nxest, int nyest, double eta, double tol, int maxit,
       int nmax, int km1, int km2, int ib1, int ib3, int nc, int intest,
       int nrest, int* nx0, double* tx, int* ny0, double* ty, double* c,
       double* fp, double* fp0, double* fpint, double* coord, double* f,
       double* ff, double* a, double* q, double* bx, double* by, double* spx,
       double* spy, double* h, int* index, int* nummer, double* wrk,
       int lwrk, int* ier)
{
    (void)km1;      // unused
    (void)km2;      // unused
    (void)ib1;      // unused
    (void)ib3;      // unused
    (void)intest;   // unused
    (void)nrest;    // unused
    double acc, arg, cos, dmax, fac1, fac2, fpmax, fpms, f1, f2, f3, hxi, p, pinv;
    double piv, p1, p2, p3, sigma, sin, sq, store, wi, x0, x1, y0, y1, zi, eps;
    int i, iband, iband1, iband3, iband4, ibb, ichang, ich1, ich3, ii;
    int in, irot, iter, i1, i2, i3, j, jrot, jxy, j1, kx, kx1, kx2, ky, ky1, ky2, l;
    int la, lf, lh, lwest, lx, ly, l1, l2, n, ncof, nk1x, nk1y, nminx, nminy, nreg;
    int nrint, num, num1, nx, nxe, nxx, ny, nye, nyy, n1, rank;
    // ..local arrays..
    double hx[6], hy[6];

    double one = 1.0, con1 = 0.1, con9 = 0.9, con4 = 0.04, ten = 10.0;

    /////////////////////////////////////////////////////////////////////////////
    // part 1: determination of the number of knots and their position.        //
    // ****************************************************************        //
    // given a set of knots we compute the least-squares spline sinf(x,y),     //
    // and the corresponding weighted sum of squared residuals fp=f(p=inf).    //
    // if iopt=-1  sinf(x,y) is the requested approximation.                   //
    // if iopt=0 or iopt=1 we check whether we can accept the knots:           //
    //   if fp <=s we will continue with the current set of knots.             //
    //   if fp > s we will increase the number of knots and compute the        //
    //      corresponding least-squares spline until finally  fp<=s.           //
    // the initial choice of knots depends on the value of s and iopt.         //
    //   if iopt=0 we first compute the least-squares polynomial of degree     //
    //     kx in x and ky in y; nx=nminx=2*kx+2 and ny=nminy=2*ky+2.           //
    //     fp0=f(0) denotes the corresponding weighted sum of squared          //
    //     residuals                                                           //
    //   if iopt=1 we start with the knots found at the last call of the       //
    //     routine, except for the case that s>=fp0; then we can compute       //
    //     the least-squares polynomial directly.                              //
    // eventually the independent variables x and y (and the corresponding     //
    // parameters) will be switched if this can reduce the bandwidth of the    //
    // system to be solved.                                                    //
    /////////////////////////////////////////////////////////////////////////////

    // ichang denotes whether(1) or not(-1) the directions have been inter-
    // changed.
    ichang = -1;
    x0 = xb;
    x1 = xe;
    y0 = yb;
    y1 = ye;
    kx = kxx;
    ky = kyy;
    kx1 = kx + 1;
    ky1 = ky + 1;
    nxe = nxest;
    nye = nyest;
    eps = sqrt(eta);

    if (iopt < 0) {
        nx = *nx0;
        ny = *ny0;
    } else if (iopt == 0) {
        // calculation of acc, the absolute tolerance for the root of f(p)=s.
        acc = tol * s;
        // initialization for the least-squares polynomial.
        nminx = 2 * kx1;
        nminy = 2 * ky1;
        nx = nminx;
        ny = nminy;
        *ier = -2;
    } else {
        acc = tol * s;
        if (*fp0 > s) {
            nx = *nx0;
            ny = *ny0;
        } else {
            nminx = 2 * kx1;
            nminy = 2 * ky1;
            nx = nminx;
            ny = nminy;
            *ier = -2;
        }
    }

    // main loop for the different sets of knots. m is a save upper bound
    // for the number of trials.
    for (iter = 1; iter <= m; iter++) {
        // find the position of the additional knots which are needed for the
        // b-spline representation of s(x,y).
        l = nx;
        for (i = 1; i <= kx1; i++) {
            tx[i - 1] = x0;
            tx[l - 1] = x1;
            l--;
        }
        l = ny;
        for (i = 1; i <= ky1; i++) {
            ty[i - 1] = y0;
            ty[l - 1] = y1;
            l--;
        }

        // find nrint, the total number of knot intervals and nreg, the number
        // of panels in which the approximation domain is subdivided by the
        // intersection of knots.
        nxx = nx - 2 * kx1 + 1;
        nyy = ny - 2 * ky1 + 1;
        nrint = nxx + nyy;
        nreg = nxx * nyy;

        // find the bandwidth of the observation matrix a.
        // if necessary, interchange the variables x and y, in order to obtain
        // a minimal bandwidth.
        iband1 = kx * (ny - ky1) + ky;
        l = ky * (nx - kx1) + kx;

        if (iband1 > l) {
            iband1 = l;
            ichang = -ichang;
            for (i = 1; i <= m; i++) {
                store = x[i - 1];
                x[i - 1] = y[i - 1];
                y[i - 1] = store;
            }
            store = x0;
            x0 = y0;
            y0 = store;
            store = x1;
            x1 = y1;
            y1 = store;
            n = (nx < ny) ? nx : ny;
            for (i = 1; i <= n; i++) {
                store = tx[i - 1];
                tx[i - 1] = ty[i - 1];
                ty[i - 1] = store;
            }
            n1 = n + 1;
            if (nx < ny) {
                for (i = n1; i <= ny; i++) {
                    tx[i - 1] = ty[i - 1];
                }
            } else if (nx > ny) {
                for (i = n1; i <= nx; i++) {
                    ty[i - 1] = tx[i - 1];
                }
            }
            l = nx;
            nx = ny;
            ny = l;
            l = nxe;
            nxe = nye;
            nye = l;
            l = nxx;
            nxx = nyy;
            nyy = l;
            l = kx;
            kx = ky;
            ky = l;
            kx1 = kx + 1;
            ky1 = ky + 1;
        }

        iband = iband1 + 1;

        // arrange the data points according to the panel they belong to.
        fporde(x, y, m, kx, ky, tx, nx, ty, ny, nummer, index, nreg);

        // find ncof, the number of b-spline coefficients.
        nk1x = nx - kx1;
        nk1y = ny - ky1;
        ncof = nk1x * nk1y;

        // initialize the observation matrix a.
        for (i = 1; i <= ncof; i++) {
            f[i - 1] = 0.0;
            for (j = 1; j <= iband; j++) {
                a[(i - 1) + (j - 1) * nc] = 0.0;
            }
        }

        // initialize the sum of squared residuals.
        *fp = 0.0;

        // fetch the data points in the new order. main loop for the
        // different panels.
        for (num = 1; num <= nreg; num++) {
            // fix certain constants for the current panel; jrot records the column
            // number of the first non-zero element in a row of the observation
            // matrix according to a data point of the panel.
            num1 = num - 1;
            lx = num1 / nyy;
            l1 = lx + kx1;
            ly = num1 - lx * nyy;
            l2 = ly + ky1;
            jrot = lx * nk1y + ly;

            // test whether there are still data points in the panel.
            in = index[num - 1];
            while (in != 0) {
                // fetch a new data point.
                wi = w[in - 1];
                zi = z[in - 1] * wi;

                // evaluate for the x-direction, the (kx+1) non-zero b-splines at x(in).
                fpbspl(tx, nx, kx, x[in - 1], l1, hx);

                // evaluate for the y-direction, the (ky+1) non-zero b-splines at y(in).
                fpbspl(ty, ny, ky, y[in - 1], l2, hy);

                // store the value of these b-splines in spx and spy respectively.
                for (i = 1; i <= kx1; i++) {
                    spx[(in - 1) + (i - 1) * m] = hx[i - 1];
                }
                for (i = 1; i <= ky1; i++) {
                    spy[(in - 1) + (i - 1) * m] = hy[i - 1];
                }

                // initialize the new row of observation matrix.
                for (i = 1; i <= iband; i++) {
                    h[i - 1] = 0.0;
                }

                // calculate the non-zero elements of the new row by making the cross
                // products of the non-zero b-splines in x- and y-direction.
                i1 = 0;
                for (i = 1; i <= kx1; i++) {
                    hxi = hx[i - 1];
                    j1 = i1;
                    for (j = 1; j <= ky1; j++) {
                        j1++;
                        h[j1 - 1] = hxi * hy[j - 1] * wi;
                    }
                    i1 = i1 + nk1y;
                }

                // rotate the row into triangle by givens transformations.
                irot = jrot;
                for (i = 1; i <= iband; i++) {
                    irot++;
                    piv = h[i - 1];
                    if (piv == 0.0) { continue; }

                    // calculate the parameters of the givens transformation.
                    fpgivs(piv, &a[irot - 1], &cos, &sin);

                    // apply that transformation to the right hand side.
                    fprota(cos, sin, &zi, &f[irot - 1]);

                    if (i == iband) {
                        break;
                    }

                    // apply that transformation to the left hand side.
                    i2 = 1;
                    i3 = i + 1;
                    for (j = i3; j <= iband; j++) {
                        i2++;
                        fprota(cos, sin, &h[j - 1], &a[(irot - 1) + (i2 - 1) * nc]);
                    }
                }

                // add the contribution of the row to the sum of squares of residual
                // right hand sides.
                *fp = *fp + zi * zi;

                // find the number of the next data point in the panel.
                in = nummer[in - 1];
            }
        }

        // find dmax, the maximum value for the diagonal elements in the reduced
        // triangle.
        dmax = 0.0;
        for (i = 1; i <= ncof; i++) {
            if (a[i - 1] > dmax) {
                dmax = a[(i - 1)];
            }
        }

        // check whether the observation matrix is rank deficient.
        sigma = eps * dmax;
        int rank_deficient = 0;
        for (i = 1; i <= ncof; i++) {
            if (a[(i - 1)] <= sigma) {
                rank_deficient = 1;
                break;
            }
        }

        if (rank_deficient) {
            // in case of rank deficiency, find the minimum norm solution.
            // check whether there is sufficient working space
            lwest = ncof * iband + ncof + iband;
            if (lwrk < lwest) {
                *ier = lwest;
                goto interpolating_solution;
            }
            for (i = 1; i <= ncof; i++) {
                ff[i - 1] = f[i - 1];
                for (j = 1; j <= iband; j++) {
                    q[(i - 1) + (j - 1) * nc] = a[(i - 1) + (j - 1) * nc];
                }
            }
            lf = 1;
            lh = lf + ncof;
            la = lh + iband;
            fprank(q, ff, ncof, iband, nc, sigma, c, &sq, &rank, &wrk[la - 1], &wrk[lf - 1], &wrk[lh - 1]);
            for (i = 1; i <= ncof; i++) {
                q[(i - 1)] = q[(i - 1)] / dmax;
            }
            // add to the sum of squared residuals, the contribution of reducing
            // the rank.
            *fp = *fp + sq;
        } else {
            // backward substitution in case of full rank.
            fpback(a, f, ncof, iband, c, nc);
            rank = ncof;
            for (i = 1; i <= ncof; i++) {
                q[i - 1] = a[i - 1] / dmax;
            }
        }

        if (*ier == -2) {
            *fp0 = *fp;
        }

        // test whether the least-squares spline is an acceptable solution.
        if (iopt < 0) {
            if (ncof != rank) {
               *ier = -rank;
            }
            goto interpolating_solution;
        }
        fpms = *fp - s;
        if (fabs(fpms) <= acc) {
            if (*fp <= 0.0) {
                *ier = -1;
                *fp = 0.0;
            }
            if (ncof != rank) {
               *ier = -rank;
            }
            goto interpolating_solution;
        }

        // test whether we can accept the choice of knots.
        if (fpms < 0.0) {
            break;
        }

        // test whether we cannot further increase the number of knots.
        if (ncof > m) {
            *ier = 4;
            goto interpolating_solution;
        }
        *ier = 0;

        // search where to add a new knot.
        // find for each interval the sum of squared residuals fpint for the
        // data points having the coordinate belonging to that knot interval.
        // calculate also coord which is the same sum, weighted by the position
        // of the data points considered.
        for (i = 1; i <= nrint; i++) {
            fpint[i - 1] = 0.0;
            coord[i - 1] = 0.0;
        }

        for (num = 1; num <= nreg; num++) {
            num1 = num - 1;
            lx = num1 / nyy;
            l1 = lx + 1;
            ly = num1 - lx * nyy;
            l2 = ly + 1 + nxx;
            jrot = lx * nk1y + ly;
            in = index[num - 1];

            while (in != 0) {
                store = 0.0;
                i1 = jrot;
                for (i = 1; i <= kx1; i++) {
                    hxi = spx[(in - 1) + (i - 1) * m];
                    j1 = i1;
                    for (j = 1; j <= ky1; j++) {
                        j1++;
                        store = store + hxi * spy[(in - 1) + (j - 1) * m] * c[j1 - 1];
                    }
                    i1 = i1 + nk1y;
                }
                store = pow(w[in - 1] * (z[in - 1] - store), 2);
                fpint[l1 - 1] = fpint[l1 - 1] + store;
                coord[l1 - 1] = coord[l1 - 1] + store * x[in - 1];
                fpint[l2 - 1] = fpint[l2 - 1] + store;
                coord[l2 - 1] = coord[l2 - 1] + store * y[in - 1];
                in = nummer[in - 1];
            }
        }

        // find the interval for which fpint is maximal on the condition that
        // there still can be added a knot.
        while (1) {
            l = 0;
            fpmax = 0.0;
            l1 = 1;
            l2 = nrint;
            if (nx == nxe) {
                l1 = nxx + 1;
            }
            if (ny == nye) {
                l2 = nxx;
            }
            if (l1 > l2) {
                *ier = 1;
                goto interpolating_solution;
            }
            for (i = l1; i <= l2; i++) {
                if (fpmax < fpint[i - 1]) {
                    l = i;
                    fpmax = fpint[i - 1];
                }
            }

            // test whether we cannot further increase the number of knots.
            if (l == 0) {
                *ier = 5;
                goto interpolating_solution;
            }

            // calculate the position of the new knot.
            arg = coord[l - 1] / fpint[l - 1];

            // test in what direction the new knot is going to be added.
            if (l <= nxx) {
                // addition in the x-direction.
                jxy = l + kx1;
                fpint[l - 1] = 0.0;
                fac1 = tx[jxy - 1] - arg;
                fac2 = arg - tx[jxy - 2];
                if ((fac1 > (ten * fac2)) || (fac2 > (ten * fac1))) { continue; }
                j = nx;
                for (i = jxy; i <= nx; i++) {
                    tx[j] = tx[j - 1];
                    j--;
                }
                tx[jxy - 1] = arg;
                nx++;
                break;
            } else {
                // addition in the y-direction.
                jxy = l + ky1 - nxx;
                fpint[l - 1] = 0.0;
                fac1 = ty[jxy - 1] - arg;
                fac2 = arg - ty[jxy - 2];
                if ((fac1 > (ten * fac2)) || (fac2 > (ten * fac1))) { continue; }
                j = ny;
                for (i = jxy; i <= ny; i++) {
                    ty[j] = ty[j - 1];
                    j--;
                }
                ty[jxy - 1] = arg;
                ny++;
                break;
            }
        }
    }

    // test whether the least-squares polynomial is a solution of our
    // approximation problem.
    if (*ier == -2) {
        goto interpolating_solution;
    }

    /////////////////////////////////////////////////////////////////////////////
    // part 2: determination of the smoothing spline sp(x,y)                   //
    // *****************************************************                   //
    // we have determined the number of knots and their position. we now       //
    // compute the b-spline coefficients of the smoothing spline sp(x,y).      //
    // the observation matrix a is extended by the rows of a matrix,           //
    // expressing that sp(x,y) must be a polynomial of degree kx in x and      //
    // ky in y. the corresponding weights of these additional rows are set     //
    // to 1./p.  iteratively we than have to determine the value of p          //
    // such that f(p)=sum((w(i)*(z(i)-sp(x(i),y(i))))**2) be = s.              //
    // we already know that the least-squares polynomial corresponds to        //
    // p=0  and that the least-squares spline corresponds to p=infinity.       //
    // the iteration process which is proposed here makes use of rational      //
    // interpolation. since f(p) is a convex and strictly decreasing           //
    // function of p, it can be approximated by a rational function r(p)=      //
    // (u*p+v)/(p+w). three values of p(p1,p2,p3) with corresponding values    //
    // of f(p) (f1=f(p1)-s,f2=f(p2)-s,f3=f(p3)-s) are used to calculate the    //
    // new value of p such that r(p)=s. convergence is guaranteed by taking    //
    // f1 > 0 and f3 < 0.                                                      //
    /////////////////////////////////////////////////////////////////////////////

    kx2 = kx1 + 1;
    // test whether there are interior knots in the x-direction.
    if (nk1x != kx1) {
        // evaluate the discotinuity jumps of the kx-th order derivative of
        // the b-splines at the knots tx(l),l=kx+2,...,nx-kx-1.
        fpdisc(tx, nx, kx2, bx, nmax);
    }
    ky2 = ky1 + 1;
    // test whether there are interior knots in the y-direction.
    if (nk1y != ky1) {
        // evaluate the discontinuity jumps of the ky-th order derivative of
        // the b-splines at the knots ty(l),l=ky+2,...,ny-ky-1.
        fpdisc(ty, ny, ky2, by, nmax);
    }

    // initial value for p.
    p1 = 0.0;
    f1 = *fp0 - s;
    p3 = -one;
    f3 = fpms;
    p = 0.0;
    for (i = 1; i <= ncof; i++) {
        p = p + a[i - 1];
    }
    p = ((double)ncof) / p;

    // find the bandwidth of the extended observation matrix.
    iband3 = kx1 * nk1y;
    iband4 = iband3 + 1;
    ich1 = 0;
    ich3 = 0;

    // iteration process to find the root of f(p)=s.
    for (iter = 1; iter <= maxit; iter++) {
        pinv = one / p;

        // store the triangularized observation matrix into q.
        for (i = 1; i <= ncof; i++) {
            ff[i - 1] = f[i - 1];
            for (j = 1; j <= iband; j++) {
                q[(i - 1) + (j - 1) * nc] = a[(i - 1) + (j - 1) * nc];
            }
            ibb = iband + 1;
            for (j = ibb; j <= iband4; j++) {
                q[(i - 1) + (j - 1) * nc] = 0.0;
            }
        }

        if (nk1y != ky1) {
            // extend the observation matrix with the rows of a matrix, expressing
            // that for x=cst. sp(x,y) must be a polynomial in y of degree ky.
            for (i = ky2; i <= nk1y; i++) {
                ii = i - ky1;
                for (j = 1; j <= nk1x; j++) {
                    // initialize the new row.
                    for (l = 0; l < iband; l++) {
                        h[l] = 0.0;
                    }
                    // fill in the non-zero elements of the row. jrot records the column
                    // number of the first non-zero element in the row.
                    for (l = 0; l < ky2; l++) {
                        h[l] = by[(ii - 1) + l * nmax] * pinv;
                    }
                    zi = 0.0;
                    jrot = (j - 1) * nk1y + ii;

                    // rotate the new row into triangle by givens transformations without
                    // square roots.
                    for (irot = jrot; irot <= ncof; irot++) {
                        piv = h[0];
                        i2 = (iband1 < (ncof - irot)) ? iband1 : (ncof - irot);
                        if (piv != 0.0) {
                            // calculate the parameters of the givens transformation.
                            fpgivs(piv, &q[irot - 1], &cos, &sin);
                            // apply that givens transformation to the right hand side.
                            fprota(cos, sin, &zi, &ff[irot - 1]);
                            if (i2 == 0) {
                                break;
                            }
                            // apply that givens transformation to the left hand side.
                            for (l = 1; l <= i2; l++) {
                                fprota(cos, sin, &h[l], &q[(irot - 1) + l * nc]);
                            }
                        } else {
                            if (i2 <= 0) { break; }
                        }
                        for (l = 0; l < i2; l++) {
                            h[l] = h[l + 1];
                        }
                        h[i2] = 0.0;
                    }
                }
            }
        }

        if (nk1x != kx1) {
            // extend the observation matrix with the rows of a matrix expressing
            // that for y=cst. sp(x,y) must be a polynomial in x of degree kx.
            for (i = kx2; i <= nk1x; i++) {
                ii = i - kx1;
                for (j = 1; j <= nk1y; j++) {
                    // initialize the new row
                    for (l = 0; l < iband4; l++) {
                        h[l] = 0.0;
                    }
                    // fill in the non-zero elements of the row. jrot records the column
                    // number of the first non-zero element in the row.
                    j1 = 1;
                    for (l = 1; l <= kx2; l++) {
                        h[j1 - 1] = bx[(ii - 1) + (l - 1) * nmax] * pinv;
                        j1 = j1 + nk1y;
                    }
                    zi = 0.0;
                    jrot = (i - kx2) * nk1y + j;

                    // rotate the new row into triangle by givens transformations.
                    for (irot = jrot; irot <= ncof; irot++) {
                        piv = h[0];
                        i2 = (iband3 < (ncof - irot)) ? iband3 : (ncof - irot);
                        if (piv != 0.0) {
                            // calculate the parameters of the givens transformation.
                            fpgivs(piv, &q[irot - 1], &cos, &sin);
                            // apply that givens transformation to the right hand side.
                            fprota(cos, sin, &zi, &ff[irot - 1]);
                            if (i2 == 0) { break; }
                            // apply that givens transformation to the left hand side.
                            for (l = 1; l <= i2; l++) {
                                fprota(cos, sin, &h[l], &q[(irot - 1) + l * nc]);
                            }
                        } else {
                            if (i2 <= 0) { break; }
                        }
                        for (l = 0; l < i2; l++) {
                            h[l] = h[l + 1];
                        }
                        h[i2] = 0.0;
                    }
                }
            }
        }

        // find dmax, the maximum value for the diagonal elements in the
        // reduced triangle.
        dmax = 0.0;
        for (i = 1; i <= ncof; i++) {
            if (q[i - 1] > dmax) {
                dmax = q[i - 1];
            }
        }

        int rank_deficient = 0;
        // check whether the matrix is rank deficient.
        sigma = eps * dmax;
        rank_deficient = 0;
        for (i = 1; i <= ncof; i++) {
            if (q[i - 1] <= sigma) {
                rank_deficient = 1;
                break;
            }
        }

        if (rank_deficient) {
            // in case of rank deficiency, find the minimum norm solution.
            lwest = ncof * iband4 + ncof + iband4;
            if (lwrk < lwest) {
                *ier = lwest;
                goto interpolating_solution;
            }
            lf = 1;
            lh = lf + ncof;
            la = lh + iband4;
            fprank(q, ff, ncof, iband4, nc, sigma, c, &sq, &rank, &wrk[la - 1], &wrk[lf - 1], &wrk[lh - 1]);
        } else {
            // backward substitution in case of full rank.
            fpback(q, ff, ncof, iband4, c, nc);
            rank = ncof;
        }

        for (i = 0; i < ncof; i++) {
            q[i] = q[i] / dmax;
        }

        // compute f(p).
        *fp = 0.0;
        for (num = 1; num <= nreg; num++) {
            num1 = num - 1;
            lx = num1 / nyy;
            ly = num1 - lx * nyy;
            jrot = lx * nk1y + ly;
            in = index[num - 1];

            while (in != 0) {
                store = 0.0;
                i1 = jrot;
                for (i = 1; i <= kx1; i++) {
                    hxi = spx[(in - 1) + (i - 1) * m];
                    j1 = i1;
                    for (j = 1; j <= ky1; j++) {
                        j1++;
                        store = store + hxi * spy[(in - 1) + (j - 1) * m] * c[j1 - 1];
                    }
                    i1 = i1 + nk1y;
                }
                *fp = *fp + pow(w[in - 1] * (z[in - 1] - store), 2);
                in = nummer[in - 1];
            }
        }

        // test whether the approximation sp(x,y) is an acceptable solution.
        fpms = *fp - s;
        if (fabs(fpms) <= acc) {
            if (ncof != rank) {
               *ier = -rank;
            }
            goto interpolating_solution;
        }

        // test whether the maximum allowable number of iterations has been
        // reached.
        if (iter == maxit) {
            *ier = 3;
            goto interpolating_solution;
        }

        // carry out one more step of the iteration process.
        p2 = p;
        f2 = fpms;
        if (ich3 == 0) {
            if ((f2 - f3) > acc) {
                if (f2 < 0.0) {
                    ich3 = 1;
                }
            } else {
                // our initial choice of p is too large.
                p3 = p2;
                f3 = f2;
                p = p * con4;
                if (p <= p1) {
                    p = p1 * con9 + p2 * con1;
                }
                continue;
            }
        }

        if (ich1 == 0) {
            if ((f1 - f2) > acc) {
                if (f2 > 0.0) {
                    ich1 = 1;
                }
            } else {
                // our initial choice of p is too small
                p1 = p2;
                f1 = f2;
                p = p / con4;
                if (p3 >= 0.0 && p >= p3) {
                    p = p2 * con1 + p3 * con9;
                }
                continue;
            }
        }

        // test whether the iteration process proceeds as theoretically
        // expected.
        if (f2 >= f1 || f2 <= f3) {
            *ier = 2;
            goto interpolating_solution;
        }

        // find the new value of p.
        p = fprati(&p1, &f1, &p2, &f2, &p3, &f3);
    }

interpolating_solution:
    // test whether x and y are in the original order.
    if (ichang >= 0) {
        // if not, interchange x and y once more.
        l1 = 1;
        for (i = 1; i <= nk1x; i++) {
            l2 = i;
            for (j = 1; j <= nk1y; j++) {
                f[l2 - 1] = c[l1 - 1];
                l1++;
                l2 = l2 + nk1x;
            }
        }
        for (i = 1; i <= ncof; i++) {
            c[i - 1] = f[i - 1];
        }
        for (i = 1; i <= m; i++) {
            store = x[i - 1];
            x[i - 1] = y[i - 1];
            y[i - 1] = store;
        }
        n = (nx < ny) ? nx : ny;
        for (i = 1; i <= n; i++) {
            store = tx[i - 1];
            tx[i - 1] = ty[i - 1];
            ty[i - 1] = store;
        }
        n1 = n + 1;
        if (nx < ny) {
            for (i = n1; i <= ny; i++) {
                tx[i - 1] = ty[i - 1];
            }
        } else if (nx > ny) {
            for (i = n1; i <= nx; i++) {
                ty[i - 1] = tx[i - 1];
            }
        }
        l = nx;
        nx = ny;
        ny = l;
    }

    if (iopt >= 0) {
        *nx0 = nx;
        *ny0 = ny;
    }

    return;
}


static void
fpsysy(double* restrict a, const int n, double* restrict g)
{
    g[0] /= a[0];
    if (n == 1) { return; }
    // decomposition of the symmetric matrix (a) = (l) * (d) *(l)'
    // with (l) a unit lower triangular matrix and (d) a diagonal
    // matrix
    for (int k = 1; k < n; k++) { a[k] /= a[0]; }
    for (int i = 2; i <= n; i++) {
        for (int k = i; k <= n; k++) {
            double fac = a[k - 1 + (i - 1)*6];
            for (int j = 1; j < i; j++) {
                fac -= a[j - 1 + (j - 1)*6] * a[k - 1 + (j - 1)*6] * a[i - 1 + (j - 1)*6];
            }
            a[k - 1 + (i - 1)*6] = fac;
            if (k > i) { a[k - 1 + (i - 1)*6] = fac / a[i - 1 + (i - 1)*6]; }
        }
    }
    // solve the system (l)*(d)*(l)'*(b) = (g).
    // first step : solve (l)*(d)*(c) = (g).
    for (int i = 2; i <= n; i++) {
        double fac = g[i - 1];
        for (int j = 1; j < i; j++) {
            fac -= g[j - 1] * a[j - 1 + (j - 1)*6] * a[i - 1 + (j - 1)*6];
        }
        g[i - 1] = fac / a[i - 1 + (i - 1)*6];
    }
    // second step : solve (l)'*(b) = (c)
    int i = n;
    for (int j = 2; j <= n; j++) {
        i--;
        double fac = g[i - 1];
        for (int k = i + 1; k <= n; k++) {
            fac -= g[k - 1] * a[k - 1 + (i - 1)*6];
        }
        g[i - 1] = fac;
    }
    return;
}


void
insert(const int iopt, const double* t, const int n, const double* c, const int k, const double x,
       double* tt, int* nn, double* cc, const int nest, int* ier)
{
    // Faithful translation of Fortran insert.f
    *ier = 10;
    if (nest <= n) { return; }
    int k1 = k + 1;
    int nk = n - k;
    if (x < t[k1 - 1] || x > t[nk - 1]) { return; }

    // Search for knot interval t[l] <= x < t[l+1]
    int l = k1;
    while (x >= t[l]) {
        l++;
        if (l == nk) {
            // if no interval found above, then reverse the search and
            // look for knot interval t[l] < x <= t[l+1]
            l = nk - 1;
            while (x <= t[l - 1]) {
                l--;
                if (l == k) { return; }
            }
            break;
        }
    }
    if (t[l - 1] >= t[l]) { return; }
    if (iopt != 0) {
        int kk = 2 * k;
        if ((l <= kk) && (l >= (n - kk))) {
            return;
        }
    }
    *ier = 0;
    // insert the new knot.
    fpinst(iopt, t, n, c, k, x, l, tt, nn, cc, nest);

    return;
}


void
parcur(const int iopt, const int ipar, const int idim, const int m, double *u, const int mx, const double *x,
       const double *w, double *ub, double *ue, const int k, const double s, const int nest, int *n, double *t,
       const int nc, double *c, double *fp, double *wrk, const int lwrk, int *iwrk, int *ier)
{
    double tol, dist;
    int i, ia, ib, ifp, ig, iq, iz, i1, i2, j, k1, k2, lwest, maxit, nmin, ncc;

    // we set up the parameters tol and maxit
    maxit = 20;
    tol = 0.1e-02;

    // before starting computations a data check is made. if the input data
    // are invalid, control is immediately repassed to the calling program.
    *ier = 10;
    if ((iopt < -1) || (iopt > 1)) { return; }
    if ((ipar < 0) || (ipar > 1)) { return; }
    if ((idim <= 0) || (idim > 10)) { return; }
    if ((k <= 0) || (k > 5)) { return; }
    k1 = k + 1;
    k2 = k1 + 1;
    nmin = 2 * k1;
    if ((m < k1) || (nest < nmin)) { return; }
    ncc = nest * idim;
    if ((mx < m * idim) || (nc < ncc)) { return; }
    lwest = m * k1 + nest * (6 + idim + 3 * k);
    if (lwrk < lwest) { return; }

    if ((ipar == 0) && (iopt <= 0)) {
        i1 = 0;
        i2 = idim;
        u[0] = 0.0;
        for (i = 2; i <= m; i++) {
            dist = 0.0;
            for (j = 1; j <= idim; j++) {
                i1++;
                i2++;
                dist += pow(x[i2 - 1] - x[i1 - 1], 2);
            }
            u[i - 1] = u[i - 2] + sqrt(dist);
        }
        if (u[m - 1] <= 0.0) { return; }
        for (i = 2; i <= m; i++) {
            u[i - 1] /= u[m - 1];
        }
        *ub = 0.0;
        *ue = 1.0;
        u[m - 1] = *ue;
    }

    if (((*ub) > u[0]) || ((*ue) < u[m - 1]) || (w[0] <= 0.0)) { return; }
    for (i = 1; i < m; i++) {
        if ((u[i - 1] >= u[i]) || (w[i] <= 0.0)) { return; }
    }

    if (iopt >= 0) {
        if (s < 0.0) { return; }
        if ((s == 0.0) && (nest < (m + k1))) { return; }
        *ier = 0;
    } else {
        if (((*n) < nmin) || ((*n) > nest)) { return; }
        j = *n;
        for (i = 1; i <= k1; i++) {
            t[i - 1] = *ub;
            t[j - 1] = *ue;
            j--;
        }
        fpchec(u, m, t, *n, k, ier);
        if (*ier != 0) { return; }
    }

    // we partition the working space and determine the spline curve.
    ifp = 0;
    iz = ifp + nest;
    ia = iz + ncc;
    ib = ia + nest * k1;
    ig = ib + nest * k2;
    iq = ig + nest * k2;
    fppara(iopt, idim, m, u, mx, x, w, *ub, *ue, k, s, nest, tol, maxit, k1, k2,
           n, t, ncc, c, fp, &wrk[ifp], &wrk[iz], &wrk[ia], &wrk[ib], &wrk[ig], &wrk[iq],
           iwrk, ier);
    return;
}


void
parder(const double *tx, int nx, const double *ty, int ny, double *c,
            int kx, int ky, int nux, int nuy, const double *x, int mx,
            const double *y, int my, double *z, double *wrk, int lwrk,
            int *iwrk, int kwrk, int *ier) {

    int kx1 = kx + 1;
    int ky1 = ky + 1;
    int nkx1 = nx - kx1;
    int nky1 = ny - ky1;
    int nc = nkx1 * nky1;
    int nxx, nyy, kkx, kky, lx, ly, l1, l2, m0, m1;
    double ak, fac;

    // before starting computations a data check is made. if the input data
    // are invalid control is immediately repassed to the calling program.
    *ier = 10;
    if ((nux < 0) || (nux >= kx)) { return; }
    if ((nuy < 0) || (nuy >= ky)) { return; }
    if (lwrk < (nc + (kx1 - nux) * mx + (ky1 - nuy) * my)) { return; }
    if (kwrk < (mx + my)) { return; }

    if (mx < 1) { return; }
    if (mx != 1) {
        for (int i = 1; i < mx; i++) {
            if (x[i] < x[i - 1]) { return; }
        }
    }

    if (my < 1) { return; }
    if (my != 1) {
        for (int i = 1; i < my; i++) {
            if (y[i] < y[i - 1]) { return; }
        }
    }

    *ier = 0;
    nxx = nkx1;
    nyy = nky1;
    kkx = kx;
    kky = ky;

    // the partial derivative of order (nux,nuy) of a bivariate spline of
    // degrees kx,ky is a bivariate spline of degrees kx-nux,ky-nuy.
    // we calculate the b-spline coefficients of this spline
    for (int i = 0; i < nc; i++) { wrk[i] = c[i]; }

    if (nux != 0) {
        lx = 1;
        for (int j = 1; j <= nux; j++) {
            ak = kkx;
            nxx--;
            l1 = lx;
            m0 = 1;
            for (int i = 1; i <= nxx; i++) {
                l1++;
                l2 = l1 + kkx;
                fac = tx[l2 - 1] - tx[l1 - 1];
                if (fac <= 0.0) { continue; }
                for (int m = 1; m <= nyy; m++) {
                    m1 = m0 + nyy;
                    wrk[m0 - 1] = (wrk[m1 - 1] - wrk[m0 - 1]) * ak / fac;
                    m0++;
                }
            }
            lx++;
            kkx--;
        }
    }

    if (nuy != 0) {
        ly = 1;
        for (int j = 1; j <= nuy; j++) {
            ak = kky;
            nyy--;
            l1 = ly;
            for (int i = 1; i <= nyy; i++) {
                l1++;
                l2 = l1 + kky;
                fac = ty[l2 - 1] - ty[l1 - 1];
                if (fac <= 0.0) { continue; }
                m0 = i;
                for (int m = 1; m <= nxx; m++) {
                    m1 = m0 + 1;
                    wrk[m0 - 1] = (wrk[m1 - 1] - wrk[m0 - 1]) * ak / fac;
                    m0 += nky1;
                }
            }
            ly++;
            kky--;
        }

        m0 = nyy;
        m1 = nky1;
        for (int m = 2; m <= nxx; m++) {
            for (int i = 1; i <= nyy; i++) {
                m0++;
                m1++;
                wrk[m0 - 1] = wrk[m1 - 1];
            }
            m1 += nuy;
        }
    }

    // Partition the workspace and evaluate the partial derivative
    int iwx = nxx * nyy;
    int iwy = iwx + mx * (kx1 - nux);
    fpbisp(&tx[nux], nx - 2*nux, &ty[nuy], ny - 2*nuy, wrk, kkx, kky, x, mx, y, my, z, &wrk[iwx], &wrk[iwy], iwrk, &iwrk[mx]);
}


void
pardeu(const double *tx, int nx, const double *ty, int ny, double *c,
       int kx, int ky, int nux, int nuy, const double *x, const double *y,
       double *z, int m, double *wrk, int lwrk, int *iwrk, int kwrk, int *ier)
{
    // Initialize variables
    int kx1 = kx + 1;
    int ky1 = ky + 1;
    int nkx1 = nx - kx1;
    int nky1 = ny - ky1;
    int nc = nkx1 * nky1;
    int nxx, nyy, kkx, kky, lx, ly, l1, l2, m0, m1;
    double ak, fac;

    // Validate input data
    *ier = 10;
    if ((nux < 0) || (nux >= kx)) { return; }
    if ((nuy < 0) || (nuy >= ky)) { return; }
    if (lwrk < (nc + (kx1 - nux) * m + (ky1 - nuy) * m)) { return; }
    if (kwrk < (m + m)) { return; }
    if (m < 1) { return; }

    *ier = 0;
    nxx = nkx1;
    nyy = nky1;
    kkx = kx;
    kky = ky;

    // Copy coefficients to workspace
    for (int i = 1; i <= nc; i++) {
        wrk[i - 1] = c[i - 1];
    }

    // Compute partial derivative with respect to x
    if (nux != 0) {
        lx = 1;
        for (int j = 1; j <= nux; j++) {
            ak = kkx;
            nxx--;
            l1 = lx;
            m0 = 1;
            for (int i = 1; i <= nxx; i++) {
                l1++;
                l2 = l1 + kkx;
                fac = tx[l2 - 1] - tx[l1 - 1];
                if (fac <= 0.0) { continue; }
                for (int mm = 1; mm <= nyy; mm++) {
                    m1 = m0 + nyy;
                    wrk[m0 - 1] = (wrk[m1 - 1] - wrk[m0 - 1]) * ak / fac;
                    m0++;
                }
            }
            lx++;
            kkx--;
        }
    }

    // Compute partial derivative with respect to y
    if (nuy != 0) {
        ly = 1;
        for (int j = 1; j <= nuy; j++) {
            ak = kky;
            nyy--;
            l1 = ly;
            for (int i = 1; i <= nyy; i++) {
                l1++;
                l2 = l1 + kky;
                fac = ty[l2 - 1] - ty[l1 - 1];
                if (fac <= 0.0) { continue; }
                m0 = i;
                for (int mm = 1; mm <= nxx; mm++) {
                    m1 = m0 + 1;
                    wrk[m0 - 1] = (wrk[m1 - 1] - wrk[m0 - 1]) * ak / fac;
                    m0 += nky1;
                }
            }
            ly++;
            kky--;
        }

        m0 = nyy;
        m1 = nky1;
        for (int mm = 2; mm <= nxx; mm++) {
            for (int i = 1; i <= nyy; i++) {
                m0++;
                m1++;
                wrk[m0 - 1] = wrk[m1 - 1];
            }
            m1 += nuy;
        }
    }

    // Partition the workspace and evaluate the partial derivative
    int iwx = nxx * nyy;
    int iwy = iwx + m * (kx1 - nux);
    for (int i = 1; i <= m; i++) {
        fpbisp(&tx[nux], nx - 2 * nux, &ty[nuy], ny - 2 * nuy, wrk, kkx, kky, &x[i - 1], 1, &y[i - 1], 1, &z[i - 1], &wrk[iwx], &wrk[iwy], iwrk, &iwrk[1]);
    }
}


void
pardtc(const double* tx, const int nx, const double* ty, const int ny, const double* c,
       const int kx, const int ky, const int nux, const int nuy, double* newc, int* ier)
{
    // before starting computations a data check is made. if the input data
    // are invalid control is immediately repassed to the calling program.
    *ier = 10;
    if ((nux < 0) || (nux >= kx)) { return; }
    if ((nuy < 0) || (nuy >= ky)) { return; }
    int kx1 = kx + 1;
    int ky1 = ky + 1;
    int nkx1 = nx - kx1;
    int nky1 = ny - ky1;
    int nc = nkx1 * nky1;
    *ier = 0;
    int nxx = nkx1;
    int nyy = nky1;
    int newkx = kx;
    int newky = ky;
    // the partial derivative of order (nux,nuy) of a bivariate spline of
    // degrees kx,ky is a bivariate spline of degrees kx-nux,ky-nuy.
    // we calculate the b-spline coefficients of this spline
    // that is to say newkx = kx - nux, newky = ky - nuy
    for (int i = 0; i < nc; i++) { newc[i] = c[i]; }

    if (nux != 0) {
        int lx = 1;
        for (int j = 1; j <= nux; j++) {
            double ak = newkx;
            nxx--;
            int l1 = lx;
            int m0 = 1;
            for (int i = 1; i <= nxx; i++) {
                l1++;
                int l2 = l1 + newkx;
                double fac = tx[l2 - 1] - tx[l1 - 1];
                if (fac <= 0.0) { continue; }
                for (int m = 1; m <= nyy; m++) {
                    int m1 = m0 + nyy;
                    newc[m0 - 1] = (newc[m1 - 1] - newc[m0 - 1]) * ak / fac;
                    m0++;
                }
            }
            lx++;
            newkx--;
        }
    }

    if (nuy != 0) {
        int ly = 1;
        for (int j = 1; j <= nuy; j++) {
            double ak = newky;
            nyy--;
            int l1 = ly;
            for (int i = 1; i <= nyy; i++) {
                l1++;
                int l2 = l1 + newky;
                double fac = ty[l2 - 1] - ty[l1 - 1];
                if (fac <= 0.0) { continue; }
                int m0 = i;
                for (int m = 1; m <= nxx; m++) {
                    int m1 = m0 + 1;
                    newc[m0 - 1] = (newc[m1 - 1] - newc[m0 - 1]) * ak / fac;
                    m0 += nky1;
                }
            }
            ly++;
            newky--;
        }

        int m0 = nyy;
        int m1 = nky1;
        for (int m = 2; m <= nxx; m++) {
            for (int i = 1; i <= nyy; i++) {
                m0++;
                m1++;
                newc[m0 - 1] = newc[m1 - 1];
            }
            m1 += nuy;
        }
    }

    return;
}


void
percur(const int iopt, const int m, const double *x, const double *y, const double *w, const int k,
       const double s, const int nest, int *n, double *t, double *c, double *fp,
       double *wrk, const int lwrk, int *iwrk, int *ier)
{
    int i, ia1, ia2, ib, ifp, ig1, ig2, iq, iz, i1, i2, j1, j2, k1, k2, lwest, maxit, m1, nmin;
    double per, tol;

    // we set up the parameters tol and maxit
    maxit = 20;
    tol = 1.0e-3;

    // before starting computations a data check is made. if the input data
    // are invalid, control is immediately repassed to the calling program.
    *ier = 10;
    if ((k <= 0) || (k > 5)) { return; }
    k1 = k + 1;
    k2 = k1 + 1;
    if ((iopt < -1) || (iopt > 1)) { return; }
    nmin = 2 * k1;
    if ((m < 2) || (nest < nmin)) { return; }
    lwest = m * k1 + nest * (8 + 5 * k);
    if (lwrk < lwest) { return; }
    m1 = m - 1;
    for (i = 0; i < m1; i++) {
        if ((x[i] >= x[i + 1]) || (w[i] <= 0.0)) { return; }
    }
    if (iopt >= 0) {
        if (s < 0.0) { return; }
        if ((s == 0.0) && (nest < (m + 2 * k))) { return; }
        *ier = 0;
    } else {
        if ((*n < nmin) || (*n > nest)) { return; }
        per = x[m - 1] - x[0];
        j1 = k1;
        t[j1 - 1] = x[0];
        i1 = *n - k;
        t[i1 - 1] = x[m - 1];
        j2 = j1;
        i2 = i1;
        for (i = 1; i <= k; i++) {
            i1++;
            i2--;
            j1++;
            j2--;
            t[j2 - 1] = t[i2 - 1] - per;
            t[i1 - 1] = t[j1 - 1] + per;
        }
        fpchep(x, m, t, *n, k, ier);
        if (*ier != 0) { return; }
    }

    // we partition the working space and determine the spline approximation.
    ifp = 0;
    iz = ifp + nest;
    ia1 = iz + nest;
    ia2 = ia1 + nest * k1;
    ib = ia2 + nest * k;
    ig1 = ib + nest * k2;
    ig2 = ig1 + nest * k2;
    iq = ig2 + nest * k1;

    fpperi(iopt, x, y, w, m, k, s, nest, tol, maxit, k1, k2, n, t, c, fp, &wrk[ifp], &wrk[iz],
           &wrk[ia1], &wrk[ia2], &wrk[ib], &wrk[ig1], &wrk[ig2], &wrk[iq], iwrk, ier);

    return;
}


void
regrid(const int iopt, const int mx, const double *x, const int my, const double *y, const double *z,
       const double xb, const double xe, const double yb, const double ye, const int kx, const int ky, const double s,
       const int nxest, const int nyest, const int maxit, int *nx, double *tx, int *ny, double *ty, double *c, double *fp,
       double *wrk, const int lwrk, int *iwrk, const int kwrk, int *ier)
{
    // regrid determines a smooth bivariate spline approximation for gridded data
    const double tol = 0.001;

    // data check
    *ier = 10;
    if ((kx <= 0) || (kx > 5)) { return; }
    int kx1 = kx + 1;
    int kx2 = kx1 + 1;
    if ((ky <= 0) || (ky > 5)) { return; }
    int ky1 = ky + 1;
    int ky2 = ky1 + 1;
    if ((iopt < -1) || (iopt > 1)) { return; }
    int nminx = 2 * kx1;
    if ((mx < kx1) || (nxest < nminx)) { return; }
    int nminy = 2 * ky1;
    if ((my < ky1) || (nyest < nminy)) { return; }
    int mz = mx * my;
    int nc = (nxest - kx1) * (nyest - ky1);
    int max_val = (nxest > my) ? nxest : my;
    int lwest = 4 + nxest * (my + 2 * kx2 + 1) + nyest * (2 * ky2 + 1) + mx * kx1 + my * ky1 + max_val;
    int kwest = 3 + mx + my + nxest + nyest;
    if ((lwrk < lwest) || (kwrk < kwest)) { return; }
    if ((xb > x[0]) || (xe < x[mx - 1])) { return; }
    for (int i = 1; i < mx; i++) {
        if (x[i - 1] >= x[i]) { return; }
    }
    if ((yb > y[0]) || (ye < y[my - 1])) { return; }
    for (int i = 1; i < my; i++) {
        if (y[i - 1] >= y[i]) { return; }
    }

    if (iopt < 0) {
        // check knots for least-squares case
        if ((*nx < nminx) || (*nx > nxest)) { return; }
        int j = *nx - 1;
        for (int i = 0; i < kx1; i++) {
            tx[i] = xb;
            tx[j] = xe;
            j--;
        }
        fpchec(x, mx, tx, *nx, kx, ier);
        if (*ier != 0) { return; }

        if ((*ny < nminy) || (*ny > nyest)) {
            *ier = 10;
            return;
        }
        j = *ny - 1;
        for (int i = 0; i < ky1; i++) {
            ty[i] = yb;
            ty[j] = ye;
            j--;
        }
        fpchec(y, my, ty, *ny, ky, ier);
        if (*ier != 0) {
            *ier = 10;
            return;
        }
    } else {
        // smoothing case
        if (s < 0.0) { return; }
        if ((s == 0.0) && ((nxest < (mx + kx1)) || (nyest < (my + ky1)))) { return; }
    }

    *ier = 0;

    // partition working space
    int lfpx = 4;
    int lfpy = lfpx + nxest;
    int lww = lfpy + nyest;
    int jwrk = lwrk - 4 - nxest - nyest;
    int knrx = 3;
    int knry = knrx + mx;
    int kndx = knry + my;
    int kndy = kndx + nxest;

    // call fpregr
    fpregr(iopt, x, mx, y, my, z, mz, xb, xe, yb, ye, kx, ky, s, nxest, nyest, tol, maxit, nc,
           nx, tx, ny, ty, c, fp, &wrk[0], &wrk[1], &wrk[2], &wrk[3], &wrk[lfpx], &wrk[lfpy],
           &iwrk[0], &iwrk[1], &iwrk[2], &iwrk[knrx], &iwrk[knry], &iwrk[kndx], &iwrk[kndy],
           &wrk[lww], jwrk, ier);
}


void
spalde(const double *t, const int n, const double *c, const int nc, const int k1, const double x,
       double* d, int* ier)
{
    (void)nc; // Unused.
    *ier = 10;
    int nk1 = n - k1;
    if ((x < t[k1 - 1]) || (x > t[nk1])) { return; }
    int l = k1;
    while (!((x < t[l]) || (l == nk1))) {
        l++;
    }
    if (t[l - 1] >= t[l]) { return; }
    *ier = 0;
    fpader(t, n, c, k1, x, l, d);

    return;
}


void
spgrid(const int *iopt, const int *ider, const int mu, const double *u, const int mv,
       const double *v, const double *r, const double r0, const double r1, const double s,
       const int nuest, const int nvest, int *nu, double *tu, int *nv, double *tv, double *c,
       double *fp, double *wrk, const int lwrk, int *iwrk, const int kwrk, int *ier)
{
    double per, tol, uu, ve, rmax, rmin, rn, rb = 0.0, re = 0.0;
    int i, i1, i2, j, jwrk, j1, j2, kndu, kndv, knru, knrv, kwest, l;
    int ldr, lfpu, lfpv, lwest, lww, m, maxit, mumin, muu, nc;

    per = PI + PI;
    ve = v[1 - 1] + per;
    maxit = 20;
    tol = 1e-03;
    *ier = 10;

    if ((iopt[0] < -1) || (iopt[0] > 1)) { return; }
    if ((iopt[1] < 0 ) || (iopt[1] > 1)) { return; }
    if ((iopt[2] < 0 ) || (iopt[2] > 1)) { return; }
    if ((ider[0] < -1) || (ider[0] > 1)) { return; }
    if ((ider[1] < 0 ) || (ider[1] > 1)) { return; }
    if ((ider[1] == 1) && (iopt[1] == 0)) { return; }
    if ((ider[2] < -1) || (ider[2] > 1)) { return; }
    if ((ider[3] < 0 ) || (ider[3] > 1)) { return; }
    if ((ider[3] == 1) && (iopt[2] == 0)) { return; }

    mumin = 4;
    if (ider[0] >= 0) { mumin = mumin - 1; }
    if ((iopt[1] == 1) && (ider[1] == 1)) { mumin = mumin - 1; }
    if (ider[2] >= 0) { mumin = mumin - 1; }
    if (iopt[2] == 1 && ider[3] == 1) { mumin = mumin - 1; }
    if (mumin == 0) { mumin = 1; }
    if ((mu < mumin) || (mv < 4)) { return; }
    if ((nuest < 8) || (nvest < 8)) { return; }

    m = mu * mv;
    nc = (nuest - 4) * (nvest - 4);

    lwest = 12 + nuest * (mv + nvest + 3) + 24 * nvest + 4 * mu + 8 * mv +
            ((nuest > (mv + nvest)) ? nuest : (mv + nvest));
    kwest = 5 + mu + mv + nuest + nvest;
    if (lwrk < lwest || kwrk < kwest) { return; }

    if (u[0] <= 0.0 || u[mu - 1] >= PI) { return; }
    if (mu != 1) {
        for (i = 1; i < mu; ++i) {
            if (u[i - 1] >= u[i]) { return; }
        }
    }
    if (v[0] < -PI || v[0] >= PI) { return; }
    if (v[mv - 1] >= v[0] + per) { return; }
    for (i = 1; i < mv; ++i) {
        if (v[i - 1] >= v[i]) { return; }
    }

    // Static conditions used in goto 140 statements
    int condition1 = s < 0.0;
    int condition2 = ((s == 0.0) && (nuest < (mu + 6 + iopt[1] + iopt[2]))) || (nvest < (mv + 7));

    if (iopt[0] > 0) {
        if (condition1 || condition2) { return; }
    } else {
        // if not given, we compute an estimate for r0.
        rn = (double)mv;
        if (ider[0] >= 0) {
            rb = r0;
        } else {
            rb = 0.0;
            for (i = 1; i <= mv; i++) {
                rb = rb + r[i - 1];
            }
            rb = rb / rn;
        }

        // if not given, we compute an estimate for r1.
        if (ider[2] >= 0) {
            re = r1;
        } else {
            re = 0.0;
            j = m;
            for (i = 1; i <= mv; i++) {
                re = re + r[j - 1];
                j--;
            }
            re = re / rn;
        }

        // we determine the range of r-values.
        rmin = rb;
        rmax = re;
        for (i = 1; i <= m; i++) {
            if (r[i - 1] < rmin) { rmin = r[i - 1]; }
            if (r[i - 1] > rmax) { rmax = r[i - 1]; }
        }
        wrk[4] = rb;
        wrk[5] = 0.0;
        wrk[6] = 0.0;
        wrk[7] = re;
        wrk[8] = 0.0;
        wrk[9] = 0.0;
        wrk[10] = rmax - rmin;
        wrk[11] = wrk[10];
        iwrk[3] = mu;
        iwrk[4] = mu;
        if (iopt[0] == 0) {
            if (condition1 || condition2) { return; }
        } else {
            if ((*nu < 8 ) || (*nu > nuest)) { return; }
            if ((*nv < 11) || (*nv > nvest)) { return; }
            j = *nu;
            for (i = 1; i <= 4; i++) {
                tu[i - 1] = 0.0;
                tu[j - 1] = PI;
                j--;
            }
            l = 13;
            wrk[l - 1] = 0.0;
            if (iopt[1] != 0) {
                l++;
                uu = u[0];
                if (uu > tu[4]) { uu = tu[4]; }
                wrk[l - 1] = uu * 0.5;
            }
            for (i = 1; i <= mu; i++) {
                l++;
                wrk[l - 1] = u[i - 1];
            }
            if (iopt[2] != 0) {
                l++;
                uu = u[mu - 1];
                if (uu < tu[*nu - 5]) { uu = tu[*nu - 5]; }
                wrk[l - 1] = uu + (PI - uu) * 0.5;
            }
            l++;
            wrk[l - 1] = PI;
            muu = l - 12;
            fpchec(&wrk[13], muu, tu, *nu, 3, ier);
            if (*ier != 0) { return; }
            j1 = 4;
            tv[j1 - 1] = v[0];
            i1 = *nv - 3;
            tv[i1 - 1] = ve;
            j2 = j1;
            i2 = i1;
            for (i = 1; i <= 3; i++) {
                i1++;
                i2--;
                j1++;
                j2--;
                tv[j2 - 1] = tv[i2 - 1] - per;
                tv[i1 - 1] = tv[j1 - 1] + per;
            }
            l = 13;
            for (i = 1; i <= mv; i++) {
                wrk[l - 1] = v[i - 1];
                l++;
            }
            wrk[l - 1] = ve;
            fpchep(&wrk[12], mv + 1, tv, *nv, 3, ier);
            if (*ier != 0) { return; }
        }
    }

    // we partition the working space and determine the spline approximation
    ldr = 5;
    lfpu = 13;
    lfpv = lfpu + nuest;
    lww = lfpv + nvest;
    jwrk = lwrk - 12 - nuest - nvest;
    knru = 6;
    knrv = knru + mu;
    kndu = knrv + mv;
    kndv = kndu + nuest;

    fpspgr(iopt, ider, u, mu, v, mv, r, m, rb, re, s, nuest, nvest, tol, maxit,
           nc, nu, tu, nv, tv, c, fp, &wrk[0], &wrk[1], &wrk[2], &wrk[3],
           &wrk[lfpu - 1], &wrk[lfpv - 1], &wrk[ldr - 1], &wrk[10],
           &iwrk[0], &iwrk[1], &iwrk[2], &iwrk[3], &iwrk[4], &iwrk[knru - 1],
           &iwrk[knrv - 1], &iwrk[kndu - 1], &iwrk[kndv - 1], &wrk[lww - 1], jwrk, ier);
}


void
sphere(const int iopt, const int m, const double *teta, const double *phi, const double *r,
       const double *w, const double s, const int ntest, const int npest, const double eps,
       int *nt, double *tt, int *np, double *tp, double *c, double *fp,
       double *wrk1, const int lwrk1, double *wrk2, const int lwrk2, int *iwrk, const int kwrk,
       int *ier)
{
    double tol, pi2;
    int i, ib1, ib3, ki, kn, kwest, la, lbt, lcc, lcs, lro, j;
    int lbp, lco, lf, lff, lfp, lh, lq, lst, lsp, lwest, maxit, ncest, ncc, ntt;
    int npp, nreg, nrint, ncof, nt4, np4;


    // we set up the parameters tol and maxit.
    maxit = 20;
    tol = 0.001;
    // before starting computations a data check is made. if the input data
    // are invalid,control is immediately repassed to the calling program.
    *ier = 10;
    if ((eps <= 0.0) || (eps >= 1.0)) { return; }
    if ((iopt < -1) || (iopt > 1)) { return; }
    if (m < 2) { return; }
    if ((ntest < 8) || (npest < 8)) { return; }
    nt4 = ntest - 4;
    np4 = npest - 4;
    ncest = nt4 * np4;
    ntt = ntest - 7;
    npp = npest - 7;
    ncc = 6 + npp * (ntt - 1);
    nrint = ntt + npp;
    nreg = ntt * npp;
    ncof = 6 + 3 * npp;
    ib1 = 4 * npp;
    ib3 = ib1 + 3;
    if (ncof > ib1) { ib1 = ncof; }
    if (ncof > ib3) { ib3 = ncof; }
    lwest = 185 + 52 * npp + 10 * ntt + 14 * ntt * npp + 8 * (m + (ntt - 1) * npp * npp);
    kwest = m + nreg;
    if ((lwrk1 < lwest) || (kwrk < kwest)) { return; }
    if (iopt <= 0) {
        pi2 = PI + PI;
        for (i = 1; i <= m; i++) {
            if (w[i - 1] <= 0.0) { return; }
            if ((teta[i - 1] < 0.0) || (teta[i - 1] > PI)) { return; }
            if ((phi[i - 1] < 0.0) || (phi[i - 1] > pi2)) { return; }
        }
        if (iopt == 0) {
            if (s < 0.0) { return; }
        } else {
            ntt = *nt - 8;
            if ((ntt < 0) || (*nt > ntest)) { return; }
            if (ntt != 0) {
                tt[3] = 0.0;
                for (i = 1; i <= ntt; i++) {
                    j = i + 4;
                    if ((tt[j - 1] <= tt[j - 2]) || (tt[j - 1] >= PI)) { return; }
                }
            }
            npp = *np - 8;
            if ((npp < 1) || (*np > npest)) { return; }
            tp[3] = 0.0;
            for (i = 1; i <= npp; i++) {
                j = i + 4;
                if ((tp[j - 1] <= tp[j - 2]) || (tp[j - 1] >= pi2)) { return; }
            }
        }
    } else {
        if (s < 0.0) { return; }
    }
    *ier = 0;
    // we partition the working space and determine the spline approximation
    kn = 1;
    ki = kn + m;
    lq = 2;
    la = lq + ncc * ib3;
    lf = la + ncc * ib1;
    lff = lf + ncc;
    lfp = lff + ncest;
    lco = lfp + nrint;
    lh = lco + nrint;
    lbt = lh + ib3;
    lbp = lbt + 5 * ntest;
    lro = lbp + 5 * npest;
    lcc = lro + npest;
    lcs = lcc + npest;
    lst = lcs + npest;
    lsp = lst + m * 4;
    fpsphe(iopt, m, teta, phi, r, w, s, ntest, npest, eps, tol, maxit,
           ib1, ib3, ncest, ncc, nrint, nreg, nt, tt, np, tp, c, fp, &wrk1[0], &wrk1[lfp - 1],
           &wrk1[lco - 1], &wrk1[lf - 1], &wrk1[lff - 1], &wrk1[lro - 1], &wrk1[lcc - 1], &wrk1[lcs - 1],
           &wrk1[la - 1], &wrk1[lq - 1], &wrk1[lbt - 1], &wrk1[lbp - 1], &wrk1[lst - 1], &wrk1[lsp - 1],
           &wrk1[lh - 1], &iwrk[ki - 1], &iwrk[kn - 1], wrk2, lwrk2, ier);
}


void
splder(const double *t, const int n, const double *c, const int nc, const int k, const int nu,
       const double *x, double *y, const int m, const int e, double *wrk, int *ier)
{
    (void)nc; // nc is not used
    int l, kk, i, j, l1, l2, nk2, k2, ll;
    double ak, arg, fac, sp;
    double h[6];
    // before starting computations a data check is made. if the input data
    // are invalid control is immediately repassed to the calling program.
    *ier = 10;
    if ((nu < 0) || (nu > k)) { return; }
    if (m < 1) { return; }

    *ier = 0;
    // fetch tb and te, the boundaries of the approximation interval.
    int k1 = k + 1;
    int k3 = k1 + 1;
    int nk1 = n - k1;
    double tb = t[k1 - 1];
    double te = t[nk1];

    // the derivative of order nu of a spline of degree k is a spline of
    // degree k-nu,the b-spline coefficients wrk(i) of which can be found
    // using the recurrence scheme of de boor.
    l = 1;
    kk = k;

    for (i = 1; i <= nk1; i++) {
        wrk[i - 1] = c[i - 1];
    }
    if (nu != 0) {
        nk2 = nk1;
        for (j = 1; j <= nu; j++) {
            ak = (double)kk;
            nk2 = nk2 - 1;
            l1 = l;
            for (i = 1; i <= nk2; i++) {
                l1++;
                l2 = l1 + kk;
                fac = t[l2 - 1] - t[l1 - 1];
                if (fac <= 0.0) { continue; }
                wrk[i - 1] = ak * (wrk[i] - wrk[i - 1]) / fac;
            }
            l++;
            kk--;
        }
        if (kk == 0) {
            // if nu=k the derivative is a piecewise constant function
            j = 1;
            for (i = 1; i <= m; i++) {
                arg = x[i - 1];
                // check if arg is in the support
                if ((arg < tb) || (arg > te)) {
                    if (e == 0) {
                        ; // go to 65
                    } else if (e == 1) {
                        y[i - 1] = 0.0;
                        continue;
                    } else if (e == 2) {
                        *ier = 1;
                        return;
                    }
                }
                // search for knot interval t(l) <= arg < t(l+1)
                while (!((arg >= t[l - 1]) || (l + 1 == k3))) {
                    l1 = l;
                    l--;
                    j--;
                }
                while (!((arg < t[l]) || (l == nk1))) {
                    l++;
                    j++;
                }
                y[i - 1] = wrk[j - 1];
            }
            return;
        }
    }

    l = k1;
    l1 = l + 1;
    k2 = k1 - nu;
    // main loop for the different points.
    for (i = 1; i <= m; i++) {
        // fetch a new x-value arg.
        arg = x[i - 1];
        // check if arg is in the support
        if ((arg < tb) || (arg > te)) {
            if (e == 0) {
                ; // go to 135
            } else if (e == 1) {
                y[i - 1] = 0.0;
                continue;
            } else if (e == 2) {
                *ier = 1;
                return;
            }
        }
        // search for knot interval t(l) <= arg < t(l+1)
        while (!(arg >= t[l - 1] || l1 == k3)) {
            l1 = l;
            l--;
        }
        while (!(arg < t[l1 - 1] || l == nk1)) {
            l = l1;
            l1 = l + 1;
        }
        // evaluate the non-zero b-splines of degree k-nu at arg.
        for (int hi = 0; hi < 6; hi++) { h[hi] = 0.0; }
        fpbspl(t, n, kk, arg, l, h);
        // find the value of the derivative at x=arg.
        sp = 0.0;
        ll = l - k1;
        for (j = 1; j <= k2; j++) {
            ll++;
            sp = sp + wrk[ll - 1] * h[j - 1];
        }
        y[i - 1] = sp;
    }
}


void
splev(const double *t, const int n, const double *c, const int nc, const int k,
      const double *x, double *y, const int m, const int e, int *ier)
{
    (void)nc; // nc is not used
    int k1, l, ll, l1, nk1, k2;
    double arg, sp, tb, te;
    double h[20];
    // before starting computations a data check is made. if the input data
    // are invalid control is immediately repassed to the calling program.
    *ier = 10;
    if (m < 1) { return; }

    *ier = 0;
    // fetch tb and te, the boundaries of the approximation interval.
    k1 = k + 1;
    k2 = k1 + 1;
    nk1 = n - k1;
    tb = t[k];
    te = t[nk1];
    l = k1;
    l1 = l + 1;

    // main loop for the different points.
    for (int i = 1; i <= m; i++) {
        // fetch a new x-value arg.
        arg = x[i - 1];
        // check if arg is in the support
        if ((arg < tb) || (arg > te)) {
            if (e == 0) {
                ; // go to 135
            } else if (e == 1) {
                y[i - 1] = 0.0;
                continue;
            } else if (e == 2) {
                *ier = 1;
                return;
            } else if (e == 3) {
                if (arg < tb) {
                    arg = tb;
                } else {
                    arg = te;
                }
            }
        }
        // search for knot interval t(l) <= arg < t(l+1)
        while (!((arg >= t[l - 1]) || (l1 == k2))) {
            l1 = l;
            l--;
        }
        while (!((arg < t[l1 - 1]) || (l == nk1))) {
            l = l1;
            l1 = l + 1;
        }
        // evaluate the non-zero b-splines at arg.
        fpbspl(t, n, k, arg, l, h);
        // find the value of s(x) at x=arg.
        sp = 0.0;
        ll = l - k1;
        for (int j = 1; j <= k1; j++) {
            ll++;
            sp = sp + c[ll - 1] * h[j - 1];
        }
        y[i - 1] = sp;
    }
}

double
splint(const double *t, const int n, const double *c, const int nc, const int k, const double a, const double b, double* wrk)
{
    (void)nc;  // Unused
    // calculate the integrals wrk(i) of the normalized b-splines
    // ni,k+1(x), i=1,2,...,n-k-1
    fpintb(t, n, wrk, n-k-1, a, b);
    // calculate the integral of s(x).
    double res = 0.0;
    for (int i = 0; i < n - k - 1; i++) {
        res += c[i] * wrk[i];
    }
    return res;
}


void
sproot(const double *t, const int n, const double *c, const int nc, double *zero, const int mest, int *m, int *ier)
{
    (void)nc; // Unused
    // subroutine sproot finds the zeros of a cubic spline s(x),which is
    // given in its normalized b-spline representation.
    int i, j, l;
    double ah, a0, a1, a2, a3, bh, b0, b1, c1, c2, c3, c4, c5, d4, d5, h1, h2;
    double t1, t2, t3, t4, t5, zz;
    int z0, z1, z2, z3, z4, nz0, nz1, nz2, nz3, nz4;
    double y[3];

    // before starting computations a data check is made. if the input data
    // are invalid, control is immediately repassed to the calling program.
    *ier = 10;
    if (n < 8) { return; }
    j = n;
    for (i = 0; i < 3; i++) {
        if (t[i] > t[i + 1]) { return; }
        if (t[j - 1] < t[j - 2]) { return; }
        j--;
    }
    for (i = 3; i < n - 4; i++) {
        if (t[i] >= t[i + 1]) { return; }
    }

    // the problem considered reduces to finding the zeros of the cubic
    // polynomials pl(x) which define the cubic spline in each knot
    // interval t(l)<=x<=t(l+1). a zero of pl(x) is also a zero of s(x) on
    // the condition that it belongs to the knot interval.
    // the cubic polynomial pl(x) is determined by computing s(t(l)),
    // s'(t(l)),s(t(l+1)) and s'(t(l+1)). in fact we only have to compute
    // s(t(l+1)) and s'(t(l+1)); because of the continuity conditions of
    // splines and their derivatives, the value of s(t(l)) and s'(t(l))
    // is already known from the foregoing knot interval.
    *ier = 0;

    // evaluate some constants for the first knot interval
    h1 = t[3] - t[2];
    h2 = t[4] - t[3];
    t1 = t[3] - t[1];
    t2 = t[4] - t[2];
    t3 = t[5] - t[3];
    t4 = t[4] - t[1];
    t5 = t[5] - t[2];

    // calculate a0 = s(t(4)) and ah = s'(t(4)).
    c1 = c[0];
    c2 = c[1];
    c3 = c[2];
    c4 = (c2 - c1) / t4;
    c5 = (c3 - c2) / t5;
    d4 = (h2 * c1 + t1 * c2) / t4;
    d5 = (t3 * c2 + h1 * c3) / t5;
    a0 = (h2 * d4 + h1 * d5) / t2;
    ah = 3.0 * (h2 * c4 + h1 * c5) / t2;
    z1 = (ah >= 0.0);
    nz1 = !z1;
    *m = 0;

    // main loop for the different knot intervals.
    for (l = 3; l < n - 4; l++) {
        // evaluate some constants for the knot interval t(l) <= x <= t(l+1).
        h1 = h2;
        h2 = t[l + 2] - t[l + 1];
        t1 = t2;
        t2 = t3;
        t3 = t[l + 3] - t[l + 1];
        t4 = t5;
        t5 = t[l + 3] - t[l];

        // find a0 = s(t(l)), ah = s'(t(l)), b0 = s(t(l+1)) and bh = s'(t(l+1)).
        c1 = c2;
        c2 = c3;
        c3 = c[l];
        c4 = c5;
        c5 = (c3 - c2) / t5;
        d4 = (h2 * c1 + t1 * c2) / t4;
        d5 = (h1 * c3 + t3 * c2) / t5;
        b0 = (h2 * d4 + h1 * d5) / t2;
        bh = 3.0 * (h2 * c4 + h1 * c5) / t2;

        // calculate the coefficients a0,a1,a2 and a3 of the cubic polynomial
        // pl(x) = ql(y) = a0+a1*y+a2*y**2+a3*y**3 ; y = (x-t(l))/(t(l+1)-t(l)).
        a1 = ah * h1;
        b1 = bh * h1;
        a2 = 3.0 * (b0 - a0) - b1 - 2.0 * a1;
        a3 = 2.0 * (a0 - b0) + b1 + a1;

        // test whether or not pl(x) could have a zero in the range
        // t(l) <= x <= t(l+1).
        z3 = (b1 >= 0.0);
        nz3 = !z3;
        if (a0 * b0 <= 0.0) {
            // find the zeros of ql(y).
            fpcuro(a3, a2, a1, a0, y, &j);
            if (j != 0) {
                // find which zeros of pl(x) are zeros of s(x).
                for (i = 0; i < j; i++) {
                    if ((y[i] < 0.0) || (y[i] > 1.0)) { continue; }
                    // test whether the number of zeros of s(x) exceeds mest.
                    if (*m >= mest) { *ier = 1; return; }
                    (*m)++;
                    zero[*m - 1] = t[l] + h1 * y[i];
                }
            }
        } else {
            z0 = (a0 >= 0.0);
            nz0 = !z0;
            z2 = (a2 >= 0.0);
            nz2 = !z2;
            z4 = ((3.0 * a3 + a2) >= 0.0);
            nz4 = !z4;
            if (!(((z0) && ((nz1 && (z3 || (z2 && nz4))) || (nz2 && z3 && z4))) ||
                  ((nz0) && ((z1 && (nz3 || (nz2 && z4))) || (z2 && nz3 && nz4))))) {
                a0 = b0;
                ah = bh;
                z1 = z3;
                nz1 = nz3;
                continue;
            }
            // find the zeros of ql(y).
            fpcuro(a3, a2, a1, a0, y, &j);
            if (j != 0) {
                // find which zeros of pl(x) are zeros of s(x).
                for (i = 0; i < j; i++) {
                    if ((y[i] < 0.0) || (y[i] > 1.0)) { continue; }
                    if (*m >= mest) { *ier = 1; return; }
                    (*m)++;
                    zero[*m - 1] = t[l] + h1 * y[i];
                }
            }
        }
        a0 = b0;
        ah = bh;
        z1 = z3;
        nz1 = nz3;
    }

    // the zeros of s(x) are arranged in increasing order.
    if (*m < 2) { return; }
    for (i = 1; i < *m; i++) {
        j = i;
        while (j > 0) {
            if (zero[j] >= zero[j - 1]) { break; }
            zz = zero[j];
            zero[j] = zero[j - 1];
            zero[j - 1] = zz;
            j--;
        }
    }
    j = *m;
    *m = 1;
    for (i = 1; i < j; i++) {
        if (zero[i] == zero[*m - 1]) { continue; }
        (*m)++;
        zero[*m - 1] = zero[i];
    }
    return;
}


void
surfit(int iopt, int m, double* x, double* y, double* z, double* w,
       double xb, double xe, double yb, double ye, int kx, int ky, double s,
       int nxest, int nyest, int nmax, double eps, int* nx, double* tx,
       int* ny, double* ty, double* c, double* fp, double* wrk1, int lwrk1,
       double* wrk2, int lwrk2, int* iwrk, int kwrk, int* ier)
{
    // ..local scalars..
    double tol;
    int i, ib1, ib3, jb1, ki, kmax, km1, km2, kn, kwest, kx1, ky1, la, lbx;
    int lby, lco, lf, lff, lfp, lh, lq, lsx, lsy, lwest, maxit, ncest, nest, nek;
    int nminx, nminy, nmx, nmy, nreg, nrint, nxk, nyk;

    // we set up the parameters tol and maxit.
    maxit = 20;
    tol = 0.001;

    // before starting computations a data check is made. if the input data
    // are invalid,control is immediately repassed to the calling program.
    *ier = 10;

    if ((eps <= 0.0) || (eps >= 1.0)) { return; }
    if ((kx <= 0) || (kx > 5)) { return; }
    kx1 = kx + 1;
    if ((ky <= 0) || (ky > 5)) { return; }
    ky1 = ky + 1;
    kmax = (kx > ky) ? kx : ky;
    km1 = kmax + 1;
    km2 = km1 + 1;
    if ((iopt < -1) || (iopt > 1)) { return; }
    if (m < (kx1 * ky1)) { return; }
    nminx = 2 * kx1;
    if ((nxest < nminx) || (nxest > nmax)) { return; }
    nminy = 2 * ky1;
    if ((nyest < nminy) || (nyest > nmax)) { return; }
    nest = (nxest > nyest) ? nxest : nyest;
    nxk = nxest - kx1;
    nyk = nyest - ky1;
    ncest = nxk * nyk;
    nmx = nxest - nminx + 1;
    nmy = nyest - nminy + 1;
    nrint = nmx + nmy;
    nreg = nmx * nmy;
    ib1 = kx * nyk + ky1;
    jb1 = ky * nxk + kx1;
    ib3 = kx1 * nyk + 1;
    if (ib1 > jb1) {
        ib1 = jb1;
        ib3 = ky1 * nxk + 1;
    }
    lwest = ncest * (2 + ib1 + ib3) + 2 * (nrint + nest * km2 + m * km1) + ib3;
    kwest = m + nreg;
    if ((lwrk1 < lwest) || (kwrk < kwest)) { return; }
    if ((xb >= xe) || (yb >= ye)) { return; }
    for (i = 0; i < m; i++) {
        if (w[i] <= 0.0) { return; }
        if ((x[i] < xb) || (x[i] > xe)) { return; }
        if ((y[i] < yb) || (y[i] > ye)) { return; }
    }

    if (iopt < 0) {
        if ((*nx < nminx) || (*nx > nxest)) { return; }
        nxk = *nx - kx1;
        tx[kx1 - 1] = xb;
        tx[nxk] = xe;
        for (i = kx1; i <= nxk; i++) {
            if (tx[i] <= tx[i - 1]) { return; }
        }
        if ((*ny < nminy) || (*ny > nyest)) { return; }
        nyk = *ny - ky1;
        ty[ky1 - 1] = yb;
        ty[nyk] = ye;
        for (i = ky1; i <= nyk; i++) {
            if (ty[i] <= ty[i - 1]) { return; }
        }
    } else {
        if (s < 0.0) { return; }
    }

    *ier = 0;

    // we partition the working space and determine the spline approximation
    kn = 0;
    ki = kn + m;
    lq = 1;
    la = lq + ncest * ib3;
    lf = la + ncest * ib1;
    lff = lf + ncest;
    lfp = lff + ncest;
    lco = lfp + nrint;
    lh = lco + nrint;
    lbx = lh + ib3;
    nek = nest * km2;
    lby = lbx + nek;
    lsx = lby + nek;
    lsy = lsx + m * km1;

    fpsurf(iopt, m, x, y, z, w, xb, xe, yb, ye, kx, ky, s, nxest, nyest,
           eps, tol, maxit, nest, km1, km2, ib1, ib3, ncest, nrint, nreg,
           nx, tx, ny, ty, c, fp, &wrk1[0], &wrk1[lfp], &wrk1[lco],
           &wrk1[lf], &wrk1[lff], &wrk1[la], &wrk1[lq],
           &wrk1[lbx], &wrk1[lby], &wrk1[lsx], &wrk1[lsy],
           &wrk1[lh], &iwrk[ki], &iwrk[kn], wrk2, lwrk2, ier);
}
