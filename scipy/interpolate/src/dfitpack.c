#include <math.h>


void
fprota(double c, double s, double *a, double *b)
{
    // fprota applies a givens rotation to the pair (a,b).
    double temp = c * (*a) - s * (*b);
    *b          = c * (*b) + s * (*a);
    *a          = temp;
}


void
fpbspl(const double* t, const int n, const int k, const double x, const int l, double* h)
{
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
    for (int j = 0; j < k; j++) {
        for (int i = 0; i <= j; i++)
        {
            hh[i] = h[i];
        }

        h[0] = 0.0;
        for (int i = 0; i <= j; i++)
        {
            li = l + i;
            lj  = li - j;
            if (t[li-1] != t[lj-1])
            {
                f = hh[i] / (t[li] - t[lj]);
                h[i] = h[i] + f * (t[li] - x);
                h[i + 1] = f * (x - t[lj]);
            }
            else {
                h[i + 1] = 0.0;
            }
        }
    }
}


void
fpback(double* a, double* z, const int n, const int k, double* c, const int nest)
{
    c[n - 1] = z[n - 1] / a[(n - 1) + 0 * nest];
    if (n == 1) { return; }

    int i = n - 2;
    for (int j = 1; j < n; j++) {
        double store = z[i];
        int i1 = k - 1;
        if (j <= k - 1) {
            i1 = j;
        }
        int m = i;
        for (int l = 0; l < i1; l++) {
            m++;
            store = store - c[m] * a[i + (l + 1) * nest];
        }
        c[i] = store / a[i + 0 * nest];
        i--;
    }
}


void
fpbisp(const double *tx, const int nx, double *ty, const int ny, const double *c,
    const int kx, const int ky, const double *x, const int mx, const double *y, const int my,
    double *z, double *wx, double *wy, int *lx, int *ly)
{
    int kx1, ky1, l, l1, l2, nkx1, nky1, i1, j1, m;
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
            wx[(i - 1) + (j - 1)*kx1] = h[j - 1];
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
            wy[(i - 1) + (j - 1)*ky1] = h[j - 1];
        }
    }

    m = 0;
    for (int i = 1; i <= mx; i++) {

        l = lx[i - 1] * nky1;
        for (int i1 = 1; i1 <= kx1; i1++) {
            h[i1 - 1] = wx[(i - 1) + (i1 - 1)*kx1];
        }
        for (int j = 1; j <= my; j++) {
            l1 = l + ly[j - 1];
            sp = 0.0;
            for (i1 = 1; i1 <= kx1; i1++) {
                l2 = l1;
                for (int j1 = 1; j1 <= ky1; j1++) {
                    l2++;
                    sp = sp + c[l2 - 1] * h[i1 - 1] * wy[(j - 1) * ky1 + (j1 - 1)];
                }
                l1 += nky1;
            }
            z[m] = sp;
            m++;
        }
    }
}


void
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


void
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


void
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


void
fprank(const double* a, const double* f, const int n, const int m, const int na,
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
            store += a[(ij - 1) + j * na] * ff[kk - 1];
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


double
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
        p = -((*p1)*(*p2)*h3 + (*p2)*(*p3)*h1 + (*p3)*(*p1)*h2) / ((*p1)*h3 + (*p2)*h1 + (*p3)*h2);
    } else {
        // value of p in case p3 == inf.
        p = ((*p1) * (*f2) - (*p2) * (*f1)) / ((*f2) - (*f1));
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


void
fpsurf(int iopt, int m, double* x, double* y, double* z, double* w,
       double xb, double xe, double yb, double ye, int kxx, int kyy,
       double s, int nxest, int nyest, double eta, double tol, int maxit,
       int nmax, int km1, int km2, int ib1, int ib3, int nc, int intest,
       int nrest, int nx0, double* tx, int ny0, double* ty, double* c,
       double* fp, double* fp0, double* fpint, double* coord, double* f,
       double* ff, double* a, double* q, double* bx, double* by, double* spx,
       double* spy, double* h, int* index, int* nummer, double* wrk,
       int lwrk, int* ier)
{
    // ..local scalars..
    double acc, arg, cos, dmax, fac1, fac2, fpmax, fpms, f1, f2, f3, hxi, p, pinv;
    double piv, p1, p2, p3, rn, sigma, sin, sq, store, wi, x0, x1, y0, y1, zi, eps;
    int i, iband, iband1, iband3, iband4, ibb, ichang, ich1, ich3, ii;
    int in, irot, iter, i1, i2, i3, j, jrot, jxy, j1, kx, kx1, kx2, ky, ky1, ky2, l;
    int la, lf, lh, lwest, lx, ly, l1, l2, n, ncof, nk1x, nk1y, nminx, nminy, nreg;
    int nrint, num, num1, nx, nxe, nxx, ny, nye, nyy, n1, rank;
    // ..local arrays..
    double hx[6], hy[6];

    double one = 1.0, con1 = 0.1, con9 = 0.9, con4 = 0.04, half = 0.5, ten = 10.0;

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
        nx = nx0;
        ny = ny0;
    } else {
        // calculation of acc, the absolute tolerance for the root of f(p)=s.
        acc = tol * s;
        if (iopt == 0) {
            // initialization for the least-squares polynomial.
            nminx = 2 * kx1;
            nminy = 2 * ky1;
            nx = nminx;
            ny = nminy;
            *ier = -2;
        } else {
            if (*fp0 > s) {
                nx = nx0;
                ny = ny0;
            } else {
                nminx = 2 * kx1;
                nminy = 2 * ky1;
                nx = nminx;
                ny = nminy;
                *ier = -2;
            }
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
                    fpgivs(piv, &a[(irot - 1) + 0 * nc], &cos, &sin);

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
            if (a[(i - 1)] > dmax) {
                dmax = a[(i - 1)];
            }
        }

        // check whether the observation matrix is rank deficient.
        sigma = eps * dmax;
        int rank_deficient = 0;
        for (i = 1; i <= ncof; i++) {
            if (a[(i - 1)] <= sigma) { rank_deficient = 1; }
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
                q[(i - 1)] = a[(i - 1)] / dmax;
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
                store = (w[in - 1] * (z[in - 1] - store)) * (w[in - 1] * (z[in - 1] - store));
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
                if (fac1 > (ten * fac2) || fac2 > (ten * fac1)) { continue; }
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
                if (fac1 > (ten * fac2) || fac2 > (ten * fac1)) { continue; }
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
    p = (double)ncof / p;

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
                    for (l = 1; l <= iband; l++) {
                        h[l - 1] = 0.0;
                    }
                    // fill in the non-zero elements of the row. jrot records the column
                    // number of the first non-zero element in the row.
                    for (l = 1; l <= ky2; l++) {
                        h[l - 1] = by[(ii - 1) + (l - 1) * nmax] * pinv;
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
                            fpgivs(piv, &q[(irot - 1)], &cos, &sin);
                            // apply that givens transformation to the right hand side.
                            fprota(cos, sin, &zi, &ff[irot - 1]);
                            if (i2 == 0) {
                                break;
                            }
                            // apply that givens transformation to the left hand side.
                            for (l = 1; l <= i2; l++) {
                                l1 = l + 1;
                                fprota(cos, sin, &h[l1 - 1], &q[(irot - 1) + l1 * nc]);
                            }
                        }
                        if (i2 <= 0) {
                            break;
                        }
                        for (l = 1; l <= i2; l++) {
                            h[l - 1] = h[l];
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
                    for (l = 1; l <= iband4; l++) {
                        h[l - 1] = 0.0;
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
                            fpgivs(piv, &q[(irot - 1)], &cos, &sin);
                            // apply that givens transformation to the right hand side.
                            fprota(cos, sin, &zi, &ff[irot - 1]);
                            if (i2 == 0) {
                                break;
                            }
                            // apply that givens transformation to the left hand side.
                            for (l = 1; l <= i2; l++) {
                                l1 = l + 1;
                                fprota(cos, sin, &h[l1 - 1], &q[(irot - 1) + l1 * nc]);
                            }
                        }
                        if (i2 <= 0) {
                            break;
                        }
                        for (l = 1; l <= i2; l++) {
                            h[l - 1] = h[l];
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
            if (q[(i - 1)] > dmax) {
                dmax = q[(i - 1)];
            }
        }

        int rank_deficient = 0;
        // check whether the matrix is rank deficient.
        sigma = eps * dmax;
        for (i = 1; i <= ncof; i++) {
            if (q[(i - 1)] <= sigma) {
                rank_deficient = 1;
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

        for (i = 1; i <= ncof; i++) {
            q[(i - 1)] = q[(i - 1)] / dmax;
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
                *fp = *fp + (w[in - 1] * (z[in - 1] - store)) * (w[in - 1] * (z[in - 1] - store));
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
        nx0 = nx;
        ny0 = ny;
    }

    return;
}


void
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
            i++;
            j++;
            f[i - 1] = c[j - 1];
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
    for (i = 1; i <= m; i++) {
        if (w[i - 1] <= 0.0) {
            *ier = 10;
            return;
        }
        if ((x[i - 1] < xb) || (x[i - 1] > xe)) { return; }
        if ((y[i - 1] < yb) || (y[i - 1] > ye)) { return; }
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
    kn = 1;
    ki = kn + m;
    lq = 2;
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
           *nx, tx, *ny, ty, c, fp, &wrk1[0], &wrk1[lfp - 1], &wrk1[lco - 1],
           &wrk1[lf - 1], &wrk1[lff - 1], &wrk1[la - 1], &wrk1[lq - 1],
           &wrk1[lbx - 1], &wrk1[lby - 1], &wrk1[lsx - 1], &wrk1[lsy - 1],
           &wrk1[lh - 1], &iwrk[ki - 1], &iwrk[kn - 1], wrk2, lwrk2, ier);
}


void parder(const double *tx, int nx, const double *ty, int ny, double *c,
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
        for (int j = 0; j < nux; j++) {
            ak = kkx;
            nxx--;
            l1 = lx;
            m0 = 1;
            for (int i = 0; i < nxx; i++) {
                l1++;
                l2 = l1 + kkx;
                fac = tx[l2 - 1] - tx[l1 - 1];
                if (fac <= 0.0) { continue; }
                for (int m = 0; m < nyy; m++) {
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
        for (int j = 0; j < nuy; j++) {
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
                wrk[m0 - 1] = wrk[m1 - 1];
                m0++;
                m1++;
            }
            m1 += nuy;
        }
    }

    // Partition the workspace and evaluate the partial derivative
    int iwx = nxx * nyy;
    int iwy = iwx + mx * (kx1 - nux);
    fpbisp(&tx[nux], nx - 2*nux, &ty[nuy], ny - 2*nuy, wrk, kkx, kky, x, mx, y, my, z, &wrk[iwx], &wrk[iwy], iwrk, &iwrk[mx]);
}


void pardeu(const double *tx, int nx, const double *ty, int ny, double *c,
            int kx, int ky, int nux, int nuy, const double *x, const double *y,
            double *z, int m, double *wrk, int lwrk, int *iwrk, int kwrk, int *ier) {
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
                wrk[m0 - 1] = wrk[m1 - 1];
                m0++;
                m1++;
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


void bispeu(const double *tx, int nx, const double *ty, int ny, const double *c,
            int kx, int ky, const double *x, const double *y, double *z, int m,
            double *wrk, int lwrk, int *ier)
{

    int lwest = kx + ky + 2;

    // before starting computations a data check is made. if the input data
    // are invalid control is immediately repassed to the calling program.
    *ier = 10;
    if (lwrk < lwest) { return; }
    if (m < 1) { return; }

    *ier = 0;
    for (int i = 1; i <= m; i++) {
        fpbisp(tx, nx, ty, ny, c, kx, ky, &x[i - 1], 1, &y[i - 1], 1, &z[i - 1], wrk, &wrk[kx + 1], wrk, &wrk[1]);
    }
}


void bispev(const double *tx, int nx, const double *ty, int ny, const double *c,
            int kx, int ky, const double *x, int mx, const double *y, int my,
            double *z, double *wrk, int lwrk, int *iwrk, int kwrk, int *ier)
{
    int i, iw, lwest;

    // before starting computations a data check is made. if the input data
    // are invalid control is immediately repassed to the calling program.
    *ier = 10;
    lwest = (kx + 1) * mx + (ky + 1) * my;
    if (lwrk < lwest || kwrk < (mx + my) || mx < 1 || my < 1) {
        return;
    }

    if (mx > 1) {
        for (i = 2; i <= mx; i++) {
            if (x[i - 1] < x[i - 2]) {
                return;
            }
        }
    }

    if (my > 1) {
        for (i = 2; i <= my; i++) {
            if (y[i - 1] < y[i - 2]) {
                return;
            }
        }
    }

    *ier = 0;
    // Partition the working space and evaluate the bivariate spline
    iw = mx * (kx + 1) + 1;
    fpbisp(tx, nx, ty, ny, c, kx, ky, x, mx, y, my, z, wrk, &wrk[iw - 1], iwrk, &iwrk[mx - 1]);
}


void fpsphe(
    const int iopt, const int m, const double *teta, const double *phi, const double *r, const double *w,
    const double s, const int ntest, const int npest, const double eta, const double tol, const int maxit,
    const int ib1, const int ib3, const int nc, const int ncc, const int intest, const int nrest,
    int *nt, double *tt, int *np, double *tp, double *c, double *fp, double *sup, double *fpint, double *coord,
    double *f, double *ff, double *row, double *coco, double *cosi, double *a, double *q, double *bt, double *bp,
    double *spt, double *spp, double *h, int *index, int *nummer, double *wrk, const int lwrk, int *ier)
{
    double aa, acc, arg, cn, co, c1, dmax, d1, d2, eps, facc, facs, fac1, fac2, fn;
    double fpmax, fpms, f1, f2, f3, hti, htj, p, pi, pinv, piv, pi2, p1, p2, p3, ri, si;
    double sigma, sq, store, wi, rn, one, con1, con9, con4, half, ten;
    int i, iband, iband1, ii, irot, j, jlt, jrot, j1, j2, l, la, lf, lh, ll, lp, lt, lwest, l1, l2, l3, l4;
    int ncof, ncoff, npp, np4, nreg, nrint, nrr, nr1, ntt, nt4, num, num1, rank, in, iter, i1, i2, i3;
    double ht[4], hp[4];

    // set constants
    one = 1.0;
    con1 = 0.1;
    con9 = 0.9;
    con4 = 0.04;
    half = 0.5;
    ten = 10.0;
    pi = 3.1415926535897932384626433832795;
    pi2 = pi + pi;
    eps = sqrt(eta);


    if (iopt >= 0) {
        acc = tol * s;
        if ((iopt == 0) || ((iopt > 0) && (s >= *sup))) {
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
            fac1 = pi * (one + half);
            fac2 = (one + one) / (pi * pi * pi);
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
                } else if (f1 != 0.0) {
                    fpgivs(f1, &d2, &co, &si);
                    fprota(co, si, &ri, &c1);
                }
                *sup += ri * ri;
            }
            if (d2 != 0.0) c1 = c1 / d2;
            if (d1 != 0.0) cn = (cn - aa * c1) / d1;
            // find the b-spline representation of this least-squares polynomial
            *nt = 8;
            *np = 8;
            for (i = 1; i <= 4; i++) {
                c[i - 1] = c1;
                c[i + 3] = c1;
                c[i + 7] = cn;
                c[i + 11] = cn;
                tt[i - 1] = 0.0;
                tt[i + 3] = pi;
                tp[i - 1] = 0.0;
                tp[i + 3] = pi2;
            }
            *fp = *sup;
            // test whether the least-squares polynomial is an acceptable solution
            fpms = *sup - s;
            if (fpms < acc){
                *ier = -2;
                return;
            }
        }
        // test whether we cannot further increase the number of knots.
        if ((npest < 11) || (ntest < 9)) {
            *ier = 1;
            return;
        }
        *np = 11;
        tp[4] = pi * half;
        tp[5] = pi;
        tp[6] = tp[4] + pi;
        *nt = 9;
        tt[4] = tp[4];
        // 60
    } else if ((iopt > 0) && (s < *sup)) {
        if (*np < 11) {
            // test whether we cannot further increase the number of knots.
            if ((npest < 11) || (ntest < 9)) {
                *ier = 1;
                return;
            }
            *np = 11;
            tp[4] = pi * half;
            tp[5] = pi;
            tp[6] = tp[4] + pi;
            *nt = 9;
            tt[4] = tp[4];
            // 60
        }
    }
    // 70

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
            l1 = l1 + 1;
            l2 = l2 - 1;
            l3 = l3 + 1;
            l4 = l4 - 1;
            tp[l2 - 1] = tp[l4 - 1] - pi2;
            tp[l3 - 1] = tp[l1 - 1] + pi2;
        }
        l = *nt;
        for (i = 1; i <= 4; i++) {
            tt[i - 1] = 0.0;
            tt[l - 1] = pi;
            l = l - 1;
        }
        // find nrint, the total number of knot intervals and nreg, the number
        // of panels in which the approximation domain is subdivided by the
        // intersection of knots.
        ntt = *nt - 7;
        npp = *np - 7;
        nrr = npp / 2;
        nr1 = nrr + 1;
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
                a[(i - 1) + (j - 1)*ib1] = 0.0;
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
                    fprota(co, si, &row[l - 1], &a[(j - 1) +  (i2 - 1)*ib1]);
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
                a[(i - 1) + (j - 1)*ib1] = 0.0;
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
                    spp[(in - 1) + (i - 1)*4] = hp[i - 1];
                    spt[(in - 1) + (i - 1)*4] = ht[i - 1];
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
                            j1 = j1 + 2;
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
                        fprota(co, si, &h[j - 1], &a[(irot - 1) + (i2 - 1)*ib1]);
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
                    q[i + j*ib3] = a[i + j*ib1];
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
            if (ncof < rank) {
                *ier = -rank;
            }
            return;
        }
        fpms = *fp - s;
        if (fabs(fpms) < acc) {
            if (*fp <= 0.0) {
                *ier = -1;
                *fp = 0.0;
            }
            if (ncof < rank) {
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
                    hti = spt[(in - 1) + (i - 1)*4];
                    j1 = i1;
                    for (j = 1; j <= 4; j++) {
                        j1++;
                        store += hti * spp[(in - 1) + (j - 1) * 4] * c[j1 - 1];
                    }
                    i1 = i1 + np4;
                }
                double wval = w[in - 1];
                double rval = r[in - 1];
                store = (wval * (rval - store)) * (wval * (rval - store));
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
                if (arg >= pi) {
                    arg = arg - pi;
                    l4 = l4 - nrr;
                }
                fpint[l - 1] = 0.0;
                fac1 = tp[l4 - 1] - arg;
                fac2 = arg - tp[l4 - 2];
            } else {
                // addition in the teta-direction
                l4 = l + 4;
                fpint[l - 1] = 0.0;
                fac1 = tt[l4 - 1] - arg;
                fac2 = arg - tt[l4 - 1 - 1];
            }


            if ((fac1 > (ten * fac2)) || (fac2 > (ten * fac1))) { continue; }

            if (l > ntt) {
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
                    tp[j - 1] = tp[i - 1] + pi;
                }
                break;
            }   else {
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
    fpdisc(tt, nt, 5, bt, ntest);
    // evaluate the discontinuity jumps of the 3-th order derivative of
    // the b-splines at the knots tp(l),l=5,...,np-4.
    fpdisc(tp, np, 5, bp, npest);

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
                q[(i - 1) + (j - 1)*ib3] = 0.0;
            }
            for (j = 1; j <= iband; j++) {
                q[(i - 1) + (j - 1)*ib3] = a[(i - 1) + (j - 1)*ib1];
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
                    if (piv != 0.0) {
                        // calculate the parameters of the givens transformation.
                        fpgivs(piv, &q[irot - 1], &co, &si);
                        // apply that givens transformation to the right hand side.
                        fprota(co, si, &ri, &ff[irot - 1]);
                        if (i2 == 0) { break; }
                        // apply that givens transformation to the left hand side.
                        for (l = 1; l <= i2; l++) {
                            l1 = l + 1;
                            fprota(co, si, &h[l1 - 1], &q[(irot - 1) + l*ib3]);
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
                    j1 = j1 + ij;
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
                    if (piv != 0.0) {
                        // calculate the parameters of the givens transformation.
                        fpgivs(piv, &q[irot - 1], &co, &si);
                        // apply that givens transformation to the right hand side.
                        fprota(co, si, &ri, &ff[irot - 1]);
                        if (i2 == 0) { break; }
                        // apply that givens transformation to the left hand side.
                        for (l = 1; l <= i2; l++) {
                            l1 = l + 1;
                            fprota(co, si, &h[l1 - 1], &q[(irot - 1) + l*ib3]);
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
                    hti = spt[(in - 1) + (i - 1)*4];
                    j1 = i1;
                    for (j = 1; j <= 4; j++) {
                        j1++;
                        store += hti * spp[(in - 1) + (j - 1)*4] * c[j1 - 1];
                    }
                    i1 += np4;
                }
                *fp += (w[in - 1] * (r[in - 1] - store)) * (w[in - 1] * (r[in - 1] - store));
                in = nummer[in - 1];
            }
        }
        // test whether the approximation sp(teta,phi) is an acceptable solution
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


void
sphere(const int iopt, const int m, const double *teta, const double *phi, const double *r,
       const double *w, const double s, const int ntest, const int npest, const double eps,
       int *nt, double *tt, int *np, double *tp, double *c, double *fp,
       double *wrk1, const int lwrk1, double *wrk2, const int lwrk2, int *iwrk, const int kwrk, int *ier)
{
    double tol, pi, pi2, one;
    int i, ib1, ib3, ki, kn, kwest, la, lbt, lcc, lcs, lro, j;
    int lbp, lco, lf, lff, lfp, lh, lq, lst, lsp, lwest, maxit, ncest, ncc, ntt;
    int npp, nreg, nrint, ncof, nt4, np4;

    // set constants
    one = 1.0;
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
        pi = atan(one) * 4;
        pi2 = pi + pi;
        for (i = 1; i <= m; i++) {
            if (w[i - 1] <= 0.0) { return; }
            if ((teta[i - 1] < 0.0) || (teta[i - 1] > pi)) { return; }
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
                    if ((tt[j - 1] <= tt[j - 2]) || (tt[j - 1] >= pi)) { return; }
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
