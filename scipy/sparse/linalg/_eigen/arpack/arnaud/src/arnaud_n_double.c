#include "arnaud_n_double.h"
#include <float.h>

typedef int ARNAUD_compare_cfunc(const double, const double, const double, const double);

static int sortc_LM(const double, const double, const double, const double);
static int sortc_SM(const double, const double, const double, const double);
static int sortc_LR(const double, const double, const double, const double);
static int sortc_SR(const double, const double, const double, const double);
static int sortc_LI(const double, const double, const double, const double);
static int sortc_SI(const double, const double, const double, const double);

static const double unfl = DBL_MIN;    // 2.2250738585072014e-308
// static const double ovfl = DBL_MAX; // 1.0 / 2.2250738585072014e-308;
static const double ulp = DBL_EPSILON; // 2.220446049250313e-16;

static void dnaup2(struct ARNAUD_state_d*, double*, double*, int, double*, int, double*, double*, double*, double*, int, double*, int*, double*);
static void dnconv(int n, double* ritzr, double* ritzi, double* bounds, const double tol, int* nconv);
static void dneigh(double*,int,double*,int,double*,double*,double*,double*,int,double*,int*);
static void dnaitr(struct ARNAUD_state_d*,int,int,double*,double*,double*,int,double*,int,int*,double*);
static void dnapps(int,int*,int,double*,double*,double*,int,double*,int,double*,double*,int,double*,double*);
static void dngets(struct ARNAUD_state_d*,int*,int*,double*,double*,double*);
static void dsortc(const enum ARNAUD_which w, const int apply, const int n, double* xreal, double* ximag, double* y);
static void dgetv0(struct ARNAUD_state_d *V, int initv, int n, int j, double* v, int ldv, double* resid, double* rnorm, int* ipntr, double* workd);

enum ARNAUD_neupd_type {
    REGULAR = 0,
    SHIFTI,
    REALPART,
    IMAGPART
};


void
ARNAUD_dneupd(struct ARNAUD_state_d *V, int rvec, int howmny, int* select,
       double* dr, double* di, double* z, int ldz, double sigmar, double sigmai,
       double* workev, double* resid, double* v, int ldv, int* ipntr, double* workd,
       double* workl)
{
    const double eps23 = pow(ulp, 2.0 / 3.0);
    int ibd, iconj, ih, iheigr, iheigi, ihbds, iuptri, invsub, iri, irr, j, jj;
    int bounds, k, ldh, ldq, np, numcnv, reord, ritzr, ritzi;
    int iwork[1] = { 0 };
    int ierr = 0, int1 = 1, tmp_int = 0, nconv2 = 0, outncv;
    double conds, rnorm, sep, temp, temp1, dbl0 = 0.0, dbl1 = 1.0, dblm1 = -1.0;
    double vl[1] = { 0.0 };
    enum ARNAUD_neupd_type TYP;

    if (V->nconv <= 0) {
        ierr = -14;
    } else if (V->n <= 0) {
        ierr = -1;
    } else if (V->nev <= 0) {
        ierr = -2;
    } else if ((V->ncv <= V->nev + 1) || (V->ncv > V->n)) {
        ierr = -3;
    } else if ((V->which > 5) || (V->which < 0)) {
        ierr = -5;
    } else if ((V->bmat != 0) && (V->bmat != 1)) {
        ierr = -6;
    } else if ((rvec) && ((howmny < 0) || (howmny > 2))) {
        ierr = -13;
    } else if (howmny == 2) {
        ierr = -12;  // NotImplementedError
    }

    if ((V->mode == 1) || (V->mode == 2)) {
        TYP = REGULAR;
    } else if ((V->mode == 3) && (sigmai == 0.0)) {
        TYP = SHIFTI;
    } else if (V->mode == 3) {
        TYP = REALPART;
    } else if (V->mode == 4) {
        TYP = IMAGPART;
    } else {
        ierr = -10;
    }

    if ((V->mode == 1) && (V->bmat)) { ierr = -11; }

    if (ierr != 0) {
        V->info = ierr;
        return;
    }

    //  Pointer into WORKL for address of H, RITZ, BOUNDS, Q
    //  etc... and the remaining workspace.
    //  Also update pointer to be used on output.
    //  Memory is laid out as follows:
    //  workl(1:ncv*ncv) := generated Hessenberg matrix
    //  workl(ncv*ncv+1:ncv*ncv+2*ncv) := real and imaginary parts of ritz values
    //  workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := error bounds

    //  The following is used and set by DNEUPD .
    //  workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed real part of the Ritz values.
    //  workl(ncv*ncv+4*ncv+1:ncv*ncv+5*ncv) := The untransformed imaginary part of the Ritz values.
    //  workl(ncv*ncv+5*ncv+1:ncv*ncv+6*ncv) := The untransformed error bounds of the Ritz values
    //  workl(ncv*ncv+6*ncv+1:2*ncv*ncv+6*ncv) := Holds the upper quasi-triangular matrix for H
    //  workl(2*ncv*ncv+6*ncv+1: 3*ncv*ncv+6*ncv) := Holds the associated matrix representation of the invariant subspace for H.
    //  GRAND total of NCV * ( 3 * NCV + 6 ) locations.

    ih     = ipntr[4];
    ritzr  = ipntr[5];
    ritzi  = ipntr[6];
    bounds = ipntr[7];
    ldh = V->ncv;
    ldq = V->ncv;
    iheigr = bounds + ldh;
    iheigi = iheigr + ldh;
    ihbds  = iheigi + ldh;
    iuptri = ihbds  + ldh;
    invsub = iuptri + ldh*V->ncv;
    ipntr[8]  = iheigr;
    ipntr[9]  = iheigi;
    ipntr[10] = ihbds;
    ipntr[11] = iuptri;
    ipntr[12] = invsub;

    //  irr points to the REAL part of the Ritz
    //      values computed by _neigh before
    //      exiting _naup2.
    //  iri points to the IMAGINARY part of the
    //      Ritz values computed by _neigh
    //      before exiting _naup2.
    //  ibd points to the Ritz estimates
    //      computed by _neigh before exiting
    //      _naup2.

    irr = ipntr[13] + (V->ncv)*(V->ncv);
    iri = irr + V->ncv;
    ibd = iri + V->ncv;

    //  RNORM is B-norm of the RESID(1:N).
    rnorm = workl[ih+2];
    workl[ih+2] = 0.0;

    if (rvec) {
        reord = 0;

        //  Use the temporary bounds array to store indices
        //  These will be used to mark the select array later

        for (j = 0; j < V->ncv; j++)
        {
            workl[bounds + j] = j;
            select[j] = 0;
        }
        // 10

        //  Select the wanted Ritz values.
        //  Sort the Ritz values so that the
        //  wanted ones appear at the tailing
        //  NEV positions of workl(irr) and
        //  workl(iri).  Move the corresponding
        //  error estimates in workl(bound)
        //  accordingly.

        np = V->ncv - V->nev;
        dngets(V, &V->nev, &np, &workl[irr], &workl[iri], &workl[bounds]);

        //  Record indices of the converged wanted Ritz values
        //  Mark the select array for possible reordering

        numcnv = 0;
        for (j = 1; j <= V->ncv; j++)
        {
            temp1 = fmax(eps23, hypot(workl[irr + V->ncv - j], workl[iri + V->ncv - j]));

            jj = (int)workl[bounds + V->ncv - j];

            if ((numcnv < V->nconv) && (workl[ibd + jj] <= V->tol*temp1))
            {
                select[jj] = 1;
                numcnv += 1;
                if (jj > V->nconv - 1) { reord = 1; }
            }
        }
        // 11

        //  Check the count (numcnv) of converged Ritz values with
        //  the number (nconv) reported by dnaupd.  If these two
        //  are different then there has probably been an error
        //  caused by incorrect passing of the dnaupd data.

        if (numcnv != V->nconv)
        {
            V->info = -15;
            return;
        }

        //  Call LAPACK routine dlahqr  to compute the real Schur form
        //  of the upper Hessenberg matrix returned by DNAUPD .
        //  Make a copy of the upper Hessenberg matrix.
        //  Initialize the Schur vector matrix Q to the identity.

        tmp_int = ldh*V->ncv;
        dcopy_(&tmp_int, &workl[ih], &int1, &workl[iuptri], &int1);
        dlaset_("A", &V->ncv, &V->ncv, &dbl0, &dbl1, &workl[invsub], &ldq);
        dlahqr_(&int1, &int1, &V->ncv, &int1, &V->ncv, &workl[iuptri], &ldh,
                &workl[iheigr], &workl[iheigi], &int1, &V->ncv, &workl[invsub],
                &ldq, &ierr);
        dcopy_(&V->ncv, &workl[invsub + V->ncv - 1], &ldq, &workl[ihbds], &int1);

        if (ierr != 0)
        {
            V->info = -8;
            return;
        }

        if (reord)
        {
            dtrsen_("N", "V", select, &V->ncv, &workl[iuptri], &ldh, &workl[invsub], &ldq,
                    &workl[iheigr], &workl[iheigi], &nconv2, &conds, &sep, &workl[ihbds],
                    &V->ncv, iwork, &int1, &ierr);

            if (nconv2 < V->nconv) { V->nconv = nconv2; }
            if (ierr == 1) {
                V->info = 1;
                return;
            }
        }

        //  Copy the last row of the Schur vector
        //  into workl(ihbds).  This will be used
        //  to compute the Ritz estimates of
        //  converged Ritz values.

        dcopy_(&V->ncv, &workl[invsub + V->ncv - 1], &ldq, &workl[ihbds], &int1);

        //  Place the computed eigenvalues of H into DR and DI
        //  if a spectral transformation was not used.

        if (TYP == REGULAR) {
            dcopy_(&V->nconv, &workl[iheigr], &int1, dr, &int1);
            dcopy_(&V->nconv, &workl[iheigi], &int1, di, &int1);
        }

        //  Compute the QR factorization of the matrix representing
        //  the wanted invariant subspace located in the first NCONV
        //  columns of workl(invsub,ldq).

        dgeqr2_(&V->ncv, &V->nconv, &workl[invsub], &ldq, workev, &workev[V->ncv], &ierr);

        //  * Postmultiply V by Q using dorm2r .
        //  * Copy the first NCONV columns of VQ into Z.
        //  * Postmultiply Z by R.
        //  The N by NCONV matrix Z is now a matrix representation
        //  of the approximate invariant subspace associated with
        //  the Ritz values in workl(iheigr) and workl(iheigi)
        //  The first NCONV columns of V are now approximate Schur
        //  vectors associated with the real upper quasi-triangular
        //  matrix of order NCONV in workl(iuptri)

        dorm2r_("R", "N", &V->n, &V->ncv, &V->nconv, &workl[invsub], &ldq, workev,
                v, &ldv, &workd[V->n], &ierr);

        dlacpy_("A", &V->n, &V->nconv, v, &ldv, z, &ldz);

        //  Perform both a column and row scaling if the
        //  diagonal element of workl(invsub,ldq) is negative
        //  I'm lazy and don't take advantage of the upper
        //  quasi-triangular form of workl(iuptri,ldq)
        //  Note that since Q is orthogonal, R is a diagonal
        //  matrix consisting of plus or minus ones

        for (j = 0; j < V->nconv; j++)
        {
            if (workl[invsub + j*ldq + j] < 0.0)
            {
                dscal_(&V->nconv, &dblm1, &workl[iuptri + j], &ldq);
                dscal_(&V->nconv, &dblm1, &workl[iuptri + j*ldq], &int1);
            }
        }
        // 20

        if (howmny == 0)
        {

            //  Compute the NCONV wanted eigenvectors of T
            //  located in workl(iuptri,ldq).

            for (j = 0; j < V->ncv; j++)
            {
                if (j < V->nconv)
                {
                    select[j] = 1;
                } else {
                    select[j] = 0;
                }
            }
            // 30

            dtrevc_("R", "S", select, &V->ncv, &workl[iuptri], &ldq, vl, &int1,
                    &workl[invsub], &ldq, &V->ncv, &outncv, workev, &ierr);

            if (ierr != 0)
            {
                V->info = -9;
                return;
            }

            //  Scale the returning eigenvectors so that their
            //  Euclidean norms are all one. LAPACK subroutine
            //  dtrevc  returns each eigenvector normalized so
            //  that the element of largest magnitude has
            //  magnitude 1;

            iconj = 0;
            for (j = 0; j < V->nconv; j++)
            {
                if (workl[iheigi + j] == 0.0)
                {

                    //  real eigenvalue case

                    temp = 1.0 / dnrm2_(&V->ncv, &workl[invsub + j*ldq], &int1);
                    dscal_(&V->ncv, &temp, &workl[invsub + j*ldq], &int1);

                } else {

                    //  Complex conjugate pair case. Note that
                    //  since the real and imaginary part of
                    //  the eigenvector are stored in consecutive
                    //  columns, we further normalize by the
                    //  square root of two.

                    if (iconj == 0)
                    {
                        temp = 1.0 / hypot(dnrm2_(&V->ncv, &workl[invsub + j*ldq], &int1),
                                           dnrm2_(&V->ncv, &workl[invsub + (j+1)*ldq], &int1));
                        dscal_(&V->ncv, &temp, &workl[invsub + j*ldq], &int1);
                        dscal_(&V->ncv, &temp, &workl[invsub + (j+1)*ldq], &int1);
                        iconj = 1;
                    } else {
                        iconj = 0;
                    }
                }
            }
            // 40

            dgemv_("T", &V->ncv, &V->nconv, &dbl1, &workl[invsub], &ldq, &workl[ihbds], &int1, &dbl0, workev, &int1);

            iconj = 0;
            for (j = 0; j < V->nconv; j++)
            {
                if (workl[iheigi + j] != 0.0)
                {

                    //  Complex conjugate pair case. Note that
                    //  since the real and imaginary part of
                    //  the eigenvector are stored in consecutive

                    if (iconj == 0)
                    {
                        workev[j] = hypot(workev[j], workev[j+1]);
                        workev[j+1] = workev[j];
                        iconj = 1;
                    } else {
                        iconj = 0;
                    }
                }
            }
            // 45

            //  Copy Ritz estimates into workl(ihbds)

            dcopy_(&V->nconv, workev, &int1, &workl[ihbds], &int1);

            //  Compute the QR factorization of the eigenvector matrix
            //  associated with leading portion of T in the first NCONV
            //  columns of workl(invsub,ldq).

            dgeqr2_(&V->ncv, &V->nconv, &workl[invsub], &ldq, workev, &workev[V->ncv], &ierr);

            //  * Postmultiply Z by Q.
            //  * Postmultiply Z by R.
            //  The N by NCONV matrix Z is now contains the
            //  Ritz vectors associated with the Ritz values
            //  in workl(iheigr) and workl(iheigi).

            dorm2r_("R", "N", &V->n, &V->ncv, &V->nconv, &workl[invsub], &ldq,
                    workev, z, &ldz, &workd[V->n], &ierr);

            dtrmm_("R", "U", "N", "N", &V->n, &V->nconv, &dbl1, &workl[invsub], &ldq, z, &ldz);

        }

    } else {

        //  An approximate invariant subspace is not needed.
        //  Place the Ritz values computed DNAUPD  into DR and DI

        dcopy_(&V->nconv, &workl[ritzr], &int1, dr, &int1);
        dcopy_(&V->nconv, &workl[ritzi], &int1, di, &int1);
        dcopy_(&V->nconv, &workl[ritzr], &int1, &workl[iheigr], &int1);
        dcopy_(&V->nconv, &workl[ritzi], &int1, &workl[iheigi], &int1);
        dcopy_(&V->nconv, &workl[bounds], &int1, &workl[ihbds], &int1);
    }

    //  Transform the Ritz values and possibly vectors
    //  and corresponding error bounds of OP to those
    //  of A*x = lambda*B*x.

    if (TYP == REGULAR)
    {
        if (rvec)
        {
            dscal_(&V->ncv, &rnorm, &workl[ihbds], &int1);
        }
    } else {

        //    A spectral transformation was used.
        //  * Determine the Ritz estimates of the
        //    Ritz values in the original system.

        if (TYP == SHIFTI)
        {
            if (rvec)
            {
                dscal_(&V->ncv, &rnorm, &workl[ihbds], &int1);
            }

            for (k = 0; k < V->ncv; k++)
            {
                temp = hypot(workl[iheigr+k], workl[iheigi+k]);
                workl[ihbds+k] = fabs(workl[ihbds+k]) / temp / temp;
            }
            // 50

        }

        //  *  Transform the Ritz values back to the original system.
        //     For TYPE = 'SHIFTI' the transformation is
        //              lambda = 1/theta + sigma
        //     For TYPE = 'REALPT' or 'IMAGPT' the user must from
        //     Rayleigh quotients or a projection. See remark 3 above.
        //  NOTES:
        //  *The Ritz vectors are not affected by the transformation.

        if (TYP == SHIFTI)
        {
            for (k = 0; k < V->ncv; k++)
            {
                temp = hypot(workl[iheigr+k], workl[iheigi+k]);
                workl[iheigr+k] =  workl[iheigr+k] / temp / temp + sigmar;
                workl[iheigi+k] = -workl[iheigi+k] / temp / temp + sigmai;
            }
            // 80

            dcopy_(&V->nconv, &workl[iheigr], &int1, dr, &int1);
            dcopy_(&V->nconv, &workl[iheigi], &int1, di, &int1);

        } else if ((TYP == REALPART) || (TYP == IMAGPART)) {
            dcopy_(&V->nconv, &workl[iheigr], &int1, dr, &int1);
            dcopy_(&V->nconv, &workl[iheigi], &int1, di, &int1);
        }
    }

    //  Eigenvector Purification step. Formally perform
    //  one of inverse subspace iteration. Only used
    //  for MODE = 2.

    if ((rvec) && (howmny == 0) && (TYP == SHIFTI))
    {

        //  Purify the computed Ritz vectors by adding a
        //  little bit of the residual vector:
        //                       T
        //           resid(:)*( e    s ) / theta
        //                       NCV
        //  where H s = s theta. Remember that when theta
        //  has nonzero imaginary part, the corresponding
        //  Ritz vector is stored across two columns of Z.

        iconj = 0;
        for (j = 0; j < V->nconv; j++)
        {
            if ((workl[iheigi+j] == 0.0) && (workl[iheigr+j] != 0.0))
            {
                workev[j] = workl[invsub + j*ldq + V->ncv - 1] / workl[iheigr+j];
            } else if (iconj == 0) {

                temp = hypot(workl[iheigr+j], workl[iheigi+j]);
                if (temp != 0.0)
                {
                    workev[j] = (workl[invsub + j*ldq + V->ncv - 1]*workl[iheigr+j] +
                                 workl[invsub + (j+1)*ldq + V->ncv - 1]*workl[iheigi+j]
                                ) / temp / temp;
                    workev[j+1] = (workl[invsub + (j+1)*ldq + V->ncv - 1]*workl[iheigr+j] -
                                 workl[invsub + j*ldq + V->ncv - 1]*workl[iheigi+j]
                                ) / temp / temp;
                }
                iconj = 1;
            } else {
                iconj = 0;
            }
        }
        // 110

        //  Perform a rank one update to Z and
        //  purify all the Ritz vectors together.

        dger_(&V->n, &V->nconv, &dbl1, resid, &int1, workev, &int1, z, &ldz);
    }

    return;
}

void
ARNAUD_dnaupd(struct ARNAUD_state_d *V, double* resid, double* v,
              int ldv, int* ipntr, double* workd, double* workl)
{
    int bounds, ih, iq, iw, j, ldh, ldq, next, iritzi, iritzr;

    if (V->ido == ido_FIRST)
    {

        // perform basic checks
        if (V->n <= 0) {
            V->info = -1;
        } else if (V->nev <= 0) {
            V->info = -2;
        } else if ((V->ncv < V->nev + 1) || (V->ncv > V->n)) {
            V->info = -3;
        } else if (V->maxiter <= 0) {
            V->info = -4;
        } else if ((V->which < 0) || (V->which > 5)) {
            V->info = -5;
        } else if ((V->bmat != 0) && (V->bmat != 1)) {
            V->info = -6;
        } else if ((V->mode < 1) || (V->mode > 4)) {
            V->info = -10;
        } else if ((V->mode == 1) && (V->bmat == 1)) {
            V->info = -11;
        } else if ((V->shift != 0) && (V->shift != 1)) {
            V->info = -12;
        }

        if (V->info < 0) {
            V->ido = ido_DONE;
            return;
        }

        if (V->tol <= 0.0) {
            V->tol = ulp;
        }
        V->np = V->ncv - V->nev;

        for (j = 0; j < 3 * (V->ncv)*(V->ncv) + 6*(V->ncv); j++)
        {
            workl[j] = 0.0;
        }
    }

    //  Pointer into WORKL for address of H, RITZ, BOUNDS, Q
    //  etc... and the remaining workspace.
    //  Also update pointer to be used on output.
    //  Memory is laid out as follows:
    //
    //  workl[0:ncv*ncv] := generated Hessenberg matrix
    //  workl[ncv**2:ncv**2+2*ncv] := ritz.real and ritz.imag values
    //  workl[ncv**2+2*ncv:ncv*ncv+3*ncv] := error bounds
    //  workl[ncv**2+3*ncv+1:2*ncv*ncv+3*ncv] := rotation matrix Q
    //  workl[2*ncv**2+3*ncv:3*ncv*ncv+6*ncv] := workspace
    //
    //  The final workspace is needed by subroutine dneigh  called
    //  by dnaup2 . Subroutine dneigh  calls LAPACK routines for
    //  calculating eigenvalues and the last row of the eigenvector
    //  matrix.

    ldh    = V->ncv;
    ldq    = V->ncv;
    ih     = 0;
    iritzr = ih     + ldh*V->ncv;
    iritzi = iritzr + V->ncv;
    bounds = iritzi + V->ncv;
    iq     = bounds + V->ncv;
    iw     = iq     + ldq*V->ncv;
    next   = iw     + (V->ncv*V->ncv) + 3*V->ncv;

    ipntr[3] = next;
    ipntr[4] = ih;
    ipntr[5] = iritzr;
    ipntr[6] = iritzi;
    ipntr[7] = bounds;
    ipntr[13]  = iw;

    //  Carry out the Implicitly restarted Arnoldi Iteration.

    dnaup2(V, resid, v, ldv, &workl[ih], ldh, &workl[iritzr], &workl[iritzi], &workl[bounds], &workl[iq], ldq, &workl[iw], ipntr, workd);

    //  ido != DONE implies use of reverse communication
    //  to compute operations involving OP or shifts.

    if (V->ido != ido_DONE) { return; }

    V->nconv = V->np;
    // iparam(9) = nopx
    // iparam(10) = nbx
    // iparam(11) = nrorth

    if (V->info < 0) { return; }
    if (V->info == 2) { V->info = 3; }

    return;
}

void
dnaup2(struct ARNAUD_state_d *V, double* resid, double* v, int ldv,
       double* h, int ldh, double* ritzr, double* ritzi, double* bounds,
       double* q, int ldq, double* workl, int* ipntr, double* workd)
{
    enum ARNAUD_which temp_which;
    int int1 = 1, j, tmp_int;
    const double eps23 = pow(ulp, 2.0 / 3.0);
    double temp = 0.0;

    if (V->ido == ido_FIRST)
    {
        V->aup2_nev0 = V->nev;
        V->aup2_np0 = V->np;

        //  kplusp is the bound on the largest
        //         Lanczos factorization built.
        //  nconv is the current number of
        //         "converged" eigenvlues.
        //  iter is the counter on the current
        //       iteration step.

        V->aup2_kplusp = V->nev + V->np;
        V->nconv = 0;
        V->aup2_iter = 0;

        //  Set flags for computing the first NEV
        //  steps of the Arnoldi factorization.

        V->aup2_getv0 = 1;
        V->aup2_update = 0;
        V->aup2_ushift = 0;
        V->aup2_cnorm = 0;

        if (V->info != 0)
        {

            //  User provides the initial residual vector.

            V->aup2_initv = 1;
            V->info = 0;
        } else {
            V->aup2_initv = 0;
        }
    }

    //  Get a possibly random starting vector and
    //  force it into the range of the operator OP.

    if (V->aup2_getv0)
    {
        V->getv0_itry = 1;
        dgetv0(V, V->aup2_initv, V->n, 0, v, ldv, resid, &V->aup2_rnorm, ipntr, workd);
        if (V->ido != ido_DONE) { return; }
        if (V->aup2_rnorm == 0.0)
        {
            V->info = -9;
            V->ido = ido_DONE;
            return;
        }
        V->aup2_getv0 = 0;
        V->ido = ido_FIRST;
    }

    //  Back from reverse communication :
    //  continue with update step

    if (V->aup2_update) { goto LINE20; }

    //  Back from computing user specified shifts

    if (V->aup2_ushift) { goto LINE50; }

    //  Back from computing residual norm
    //  at the end of the current iteration

    if (V->aup2_cnorm) { goto LINE100; }

    //  Compute the first NEV steps of the Arnoldi factorization

    dnaitr(V, 0, V->nev, resid, &V->aup2_rnorm, v, ldv, h, ldh, ipntr, workd);

    //  ido .ne. 99 implies use of reverse communication
    //  to compute operations involving OP and possibly B

    if (V->ido != ido_DONE) { return; }

    if (V->info > 0)
    {
        V->np = V->info;
        V->iter = V->aup2_iter;
        V->info = -9999;
        V->ido = ido_DONE;
        return;
    }

    //
    //            M A I N  ARNOLDI  I T E R A T I O N  L O O P
    //            Each iteration implicitly restarts the Arnoldi
    //            factorization in place.
    //

LINE1000:
    V->aup2_iter += 1;

    //  Compute NP additional steps of the Arnoldi factorization.
    //  Adjust NP since NEV might have been updated by last call
    //  to the shift application routine dnapps .

    V->np = V->aup2_kplusp - V->nev;

    //  Compute NP additional steps of the Arnoldi factorization.

    V->ido = ido_FIRST;

LINE20:
    V->aup2_update = 1;

    dnaitr(V, V->nev, V->np, resid, &V->aup2_rnorm, v, ldv, h, ldh, ipntr, workd);

    //  ido .ne. 99 implies use of reverse communication
    //  to compute operations involving OP and possibly B

    if (V->ido != ido_DONE) { return; }

    if (V->info > 0) {
        V->np = V->info;
        V->iter = V->aup2_iter;
        V->info = -9999;
        V->ido = ido_DONE;
        return;
    }
    V->aup2_update = 0;

    //  Compute the eigenvalues and corresponding error bounds
    //  of the current upper Hessenberg matrix.

    dneigh(&V->aup2_rnorm, V->aup2_kplusp, h, ldh, ritzr, ritzi, bounds, q, ldq, workl, &V->info);

    if (V->info != 0)
    {
       V->info = -8;
       V->ido = ido_DONE;
       return;
    }

    //  Make a copy of eigenvalues and corresponding error
    //  bounds obtained from dneigh.

    tmp_int = V->aup2_kplusp * V->aup2_kplusp;
    dcopy_(&V->aup2_kplusp, ritzr, &int1, &workl[tmp_int], &int1);
    tmp_int += V->aup2_kplusp;
    dcopy_(&V->aup2_kplusp, ritzi, &int1, &workl[tmp_int], &int1);
    tmp_int += V->aup2_kplusp;
    dcopy_(&V->aup2_kplusp, bounds, &int1, &workl[tmp_int], &int1);

    //  Select the wanted Ritz values and their bounds
    //  to be used in the convergence test.
    //  The wanted part of the spectrum and corresponding
    //  error bounds are in the last NEV loc. of RITZR,
    //  RITZI and BOUNDS respectively. The variables NEV
    //  and NP may be updated if the NEV-th wanted Ritz
    //  value has a non zero imaginary part. In this case
    //  NEV is increased by one and NP decreased by one.
    //  NOTE: The last two arguments of dngets  are no
    //  longer used as of version 2.1.

    V->nev = V->aup2_nev0;
    V->np = V->aup2_np0;
    V->aup2_numcnv = V->nev;

    dngets(V, &V->nev, &V->np, ritzr, ritzi, bounds);

    if (V->nev == V->aup2_nev0 + 1) { V->aup2_numcnv = V->aup2_nev0 + 1;}

    //  Convergence test.

    dcopy_(&V->nev, &bounds[V->np], &int1, &workl[2*V->np], &int1);
    dnconv(V->nev, &ritzr[V->np], &ritzi[V->np], &workl[2*V->np], V->tol, &V->nconv);

    //  Count the number of unwanted Ritz values that have zero
    //  Ritz estimates. If any Ritz estimates are equal to zero
    //  then a leading block of H of order equal to at least
    //  the number of Ritz values with zero Ritz estimates has
    //  split off. None of these Ritz values may be removed by
    //  shifting. Decrease NP the number of shifts to apply. If
    //  no shifts may be applied, then prepare to exit

    // We are modifying V->np hence the temporary variable.
    int nptemp = V->np;

    for (j = 0; j < nptemp; j++)
    {
        if (bounds[j] == 0.0)
        {
            V->np -= 1;
            V->nev += 1;
        }
    }
    // 30

    if ((V->nconv >= V->aup2_numcnv) || (V->aup2_iter > V->maxiter) || (V->np == 0))
    {

        //  Prepare to exit. Put the converged Ritz values
        //  and corresponding bounds in RITZ(1:NCONV) and
        //  BOUNDS(1:NCONV) respectively. Then sort. Be
        //  careful when NCONV > NP

        //   Use h( 3,1 ) as storage to communicate
        //   rnorm to _neupd if needed

         h[2] = V->aup2_rnorm;

        //  To be consistent with dngets , we first do a
        //  pre-processing sort in order to keep complex
        //  conjugate pairs together.  This is similar
        //  to the pre-processing sort used in dngets
        //  except that the sort is done in the opposite
        //  order.

        // Translation note: Is this all because ARNAUD did not have complex sort?

        if (V->which == which_LM) { temp_which = which_SR; }
        if (V->which == which_SM) { temp_which = which_LR; }
        if (V->which == which_LR) { temp_which = which_SM; }
        if (V->which == which_SR) { temp_which = which_LM; }
        if (V->which == which_LI) { temp_which = which_SM; }
        if (V->which == which_SI) { temp_which = which_LM; }

        dsortc(temp_which, 1, V->aup2_kplusp, ritzr, ritzi, bounds);

        if (V->which == which_LM) { temp_which = which_SM; }
        if (V->which == which_SM) { temp_which = which_LM; }
        if (V->which == which_LR) { temp_which = which_SR; }
        if (V->which == which_SR) { temp_which = which_LR; }
        if (V->which == which_LI) { temp_which = which_SI; }
        if (V->which == which_SI) { temp_which = which_LI; }

        dsortc(temp_which, 1, V->aup2_kplusp, ritzr, ritzi, bounds);

        //  Scale the Ritz estimate of each Ritz value
        //  by 1 / max(eps23,magnitude of the Ritz value).

        for (j = 0; j < V->aup2_numcnv; j++)
        {
            temp = fmax(eps23, hypot(ritzr[j], ritzi[j]));
            bounds[j] = bounds[j] / temp;
        }
        // 35

        //  Sort the Ritz values according to the scaled Ritz
        //  estimates.  This will push all the converged ones
        //  towards the front of ritzr, ritzi, bounds
        //  (in the case when NCONV < NEV.)

        temp_which = which_LR;
        dsortc(temp_which, 1, V->aup2_numcnv, bounds, ritzr, ritzi);

        //  Scale the Ritz estimate back to its original
        //  value.

        for (j = 0; j < V->aup2_numcnv; j++)
        {
            temp = fmax(eps23, hypot(ritzr[j], ritzi[j]));
            bounds[j] = bounds[j] * temp;
        }
        // 40

        //  Sort the converged Ritz values again so that
        //  the "threshold" value appears at the front of
        //  ritzr, ritzi and bound.

        dsortc(V->which, 1, V->nconv, ritzr, ritzi, bounds);

        if ((V->aup2_iter > V->maxiter) && (V->nconv < V->aup2_numcnv))
        {

            //  Max iterations have been exceeded.

            V->info = 1;
        }

        if ((V->np == 0) && (V->nconv < V->aup2_numcnv))
        {

            //  No shifts to apply.

            V->info = 2;
        }

        V->np = V->nconv;
        V->iter = V->aup2_iter;
        V->nev = V->aup2_numcnv;
        V->ido = ido_DONE;
        return;

    } else if ((V->nconv < V->aup2_numcnv) && (V->shift)) {

        //  Do not have all the requested eigenvalues yet.
        //  To prevent possible stagnation, adjust the size
        //  of NEV.

        int nevbef = V->nev;
        V->nev += (V->nconv > (V->np / 2) ? (V->np / 2) : V->nconv);
        if ((V->nev == 1) && (V->aup2_kplusp >= 6)) {
            V->nev = V->aup2_kplusp / 2;
        } else if ((V->nev == 1) && (V->aup2_kplusp > 3)) {
            V->nev = 2;
        }

        //  SciPy Fix
        //  We must keep nev below this value, as otherwise we can get
        //  np == 0 (note that dngets below can bump nev by 1). If np == 0,
        // the next call to `dnaitr` will write out-of-bounds.

        if (V->nev > (V->aup2_kplusp - 2)) {
            V->nev = V->aup2_kplusp - 2;
        }
        //  SciPy Fix End

        V->np = V->aup2_kplusp - V->nev;

        if (nevbef < V->nev) {
            dngets(V, &V->nev, &V->np, ritzr, ritzi, bounds);
        }

    }

    if (V->shift == 0)
    {

        //  User specified shifts: reverse communication to
        //  compute the shifts. They are returned in the first
        //  2*NP locations of WORKL.

        V->aup2_ushift = 1;
        V->ido = ido_USER_SHIFT;
        return;
    }

LINE50:

    //  Back from reverse communication;
    //  User specified shifts are returned
    //  in WORKL(1:2*NP)

    V->aup2_ushift = 0;

    if (V->shift == 0)
    {

        //  Move the NP shifts from WORKL to
        //  RITZR, RITZI to free up WORKL
        //  for non-exact shift case.

        dcopy_(&V->np, workl, &int1, ritzr, &int1);
        dcopy_(&V->np, &workl[V->np], &int1, ritzi, &int1);
    }

    //  Apply the NP implicit shifts by QR bulge chasing.
    //  Each shift is applied to the whole upper Hessenberg
    //  matrix H.
    //  The first 2*N locations of WORKD are used as workspace.

    dnapps(V->n, &V->nev, V->np, ritzr, ritzi, v, ldv, h, ldh, resid, q, ldq, workl, workd);

    //  Compute the B-norm of the updated residual.
    //  Keep B*RESID in WORKD(1:N) to be used in
    //  the first step of the next call to dnaitr .

    V->aup2_cnorm = 1;
    if (V->bmat)
    {
        dcopy_(&V->n, resid, &int1, &workd[V->n], &int1);
        ipntr[0] = V->n;
        ipntr[1] = 0;
        V->ido = ido_BX;

        //  Exit in order to compute B*RESID

        return;
    } else {
        dcopy_(&V->n, resid, &int1, workd, &int1);
    }

LINE100:

    //  Back from reverse communication;
    //  WORKD(1:N) := B*RESID

    if (V->bmat)
    {
        V->aup2_rnorm = ddot_(&V->n, resid, &int1, workd, &int1);
        V->aup2_rnorm = sqrt(fabs(V->aup2_rnorm));
    } else {
        V->aup2_rnorm = dnrm2_(&V->n, resid, &int1);
    }
    V->aup2_cnorm = 0;

    goto LINE1000;

    //
    //   E N D     O F     M A I N     I T E R A T I O N     L O O P
    //

}

void
dnconv(int n, double* ritzr, double* ritzi, double* bounds, const double tol, int* nconv)
{
    const double eps23 = pow(ulp, 2.0 / 3.0);
    double temp;

    *nconv = 0;
    for (int i = 0; i < n; i++)
    {
        temp = fmax(eps23, hypot(ritzr[i], ritzi[i]));
        if (bounds[i] <= tol*temp)
        {
            *nconv += 1;
        }
    }

    return;
}

void
dneigh(double* rnorm, int n, double* h, int ldh, double* ritzr, double* ritzi,
       double* bounds, double* q, int ldq, double* workl, int* ierr)
{
    int select[1] = { 0 };
    int i, iconj, int1 = 1, j;
    double dbl1 = 1.0, dbl0 = 0.0, temp, tmp_dbl, vl[1] = { 0.0 };

    //  1. Compute the eigenvalues, the last components of the
    //     corresponding Schur vectors and the full Schur form T
    //     of the current upper Hessenberg matrix H.
    //  dlahqr returns the full Schur form of H in WORKL(1:N**2)
    //  and the last components of the Schur vectors in BOUNDS.

    dlacpy_("A", &n, &n, h, &ldh, workl, &n);
    for (j = 0; j < n-1; j++)
    {
        bounds[j] = 0.0;
    }
    bounds[n-1] = 1.0;
    dlahqr_(&int1, &int1, &n, &int1, &n, workl, &n, ritzr, ritzi, &int1, &int1, bounds, &int1, ierr);

    if (*ierr != 0) { return; }

    //  2. Compute the eigenvectors of the full Schur form T and
    //     apply the last components of the Schur vectors to get
    //     the last components of the corresponding eigenvectors.
    //  Remember that if the i-th and (i+1)-st eigenvalues are
    //  complex conjugate pairs, then the real & imaginary part
    //  of the eigenvector components are split across adjacent
    //  columns of Q.

    dtrevc_("R", "A", select, &n, workl, &n, vl, &n, q, &ldq, &n, &n, &workl[n*n], ierr);
    if (*ierr != 0) { return; }

    //  Scale the returning eigenvectors so that their
    //  euclidean norms are all one. LAPACK subroutine
    //  dtrevc returns each eigenvector normalized so
    //  that the element of largest magnitude has
    //  magnitude 1; here the magnitude of a complex
    //  number (x,y) is taken to be |x| + |y|.

    iconj = 0;
    for (i = 0; i < n; i++)
    {
        if (fabs(ritzi[i]) == 0.0)
        {

            //  Real eigenvalue case

            temp = dnrm2_(&n, &q[ldq*i], &int1);
            tmp_dbl = 1.0 / temp;
            dscal_(&n, &tmp_dbl, &q[ldq*i], &int1);

        } else {

            //  Complex conjugate pair case. Note that
            //  since the real and imaginary part of
            //  the eigenvector are stored in consecutive
            //  columns, we further normalize by the
            //  square root of two.

            if (iconj == 0)
            {
                temp = hypot(dnrm2_(&n, &q[ldq*i], &int1),
                             dnrm2_(&n, &q[ldq*(i+1)], &int1));
                tmp_dbl = 1.0 / temp;
                dscal_(&n, &tmp_dbl, &q[ldq*i], &int1);
                dscal_(&n, &tmp_dbl, &q[ldq*(i+1)], &int1);
                iconj = 1;
            } else {
                iconj = 0;
            }
        }
    }
    // 10

    dgemv_("T", &n, &n, &dbl1, q, &ldq, bounds, &int1, &dbl0, workl, &int1);

    //  Compute the Ritz estimates

    iconj = 0;
    for (i = 0; i < n; i++)
    {
        if (fabs(ritzi[i]) == 0.0)
        {

            //  Real eigenvalue case

            bounds[i] = *rnorm * fabs(workl[i]);

        } else {

            //  Complex conjugate pair case. Note that
            //  since the real and imaginary part of
            //  the eigenvector are stored in consecutive
            //  columns, we need to take the magnitude
            //  of the last components of the two vectors

            if (iconj == 0)
            {
                bounds[i] = *rnorm * hypot(workl[i], workl[i+1]);
                bounds[i+1] = bounds[i];
                iconj = 1;
            } else {
                iconj = 0;
            }
        }
    }
    // 20

    return;
}

void
dnaitr(struct ARNAUD_state_d *V, int k, int np, double* resid, double* rnorm,
       double* v, int ldv, double* h, int ldh, int* ipntr, double* workd)
{
    int i = 0, infol, ipj, irj, ivj, jj, n, tmp_int;
    double smlnum = unfl * ( V->n / ulp);
    const double sq2o2 = sqrt(2.0) / 2.0;

    int int1 = 1;
    double dbl1 = 1.0, dbl0 = 0.0, dblm1 = -1.0, temp1, tst1;

    n = V->n;  // n is constant, this is just for typing convenience
    ipj = 0;
    irj = ipj + n;
    ivj = irj + n;

    if (V->ido == ido_FIRST)
    {

        //  Initial call to this routine

        V->aitr_j = k;
        V->info = 0;
        V->aitr_step3 = 0;
        V->aitr_step4 = 0;
        V->aitr_orth1 = 0;
        V->aitr_orth2 = 0;
        V->aitr_restart = 0;
    }

    //  When in reverse communication mode one of:
    //  STEP3, STEP4, ORTH1, ORTH2, RSTART
    //  will be .true. when ....
    //  STEP3: return from computing OP*v_{j}.
    //  STEP4: return from computing B-norm of OP*v_{j}
    //  ORTH1: return from computing B-norm of r_{j+1}
    //  ORTH2: return from computing B-norm of correction to the residual vector.
    //  RSTART: return from OP computations needed by dgetv0.

    if (V->aitr_step3) { goto LINE50; }
    if (V->aitr_step4) { goto LINE60; }
    if (V->aitr_orth1) { goto LINE70; }
    if (V->aitr_orth2) { goto LINE90; }
    if (V->aitr_restart) { goto LINE30; }

    //  Else this is the first step

    //
    //         A R N O L D I     I T E R A T I O N     L O O P
    //
    //  Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1)

LINE1000:

    //  STEP 1: Check if the B norm of j-th residual
    //  vector is zero. Equivalent to determining whether
    //  an exact j-step Arnoldi factorization is present.

    V->aitr_betaj = *rnorm;

    if (*rnorm > 0.0) { goto LINE40; }

    //  Invariant subspace found, generate a new starting
    //  vector which is orthogonal to the current Arnoldi
    //  basis and continue the iteration.

    V->aitr_betaj = 0.0;
    V->getv0_itry = 1;

LINE20:
    V->aitr_restart = 1;
    V->ido = ido_FIRST;

LINE30:

    // If in reverse communication mode and aitr_restart = 1, flow returns here.

    dgetv0(V, 0, n, V->aitr_j, v, ldv, resid, rnorm, ipntr, workd);

    if (V->ido != ido_DONE) { return; }
    V->aitr_ierr = V->info;
    if (V->aitr_ierr < 0)
    {
        V->getv0_itry += 1;
        if (V->getv0_itry <= 3) { goto LINE20; }

        //  Give up after several restart attempts.
        //  Set INFO to the size of the invariant subspace
        //  which spans OP and exit.

        V->info = V->aitr_j;
        V->ido = ido_DONE;
        return;
    }

LINE40:

    //  STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm
    //  Note that p_{j} = B*r_{j-1}. In order to avoid overflow
    //  when reciprocating a small RNORM, test against lower
    //  machine bound.

    dcopy_(&n, resid, &int1, &v[ldv*(V->aitr_j)], &int1);
    if (*rnorm >= unfl)
    {
        temp1 = 1.0 / *rnorm;
        dscal_(&n, &temp1, &v[ldv*(V->aitr_j)], &int1);
        dscal_(&n, &temp1, &workd[ipj], &int1);
    } else {
        dlascl_("G", &i, &i, rnorm, &dbl1, &n, &int1, &v[ldv*(V->aitr_j)], &n, &infol);
        dlascl_("G", &i, &i, rnorm, &dbl1, &n, &int1, &workd[ipj], &n, &infol);
    }

    //  STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j}
    //  Note that this is not quite yet r_{j}. See STEP 4

    V->aitr_step3 = 1;
    dcopy_(&n, &v[ldv*(V->aitr_j)], &int1, &workd[ivj], &int1);
    ipntr[0] = ivj;
    ipntr[1] = irj;
    ipntr[2] = ipj;
    V->ido = ido_OPX;

    //  Exit in order to compute OP*v_{j}

    return;

LINE50:

    //  Back from reverse communication;
    //  WORKD(IRJ:IRJ+N-1) := OP*v_{j}
    //  if step3 = .true.

    V->aitr_step3 = 0;

    //  Put another copy of OP*v_{j} into RESID.

    dcopy_(&n, &workd[irj], &int1, resid, &int1);

    //  STEP 4:  Finish extending the Arnoldi
    //           factorization to length j.

    if (V->bmat)
    {
        V->aitr_step4 = 1;
        ipntr[0] = irj;
        ipntr[1] = ipj;
        V->ido = ido_BX;

        //  Exit in order to compute B*OP*v_{j}

        return;
    } else {
        dcopy_(&n, resid, &int1, &workd[ipj], &int1);
    }

LINE60:

    //  Back from reverse communication;
    //  WORKD(IPJ:IPJ+N-1) := B*OP*v_{j}
    //  if step4 = .true.

    V->aitr_step4 = 0;

    //  The following is needed for STEP 5.
    //  Compute the B-norm of OP*v_{j}.

    if (V->bmat)
    {
        V->aitr_wnorm = ddot_(&n, resid, &int1, &workd[ipj], &int1);
        V->aitr_wnorm = sqrt(fabs(V->aitr_wnorm));
    } else {
        V->aitr_wnorm = dnrm2_(&n, resid, &int1);
    }

    //  Compute the j-th residual corresponding
    //  to the j step factorization.
    //  Use Classical Gram Schmidt and compute:
    //  w_{j} <-  V_{j}^T * B * OP * v_{j}
    //  r_{j} <-  OP*v_{j} - V_{j} * w_{j}

    //  Compute the j Fourier coefficients w_{j}
    //  WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.
    tmp_int = V->aitr_j + 1;
    dgemv_("T", &n, &tmp_int, &dbl1, v, &ldv, &workd[ipj], &int1, &dbl0, &h[ldh*(V->aitr_j)], &int1);

    //  Orthogonalize r_{j} against V_{j}.
    //  RESID contains OP*v_{j}. See STEP 3.

    dgemv_("N", &n, &tmp_int, &dblm1, v, &ldv, &h[ldh*(V->aitr_j)], &int1, &dbl1, resid, &int1);

    if (V->aitr_j > 0) { h[V->aitr_j + ldh*(V->aitr_j-1)] = V->aitr_betaj; }

    V->aitr_orth1 = 1;
    if (V->bmat)
    {
        dcopy_(&n, resid, &int1, &workd[irj], &int1);
        ipntr[0] = irj;
        ipntr[1] = ipj;
        V->ido = ido_BX;

        //  Exit in order to compute B*r_{j}

        return;
    } else {
        dcopy_(&n, resid, &int1, &workd[ipj], &int1);
    }

LINE70:

    //  Back from reverse communication if ORTH1 = .true.
    //  WORKD(IPJ:IPJ+N-1) := B*r_{j}.

    V->aitr_orth1 = 0;

    //  Compute the B-norm of r_{j}.

    if (V->bmat)
    {
        *rnorm = ddot_(&n, resid, &int1, &workd[ipj], &int1);
        *rnorm = sqrt(fabs(*rnorm));
    } else {
        *rnorm = dnrm2_(&n, resid, &int1);
    }

    //  STEP 5: Re-orthogonalization / Iterative refinement phase
    //  Maximum NITER_ITREF tries.
    //
    //           s      = V_{j}^T * B * r_{j}
    //           r_{j}  = r_{j} - V_{j}*s
    //           alphaj = alphaj + s_{j}
    //
    //  The stopping criteria used for iterative refinement is
    //  discussed in Parlett's book SEP, page 107 and in Gragg &
    //  Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.
    //  Determine if we need to correct the residual. The goal is
    //  to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||
    //  The following test determines whether the sine of the
    //  angle between  OP*x and the computed residual is less
    //  than or equal to 0.7071.

    if (*rnorm > sq2o2*V->aitr_wnorm) { goto LINE100; }
    V->aitr_iter = 0;

    //  Enter the Iterative refinement phase. If further
    //  refinement is necessary, loop back here. The loop
    //  variable is ITER. Perform a step of Classical
    //  Gram-Schmidt using all the Arnoldi vectors V_{j}

LINE80:

    //  Compute V_{j}^T * B * r_{j}.
    //  WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1).
    tmp_int = V->aitr_j + 1;
    dgemv_("T", &n, &tmp_int, &dbl1, v, &ldv, &workd[ipj], &int1, &dbl0, &workd[irj], &int1);

    //  Compute the correction to the residual:
    //  r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1).
    //  The correction to H is v(:,1:J)*H(1:J,1:J)
    //  + v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j.

    dgemv_("N", &n, &tmp_int, &dblm1, v, &ldv, &workd[irj], &int1, &dbl1, resid, &int1);
    daxpy_(&tmp_int, &dbl1, &workd[irj], &int1, &h[ldh*(V->aitr_j)], &int1);

    V->aitr_orth2 = 1;

    if (V->bmat)
    {
        dcopy_(&n, resid, &int1, &workd[irj], &int1);
        ipntr[0] = irj;
        ipntr[1] = ipj;
        V->ido = ido_BX;

        //  Exit in order to compute B*r_{j}.
        //  r_{j} is the corrected residual.

        return;
    } else {
        dcopy_(&n, resid, &int1, &workd[ipj], &int1);
    }

LINE90:

    //  Back from reverse communication if ORTH2 = .true.

    //  Compute the B-norm of the corrected residual r_{j}.

    if (V->bmat)
    {
        V->aitr_rnorm1 = ddot_(&n, resid, &int1, &workd[ipj], &int1);
        V->aitr_rnorm1 = sqrt(fabs(V->aitr_rnorm1));
    } else {
        V->aitr_rnorm1 = dnrm2_(&n, resid, &int1);
    }

    //  Determine if we need to perform another
    //  step of re-orthogonalization.

    if (V->aitr_rnorm1 > sq2o2*(*rnorm))
    {

        //  No need for further refinement.
        //  The cosine of the angle between the
        //  corrected residual vector and the old
        //  residual vector is greater than 0.717
        //  In other words the corrected residual
        //  and the old residual vector share an
        //  angle of less than arcCOS(0.717)

        *rnorm = V->aitr_rnorm1;

    } else {

        //  Another step of iterative refinement step
        //  is required.

        *rnorm = V->aitr_rnorm1;
        V->aitr_iter += 1;
        if (V->aitr_iter < 2) { goto LINE80; }

        //  Otherwise RESID is numerically in the span of V

        for (jj = 0; jj < n; jj++)
        {
            resid[jj] = 0.0;
        }
        *rnorm = 0.0;
    }

    // Branch here directly if iterative refinement
    // wasn't necessary or after at most NITER_REF
    // steps of iterative refinement.

LINE100:

    V->aitr_restart = 0;
    V->aitr_orth2 = 0;

    //  STEP 6: Update  j = j+1;  Continue

    V->aitr_j += 1;
    if (V->aitr_j >= k + np)
    {
        V->ido = ido_DONE;
        for (i = (k > 0 ? k-1 : k); i < k + np - 1; i++)
        {

            //  Check for splitting and deflation.
            //  Use a standard test as in the QR algorithm
            //  REFERENCE: LAPACK subroutine dlahqr

            tst1 = fabs(h[i + ldh*i]) + fabs(h[i+1 + ldh*(i+1)]);
            if (tst1 == 0.0)
            {
                tmp_int = k + np;
                tst1 = dlanhs_("1", &tmp_int, h, &ldh, &workd[n]);
            }
            if (fabs(h[i+1 + ldh*i]) <= fmax(ulp*tst1, smlnum))
            {
                h[i+1 + ldh*i] = 0.0;
            }
        }
        // 110
        return;
    }
    goto LINE1000;

}


void
dnapps(int n, int* kev, int np, double* shiftr, double* shifti, double* v,
       int ldv, double* h, int ldh, double* resid, double* q, int ldq, double* workl,
       double* workd)
{
    int cconj;
    int i, ir, j, jj, int1 = 1, istart, iend = 0, nr, tmp_int;
    int kplusp = *kev + np;
    double smlnum = unfl * ( n / ulp);
    double c, f, g, h11, h21, h12, h22, h32, s, sigmar, sigmai, r, t, tau, tst1;
    double dbl1 = 1.0, dbl0 = 0.0, dblm1 = -1.0;
    double u[3] = { 0.0 };

    //  Initialize Q to the identity to accumulate
    //  the rotations and reflections
    dlaset_("A", &kplusp, &kplusp, &dbl0, &dbl1, q, &ldq);

    //  Quick return if there are no shifts to apply

    if (np == 0) { return; }

    //  Chase the bulge with the application of each
    //  implicit shift. Each shift is applied to the
    //  whole matrix including each block.

    cconj = 0;

    // Loop over the shifts

    for (jj = 0; jj < np; jj++)
    {
        sigmar = shiftr[jj];
        sigmai = shifti[jj];

        if (cconj)
        {

            // Skip flag is on; turn it off and proceed to the next shift.

            cconj = 0;
            continue;

        } else if ((jj < np - 1) && fabs(sigmai) != 0.0) {

            // This shift has nonzero imaginary part, so we will apply
            // together with the next one; turn on the skip flag.

            cconj = 1;

        } else if ((jj == np - 1) && (fabs(sigmai) != 0.0)) {

            // We have one block left but the shift has nonzero imaginary part.
            // Don't apply it and reduce the number of shifts by incrementing
            // kev by one.

            *kev += 1;
            continue;
        }

        // if sigmai = 0 then
        //    Apply the jj-th shift ...
        // else
        //    Apply the jj-th and (jj+1)-th together ...
        //    (Note that jj < np at this point in the code)
        // end
        // to the current block of H

        istart = 0;
        while (istart < kplusp - 1)
        {
            for (iend = istart; iend < kplusp - 1; iend++)
            {
                tst1 = fabs(h[iend + (iend * ldh)]) + fabs(h[iend + 1 + (iend + 1) * ldh]);
                if (tst1 == 0.0)
                {
                    tmp_int = kplusp - jj;
                    tst1 = dlanhs_("1", &tmp_int, h, &ldh, workl);
                }
                if (fabs(h[iend+1 + (iend * ldh)]) <= fmax(smlnum, ulp * tst1))
                {
                    break;
                }
            }
            if (istart == iend)
            {
                istart += 1;
                continue;
            } else if  ((istart + 1 == iend) && fabs(sigmai) > 0.0) {
                istart += 2;
                continue;
            } else {
                h[iend+1 + (iend * ldh)] = 0.0;
            }

            // We have a block [istart, iend] inclusive.
            h11 = h[istart + istart * ldh];
            h21 = h[istart + 1 + istart * ldh];

            if (fabs(sigmai) == 0.0)
            {

                f = h11 - sigmar;
                g = h21;
                for (i = istart; i < iend; i++)
                {
                    dlartgp_(&f, &g, &c, &s, &r);
                    if (i > istart)
                    {
                        h[i + (i - 1) * ldh] = r;
                        h[i + 1 + (i - 1) * ldh] = 0.0;
                    }
                    tmp_int = kplusp - i;
                    drot_(&tmp_int, &h[i + ldh*i], &ldh, &h[i + 1 + ldh*i], &ldh, &c, &s);
                    tmp_int = (i+2 > iend ? iend : i + 2) + 1;
                    drot_(&tmp_int, &h[ldh*i], &int1, &h[ldh*(i+1)], &int1, &c, &s);
                    tmp_int = (i+jj+2 > kplusp ? kplusp : i + jj + 2);
                    drot_(&tmp_int, &q[ldq*i], &int1, &q[ldq*(i+1)], &int1, &c, &s);

                    if (i < iend - 1)
                    {
                        f = h[i+1 + i * ldh];
                        g = h[i+2 + i * ldh];
                    }
                }
            } else {

                h12 = h[istart + ldh*(istart + 1)];
                h22 = h[istart + 1 + ldh*(istart + 1)];
                h32 = h[istart + 2 + ldh*(istart + 1)];

                s = 2.0*sigmar;
                t = hypot(sigmar, sigmai);
                u[0] = (h11*(h11 - s) + t*t) / h21 + h12;
                u[1] = h11 + h22 - s;
                u[2] = h32;

                for (i = istart; i < iend; i++)
                {
                    nr = iend - i + 1;
                    nr = (nr > 3? 3 : nr);
                    dlarfg_(&nr, &u[0], &u[1], &int1, &tau);
                    if (i > istart)
                    {
                        h[i + (i - 1) * ldh] = u[0];
                        h[i + 1 + (i - 1) * ldh] = 0.0;
                        if (i < iend - 1) { h[i + 2 + (i - 1) * ldh] = 0.0; }
                    }
                    u[0] = 1.0;

                    tmp_int = kplusp - i;
                    dlarf_("L", &nr, &tmp_int, u, &int1, &tau, &h[i + ldh*i], &ldh, workl);
                    ir = (i + 3 > iend ? iend : i + 3) + 1;
                    dlarf_("R", &ir, &nr, u, &int1, &tau, &h[ldh*i], &ldh, workl);
                    dlarf_("R", &kplusp, &nr, u, &int1, &tau, &q[ldq*i], &ldq, workl);
                    if (i < iend - 1)
                    {
                        u[0] = h[i+1 + i * ldh];
                        u[1] = h[i+2 + i * ldh];
                        if (i < iend-2) { u[2] = h[i+3 + i * ldh]; }
                    }
                }
            }
            istart = iend + 1;
        }
    }
    //  Perform a similarity transformation that makes
    //  sure that H will have non negative sub diagonals

    for (j = 0; j < *kev; j++)
    {
        if (h[j+1 + ldh*j] < 0.0)
        {
            tmp_int = kplusp - j;
            dscal_(&tmp_int, &dblm1, &h[j+1 + ldh*j], &ldh);
            tmp_int = (j+3 > kplusp ? kplusp : j+3);
            dscal_(&tmp_int, &dblm1, &h[ldh*(j+1)], &int1);
            tmp_int = (j+np+2 > kplusp ? kplusp : j+np+2);
            dscal_(&tmp_int, &dblm1, &q[ldq*(j+1)], &int1);
        }
    }
    // 120

    for (i = 0; i < *kev; i++)
    {

        //  Final check for splitting and deflation.
        //  Use a standard test as in the QR algorithm
        //  REFERENCE: LAPACK subroutine dlahqr

        tst1 = fabs(h[i + ldh*i]) + fabs(h[i+1 + ldh*(i+1)]);
        if (tst1 == 0.0)
        {
            tst1 = dlanhs_("1", kev, h, &ldh, workl);
        }
        if (h[i+1 + ldh*i] <= fmax(ulp*tst1, smlnum))
        {
            h[i+1 + ldh*i] = 0.0;
        }
    }
    // 130

    //  Compute the (kev+1)-st column of (V*Q) and
    //  temporarily store the result in WORKD(N+1:2*N).
    //  This is needed in the residual update since we
    //  cannot GUARANTEE that the corresponding entry
    //  of H would be zero as in exact arithmetic.

    if (h[*kev + ldh*(*kev-1)] > 0.0)
    {
        dgemv_("N", &n, &kplusp, &dbl1, v, &ldv, &q[(*kev)*ldq], &int1, &dbl0, &workd[n], &int1);
    }

    //  Compute column 1 to kev of (V*Q) in backward order
    //  taking advantage of the upper Hessenberg structure of Q.

    for (i = 0; i < *kev; i++)
    {
        tmp_int = kplusp - i;
        dgemv_("N", &n, &tmp_int, &dbl1, v, &ldv, &q[(*kev-i-1)*ldq], &int1, &dbl0, workd, &int1);
        dcopy_(&n, workd, &int1, &v[(kplusp-i-1)*ldv], &int1);
    }

    //   Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev).

    for (i = 0; i < *kev; i++)
    {
        dcopy_(&n, &v[(kplusp-*kev+i)*ldv], &int1, &v[i*ldv], &int1);
    }

    //  Copy the (kev+1)-st column of (V*Q) in the appropriate place

    if (h[*kev + ldh*(*kev-1)] > 0.0){
        dcopy_(&n, &workd[n], &int1, &v[ldv*(*kev)], &int1);
    }

    //  Update the residual vector:
    //     r <- sigmak*r + betak*v(:,kev+1)
    //  where
    //     sigmak = (e_{kplusp}'*Q)*e_{kev}
    //     betak = e_{kev+1}'*H*e_{kev}

    dscal_(&n, &q[kplusp-1 + ldq*(*kev-1)], resid, &int1);

    if (h[*kev + ldh*(*kev-1)] > 0.0)
    {
        daxpy_(&n, &h[*kev + ldh*(*kev-1)], &v[ldv*(*kev)], &int1, resid, &int1);
    }

    return;

}


void
dngets(struct ARNAUD_state_d *V, int* kev, int* np,
       double* ritzr, double* ritzi, double* bounds)
{

    //  LM, SM, LR, SR, LI, SI case.
    //  Sort the eigenvalues of H into the desired order
    //  and apply the resulting order to BOUNDS.
    //  The eigenvalues are sorted so that the wanted part
    //  are always in the last KEV locations.
    //  We first do a pre-processing sort in order to keep
    //  complex conjugate pairs together

    switch (V->which)
    {
        case which_LM:
            dsortc(which_LR, 1, *kev + *np, ritzr, ritzi, bounds);
            break;
        case which_SM:
            dsortc(which_SR, 1, *kev + *np, ritzr, ritzi, bounds);
            break;
        case which_LR:
            dsortc(which_LM, 1, *kev + *np, ritzr, ritzi, bounds);
            break;
        case which_SR:
            dsortc(which_SM, 1, *kev + *np, ritzr, ritzi, bounds);
            break;
        case which_LI:
            dsortc(which_LM, 1, *kev + *np, ritzr, ritzi, bounds);
            break;
        case which_SI:
            dsortc(which_SM, 1, *kev + *np, ritzr, ritzi, bounds);
            break;
        default:
            dsortc(which_LR, 1, *kev + *np, ritzr, ritzi, bounds);
            break;
    }
    dsortc(V->which, 1, *kev + *np, ritzr, ritzi, bounds);

    //  Increase KEV by one if the ( ritzr(np),ritzi(np) )
    //  = ( ritzr(np+1),-ritzi(np+1) ) and ritz(np) .ne. zero
    //  Accordingly decrease NP by one. In other words keep
    //  complex conjugate pairs together.

    if ((ritzr[*np] - ritzr[*np-1] == 0.0) && (ritzi[*np] + ritzi[*np-1] == 0.0))
    {
        *np -= 1;
        *kev += 1;
    }

    if (V->shift == 1)
    {

        //  Sort the unwanted Ritz values used as shifts so that
        //  the ones with largest Ritz estimates are first
        //  This will tend to minimize the effects of the
        //  forward instability of the iteration when they shifts
        //  are applied in subroutine dnapps.
        //  Be careful and use 'SR' since we want to sort BOUNDS!

        dsortc(which_SR, 1, *np, bounds, ritzr, ritzi);
    }

    return;
}

void
dgetv0(struct ARNAUD_state_d *V, int initv, int n, int j,
       double* v, int ldv, double* resid, double* rnorm, int* ipntr, double* workd)
{
    int jj, int1 = 1;
    const double sq2o2 = sqrt(2.0) / 2.0;
    double dbl1 = 1.0, dbl0 = 0.0, dblm1 = -1.0;

    if (V->ido == ido_FIRST)
    {
        V->info = 0;
        V->getv0_iter = 0;
        V->getv0_first = 0;
        V->getv0_orth = 0;

        //  Possibly generate a random starting vector in RESID
        //  Skip if this the return of ido_RANDOM.

        if (!(initv))
        {
            // Request a random vector from the user into resid
            V->ido = ido_RANDOM;
            return;
        } else {
            V->ido = ido_RANDOM;
        }
    }

    // Back from random vector generation
    if (V->ido == ido_RANDOM)
    {
        //  Force the starting vector into the range of OP to handle
        //  the generalized problem when B is possibly (singular).

        if (V->getv0_itry == 1)
        {
            ipntr[0] = 0;
            ipntr[1] = n;
            dcopy_(&n, resid, &int1, workd, &int1);
            V->ido = ido_RANDOM_OPX;
            return;
        } else if ((V->getv0_itry > 1) && (V->bmat == 1))
        {
            dcopy_(&n, resid, &int1, &workd[n], &int1);
        }
    }

    //  Back from computing OP*(initial-vector)

    if (V->getv0_first) { goto LINE20; }

    //  Back from computing OP*(orthogonalized-vector)

    if (V->getv0_orth) { goto LINE40; }

    //  Starting vector is now in the range of OP; r = OP*r;
    //  Compute B-norm of starting vector.

    V->getv0_first = 1;
    if (V->getv0_itry == 1)
    {
        dcopy_(&n, &workd[n], &int1, resid, &int1);
    }
    if (V->bmat)
    {
        ipntr[0] = n;
        ipntr[1] = 0;
        V->ido = ido_BX;
        return;
    } else {
        dcopy_(&n, resid, &int1, workd, &int1);
    }

LINE20:

    V->getv0_first = 0;
    if (V->bmat)
    {
        V->getv0_rnorm0 = ddot_(&n, resid, &int1, workd, &int1);
        V->getv0_rnorm0 = sqrt(fabs(V->getv0_rnorm0));
    } else {
        V->getv0_rnorm0 = dnrm2_(&n, resid, &int1);
    }
    *rnorm = V->getv0_rnorm0;

    //  Exit if this is the very first Arnoldi step

    if (j == 0)
    {
        V->ido = ido_DONE;
        return;
    }

    //  Otherwise need to B-orthogonalize the starting vector against
    //  the current Arnoldi basis using Gram-Schmidt with iter. ref.
    //  This is the case where an invariant subspace is encountered
    //  in the middle of the Arnoldi factorization.
    //
    //        s = V^{T}*B*r;   r = r - V*s;
    //
    //  Stopping criteria used for iter. ref. is discussed in
    //  Parlett's book, page 107 and in Gragg & Reichel TOMS paper.

    V->getv0_orth = 1;

LINE30:

    dgemv_("T", &n, &j, &dbl1, v, &ldv, workd, &int1, &dbl0, &workd[n], &int1);
    dgemv_("N", &n, &j, &dblm1, v, &ldv, &workd[n], &int1, &dbl1, resid, &int1);

    //  Compute the B-norm of the orthogonalized starting vector

    if (V->bmat)
    {
        dcopy_(&n, resid, &int1, &workd[n], &int1);
        ipntr[0] = n;
        ipntr[1] = 0;
        V->ido = ido_BX;
        return;
    } else {
        dcopy_(&n, resid, &int1, workd, &int1);
    }

LINE40:
    if (V->bmat)
    {
        *rnorm = ddot_(&n, resid, &int1, workd, &int1);
        *rnorm = sqrt(fabs(*rnorm));
    } else {
        *rnorm = dnrm2_(&n, resid, &int1);
    }

    //  Check for further orthogonalization.

    if (*rnorm > sq2o2*V->getv0_rnorm0)
    {
        V->ido = ido_DONE;
        return;
    }

    V->getv0_iter += 1;
    if (V->getv0_iter < 5)
    {

        //  Perform iterative refinement step

        V->getv0_rnorm0 = *rnorm;
        goto LINE30;
    } else {

        //  Iterative refinement step "failed"

        for (jj = 0; jj < n; jj++) { resid[jj] = 0.0; }
        *rnorm = 0.0;
        V->info = -1;
    }

    V->ido = ido_DONE;

    return;
}


static void
dsortc(const enum ARNAUD_which w, const int apply, const int n, double* xreal, double* ximag, double* y)
{
    int i, gap, pos;
    double temp;
    ARNAUD_compare_cfunc *f;

    switch (w)
    {
        case which_LM:
            f = sortc_LM;
            break;
        case which_SM:
            f = sortc_SM;
            break;
        case which_LR:
            f = sortc_LR;
            break;
        case which_LI:
            f = sortc_LI;
            break;
        case which_SR:
            f = sortc_SR;
            break;
        case which_SI:
            f = sortc_SI;
            break;
        default:
            f = sortc_LM;
            break;
    }

    gap = n / 2;

    while (gap != 0)
    {
        for (i = gap; i < n; i++)
        {
            pos = i - gap;
            while ((pos >= 0) && (f(xreal[pos], ximag[pos], xreal[pos+gap], ximag[pos+gap])))
            {
                temp = xreal[pos];
                xreal[pos] = xreal[pos+gap];
                xreal[pos+gap] = temp;
                temp = ximag[pos];
                ximag[pos] = ximag[pos+gap];
                ximag[pos+gap] = temp;

                if (apply)
                {
                    temp = y[pos];
                    y[pos] = y[pos+gap];
                    y[pos+gap] = temp;
                }
                pos -= gap;
            }
        }
        gap = gap / 2;
    }
}


// The void casts are to avoid compiler warnings for unused parameters
int
sortc_LM(const double xre, const double xim, const double xreigap, const double ximigap)
{
    return (hypot(xre, xim) > hypot(xreigap, ximigap));
}

int
sortc_SM(const double xre, const double xim, const double xreigap, const double ximigap)
{
    return (hypot(xre, xim) < hypot(xreigap, ximigap));
}

int
sortc_LR(const double xre, const double xim, const double xreigap, const double ximigap)
{
    (void)xim; (void)ximigap;
    return (xre > xreigap);
}

int
sortc_SR(const double xre, const double xim, const double xreigap, const double ximigap)
{
    (void)xim; (void)ximigap;
    return (xre < xreigap);
}

int
sortc_LI(const double xre, const double xim, const double xreigap, const double ximigap)
{
    (void)xre; (void)xreigap;
    return (fabs(xim) > fabs(ximigap));
}

int
sortc_SI(const double xre, const double xim, const double xreigap, const double ximigap)
{
    (void)xre; (void)xreigap;
    return (fabs(xim) < fabs(ximigap));
}
