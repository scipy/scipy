#include "arnaud_s_double.h"
#include <float.h>

typedef int ARNAUD_compare_rfunc(const double, const double);

static int sortr_LM(const double, const double);
static int sortr_SM(const double, const double);
static int sortr_LA(const double, const double);
static int sortr_SA(const double, const double);

static const double unfl = DBL_MIN;    // 2.2250738585072014e-308;
static const double ulp = DBL_EPSILON; // 2.220446049250313e-16;

static void dsaup2(struct ARNAUD_state_d*, double*, double*, int, double*, int, double*, double*, double*, int, double*, int*, double*);
static void dsconv(int, double*, double*, double, int*);
static void dseigt(double, int, double*, int, double*, double*, double*, int*);
static void dsaitr(struct ARNAUD_state_d*, int, int, double*, double*, double*, int, double*, int, int*, double*);
static void dsapps(int, int*, int, double*, double*, int, double*, int, double*, double* , int, double*);
static void dsgets(struct ARNAUD_state_d*, int*, int*, double*, double*, double*);
static void dgetv0(struct ARNAUD_state_d *, int, int, int, double*, int, double*, double*, int*, double*);
static void dsortr(const enum ARNAUD_which w, const int apply, const int n, double* x1, double* x2);
static void dsesrt(const enum ARNAUD_which w, const int apply, const int n, double* x, int na, double* a, const int lda);
static void dstqrb(int n, double* d, double* e, double* z, double* work, int* info);

enum ARNAUD_seupd_type {
    REGULAR,
    SHIFTI,
    BUCKLE,
    CAYLEY
};


void
ARNAUD_dseupd(struct ARNAUD_state_d *V, int rvec, int howmny, int* select,
       double* d, double* z, int ldz, double sigma, double* resid, double* v,
       int ldv, int* ipntr, double* workd, double* workl)
{
    const double eps23 = pow(ulp, 2.0 / 3.0);
    int j, jj, k;
    int ibd, ih, ihb, ihd, iq, irz, iw, ldh, ldq, ritz, bounds, next, np;
    int ierr = 0, int1 = 1, tmp_int = 0, numcnv, reord;
    double bnorm2, rnorm, temp, temp1, dbl1 = 1.0;

    if (V->nconv == 0) { return; }

    ierr = 0;
    enum ARNAUD_seupd_type TYP = REGULAR;

    if (V->nconv <= 0) {
        ierr = -14;
    } else if (V->n <= 0) {
        ierr = -1;
    } else if (V->nev <= 0) {
        ierr = -2;
    } else if ((V->ncv <= V->nev) || (V->ncv > V->n)) {
        ierr = -3;
    } else if ((V->which != 0) &&
               (V->which != 1) &&
               (V->which != 6) &&
               (V->which != 7) &&
               (V->which != 8)) {
        ierr = -5;
    } else if ((V->bmat != 0) && (V->bmat != 1)) {
        ierr = -6;
    } else if ((rvec) && ((howmny < 0) || (howmny > 2))) {
        ierr = -15;
    } else if ((rvec) && (howmny == 2)) {
        ierr = -16;  // NotImplementedError
    }

    if ((V->mode == 1) || (V->mode == 2)) {
        TYP = REGULAR;
    } else if (V->mode == 3) {
        TYP = SHIFTI;
    } else if (V->mode == 4) {
        TYP = BUCKLE;
    } else if (V->mode == 5) {
        TYP = CAYLEY;
    } else {
        ierr = -10;
    }

    if ((V->mode == 1) && (V->bmat)) { ierr = -11; }
    if ((V->nev == 1) && (V->which == which_BE)) { ierr = -12; }

    if (ierr != 0) {
        V->info = ierr;
        return;
    }


    // Pointer into WORKL for address of H, RITZ, BOUNDS, Q
    // etc... and the remaining workspace.
    // Also update pointer to be used on output.
    // Memory is laid out as follows:
    // workl(1:2*ncv) := generated tridiagonal matrix H
    //       The subdiagonal is stored in workl(2:ncv).
    //       The dead spot is workl(1) but upon exiting
    //       dsaupd  stores the B-norm of the last residual
    //       vector in workl(1). We use this !!!
    // workl(2*ncv+1:2*ncv+ncv) := ritz values
    //       The wanted values are in the first NCONV spots.
    // workl(3*ncv+1:3*ncv+ncv) := computed Ritz estimates
    //       The wanted values are in the first NCONV spots.
    // NOTE: workl(1:4*ncv) is set by dsaupd  and is not
    //       modified by dseupd .
    // ----------------------------------------------------
    // The following is used and set by dseupd .
    // workl(4*ncv+1:4*ncv+ncv) := used as workspace during
    //       computation of the eigenvectors of H. Stores
    //       the diagonal of H. Upon EXIT contains the NCV
    //       Ritz values of the original system. The first
    //       NCONV spots have the wanted values. If MODE =
    //       1 or 2 then will equal workl(2*ncv+1:3*ncv).
    // workl(5*ncv+1:5*ncv+ncv) := used as workspace during
    //       computation of the eigenvectors of H. Stores
    //       the subdiagonal of H. Upon EXIT contains the
    //       NCV corresponding Ritz estimates of the
    //       original system. The first NCONV spots have the
    //       wanted values. If MODE = 1,2 then will equal
    //       workl(3*ncv+1:4*ncv).
    // workl(6*ncv+1:6*ncv+ncv*ncv) := orthogonal Q that is
    //       the eigenvector matrix for H as returned by
    //       dsteqr . Not referenced if RVEC = .False.
    //       Ordering follows that of workl(4*ncv+1:5*ncv)
    // workl(6*ncv+ncv*ncv+1:6*ncv+ncv*ncv+2*ncv) :=
    //       Workspace. Needed by dsteqr  and by dseupd .
    // GRAND total of NCV*(NCV+8) locations.

    ih = ipntr[4];
    ritz = ipntr[5];
    bounds = ipntr[6];
    ldh = V->ncv;
    ldq = V->ncv;
    ihd = bounds + ldh;
    ihb = ihd    + ldh;
    iq  = ihb    + ldh;
    iw  = iq     + ldh*V->ncv;
    next = iw    + 2*V->ncv;
    ipntr[3] = next;
    ipntr[7] = ihd;
    ipntr[8] = ihb;
    ipntr[9] = iq;

    // irz points to the Ritz values computed
    //     by _seigt before exiting _saup2.
    // ibd points to the Ritz estimates
    //     computed by _seigt before exiting
    //     _saup2.

    irz = ipntr[10] + V->ncv;
    ibd = irz + V->ncv;

    // RNORM is B-norm of the RESID(1:N).
    // BNORM2 is the 2 norm of B*RESID(1:N).
    // Upon exit of dsaupd  WORKD(1:N) has
    // B*RESID(1:N).

    rnorm = workl[ih];
    if (V->bmat)
    {
        bnorm2 = dnrm2_(&V->n, workd, &int1);
    } else {
        bnorm2 = rnorm;
    }

    if (rvec) {
        reord = 0;

        // Use the temporary bounds array to store indices
        // These will be used to mark the select array later

        for (j = 0; j < V->ncv; j++)
        {
            workl[bounds + j] = j;
            select[j] = 0;
        }
        // 10

        // Select the wanted Ritz values.
        // Sort the Ritz values so that the
        // wanted ones appear at the tailing
        // NEV positions of workl(irr) and
        // workl(iri).  Move the corresponding
        // error estimates in workl(bound)
        // accordingly.

        np = V->ncv - V->nev;
        V->shift = 0;
        dsgets(V, &V->nev, &np, &workl[irz], &workl[bounds], workl);

        //  Record indices of the converged wanted Ritz values
        //  Mark the select array for possible reordering

        numcnv = 0;
        for (j = 1; j <= V->ncv; j++)
        {
            temp1 = fmax(eps23, fabs(workl[irz + V->ncv - j]));

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
        //  the number (nconv) reported by saupd.  If these two
        //  are different then there has probably been an error
        //  caused by incorrect passing of the saupd data.

        if (numcnv != V->nconv)
        {
            V->info = -17;
            return;
        }

        // Call LAPACK routine _steqr to compute the eigenvalues and
        // eigenvectors of the final symmetric tridiagonal matrix H.
        // Initialize the eigenvector matrix Q to the identity.

        tmp_int = V->ncv - 1;
        dcopy_(&tmp_int, &workl[ih+1], &int1, &workl[ihb], &int1);
        dcopy_(&V->ncv, &workl[ih+ldh], &int1, &workl[ihd], &int1);

        dsteqr_("I", &V->ncv, &workl[ihd], &workl[ihb], &workl[iq], &ldq, &workl[iw], &ierr);

        if (ierr != 0)
        {
            V->info = -8;
            return;
        }

        if (reord)
        {

            // Reordered the eigenvalues and eigenvectors
            // computed by _steqr so that the "converged"
            // eigenvalues appear in the first NCONV
            // positions of workl(ihd), and the associated
            // eigenvectors appear in the first NCONV
            // columns.

            int leftptr = 0;
            int rightptr = V->ncv - 1;

            if (V->ncv > 1)
            {
                do
                {
                    if (select[leftptr])
                    {

                        // Search, from the left, for the first non-converged Ritz value.

                        leftptr += 1;

                    } else if (!(select[rightptr])) {

                        // Search, from the right, the first converged Ritz value

                        rightptr -= 1;

                    } else {

                        // Swap the Ritz value on the left that has not
                        // converged with the Ritz value on the right
                        // that has converged.  Swap the associated
                        // eigenvector of the tridiagonal matrix H as
                        // well.

                        temp = workl[ihd + leftptr];
                        workl[ihd + leftptr] = workl[ihd + rightptr];
                        workl[ihd + rightptr] = temp;

                        dcopy_(&V->ncv, &workl[iq + V->ncv*leftptr], &int1, &workl[iw], &int1);
                        dcopy_(&V->ncv, &workl[iq + V->ncv*rightptr], &int1, &workl[iq + V->ncv*leftptr], &int1);
                        dcopy_(&V->ncv, &workl[iw], &int1, &workl[iq + V->ncv*rightptr], &int1);

                        leftptr += 1;
                        rightptr -= 1;
                    }
                } while (leftptr < rightptr);
            }
        }

        // Load the converged Ritz values into D.

        dcopy_(&V->nconv, &workl[ihd], &int1, d, &int1);

    } else {

        // Ritz vectors not required. Load Ritz values into D.

        dcopy_(&V->nconv, &workl[ritz], &int1, d, &int1);
        dcopy_(&V->ncv, &workl[ritz], &int1, &workl[ihd], &int1);
    }

    // Transform the Ritz values and possibly vectors and corresponding
    // Ritz estimates of OP to those of A*x=lambda*B*x. The Ritz values
    // (and corresponding data) are returned in ascending order.

    if (TYP == REGULAR)
    {

        // Ascending sort of wanted Ritz values, vectors and error
        // bounds. Not necessary if only Ritz values are desired.

        if (rvec) {
            dsesrt(which_LA, rvec, V->nconv, d, V->ncv, &workl[iq], ldq);
        } else {
            dcopy_(&V->ncv, &workl[bounds], &int1, &workl[ihb], &int1);
        }

    } else {

        // *  Make a copy of all the Ritz values.
        // *  Transform the Ritz values back to the original system.
        //    For TYPE = 'SHIFTI' the transformation is
        //             lambda = 1/theta + sigma
        //    For TYPE = 'BUCKLE' the transformation is
        //             lambda = sigma * theta / ( theta - 1 )
        //    For TYPE = 'CAYLEY' the transformation is
        //             lambda = sigma * (theta + 1) / (theta - 1 )
        //    where the theta are the Ritz values returned by dsaupd.
        // NOTES:
        // *The Ritz vectors are not affected by the transformation.
        //  They are only reordered.

        dcopy_(&V->ncv, &workl[ihd], &int1, &workl[iw], &int1);
        if (TYP == SHIFTI)
        {
            for (int k = 0; k < V->ncv; k++)
            {
                workl[ihd + k] = 1.0 / workl[ihd + k] + sigma;
            }
        } else if (TYP == BUCKLE) {
            for (int k = 0; k < V->ncv; k++) {
                workl[ihd + k] = sigma * workl[ihd + k] / (workl[ihd + k] - 1.0);
            }
        } else if (TYP == CAYLEY) {
            for (int k = 0; k < V->ncv; k++) {
                workl[ihd + k] = sigma * (workl[ihd + k] + 1.0) / (workl[ihd + k] - 1.0);
            }
        }

        // *  Store the wanted NCONV lambda values into D.
        // *  Sort the NCONV wanted lambda in WORKL(IHD:IHD+NCONV-1)
        //    into ascending order and apply sort to the NCONV theta
        //    values in the transformed system. We will need this to
        //    compute Ritz estimates in the original system.
        // *  Finally sort the lambda`s into ascending order and apply
        //    to Ritz vectors if wanted. Else just sort lambda`s into
        //    ascending order.
        // NOTES:
        // *workl(iw:iw+ncv-1) contain the theta ordered so that they
        //  match the ordering of the lambda. We`ll use them again for
        //  Ritz vector purification.

        dcopy_(&V->nconv, &workl[ihd], &int1, d, &int1);
        dsortr(which_LA, 1, V->nconv, &workl[ihd], &workl[iw]);
        if (rvec) {
            dsesrt(which_LA, rvec, V->nconv, d, V->ncv, &workl[iq], ldq);
        } else {
            dcopy_(&V->ncv, &workl[bounds], &int1, &workl[ihb], &int1);
            temp = bnorm2 / rnorm;
            dscal_(&V->ncv, &temp, &workl[ihb], &int1);
            dsortr(which_LA, 1, V->nconv, d, &workl[ihb]);
        }
    }

    // Compute the Ritz vectors. Transform the wanted
    // eigenvectors of the symmetric tridiagonal H by
    // the Lanczos basis matrix V.

    if ((rvec) && (howmny == 0))
    {

        // Compute the QR factorization of the matrix representing
        // the wanted invariant subspace located in the first NCONV
        // columns of workl(iq, ldq).

        dgeqr2_(&V->ncv, &V->nconv, &workl[iq], &ldq, &workl[iw + V->ncv], &workl[ihb], &ierr);

        // * Postmultiply V by Q.
        // * Copy the first NCONV columns of VQ into Z.
        // The N by NCONV matrix Z is now a matrix representation
        // of the approximate invariant subspace associated with
        // the Ritz values in workl(ihd).

        dorm2r_("R", "N", &V->n, &V->ncv, &V->nconv, &workl[iq], &ldq, &workl[iw + V->ncv], v, &ldv, &workd[V->n], &ierr);
        dlacpy_("A", &V->n, &V->nconv, v, &ldv, z, &ldz);

        // In order to compute the Ritz estimates for the Ritz
        // values in both systems, need the last row of the
        // eigenvector matrix. Remember, it's in factored form.

        for (int j = 0; j < V->ncv - 1; j++)
        {
            workl[ihb + j] = 0.0;
        }
        workl[ihb + V->ncv - 1] = 1.0;
        dorm2r_("L", "T", &V->ncv, &int1, &V->nconv, &workl[iq], &ldq, &workl[iw + V->ncv], &workl[ihb], &V->ncv, &temp, &ierr);

        //  Make a copy of the last row into
        //  workl(iw+ncv:iw+2*ncv), as it is needed again in
        //  the Ritz vector purification step below

        for (int j = 0; j < V->nconv; j++)
        {
            workl[iw + V->ncv + j] = workl[ihb + j];
        }
        // 67

    } else if ((rvec) && (howmny == 2)) {
        // Not yet implemented
        ;
    }

    if ((TYP == REGULAR) && (rvec))
    {
        for (j = 0; j < V->ncv; j++)
        {
            workl[ihb + j] = rnorm * fabs(workl[ihb + j]);
        }

    } else if ((TYP != REGULAR) && (rvec)) {

        // *  Determine Ritz estimates of the theta.
        //    If RVEC = .true. then compute Ritz estimates
        //               of the theta.
        //    If RVEC = .false. then copy Ritz estimates
        //              as computed by dsaupd .
        // *  Determine Ritz estimates of the lambda.


        dscal_(&V->ncv, &bnorm2, &workl[ihb], &int1);

        for (k = 0; k < V->ncv; k++)
        {
            if (TYP == SHIFTI)
            {
                workl[ihb + k] = fabs(workl[ihb + k]) / pow(workl[iw + k], 2.0);
            } else if (TYP == BUCKLE) {
                workl[ihb + k] = sigma * fabs(workl[ihb + k]) / pow(workl[iw + k] - 1.0, 2.0);
            } else if (TYP == CAYLEY) {
                workl[ihb + k] = fabs(workl[ihb + k] / workl[iw + k] * (workl[iw + k] - 1.0));
            }
        }
        // 80, 90, 100
    }

    // Ritz vector purification step. Formally perform
    // one of inverse subspace iteration. Only used
    // for MODE = 3,4,5. See reference 7

    if ((rvec) && ((TYP == SHIFTI) || (TYP == CAYLEY)))
    {
        for (k = 0; k < V->nconv; k++)
        {
            workl[iw + k] = workl[iw + V->ncv + k] / workl[iw + k];
        }
        // 110
    } else if ((rvec) && (TYP == BUCKLE))
    {
        for (k = 0; k < V->nconv; k++)
        {
            workl[iw + k] = workl[iw + V->ncv + k] / (workl[iw + k] - 1.0);
        }
        // 120
    }

    if ((rvec) && (TYP != REGULAR))
    {
        dger_(&V->n, &V->nconv, &dbl1, resid, &int1, &workl[iw], &int1, z, &ldz);
    }

    return;
}


void
ARNAUD_dsaupd(struct ARNAUD_state_d *V, double* resid, double* v, int ldv,
       int* ipntr, double* workd, double* workl)
{

    int bounds = 0, ih = 0, iq = 0, iw = 0, j = 0, ldh = 0, ldq = 0, next = 0, ritz = 0;

    if (V->ido == ido_FIRST)
    {
        if (V->n <= 0) {
            V->info = -1;
        } else if (V->nev <= 0) {
            V->info = -2;
        } else if ((V->ncv < V->nev + 1) || (V->ncv > V->n)) {
            V->info = -3;
        } else if (V->maxiter <= 0) {
            V->info = -4;
        } else if ((V->which != 0) &&
            (V->which != 1) &&
            (V->which != 6) &&
            (V->which != 7) &&
            (V->which != 8)) {
            V->info = -5;
        } else if ((V->bmat != 0) && (V->bmat != 1)) {
            V->info = -6;
        } else if ((V->mode < 1) || (V->mode > 5)) {
            V->info = -10;
        } else if ((V->mode == 1) && (V->bmat == 1)) {
            V->info = -11;
        } else if ((V->shift != 0) && (V->shift != 1)) {
            V->info = -12;
        } else if ((V->nev == 1) && (V->which == which_BE)) {
            V->info = -13;
        }

        if (V->info < 0)
        {
            V->ido = ido_DONE;
            return;
        }

        if (V->tol <= 0.0) { V->tol = ulp; }

        // NP is the number of additional steps to
        // extend the length NEV Lanczos factorization.
        V->np = V->ncv - V->nev;

        // Zero out internal workspace
        for (j = 0; j < (V->ncv)*((V->ncv) + 8); j++) { workl[j] = 0.0; }
    }

    // Pointer into WORKL for address of H, RITZ, BOUNDS, Q
    // etc... and the remaining workspace.
    // Also update pointer to be used on output.
    // Memory is laid out as follows:
    // workl(1:2*ncv) := generated tridiagonal matrix
    // workl(2*ncv+1:2*ncv+ncv) := ritz values
    // workl(3*ncv+1:3*ncv+ncv) := computed error bounds
    // workl(4*ncv+1:4*ncv+ncv*ncv) := rotation matrix Q
    // workl(4*ncv+ncv*ncv+1:7*ncv+ncv*ncv) := workspace

    ldh = V->ncv;
    ldq = V->ncv;
    ih = 0;
    ritz = ih + 2*ldh;
    bounds = ritz + V->ncv;
    iq = bounds + V->ncv;
    iw = iq + (V->ncv)*(V->ncv);
    next = iw + 3*(V->ncv);

    ipntr[3] = next;
    ipntr[4] = ih;
    ipntr[5] = ritz;
    ipntr[6] = bounds;
    ipntr[10] = iw;

    // Carry out the Implicitly restarted Lanczos Iteration.
    dsaup2(V, resid, v, ldv, &workl[ih], ldh, &workl[ritz], &workl[bounds],
           &workl[iq], ldq, &workl[iw], ipntr, workd);

     /*-------------------------------------------------*
     | ido .ne. 99 implies use of reverse communication |
     | to compute operations involving OP or shifts.    |
     *-------------------------------------------------*/

    if (V->ido != ido_DONE) { return; }

    if (V->info < 0) { return; }
    if (V->info == 2) { V->info = 3; }

    return;
}


void
dsaup2(struct ARNAUD_state_d *V, double* resid, double* v, int ldv,
       double* h, int ldh, double* ritz, double* bounds,
       double* q, int ldq, double* workl, int* ipntr, double* workd)
{
    int int1 = 1, j, tmp_int, nevd2, nevm2;
    const double eps23 = pow(ulp, 2.0 / 3.0);
    double temp = 0.0;
    // Initialize to silence the compiler warning
    enum ARNAUD_which temp_which = which_LM;

    if (V->ido == ido_FIRST)
    {
        // nev0 and np0 are integer variables hold the initial values of NEV & NP
        V->aup2_nev0 = V->nev;
        V->aup2_np0 = V->np;

        // kplusp is the bound on the largest Lanczos factorization built.
        // nconv is the current number of "converged" eigenvalues.
        // iter is the counter on the current iteration step.
        V->aup2_kplusp = V->nev + V->np;
        V->nconv = 0;
        V->aup2_iter = 0;

        // Set flags for computing the first NEV steps of the Lanczos factorization.               |

        V->aup2_getv0 = 1;
        V->aup2_update = 0;
        V->aup2_cnorm = 0;
        V->aup2_ushift = 0;

        if (V->info != 0)
        {
            // User provides the initial residual vector.
            V->aup2_initv = 1;
            V->info = 0;
        } else {
            V->aup2_initv = 0;
        }
    }

    // Get a possibly random starting vector and
    // force it into the range of the operator OP.

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

    // Back from reverse communication :
    // continue with update step

    if (V->aup2_update) { goto LINE20; }

    //  Back from computing user specified shifts

    if (V->aup2_ushift) { goto LINE50; }

    //  Back from computing residual norm
    //  at the end of the current iteration

    if (V->aup2_cnorm) { goto LINE100; }

    //  Compute the first NEV steps of the Arnoldi factorization

    dsaitr(V, 0, V->aup2_nev0, resid, &V->aup2_rnorm, v, ldv, h, ldh, ipntr, workd);

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

    // Compute NP additional steps of the Lanczos factorization.

    V->ido = ido_FIRST;

LINE20:
    V->aup2_update = 1;

    dsaitr(V, V->nev, V->np, resid, &V->aup2_rnorm, v, ldv, h, ldh, ipntr, workd);

     /*--------------------------------------------------*
     | ido .ne. 99 implies use of reverse communication  |
     | to compute operations involving OP and possibly B |
     *--------------------------------------------------*/

    if (V->ido != ido_DONE) { return; }

    if (V->info > 0)
    {
        V->np = V->info;
        V->iter = V->aup2_iter;
        V->info = -9999;
        V->ido = ido_DONE;
        return;
    }

    V->aup2_update = 0;

    // Compute the eigenvalues and corresponding error bounds
    // of the current symmetric tridiagonal matrix.
    dseigt(V->aup2_rnorm, V->aup2_kplusp, h, ldh, ritz, bounds, workl, &V->info);

    if (V->info != 0)
    {
        V->info = -8;
        V->ido = ido_DONE;
        return;
    }

    // Make a copy of eigenvalues and corresponding error
    // bounds obtained from _seigt.

    dcopy_(&V->aup2_kplusp, ritz, &int1, &workl[V->aup2_kplusp], &int1);
    dcopy_(&V->aup2_kplusp, bounds, &int1, &workl[2*V->aup2_kplusp], &int1);

    // Select the wanted Ritz values and their bounds
    // to be used in the convergence test.
    // The selection is based on the requested number of
    // eigenvalues instead of the current NEV and NP to
    // prevent possible misconvergence.
    // * Wanted Ritz values := RITZ(NP+1:NEV+NP)
    // * Shifts := RITZ(1:NP) := WORKL(1:NP)

    V->nev = V->aup2_nev0;
    V->np = V->aup2_np0;

    dsgets(V, &V->nev, &V->np, ritz, bounds, workl);

    // Convergence test

    dcopy_(&V->nev, &bounds[V->np], &int1, &workl[V->np], &int1);
    dsconv(V->nev, &ritz[V->np], &workl[V->np], V->tol, &V->nconv);

    // Count the number of unwanted Ritz values that have zero
    // Ritz estimates. If any Ritz estimates are equal to zero
    // then a leading block of H of order equal to at least
    // the number of Ritz values with zero Ritz estimates has
    // split off. None of these Ritz values may be removed by
    // shifting. Decrease NP the number of shifts to apply. If
    // no shifts may be applied, then prepare to exit

    tmp_int = V->np;
    for (j = 0; j < tmp_int; j++)
    {
        if (bounds[j] == 0.0)
        {
            V->np -= 1;
            V->nev += 1;
        }
    }
    // 30

    if ((V->nconv >= V->aup2_nev0) || (V->aup2_iter > V->maxiter) || (V->np == 0))
    {
        // Prepare to exit. Put the converged Ritz values
        // and corresponding bounds in RITZ(1:NCONV) and
        // BOUNDS(1:NCONV) respectively. Then sort. Be
        // careful when NCONV > NP since we don't want to
        // swap overlapping locations.

        if (V->which == which_BE)
        {

            // Both ends of the spectrum are requested.
            // Sort the eigenvalues into algebraically decreasing
            // order first then swap low end of the spectrum next
            // to high end in appropriate locations.
            // NOTE: when np < floor(nev/2) be careful not to swap
            // overlapping locations.

            dsortr(which_SA, 1, V->aup2_kplusp, ritz, bounds);
            nevd2 = V->aup2_nev0 / 2;
            nevm2 = V->aup2_nev0 - nevd2;
            if (V->nev > 1)
            {
                V->np = V->aup2_kplusp - V->aup2_nev0;

                tmp_int = (nevd2 < V->np ? nevd2 : V->np);
                int tmp_int2 = V->aup2_kplusp - tmp_int;

                dswap_(&tmp_int, &ritz[nevm2], &int1, &ritz[tmp_int2], &int1);
                dswap_(&tmp_int, &bounds[nevm2], &int1, &bounds[tmp_int2], &int1);
            }

        } else {

            // LM, SM, LA, SA case.
            // Sort the eigenvalues of H into the an order that
            // is opposite to WHICH, and apply the resulting
            // order to BOUNDS.  The eigenvalues are sorted so
            // that the wanted part are always within the first
            // NEV locations.

            if (V->which == which_LM) { temp_which = which_SM; }
            if (V->which == which_SM) { temp_which = which_LM; }
            if (V->which == which_LA) { temp_which = which_SA; }
            if (V->which == which_SA) { temp_which = which_LA; }

            dsortr(temp_which, 1, V->aup2_kplusp, ritz, bounds);
        }

        // Scale the Ritz estimate of each Ritz value
        // by 1 / max(eps23,magnitude of the Ritz value).

        for (j = 0; j < V->aup2_nev0; j++)
        {
            temp = fmax(eps23, fabs(ritz[j]));
            bounds[j] = bounds[j] / temp;
        }
        // 35

        // Sort the Ritz values according to the scaled Ritz
        // estimates.  This will push all the converged ones
        // towards the front of ritzr, ritzi, bounds
        // (in the case when NCONV < NEV.)

        dsortr(which_LA, 1, V->aup2_nev0, bounds, ritz);

        // Scale the Ritz estimate back to its original
        // value.

        for (j = 0; j < V->aup2_nev0; j++)
        {
            temp = fmax(eps23, fabs(ritz[j]));
            bounds[j] = bounds[j] * temp;
        }
        // 40

        // Sort the "converged" Ritz values again so that
        // the "threshold" values and their associated Ritz
        // estimates appear at the appropriate position in
        // ritz and bound.

        if (V->which == which_BE)
        {

            // Sort the "converged" Ritz values in increasing
            // order.  The "threshold" values are in the
            // middle.

            dsortr(which_LA, 1, V->nconv, ritz, bounds);
        } else {

            // In LM, SM, LA, SA case, sort the "converged"
            // Ritz values according to WHICH so that the
            // "threshold" value appears at the front of
            // ritz.

            dsortr(V->which, 1, V->nconv, ritz, bounds);
        }

        //  Use h( 1,1 ) as storage to communicate
        //  rnorm to _seupd if needed

        h[0] = V->aup2_rnorm;

        // Max iterations have been exceeded.

        if ((V->aup2_iter > V->maxiter) && (V->nconv < V->nev))
        {
            V->info = 1;
        }

        // No shifts to apply.

        if ((V->np ==  0) && (V->nconv < V->aup2_nev0))
        {
            V->info = 2;
        }

        V->np = V->nconv;
        V->nev = V->nconv;
        V->iter = V->aup2_iter;
        V->ido = ido_DONE;
        return;

    } else if ((V->nconv < V->nev) && (V->shift == 1)) {

        // Do not have all the requested eigenvalues yet.
        // To prevent possible stagnation, adjust the number
        // of Ritz values and the shifts.

        int nevbef = V->nev;
        V->nev += (V->nconv > (V->np / 2) ? (V->np / 2) : V->nconv);
        if ((V->nev == 1) && (V->aup2_kplusp >= 6))
        {
            V->nev = V->aup2_kplusp / 2;
        } else if ((V->nev == 1) && (V->aup2_kplusp > 2))
        {
            V->nev = 2;
        }
        V->np = V->aup2_kplusp - V->nev;

        // If the size of NEV was just increased resort the eigenvalues.

        if (nevbef < V->nev)
        {
            dsgets(V, &V->nev, &V->np, ritz, bounds, workl);
        }
    }

    if (V->shift == 0)
    {
        // User specified shifts: reverse communication to
        // compute the shifts. They are returned in the first
        // NP locations of WORKL.

        V->aup2_ushift = 1;
        V->ido = ido_USER_SHIFT;
        return;
    }

LINE50:

    // Back from reverse communication;
    // User specified shifts are returned
    // in WORKL(1:*NP)

    V->aup2_ushift = 0;

    // Move the NP shifts to the first NP locations of RITZ to
    // free up WORKL.  This is for the non-exact shift case;
    // in the exact shift case, dsgets already handles this.

    if (V->shift == 0) { dcopy_(&V->np, workl, &int1, ritz, &int1); }

     /*--------------------------------------------------------*
     | Apply the NP0 implicit shifts by QR bulge chasing.      |
     | Each shift is applied to the entire tridiagonal matrix. |
     | The first 2*N locations of WORKD are used as workspace. |
     | After dsapps is done, we have a Lanczos                 |
     | factorization of length NEV.                            |
     *--------------------------------------------------------*/

    dsapps(V->n, &V->nev, V->np, ritz, v, ldv, h, ldh, resid, q, ldq, workd);

    // Compute the B-norm of the updated residual.
    // Keep B*RESID in WORKD(1:N) to be used in
    // the first step of the next call to dsaitr.

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

     /*---------------------------------*
     | Back from reverse communication; |
     | WORKD(1:N) := B*RESID            |
     *---------------------------------*/

    if (V->bmat)
    {
        V->aup2_rnorm = sqrt(fabs(ddot_(&V->n, resid, &int1, workd, &int1)));
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
dsconv(int n, double *ritz, double *bounds, double tol, int* nconv)
{
    const double eps23 = pow(ulp, 2.0 / 3.0);
    double temp = 0.0;

    *nconv = 0;
    // Convergence test
    for (int i = 0; i < n; i++)
    {
        temp = fmax(eps23, fabs(ritz[i]));
        if (fabs(bounds[i]) <= tol * temp) { *nconv += 1; }
    }

    return;
}


void
dseigt(double rnorm, int n, double* h, int ldh, double* eig, double* bounds,
       double* workl, int* ierr)
{
    int int1 = 1, tmp_int;
    dcopy_(&n, &h[ldh], &int1, eig, &int1);
    tmp_int = n - 1;
    dcopy_(&tmp_int, &h[1], &int1, workl, &int1);
    dstqrb(n, eig, workl, bounds, &workl[n], ierr);
    if (*ierr != 0) { return; }
    for (int k = 0; k < n; k++) { bounds[k] = rnorm * fabs(bounds[k]); }
    return;
}


void
dsaitr(struct ARNAUD_state_d *V, int k, int np, double* resid, double* rnorm,
       double* v, int ldv, double* h, int ldh, int* ipntr, double* workd)
{
    int i, infol, ipj, irj, ivj, jj, n, tmp_int;
    const double sq2o2 = sqrt(2.0) / 2.0;

    int int1 = 1;
    double dbl1 = 1.0, dbl0 = 0.0, dblm1 = -1.0, temp1;

    n = V->n;  // n is constant, this is just for typing convenience
    ipj = 0;
    irj = ipj + n;
    ivj = irj + n;

    if (V->ido == ido_FIRST)
    {

        // Initial call to this routine

        V->aitr_j = k;
        V->info = 0;
        V->aitr_step3 = 0;
        V->aitr_step4 = 0;
        V->aitr_orth1 = 0;
        V->aitr_orth2 = 0;
        V->aitr_restart = 0;
    }

    // When in reverse communication mode one of:
    // STEP3, STEP4, ORTH1, ORTH2, RSTART
    // will be .true. when ....
    // STEP3: return from computing OP*v_{j}.
    // STEP4: return from computing B-norm of OP*v_{j}
    // ORTH1: return from computing B-norm of r_{j+1}
    // ORTH2: return from computing B-norm of
    //        correction to the residual vector.
    // RSTART: return from OP computations needed by
    //         dgetv0.

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

        // Check for exact zero. Equivalent to determining whether
        // a j-step Arnoldi factorization is present.

    if (*rnorm > 0.0) { goto LINE40; }

        // Invariant subspace found, generate a new starting
        // vector which is orthogonal to the current Arnoldi
        // basis and continue the iteration.

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

    V->aitr_step3 = 0;

    // Put another copy of OP*v_{j} into RESID.

    dcopy_(&n, &workd[irj], &int1, resid, &int1);

    // STEP 4:  Finish extending the symmetric
    //          Arnoldi to length j. If MODE = 2
    //          then B*OP = B*inv(B)*A = A and
    //          we don't need to compute B*OP.
    // NOTE: If MODE = 2 WORKD(IVJ:IVJ+N-1) is
    // assumed to have A*v_{j}.

    if (V->mode == 2) { goto LINE65; }
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

    // Back from reverse communication;
    // WORKD(IPJ:IPJ+N-1) := B*OP*v_{j}.

    V->aitr_step4 = 0;

    //  The following is needed for STEP 5.
    //  Compute the B-norm of OP*v_{j}.

LINE65:

    if (V->mode == 2)
    {

        // Note that the B-norm of OP*v_{j}
        // is the inv(B)-norm of A*v_{j}.

        V->aitr_wnorm = ddot_(&n, resid, &int1, &workd[ivj], &int1);
        V->aitr_wnorm = sqrt(fabs(V->aitr_wnorm));
    } else if (V->bmat) {
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
    if (V->mode != 2)
    {
        dgemv_("T", &n, &tmp_int, &dbl1, v, &ldv, &workd[ipj], &int1, &dbl0, &workd[irj], &int1);
    } else {
        dgemv_("T", &n, &tmp_int, &dbl1, v, &ldv, &workd[ivj], &int1, &dbl0, &workd[irj], &int1);
    }

    //  Orthogonalize r_{j} against V_{j}.
    //  RESID contains OP*v_{j}. See STEP 3.

    dgemv_("N", &n, &tmp_int, &dblm1, v, &ldv, &workd[irj], &int1, &dbl1, resid, &int1);

    // Extend H to have j rows and columns.

    h[V->aitr_j + ldh] = workd[irj + V->aitr_j];

    if ((V->aitr_j == 0) || (V->aitr_restart))
    {
        h[V->aitr_j] = 0.0;
    } else {
        h[V->aitr_j] = *rnorm;
    }

    V->aitr_orth1 = 1;
    V->aitr_iter = 0;

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

    // STEP 5: Re-orthogonalization / Iterative refinement phase
    // Maximum NITER_ITREF tries.
    //
    //          s      = V_{j}^T * B * r_{j}
    //          r_{j}  = r_{j} - V_{j}*s
    //          alphaj = alphaj + s_{j}
    //
    // The stopping criteria used for iterative refinement is
    // discussed in Parlett's book SEP, page 107 and in Gragg &
    // Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.
    // Determine if we need to correct the residual. The goal is
    // to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||

    if (*rnorm > sq2o2*V->aitr_wnorm) { goto LINE100; }

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

    if ((V->aitr_j == 0) || (V->aitr_restart))
    {
        h[V->aitr_j] = 0.0;
    }
    h[V->aitr_j + ldh] += workd[irj + V->aitr_j];

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

        *rnorm = V->aitr_rnorm1;

    } else {

        // Another step of iterative refinement step is required.

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

    // Make sure the last off-diagonal element is non negative
    // If not perform a similarity transformation on H(1:j,1:j)
    // and scale v(:,j) by -1.

    if (h[V->aitr_j] < 0.0)
    {
        h[V->aitr_j] = -h[V->aitr_j];
        if (V->aitr_j < k + np - 1)
        {
            dscal_(&n, &dblm1, &v[V->aitr_j + 1], &int1);
        } else {
            dscal_(&n, &dblm1, resid, &int1);
        }
    }

    // STEP 6: Update  j = j+1;  Continue

    V->aitr_j += 1;
    if (V->aitr_j >= k + np)
    {
        V->ido = ido_DONE;
        return;
    }

    // Loop back to extend the factorization by another step.

    goto LINE1000;

}


void
dsapps(int n, int* kev, int np, double* shift, double* v, int ldv, double* h, int ldh,
       double* resid, double* q, int ldq, double* workd)
{
    int i, iend, istart, jj, kplusp, tmp_int, int1 = 1;
    double a1, a2, a3, a4, c, f, g, r, s, sigma, tst1;
    double dbl0 = 0.0, dbl1 = 1.0, dblm1 = -1.0;

    iend = 0;
    kplusp = *kev + np;

    // Initialize Q to the identity to accumulate
    // the rotations and reflections

    dlaset_("A", &kplusp, &kplusp, &dbl0, &dbl1, q, &ldq);

    // Quick return if there are no shifts to apply

    if (np == 0) { return; }

    for (jj = 0; jj < np; jj++)
    {
        sigma = shift[jj];

        // Check for splitting and deflation. Currently we consider
        // an off-diagonal element h(i+1,1) negligible if
        //         h(i+1,1) .le. epsmch*( |h(i,2)| + |h(i+1,2)| )
        // for i=1:KEV+NP-1.
        // If above condition tests true then we set h(i+1,1) = 0.
        // Note that h(1:KEV+NP,1) are assumed to be non negative.

        istart = 0;
        while (istart < kplusp - 1)
        {

            for (iend = istart; iend < kplusp - 1; iend++)
            {
                tst1 = fabs(h[iend + ldh]) + fabs(h[iend + 1 + ldh]);
                if (h[iend + 1] <= ulp * tst1)
                {
                    h[iend + 1] = 0.0;
                    break;
                }
            }

            // Scalar block, skipping, correct the sign if necessary
            if (istart == iend)
            {
                istart += 1;
                if (h[iend] < 0.0)
                {
                    h[iend] = -h[iend];
                    dscal_(&kplusp, &dblm1, &q[ldq*(iend)], &int1);
                }
                continue;
            }

            // We have a valid block [istart, iend] inclusive
            f = h[istart + ldh] - sigma;
            g = h[istart + 1];

            for (i = istart; i < iend; i++)
            {
                // Applying the plane rotations that create and chase the bulge X
                //
                // [c,  s] [ x x        ] [c, -s]       [ x x X      ]
                // [-s, c] [ x x x      ] [s,  c]       [ x x x      ]
                //         [   x x x    ]          =>   [ X x x x    ]
                //         [     x x x  ]               [     x x x  ]
                //         [         ...]               [         ...]
                //
                // dlartgp (instead of dlartg) is used to make sure that the
                // off-diagonal elements stay non-negative, (cf. F77 code for
                // manual handling).

                //                                  a1          a2
                // [c,  s] [ k m ] [c, -s]    [ c*k + s*m,  s*k + c*m] [c, -s]
                // [-s, c] [ m n ] [s,  c]    [-s*k + c*m, -s*m + c*n] [s,  c]
                //                                  a3          a4

                dlartgp_(&f, &g, &c, &s, &r);
                if (i > istart)
                {
                    h[i] = r;
                }
                a1 = c*h[i + ldh]     + s*h[i + 1];
                a2 = c*h[i + 1]       + s*h[i + 1 + ldh];
                a4 = c*h[i + 1 + ldh] - s*h[i + 1];
                a3 = c*h[i + 1]       - s*h[i + ldh];
                h[i + ldh]     = c*a1 + s*a2;              // h[i  , i  ]
                h[i + 1 + ldh] = c*a4 - s*a3;              // h[i+1, i+1]
                h[i + 1]       = c*a3 + s*a4;              // h[i+1, i  ]

                // Accumulate the rotation also in Q
                tmp_int = (i + jj + 2 > kplusp ? kplusp : i + jj + 2);
                drot_(&tmp_int, &q[ldq*i], &int1, &q[ldq*(i+1)], &int1, &c, &s);

                if (i < iend - 1)
                {
                    // g is the bulge created by the rotation
                    f = h[i + 1];
                    g = s*h[i + 2];
                    h[i + 2] = c*h[i + 2];
                }
            }
            istart = iend + 1;
            if (h[iend] < 0.0)
            {
                h[iend] = -h[iend];
                dscal_(&kplusp, &dblm1, &q[ldq*(iend)], &int1);
            }
        }
    }
    // 90

    // All shifts have been applied. Check for
    // more possible deflation that might occur
    // after the last shift is applied.

    for (i = 0; i < kplusp - 1; i++)
    {
        tst1 = fabs(h[i + ldh]) + fabs(h[i+1 + ldh]);
        if (h[i+1] <= ulp*tst1)
        {
            h[i+1] = 0.0;
        }
    }
    // 100

    // Compute the (kev+1)-st column of (V*Q) and
    // temporarily store the result in WORKD(N+1:2*N).
    // This is not necessary if h(kev+1,1) = 0.

    if (h[*kev] > 0.0)
    {
        dgemv_("N", &n, &kplusp, &dbl1, v, &ldv, &q[ldq*(*kev)], &int1, &dbl0, &workd[n], &int1);
    }

    // Compute column 1 to kev of (V*Q) in backward order
    // taking advantage that Q is an upper triangular matrix
    // with lower bandwidth np.
    // Place results in v(:,kplusp-kev:kplusp) temporarily.

    for (i = 0; i < *kev; i++)
    {
        tmp_int = kplusp - i;
        dgemv_("N", &n, &tmp_int, &dbl1, v, &ldv, &q[ldq*(*kev-i-1)], &int1, &dbl0, workd, &int1);
        dcopy_(&n, workd, &int1, &v[ldv*(kplusp-i-1)], &int1);
    }
    // 130

    // Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev).

    for (i = 0; i < *kev; i++)
    {
        dcopy_(&n, &v[ldv*(np+i)], &int1, &v[ldv*i], &int1);
    }
    // 140

    if (h[*kev] > 0.0)
    {
        dcopy_(&n, &workd[n], &int1, &v[ldv*(*kev)], &int1);
    }

    // Update the residual vector:
    //    r <- sigmak*r + betak*v(:,kev+1)
    // where
    //    sigmak = (e_{kev+p}'*Q)*e_{kev}
    //    betak = e_{kev+1}'*H*e_{kev}

    dscal_(&n, &q[kplusp-1 + (*kev-1)*ldq], resid, &int1);
    if (h[*kev] > 0.0)
    {
        daxpy_(&n, &h[*kev], &v[ldv*(*kev)], &int1, resid, &int1);
    }

    return;
}


void
dsgets(struct ARNAUD_state_d *V, int* kev, int* np, double* ritz,
       double* bounds, double* shifts)
{
    int kevd2, tmp1, tmp2, int1 = 1;
    if (V->which == which_BE)
    {
        // Both ends of the spectrum are requested.
        // Sort the eigenvalues into algebraically increasing
        // order first then swap high end of the spectrum next
        // to low end in appropriate locations.
        // NOTE: when np < floor(kev/2) be careful not to swap
        // overlapping locations.

        dsortr(which_LA, 1, *kev + *np, ritz, bounds);
        kevd2 = *kev / 2;
        if (*kev > 1)
        {
            tmp1 = (kevd2 > *np ? *np : kevd2);
            tmp2 = (kevd2 > *np ? kevd2 : *np);
            dswap_(&tmp1, ritz, &int1, &ritz[tmp2], &int1);
            dswap_(&tmp1, bounds, &int1, &bounds[tmp2], &int1);
        }
    } else {

        // LM, SM, LA, SA case.
        // Sort the eigenvalues of H into the desired order
        // and apply the resulting order to BOUNDS.
        // The eigenvalues are sorted so that the wanted part
        // are always in the last KEV locations.

        dsortr(V->which, 1, *kev + *np, ritz, bounds);
    }

    if ((V->shift == 1) && (*np > 0))
    {

        // Sort the unwanted Ritz values used as shifts so that
        // the ones with largest Ritz estimates are first.
        // This will tend to minimize the effects of the
        // forward instability of the iteration when the shifts
        // are applied in subroutine dsapps.

        dsortr(which_SM, 1, *np, bounds, ritz);
        dcopy_(np, ritz, &int1, shifts, &int1);
    }
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
dsortr(const enum ARNAUD_which w, const int apply, const int n, double* x1, double* x2)
{
    int i, gap, pos;
    double temp;
    ARNAUD_compare_rfunc *f;

    switch (w)
    {
        case which_LM:
            f = sortr_LM;
            break;
        case which_SM:
            f = sortr_SM;
            break;
        case which_LA:
            f = sortr_LA;
            break;
        case which_SA:
            f = sortr_SA;
            break;
        default:
            f = sortr_LM;
            break;
    }

    gap = n / 2;

    while (gap != 0)
    {
        for (i = gap; i < n; i++)
        {
            pos = i - gap;
            while ((pos >= 0) && (f(x1[pos], x1[pos+gap])))
            {
                temp = x1[pos];
                x1[pos] = x1[pos+gap];
                x1[pos+gap] = temp;

                if (apply)
                {
                    temp = x2[pos];
                    x2[pos] = x2[pos+gap];
                    x2[pos+gap] = temp;
                }
                pos -= gap;
            }
        }
        gap = gap / 2;
    }
}


static void
dsesrt(const enum ARNAUD_which w, const int apply, const int n, double* x, int na, double* a, const int lda)
{
    int i, gap, pos, int1 = 1;
    double temp;
    ARNAUD_compare_rfunc *f;

    switch (w)
    {
        case which_LM:
            f = sortr_LM;
            break;
        case which_SM:
            f = sortr_SM;
            break;
        case which_LA:
            f = sortr_LA;
            break;
        case which_SA:
            f = sortr_SA;
            break;
        default:
            f = sortr_LM;
            break;
    }

    gap = n / 2;

    while (gap != 0)
    {
        for (i = gap; i < n; i++)
        {
            pos = i - gap;
            while ((pos >= 0) && (f(x[pos], x[pos + gap])))
            {
                temp = x[pos];
                x[pos] = x[pos+gap];
                x[pos+gap] = temp;

                if (apply)
                {
                    dswap_(&na, &a[lda*pos], &int1, &a[lda*(pos+gap)], &int1);
                }
                pos -= gap;
            }
        }
        gap = gap / 2;
    }
}


void
dstqrb(int n, double* d, double* e, double* z, double* work, int* info)
{
    int int1 = 1, int0 = 0;
    double eps2 = pow(ulp, 2.0);
    double safmin = unfl;
    double safmax = (1.0 / safmin);
    double ssfmax = sqrt(safmax) / 3.0;
    double ssfmin = sqrt(safmin) / eps2;

    int nmaxit, jtot, i, ii, j, k, l1, m = 0, tmp_int = 0, l, lsv, lend, lendsv, iscale;
    double anorm = 0.0, rt1 = 0.0, rt2 = 0.0, c = 0.0, s = 0.0, g = 0.0, r = 0.0, p = 0.0;
    double b, f, tst;

    *info = 0;
    if (n == 0) { return; }
    if (n == 1) { z[0] = 1.0; return; }

    // Set z as the last row of identity matrix
    for (i = 0; i < n-1; i++) { z[i] = 0.0; }
    z[n-1] = 1.0;

    nmaxit = n*30;
    jtot = 0;

    //  Determine where the matrix splits and choose QL or QR iteration
    //  for each block, according to whether top or bottom diagonal
    //  element is smaller.

    // Translation Note:
    // All indices are 1-based, since the F77 code is very complicated.
    // Instead array indices are decremented where necessary.

    l1 = 1;

    while (jtot < nmaxit)
    {
        if (l1 > n) { break; }

        if (l1 > 1) { e[l1 - 2] = 0.0; }
        if (l1 <= n - 1)
        {
            for (m = l1; m <= n - 1; m++)
            {
                tst = fabs(e[m - 1]);
                if (tst == 0.0) { break; }
                if (tst <= (sqrt(fabs(d[m-1]))*sqrt(fabs(d[m])))*ulp)
                {
                    e[m-1] = 0.0;
                    break;
                }
            }
            // 20
        } else {
            m = n;
        }
        // 30

        // m will mark the splitting point, if any.
        l = l1;
        lsv = l;
        lend = m;
        lendsv = lend;
        l1 = m + 1;

        // Scalar block, skipping
        if (lend == l) { continue; }

        // Scale submatrix in rows and columns L to LEND
        tmp_int = lend - l + 1;
        anorm = dlanst_("I", &tmp_int, &d[l-1], &e[l-1]);
        iscale = 0;

        if (anorm == 0.0) { continue; }

        if (anorm > ssfmax)
        {
            iscale = 1;
            dlascl_("G", &int0, &int0, &anorm, &ssfmax, &tmp_int, &int1, &d[l-1], &n, info);
            tmp_int -= 1;
            dlascl_("G", &int0, &int0, &anorm, &ssfmax, &tmp_int, &int1, &e[l-1], &n, info);
        } else if (anorm < ssfmin) {
            iscale = 2;
            dlascl_("G", &int0, &int0, &anorm, &ssfmin, &tmp_int, &int1, &d[l-1], &n, info);
            tmp_int -= 1;
            dlascl_("G", &int0, &int0, &anorm, &ssfmin, &tmp_int, &int1, &e[l-1], &n, info);
        }
        // Choose between QL and QR iteration

        if (fabs(d[lend-1]) < fabs(d[l-1]))
        {
            lend = lsv;
            l = lendsv;
        }
        if (lend > l)
        {
            // QL Iteration
            while (1)
            {
                // Look for small subdiagonal element.
                // 40
                if (l != lend)
                {
                    for (m = l; m < lend; m++)
                    {
                        tst = fabs(e[m-1]);
                        tst = tst*tst;
                        if (tst <= (eps2*fabs(d[m - 1]))*fabs(d[m]) + safmin) { break; }
                        if (m == lend - 1) { m = lend; break; }  // No break condition
                    }
                    // 50, 60
                } else {
                    m = lend;
                }
                if (m < lend) { e[m - 1] = 0.0; }

                p = d[l - 1];
                if (m == l)
                {
                        // 80
                        // Eigenvalue found
                        d[l - 1] = p;
                        l += 1;
                        if (l <= lend) { continue; }  // Top of QL iteration
                        break;  // go to 140 undo scaling
                }
                // If remaining matrix is 2x2, use dlaev2 to compute its eigensystem
                if (m == l + 1)
                {
                    dlaev2_(&d[l - 1], &e[l - 1], &d[l], &rt1, &rt2, &c, &s);
                    work[l - 1] = c;
                    work[n - 1 + l - 1] = s;
                    tst    = z[l];
                    z[l]   = c*tst - s*z[l-1];
                    z[l-1] = s*tst + c*z[l-1];
                    d[l-1] = rt1;
                    d[l]   = rt2;
                    e[l-1] = 0.0;
                    l += 2;
                    if (l <= lend) { continue; }  // go to 40
                    break;  // go to 140
                }

                if (jtot == nmaxit) { break; } // go to 140
                jtot += 1;

                // Form shift
                g = (d[l]- p) / (2.0 * e[l-1]);
                r = hypot(g, 1.0);
                g = d[m-1] - p + (e[l-1] / (g + copysign(r, g)));

                s = 1.0;
                c = 1.0;
                p = 0.0;

                // Inner loop
                for (i = m - 1; i >= l; i--)
                {
                    f = s * e[i-1];
                    b = c * e[i-1];
                    dlartg_(&g, &f, &c, &s, &r);
                    if (i != m - 1) { e[i] = r; }
                    g = d[i] - p;
                    r = (d[i-1] - g)*s + 2.0*c*b;
                    p = s*r;
                    d[i] = g + p;
                    g = c*r - b;
                    work[i-1] = c;
                    work[n-1+i-1] = -s;
                }
                // 70
                tmp_int = m - l + 1;
                dlasr_("R", "V", "B", &int1, &tmp_int, &work[l-1], &work[n-1+l-1], &z[l-1], &int1);

                d[l-1] = d[l-1] - p;
                e[l-1] = g;

            }
        } else {
            // QR Iteration

            // Look for small subdiagonal element.
            while (1)
            {
                if (l != lend)
                {
                    for (m = l; m > lend; m--)
                    {
                        tst = fabs(e[m-2]);
                        tst = tst*tst;
                        if (tst <= (eps2*fabs(d[m-1]))*fabs(d[m-2]) + safmin) { break; }
                        if (m == lend+1) { m = lend; break; }  // No break
                    }
                } else {
                    m = lend;
                }
                // 100, 110
                if (m > lend) { e[m-2] = 0.0; }
                p = d[l-1];
                if (m == l)
                {
                    // 130
                    // Eigenvalue found
                    d[l-1] = p;
                    l -= 1;
                    if (l >= lend) { continue; } // Top of QR iteration
                    break; // go to 140 undo scaling
                }
                // If remaining matrix is 2x2, use dlaev2 to compute its eigensystem
                if (m == l - 1)
                {
                    dlaev2_(&d[l-2], &e[l-2], &d[l-1], &rt1, &rt2, &c, &s);
                    tst    = z[l-1];
                    z[l-1] = c*tst - s*z[l-2];
                    z[l-2] = s*tst + c*z[l-2];
                    d[l-2] = rt1;
                    d[l-1] = rt2;
                    e[l-2] = 0.0;
                    l -= 2;

                    if (l >= lend) { continue; } // Top of QR iteration
                    break; // go to 140 undo scaling
                }

                if (jtot == nmaxit) { break; }  // go to 140
                jtot += 1;

                // Form the shift
                g = (d[l-2] - p) / (2.0*e[l-2]);
                r = hypot(g, 1.0);
                g = d[m-1] - p + (e[l-2] / (g + copysign(r, g)));

                s = 1.0;
                c = 1.0;
                p = 0.0;

                // Inner loop
                for (i = m; i < l; i++)
                {
                    f = s * e[i-1];
                    b = c * e[i-1];
                    dlartg_(&g, &f, &c, &s, &r);
                    if (i != m) { e[i-2] = r; }
                    g = d[i-1] - p;
                    r = (d[i] - g)*s + 2.0*c*b;
                    p = s*r;
                    d[i-1] = g + p;
                    g = c*r - b;

                    // Save rotations
                    work[i-1] = c;
                    work[n-1+i-1] = s;
                }
                // 120
                // Apply saved rotations.
                tmp_int = l - m + 1;
                dlasr_("R", "V", "F", &int1, &tmp_int, &work[m-1], &work[n-1+m-1], &z[m-1], &int1);

                d[l-1] = d[l-1] - p;
                e[l - 2] = g;

            }
        }
        // 140 Still in the outer while loop; it breaks at the top

        // Undo scaling if necessary
        tmp_int = lendsv-lsv+1;
        if (iscale == 1)
        {

            dlascl_("G", &int0, &int0, &ssfmax, &anorm, &tmp_int, &int1, &d[lsv-1], &n, info);
            tmp_int -= 1;
            dlascl_("G", &int0, &int0, &ssfmax, &anorm, &tmp_int, &int1, &e[lsv-1], &n, info);

        } else if (iscale == 2) {

            dlascl_("G", &int0, &int0, &ssfmin, &anorm, &tmp_int, &int1, &d[lsv-1], &n, info);
            tmp_int -= 1;
            dlascl_("G", &int0, &int0, &ssfmin, &anorm, &tmp_int, &int1, &e[lsv-1], &n, info);

        }

        // Check for no convergence to an eigenvalue after a total of n*maxit iterations
        if (jtot >= nmaxit)
        {
            for (i = 0; i < n-1; i++) { if (e[i] != 0.0) { *info += 1; } }
            return;  // 150
        }
    }
    // Out of the while loop

    //  Order eigenvalues and eigenvectors.
    //  Use selection sort to minimize swaps of eigenvectors.
    for (ii = 1; ii < n; ii++)
    {
        i = ii - 1;
        k = i;
        p = d[i];

        for (j = ii; j < n; j++)
        {
            if (d[j] < p)
            {
                k = j;
                p = d[j];
            }
        }
        // 170
        if (k != i)
        {
            d[k] = d[i];
            d[i] = p;
            p = z[k];
            z[k] = z[i];
            z[i] = p;
        }
    }
    // 180

    return;
}


int sortr_LM(const double x1, const double x2) { return (fabs(x1) > fabs(x2)); }
int sortr_SM(const double x1, const double x2) { return (fabs(x1) < fabs(x2)); }
int sortr_LA(const double x1, const double x2) { return (x1 > x2); }
int sortr_SA(const double x1, const double x2) { return (x1 < x2); }
