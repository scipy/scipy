#include "_arpack_n_single_complex.h"


typedef int ARPACK_compare_cfunc(const ARPACK_CPLXF_TYPE, const ARPACK_CPLXF_TYPE);
typedef int ARPACK_compare_rfunc(const float, const float);

static int sortc_LM(const ARPACK_CPLXF_TYPE, const ARPACK_CPLXF_TYPE);
static int sortc_SM(const ARPACK_CPLXF_TYPE, const ARPACK_CPLXF_TYPE);
static int sortc_LR(const ARPACK_CPLXF_TYPE, const ARPACK_CPLXF_TYPE);
static int sortc_SR(const ARPACK_CPLXF_TYPE, const ARPACK_CPLXF_TYPE);
static int sortc_LI(const ARPACK_CPLXF_TYPE, const ARPACK_CPLXF_TYPE);
static int sortc_SI(const ARPACK_CPLXF_TYPE, const ARPACK_CPLXF_TYPE);


static const float unfl = 1.1754943508222875e-38;
// static const float ovfl = 1.0 / 1.1754943508222875e-38;
static const float ulp = 1.1920928955078125e-07;

static void cnaup2(struct ARPACK_arnoldi_update_vars_s *V, ARPACK_CPLXF_TYPE* resid,
    ARPACK_CPLXF_TYPE* v, int ldv, ARPACK_CPLXF_TYPE* h, int ldh, ARPACK_CPLXF_TYPE* ritz,
    ARPACK_CPLXF_TYPE* bounds, ARPACK_CPLXF_TYPE* q, int ldq, ARPACK_CPLXF_TYPE* workl,
    int* ipntr, ARPACK_CPLXF_TYPE* workd, ARPACK_CPLXF_TYPE* rwork
);


enum ARPACK_neupd_type {
    REGULAR,
    SHIFTI,
    REALPART,
    IMAGPART
};

void
cneupd(struct ARPACK_arnoldi_update_vars_s *V, int rvec, int howmny, int* select,
       ARPACK_CPLXF_TYPE* d, ARPACK_CPLXF_TYPE* z, int ldz, float sigma,
       ARPACK_CPLXF_TYPE* workev, ARPACK_CPLXF_TYPE* resid, ARPACK_CPLXF_TYPE* v, int ldv,
       int* ipntr, ARPACK_CPLXF_TYPE* workd, ARPACK_CPLXF_TYPE* workl, float* rwork)
{
    const float eps23 = powf(ulp, 2.0 / 3.0);
    int ibd, iconj, ih, iheig, ihbds, iuptri, invsub, iq, irz, iwev, j, jj;
    int bounds, k, ldh, ldq, np, numcnv, reord, ritz, wr;
    int iwork[1] = { 0 };
    int ierr = 0, int1 = 1, tmp_int = 0, nconv2 = 0, outncv;
    float conds, sep, temp1, dbl0 = 0.0, dbl1 = 1.0, dblm1 = -1.0;
    ARPACK_CPLXF_TYPE rnorm, temp;
    ARPACK_CPLXF_TYPE cdbl0 = ARPACK_cplxf(0.0, 0.0);
    ARPACK_CPLXF_TYPE cdbl1 = ARPACK_cplxf(1.0, 0.0);
    ARPACK_CPLXF_TYPE cdblm1 = ARPACK_cplxf(-1.0, 0.0);
    ARPACK_CPLXF_TYPE vl[1] = { 0.0 };
    enum ARPACK_neupd_type TYP;

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
    } else if (V->mode == 3) {
        TYP = SHIFTI;
    } else {
        ierr = -10;
    }

    if ((V->mode == 1) && (V->bmat)) { ierr = -11; }

    if (ierr != 0) {
        V->info = ierr;
        return;
    }

     /*-------------------------------------------------------*
     | Pointer into WORKL for address of H, RITZ, WORKEV, Q   |
     | etc... and the remaining workspace.                    |
     | Also update pointer to be used on output.              |
     | Memory is laid out as follows:                         |
     | workl(1:ncv*ncv) := generated Hessenberg matrix        |
     | workl(ncv*ncv+1:ncv*ncv+ncv) := ritz values            |
     | workl(ncv*ncv+ncv+1:ncv*ncv+2*ncv) := error bounds     |
     *-------------------------------------------------------*/

     /*----------------------------------------------------------*
     | The following is used and set by ZNEUPD.                  |
     | workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := The untransformed |
     |                                      Ritz values.         |
     | workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed |
     |                                      error bounds of      |
     |                                      the Ritz values      |
     | workl(ncv*ncv+4*ncv+1:2*ncv*ncv+4*ncv) := Holds the upper |
     |                                      triangular matrix    |
     |                                      for H.               |
     | workl(2*ncv*ncv+4*ncv+1: 3*ncv*ncv+4*ncv) := Holds the    |
     |                                      associated matrix    |
     |                                      representation of    |
     |                                      the invariant        |
     |                                      subspace for H.      |
     | GRAND total of NCV * ( 3 * NCV + 4 ) locations.           |
     *----------------------------------------------------------*/

    ih     = ipntr[4];
    ritz   = ipntr[5];
    iq     = ipntr[6];
    bounds = ipntr[7];
    ldh = V->ncv;
    ldq = V->ncv;
    iheig  = bounds + ldh;
    ihbds  = iheig  + ldh;
    iuptri = ihbds  + ldh;
    invsub = iuptri + ldh*V->ncv;
    ipntr[8]  = iheig;
    ipntr[10] = ihbds;
    ipntr[11] = iuptri;
    ipntr[12] = invsub;
    wr = 0;
    iwev = wr + V->ncv;

     /*----------------------------------------*
     | irz points to the Ritz values computed  |
     |     by _neigh before exiting _naup2.    |
     | ibd points to the Ritz estimates        |
     |     computed by _neigh before exiting   |
     |     _naup2.                             |
     *----------------------------------------*/

    irz = ipntr[13] + (V->ncv)*(V->ncv);
    ibd = irz + V->ncv;

     /*-----------------------------------*
     | RNORM is B-norm of the RESID(1:N). |
     *-----------------------------------*/
    rnorm = workl[ih+1];
    workl[ih+1] = ARPACK_cplxf(0.0, 0.0);

    if (rvec) {
        reord = 0;

         /*--------------------------------------------------*
         | Use the temporary bounds array to store indices   |
         | These will be used to mark the select array later |
         *--------------------------------------------------*/
        for (j = 0; j < V->ncv; j++)
        {
            workl[bounds + j] = ARPACK_cplxf(j, 0.0);
            select[j] = 0;
        }
        // 10

         /*------------------------------------*
         | Select the wanted Ritz values.      |
         | Sort the Ritz values so that the    |
         | wanted ones appear at the tailing   |
         | NEV positions of workl(irr) and     |
         | workl(iri).  Move the corresponding |
         | error estimates in workl(ibd)       |
         | accordingly.                        |
         *------------------------------------*/

        np = V->ncv - V->nev;
        cngets(V, V->nev, &np, &workl[irz], &workl[bounds]);

         /*----------------------------------------------------*
         | Record indices of the converged wanted Ritz values  |
         | Mark the select array for possible reordering       |
         *----------------------------------------------------*/

        numcnv = 0;
        for (j = 0; j < V->ncv; j++)
        {
            temp1 = fmaxf(eps23, cabsf(workl[irz+V->ncv-j]));
            jj = (int)crealf(workl[bounds + V->ncv - j]);

            if ((numcnv < V->nconv) && (cabsf(workl[ibd + jj]) <= V->tol*temp1))
            {
                select[jj - 1] = 1;
                numcnv += 1;
                if (jj > V->nconv) { reord = 1; }
            }
        }
        // 11

         /*----------------------------------------------------------*
         | Check the count (numcnv) of converged Ritz values with    |
         | the number (nconv) reported by znaupd.  If these two      |
         | are different then there has probably been an error       |
         | caused by incorrect passing of the znaupd data.           |
         *----------------------------------------------------------*/

        if (numcnv != V->nconv)
        {
            V->info = -15;
            return;
        }

         /*------------------------------------------------------*
         | Call LAPACK routine zlahqr  to compute the Schur form |
         | of the upper Hessenberg matrix returned by ZNAUPD .   |
         | Make a copy of the upper Hessenberg matrix.           |
         | Initialize the Schur vector matrix Q to the identity. |
         *------------------------------------------------------*/
        tmp_int = ldh*V->ncv;
        zcopy_(&tmp_int, &workl[ih], &int1, &workl[iuptri], &int1);
        zlaset_("A", &V->ncv, &V->ncv, &cdbl0, &cdbl1, &workl[invsub], &ldq);
        zlahqr_(&int1, &int1, &V->ncv, &int1, &V->ncv, &workl[iuptri], &ldh,
                &workl[iheig], &int1, &V->ncv, &workl[invsub], &ldq, &ierr);
        zcopy_(&V->ncv, &workl[invsub + V->ncv - 1], &ldq, &workl[ihbds], &int1);

        if (ierr != 0)
        {
            V->info = -8;
            return;
        }

        if (reord)
        {
            ztrsen_("N", "V", select, &V->ncv, &workl[iuptri], &ldh, &workl[invsub], &ldq,
                    &workl[iheig], &nconv2, &conds, &sep, workev, &V->ncv, &ierr);

            if (nconv2 < V->nconv) { V->nconv = nconv2; }
            if (ierr == 1) {
                V->info = 1;
                return;
            }
        }

    } else {
        zcopy_(&V->nconv, &workl[ritz], &int1, d, &int1);
        zcopy_(&V->nconv, &workl[ritz], &int1, &workl[iheig], &int1);
        zcopy_(&V->nconv, &workl[bounds], &int1, &workl[ihbds], &int1);
    }

     /*-----------------------------------------------*
     | Transform the Ritz values and possibly vectors |
     | and corresponding error bounds of OP to those  |
     | of A*x = lambda*B*x.                           |
     *-----------------------------------------------*/

    if (TYP == REGULAR)
    {
        if (rvec)
        {
            zscal_(&V->ncv, &rnorm, &workl[ihbds], &int1);
        }
    } else {
         /*--------------------------------------*
         |   A spectral transformation was used. |
         | * Determine the Ritz estimates of the |
         |   Ritz values in the original system. |
         *--------------------------------------*/
        if (rvec)
        {
            zscal_(&V->ncv, &rnorm, &workl[ihbds], &int1);
        }
        for (k = 0; k < V->ncv; k++)
        {
#if defined(_MSC_VER)
            // Complex division is not supported in MSVC, multiply with reciprocal
            temp = _FCmulcr(conjf(workl[iheig + k]), 1.0 / cabsf(workl[iheig + k]));
            workl[ihbds + k] = _FCmulcc(_FCmulcc(workl[ihbds + k], temp), temp);
#else
            temp = workl[iheig + k];
            workl[ihbds + k] = workl[ihbds + k] / temp / temp;
#endif
        }
        // 50
    }

     /*----------------------------------------------------------*
     | *  Transform the Ritz values back to the original system. |
     |    For TYPE = 'SHIFTI' the transformation is              |
     |             lambda = 1/theta + sigma                      |
     | NOTES:                                                    |
     | *The Ritz vectors are not affected by the transformation. |
     *----------------------------------------------------------*/

    if (TYP == SHIFTI)
    {
        for (k = 0; k < V->nconv; k++)
        {
#if defined(_MSC_VER)
            // Complex division is not supported in MSVC
            temp = _FCmulcr(conjf(workl[iheig + k]), 1.0 / cabsf(workl[iheig + k]));
            d[k] = ARPACK_cplxf(crealf(temp) + sigma, cimagf(temp));
#else
            d[k] = 1.0 / workl[iheig + k] + sigma;
#endif
        }
        // 60
    }

     /*------------------------------------------------*
     | Eigenvector Purification step. Formally perform |
     | one of inverse subspace iteration. Only used    |
     | for MODE = 3. See reference 3.                  |
     *------------------------------------------------*/

    if ((rvec) && (howmny == 0) && (TYP == SHIFTI))
    {
         /*-----------------------------------------------*
         | Purify the computed Ritz vectors by adding a   |
         | little bit of the residual vector:             |
         |                      T                         |
         |          resid(:)*( e    s ) / theta           |
         |                      NCV                       |
         | where H s = s theta.                           |
         *-----------------------------------------------*/

        for (j = 0; j < V->nconv; j++)
        {
            if ((crealf(workl[iheig+j]) != 0.0) || (cimagf(workl[iheig+j]) != 0.0))
            {
#if defined(_MSC_VER)
                // Complex division is not supported in MSVC
                temp = _FCmulcr(conjf(workl[iheig + j]), 1.0 / cabsf(workl[iheig + j]));
                workev[j] = _FCmulcc(workl[invsub + j*ldq + V->ncv], temp);
#else
                workev[j] = workl[invsub + j*ldq + V->ncv] / workl[iheig+j];
#endif
            }
        }
        // 100

         /*--------------------------------------*
         | Perform a rank one update to Z and    |
         | purify all the Ritz vectors together. |
         *--------------------------------------*/
        zgeru_(&V->n, &V->nconv, &cdbl1, resid, &int1, workev, &int1, z, &ldz);
    }

    return;
};


void
cnaupd(struct ARPACK_arnoldi_update_vars_s *V, ARPACK_CPLXF_TYPE* resid,
       ARPACK_CPLXF_TYPE* v, int ldv, int* ipntr, ARPACK_CPLXF_TYPE* workd,
       ARPACK_CPLXF_TYPE* workl, float* rwork)
{
    int bounds, ierr, ih, iq, iw, ldh, ldq, nev0, next, iritz;

    // perform basic checks
    if (V->n <= 0) {
        ierr = -1;
    } else if (V->nev <= 0) {
        ierr = -2;
    } else if ((V->ncv < V->nev + 1) || (V->ncv > V->n)) {
        ierr = -3;
    } else if (V->maxiter <= 0) {
        ierr = -4;
    } else if ((V->which < 0) || (V->which > 5)) {
        ierr = -5;
    } else if ((V->bmat != 0) && (V->bmat != 1)) {
        ierr = -6;
    } else if ((V->mode < 1) || (V->mode > 3)) {
        ierr = -10;
    } else if ((V->mode == 1) && (V->bmat == 1)) {
        ierr = -11;
    }

    if (ierr != 0) {
        V->info = ierr;
        V->ido = 99;
        return;
    }

    if (V->tol <= 0.0) {
        V-> tol = ulp;
    }

    if ((V->ishift != 0) || (V->ishift != 1) || (V->ishift != 2))
    {
        V->ishift = 1;
    }

     /*---------------------------------------------*
     | NP is the number of additional steps to      |
     | extend the length NEV Lanczos factorization. |
     | NEV0 is the local variable designating the   |
     | size of the invariant subspace desired.      |
     *---------------------------------------------*/
    V->np = V->ncv - V->nev;
    nev0 = V->nev;

    for (int j = 0; j < 3 * (V->ncv*V->ncv) + 6*V->ncv; j++)
    {
        workl[j] = ARPACK_cplxf(0.0, 0.0);
    }

     /*------------------------------------------------------------*
     | Pointer into WORKL for address of H, RITZ, BOUNDS, Q        |
     | etc... and the remaining workspace.                         |
     | Also update pointer to be used on output.                   |
     | Memory is laid out as follows:                              |
     | workl(1:ncv*ncv) := generated Hessenberg matrix             |
     | workl(ncv*ncv+1:ncv*ncv+ncv) := the ritz values             |
     | workl(ncv*ncv+ncv+1:ncv*ncv+2*ncv)   := error bounds        |
     | workl(ncv*ncv+2*ncv+1:2*ncv*ncv+2*ncv) := rotation matrix Q |
     | workl(2*ncv*ncv+2*ncv+1:3*ncv*ncv+5*ncv) := workspace       |
     | The final workspace is needed by subroutine zneigh  called  |
     | by znaup2 . Subroutine zneigh  calls LAPACK routines for    |
     | calculating eigenvalues and the last row of the eigenvector |
     | matrix.                                                     |
     *------------------------------------------------------------*/

    ldh    = V->ncv;
    ldq    = V->ncv;
    ih     = 0;
    iritz  = ih     + ldh*V->ncv;
    bounds = iritz  + V->ncv;
    iq     = bounds + V->ncv;
    iw     = iq     + ldq*V->ncv;
    next   = iw     + (V->ncv*V->ncv) + 3*V->ncv;

    ipntr[3] = next;
    ipntr[4] = ih;
    ipntr[5] = iritz;
    ipntr[6] = iq;
    ipntr[7] = bounds;
    ipntr[13]  = iw;


    cnaup2(V, resid, v, ldv, &workl[ih], ldh, &workl[iritz], &workl[bounds],
           &workl[iq], ldq, &workl[iw], ipntr, workd, rwork);

     /*-------------------------------------------------*
     | ido .ne. 99 implies use of reverse communication |
     | to compute operations involving OP or shifts.    |
     *-------------------------------------------------*/

    if (V->ido != 99) { return; }

    V->nconv = V->np;

    if (V->info < 0) { return; }
    if (V->info == 2) { V->info = 3; }

    return;
}


void
cnaup2(struct ARPACK_arnoldi_update_vars_s *V, ARPACK_CPLXF_TYPE* resid,
       ARPACK_CPLXF_TYPE* v, int ldv, ARPACK_CPLXF_TYPE* h, int ldh,
       ARPACK_CPLXF_TYPE* ritz, ARPACK_CPLXF_TYPE* bounds,
       ARPACK_CPLXF_TYPE* q, int ldq, ARPACK_CPLXF_TYPE* workl, int* ipntr,
       ARPACK_CPLXF_TYPE* workd, ARPACK_CPLXF_TYPE* rwork)
{
    enum ARPACK_which temp_which;
    int i, int1 = 1, j;
    const float eps23 = powf(ulp, 2.0 / 3.0);
    float temp = 0.0, rtemp;

    if (V->ido == 0)
    {
        V->aup2_nev0 = V->nev;
        V->aup2_np0 = V->np;

         /*------------------------------------*
         | kplusp is the bound on the largest  |
         |        Lanczos factorization built. |
         | nconv is the current number of      |
         |        "converged" eigenvlues.      |
         | iter is the counter on the current  |
         |      iteration step.                |
         *------------------------------------*/

        V->aup2_kplusp = V->nev + V->np;
        V->nconv = 0;
        V->aup2_iter = 0;

         /*--------------------------------------*
         | Set flags for computing the first NEV |
         | steps of the Arnoldi factorization.   |
         *--------------------------------------*/
        V->aup2_getv0 = 1;
        V->aup2_update = 0;
        V->aup2_ushift = 0;
        V->aup2_cnorm = 0;

        if (V->info != 0)
        {
             /*-------------------------------------------*
             | User provides the initial residual vector. |
             *-------------------------------------------*/
            V->aup2_initv = 1;
            V->info = 0;
        } else {
            V->aup2_initv = 0;
        }
    }

     /*--------------------------------------------*
     | Get a possibly random starting vector and   |
     | force it into the range of the operator OP. |
     *--------------------------------------------*/
    if (V->aup2_getv0)
    {
        cgetv0(V, V->aup2_initv, V->n, 0, v, ldv, resid, &V->aup2_rnorm, ipntr, workd);
        if (V->ido != 99) { return; }
        if (V->aup2_rnorm == 0.0)
        {
             /*------------------------------------------------------*
             | The initial vector is zero. It is improbable to hit   |
             | all zeros randomly. Hence regenerate a vector and let |
             | the rest figure things out for zero eigenvalues.      |
             *------------------------------------------------------*/
            generate_random_vector_z(&V->n, resid);
            V->aup2_rnorm = dznrm2_(&V->n, resid, &int1);
            zdscal_(&V->n, &V->aup2_rnorm, resid, &int1);
        }
        V->aup2_getv0 = 0;
        V->ido = 0;
    }

     /*----------------------------------*
     | Back from reverse communication : |
     | continue with update step         |
     *----------------------------------*/
    if (V->aup2_update) { goto LINE20; }

     /*------------------------------------------*
     | Back from computing user specified shifts |
     *------------------------------------------*/
    if (V->aup2_ushift) { goto LINE50; }

     /*------------------------------------*
     | Back from computing residual norm   |
     | at the end of the current iteration |
     *------------------------------------*/
    if (V->aup2_cnorm) { goto LINE100; }

     /*---------------------------------------------------------*
     | Compute the first NEV steps of the Arnoldi factorization |
     *---------------------------------------------------------*/

    cnaitr(V, resid, &V->aup2_rnorm, v, ldv, h, ldh, ipntr, workd, &V->info);

     /*--------------------------------------------------*
     | ido .ne. 99 implies use of reverse communication  |
     | to compute operations involving OP and possibly B |
     *--------------------------------------------------*/
    if (V->ido != 99) { return; }

    if (V->info > 0)
    {
        V->np = V->info;
        V->maxiter = V->aup2_iter;
        V->info = -9999;
        V->ido = 99;
        return;
    }

     /*-------------------------------------------------------------*
     |                                                              |
     |           M A I N  ARNOLDI  I T E R A T I O N  L O O P       |
     |           Each iteration implicitly restarts the Arnoldi     |
     |           factorization in place.                            |
     |                                                              |
     *-------------------------------------------------------------*/

LINE1000:
    V->aup2_iter += 1;

     /*----------------------------------------------------------*
     | Compute NP additional steps of the Arnoldi factorization. |
     | Adjust NP since NEV might have been updated by last call  |
     | to the shift application routine dnapps .                 |
     *----------------------------------------------------------*/
    V->np = V->aup2_kplusp - V->nev;
    V->ido = 0;

LINE20:
    V->aup2_update = 1;

    cnaitr(V, resid, &V->aup2_rnorm, v, ldv, h, ldh, ipntr, workd);

     /*--------------------------------------------------*
     | ido .ne. 99 implies use of reverse communication  |
     | to compute operations involving OP and possibly B |
     *--------------------------------------------------*/

    if (V->ido != 99) { return; }

    if (V->info > 0) {
        V->np = V->info;
        V->maxiter = V->aup2_iter;
        V->info = -9999;
        V->ido = 99;
        return;
    }
    V->aup2_update = 0;

     /*-------------------------------------------------------*
     | Compute the eigenvalues and corresponding error bounds |
     | of the current upper Hessenberg matrix.                |
     *-------------------------------------------------------*/

    cneigh(&V->aup2_rnorm, V->aup2_kplusp, h, ldh, ritz, bounds,
           q, ldq, workl, rwork, &V->info);

    if (V->info != 0)
    {
       V->info = -8;
       V->ido = 99;
       return;
    }

     /*--------------------------------------------------*
     | Select the wanted Ritz values and their bounds    |
     | to be used in the convergence test.               |
     | The wanted part of the spectrum and corresponding |
     | error bounds are in the last NEV loc. of RITZ,    |
     | and BOUNDS respectively.                          |
     *--------------------------------------------------*/

    V->nev = V->aup2_nev0;
    V->np = V->aup2_np0;

     /*-------------------------------------------------*
     | Make a copy of Ritz values and the corresponding |
     | Ritz estimates obtained from zneigh .            |
     *-------------------------------------------------*/

    zcopy_(&V->aup2_kplusp, ritz, &int1, &workl[V->aup2_kplusp*V->aup2_kplusp], &int1);
    zcopy_(&V->aup2_kplusp, bounds, &int1, &workl[V->aup2_kplusp*V->aup2_kplusp + V->aup2_kplusp], &int1);

     /*--------------------------------------------------*
     | Select the wanted Ritz values and their bounds    |
     | to be used in the convergence test.               |
     | The wanted part of the spectrum and corresponding |
     | bounds are in the last NEV loc. of RITZ           |
     | BOUNDS respectively.                              |
     *--------------------------------------------------*/

    cngets(V, &V->nev, &V->np, ritz, bounds);

     /*-----------------------------------------------------------*
     | Convergence test: currently we use the following criteria. |
     | The relative accuracy of a Ritz value is considered        |
     | acceptable if:                                             |
     |                                                            |
     | error_bounds(i) .le. tol*max(eps23, magnitude_of_ritz(i)). |
     |                                                            |
     *-----------------------------------------------------------*/

    for (i = 0; i < V->nev; i++)
    {
        rtemp = fmaxf(eps23, hypot(crealf(ritz[V->np + i]), cimagf(ritz[V->np + i])));
        if (hypot(crealf(bounds[V->np + i]), cimagf(bounds[V->np + i])) <= V->tol*rtemp)
        {
            V->nconv += 1;
        }
    }
    // 25


     /*--------------------------------------------------------*
     | Count the number of unwanted Ritz values that have zero |
     | Ritz estimates. If any Ritz estimates are equal to zero |
     | then a leading block of H of order equal to at least    |
     | the number of Ritz values with zero Ritz estimates has  |
     | split off. None of these Ritz values may be removed by  |
     | shifting. Decrease NP the number of shifts to apply. If |
     | no shifts may be applied, then prepare to exit          |
     *--------------------------------------------------------*/

    // We are modifying V->np hence the temporary variable.
    int nptemp = V->np;

    for (j = 0; j < nptemp; j++)
    {
        if ((crealf(bounds[j]) == 0.0) && (cimagf(bounds[j]) == 0.0))
        {
            V->np -= 1;
            V->nev += 1;
        }
    }
    // 30

    if ((V->nconv >= V->aup2_nev0) || (V->aup2_iter > V->maxiter) || (V->np == 0))
    {
         /*-----------------------------------------------*
         | Prepare to exit. Put the converged Ritz values |
         | and corresponding bounds in RITZ(1:NCONV) and  |
         | BOUNDS(1:NCONV) respectively. Then sort. Be    |
         | careful when NCONV > NP                        |
         *-----------------------------------------------*/

         /*-----------------------------------------*
         |  Use h( 3,1 ) as storage to communicate  |
         |  rnorm to _neupd if needed               |
         *-----------------------------------------*/

         h[2] = ARPACK_cplxf(V->aup2_rnorm, 0.0);

         /*---------------------------------------------*
         | To be consistent with dngets , we first do a |
         | pre-processing sort in order to keep complex |
         | conjugate pairs together.  This is similar   |
         | to the pre-processing sort used in dngets    |
         | except that the sort is done in the opposite |
         | order.                                       |
         *---------------------------------------------*/

        // Translation note: Is this all because ARPACK did not have complex sort?

        if (V->which == which_LM) { temp_which = which_SM; }
        if (V->which == which_SM) { temp_which = which_LM; }
        if (V->which == which_LR) { temp_which = which_SR; }
        if (V->which == which_SR) { temp_which = which_LR; }
        if (V->which == which_LI) { temp_which = which_SI; }
        if (V->which == which_SI) { temp_which = which_LI; }

        csortc(temp_which, 1, V->aup2_kplusp, ritz, bounds);

         /*-------------------------------------------------*
         | Scale the Ritz estimate of each Ritz value       |
         | by 1 / max(eps23,magnitude of the Ritz value).   |
         *-------------------------------------------------*/

        for (j = 0; j < V->aup2_nev0; j++)
        {
            temp = fmaxf(eps23, hypot(crealf(ritz[j]), cimagf(ritz[j])));
            bounds[j] = ARPACK_cplxf(crealf(bounds[j]) / temp, cimagf(bounds[j]) / temp);
        }
        // 35

         /*---------------------------------------------------*
         | Sort the Ritz values according to the scaled Ritz  |
         | estimates.  This will push all the converged ones  |
         | towards the front of ritzr, ritzi, bounds          |
         | (in the case when NCONV < NEV.)                    |
         *---------------------------------------------------*/
        temp_which = which_LM;
        csortc(temp_which, 1, V->aup2_nev0, bounds, ritz);

         /*---------------------------------------------*
         | Scale the Ritz estimate back to its original |
         | value.                                       |
         *---------------------------------------------*/

        for (j = 0; j < V->aup2_nev0; j++)
        {
            temp = fmaxf(eps23, hypot(crealf(ritz[j]), cimagf(ritz[j])));
            bounds[j] = ARPACK_cplxf(crealf(bounds[j]) / temp, cimagf(bounds[j]) / temp);
        }
        // 40

         /*-----------------------------------------------*
         | Sort the converged Ritz values again so that   |
         | the "threshold" value appears at the front of  |
         | ritzr, ritzi and bound.                        |
         *-----------------------------------------------*/

        csortc(V->which, 1, V->nconv, ritz, bounds);

        if ((V->aup2_iter > V->maxiter) && (V->nconv < V->aup2_nev0))
        {
             /*-----------------------------------*
             | Max iterations have been exceeded. |
             *-----------------------------------*/
            V->info = 1;
        }

        if ((V->np == 0) && (V->nconv < V->aup2_nev0))
        {
             /*--------------------*
             | No shifts to apply. |
             *--------------------*/
            V->info = 2;
        }

        V->np = V->nconv;
        V->maxiter = V->aup2_iter;
        V->nev = V->nconv;
        V->ido = 99;
        return;

    } else if ((V->nconv < V->aup2_nev0) && (V->ishift)) {
         /*------------------------------------------------*
         | Do not have all the requested eigenvalues yet.  |
         | To prevent possible stagnation, adjust the size |
         | of NEV.                                         |
         *------------------------------------------------*/
        int nevbef = V->nev;
        V->nev += (V->nconv > (V->np / 2) ? (V->np / 2) : V->nconv);
        if ((V->nev == 1) && (V->aup2_kplusp >= 6)) {
            V->nev = V->aup2_kplusp / 2;
        } else if ((V->nev == 1) && (V->aup2_kplusp > 3)) {
            V->nev = 2;
        }

        V->np = V->aup2_kplusp - V->nev;

        if (nevbef < V->nev) {
            cngets(V, &V->nev, &V->np, ritz, bounds);
        }

    }

    if (V->ishift == 0)
    {
         /*------------------------------------------------------*
         | User specified shifts: pop back out to get the shifts |
         | and return them in the first 2*NP locations of WORKL. |
         *------------------------------------------------------*/
        V->aup2_ushift = 1;
        V->ido = 3;
        return;
    }

LINE50:
     /*-----------------------------------*
     | Back from reverse communication;   |
     | User specified shifts are returned |
     | in WORKL(1:2*NP)                   |
     *-----------------------------------*/

    V->aup2_ushift = 0;

    if (V->ishift != 1)
    {
         /*---------------------------------*
         | Move the NP shifts from WORKL to |
         | RITZR, RITZI to free up WORKL    |
         | for non-exact shift case.        |
         *---------------------------------*/
        zcopy_(&V->np, workl, &int1, ritz, &int1);
    }

     /*--------------------------------------------------------*
     | Apply the NP implicit shifts by QR bulge chasing.       |
     | Each shift is applied to the whole upper Hessenberg     |
     | matrix H.                                               |
     | The first 2*N locations of WORKD are used as workspace. |
     *--------------------------------------------------------*/

    cnapps(V->n, &V->nev, V->np, ritz, v, ldv, h, ldh, resid, q, ldq, workl, workd);

     /*--------------------------------------------*
     | Compute the B-norm of the updated residual. |
     | Keep B*RESID in WORKD(1:N) to be used in    |
     | the first step of the next call to dnaitr . |
     *--------------------------------------------*/

    V->aup2_cnorm = 1;
    if (V->bmat)
    {
        zcopy_(&V->n, resid, &int1, &workd[V->n], &int1);
        ipntr[0] = V->n;
        ipntr[1] = 0;
        V->ido = 2;
         /*---------------------------------*
         | Exit in order to compute B*RESID |
         *---------------------------------*/
        return;
    } else {
        zcopy_(&V->n, resid, &int1, workd, &int1);
    }

LINE100:
     /*---------------------------------*
     | Back from reverse communication; |
     | WORKD(1:N) := B*RESID            |
     *---------------------------------*/

    if (V->bmat)
    {
        V->aup2_rnorm = cabsf(cdotc_(&V->n, resid, &int1, workd, &int1));
        V->aup2_rnorm = sqrtf(V->aup2_rnorm);
    } else {
        V->aup2_rnorm = dznrm2_(&V->n, resid, &int1);
    }

    goto LINE1000;

     /*--------------------------------------------------------------*
     |                                                               |
     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
     |                                                               |
     *--------------------------------------------------------------*/

}


void
cnaitr(struct ARPACK_arnoldi_update_vars_s *V,  ARPACK_CPLXF_TYPE* resid, float* rnorm,
       ARPACK_CPLXF_TYPE* v, int ldv, ARPACK_CPLXF_TYPE* h, int ldh, int* ipntr,
       ARPACK_CPLXF_TYPE* workd)
{
    int i, infol, iter, ipj, irj, ivj, jj, n, tmp_int;
    float smlnum = unfl * ( V->n / ulp);
    // float xtemp[2] = { 0.0 };
    const float sq2o2 = sqrtf(2.0) / 2.0;

    char *MTYPE = "G", *TRANS = "T", *NORM = "1";
    int int1 = 1;
    float dbl1 = 1.0, dbl0 = 0.0, temp1, tmp_dbl, tst1;
    ARPACK_CPLXF_TYPE cdbl1 = ARPACK_cplxf(1.0, 0.0);
    ARPACK_CPLXF_TYPE cdblm1 = ARPACK_cplxf(-1.0, 0.0);
    ARPACK_CPLXF_TYPE cdbl0 = ARPACK_cplxf(0.0, 0.0);

    n = V->n;  // constant value, just for typing convenience
    ipj = 0;
    irj = ipj + n;
    ivj = irj + n;

    if (V->ido == 0)
    {
         /*-----------------------------*
         | Initial call to this routine |
         *-----------------------------*/
        V->info = 0;
        V->aitr_step3 = 0;
        V->aitr_step4 = 0;
        V->aitr_orth1 = 0;
        V->aitr_orth2 = 0;
        V->aitr_restart = 0;
        V->aitr_j = V->nev;
    }

     /*------------------------------------------------*
     | When in reverse communication mode one of:      |
     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
     | will be .true. when ....                        |
     | STEP3: return from computing OP*v_{j}.          |
     | STEP4: return from computing B-norm of OP*v_{j} |
     | ORTH1: return from computing B-norm of r_{j+1}  |
     | ORTH2: return from computing B-norm of          |
     |        correction to the residual vector.       |
     | RSTART: return from OP computations needed by   |
     |         dgetv0.                                 |
     *------------------------------------------------*/
    if (V->aitr_step3) { goto LINE50; }
    if (V->aitr_step4) { goto LINE60; }
    if (V->aitr_orth1) { goto LINE70; }
    if (V->aitr_orth2) { goto LINE90; }
    if (V->aitr_restart) { goto LINE30; }

     /*----------------------------*
     | Else this is the first step |
     *----------------------------*/

     /*-------------------------------------------------------------*
     |                                                              |
     |        A R N O L D I     I T E R A T I O N     L O O P       |
     |                                                              |
     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
     *-------------------------------------------------------------*/

LINE1000:

     /*--------------------------------------------------*
     | STEP 1: Check if the B norm of j-th residual      |
     | vector is zero. Equivalent to determining whether |
     | an exact j-step Arnoldi factorization is present. |
     *--------------------------------------------------*/

    V->aitr_betaj = *rnorm;
    if (*rnorm > 0.0) { goto LINE40; }

     /*--------------------------------------------------*
     | Invariant subspace found, generate a new starting |
     | vector which is orthogonal to the current Arnoldi |
     | basis and continue the iteration.                 |
     *--------------------------------------------------*/
    V->aitr_betaj = 0.0;
    V->getv0_itry = 1;

LINE20:
    V->aitr_restart = 1;
    V->ido = 0;

LINE30:
    dgetv0(V, 0, n, V->aitr_j, v, ldv, resid, rnorm, ipntr, workd, &V->aitr_ierr);
    if (V->ido != 99) { return; }

    if (V->aitr_ierr < 0)
    {
        V->getv0_itry += 1;
        if (V->getv0_itry <= 3) { goto LINE20; }
         /*-----------------------------------------------*
         | Give up after several restart attempts.        |
         | Set INFO to the size of the invariant subspace |
         | which spans OP and exit.                       |
         *-----------------------------------------------*/
        V->info = V->aitr_j;
        V->ido = 99;
        return;
    }

LINE40:

     /*--------------------------------------------------------*
     | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
     | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
     | when reciprocating a small RNORM, test against lower    |
     | machine bound.                                          |
     *--------------------------------------------------------*/
    zcopy_(&n, resid, &int1, &v[ldv*V->aitr_j], &int1);
    if (*rnorm >= unfl)
    {
        temp1 = 1.0 / *rnorm;
        zdscal_(&n, &temp1, &v[ldv*V->aitr_j], &int1);
        zdscal_(&n, &temp1, &workd[ipj], &int1);
    } else {
        zlascl_(MTYPE, &i, &i, rnorm, &dbl1, &n, &int1, &v[ldv*V->aitr_j], &n, &infol);
        zlascl_(MTYPE, &i, &i, rnorm, &dbl1, &n, &int1, &workd[ipj], &n, &infol);
    }

     /*-----------------------------------------------------*
     | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
     | Note that this is not quite yet r_{j}. See STEP 4    |
     *-----------------------------------------------------*/
    V->aitr_step3 = 1;
    zcopy_(&n, &v[ldv*V->aitr_j], &int1, &workd[ivj], &int1);
    ipntr[0] = ivj;
    ipntr[1] = irj;
    ipntr[2] = ipj;
    V->ido = 1;

     /*----------------------------------*
     | Exit in order to compute OP*v_{j} |
     *----------------------------------*/
    return;

LINE50:
     /*---------------------------------*
     | Back from reverse communication; |
     | WORKD(IRJ:IRJ+N-1) := OP*v_{j}   |
     | if step3 = .true.                |
     *---------------------------------*/
    V->aitr_step3 = 0;

     /*-----------------------------------------*
     | Put another copy of OP*v_{j} into RESID. |
     *-----------------------------------------*/
    zcopy_(&n, &workd[irj], &int1, resid, &int1);

     /*--------------------------------------*
     | STEP 4:  Finish extending the Arnoldi |
     |          factorization to length j.   |
     *--------------------------------------*/
    if (V->bmat)
    {
        V->aitr_step4 = 1;
        ipntr[0] = irj;
        ipntr[1] = ipj;
        V->ido = 2;
         /*------------------------------------*
         | Exit in order to compute B*OP*v_{j} |
         *------------------------------------*/
        return;
    } else {
        zcopy_(&n, resid, &int1, &workd[ipj], &int1);
    }

LINE60:
     /*---------------------------------*
     | Back from reverse communication; |
     | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j} |
     | if step4 = .true.                |
     *---------------------------------*/
    V->aitr_step4 = 0;

     /*------------------------------------*
     | The following is needed for STEP 5. |
     | Compute the B-norm of OP*v_{j}.     |
     *------------------------------------*/

    if (V->bmat)
    {
        V->aitr_wnorm = cabsf(cdotc_(&n, resid, &int1, &workd[ipj], &int1));
    } else {
        V->aitr_wnorm = dznrm2_(&n, resid, &int1);
    }

     /*----------------------------------------*
     | Compute the j-th residual corresponding |
     | to the j step factorization.            |
     | Use Classical Gram Schmidt and compute: |
     | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
     | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
     *----------------------------------------*/

     /*-----------------------------------------*
     | Compute the j Fourier coefficients w_{j} |
     | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
     *-----------------------------------------*/
    zgemv_(TRANS, &n, &V->aitr_j, &cdbl1, v, &ldv, &workd[ipj], &int1, &dbl0, &h[ldh*V->aitr_j], &int1);

     /*-------------------------------------*
     | Orthogonalize r_{j} against V_{j}.   |
     | RESID contains OP*v_{j}. See STEP 3. |
     *-------------------------------------*/
    TRANS = "N";
    zgemv_(TRANS, &n, &V->aitr_j, &cdblm1, v, &ldv, &h[ldv*V->aitr_j], &int1, &cdbl1, resid, &int1);
    if (V->aitr_j > 0) { h[V->aitr_j + ldh*(V->aitr_j-1)] = ARPACK_cplxf(V->aitr_betaj, 0.0); }

    V->aitr_orth1 = 1;
    if (V->bmat)
    {
        zcopy_(&n, resid, &int1, &workd[irj], &int1);
        ipntr[0] = irj;
        ipntr[1] = ipj;
        V->ido = 2;
         /*---------------------------------*
         | Exit in order to compute B*r_{j} |
         *---------------------------------*/
        return;
    } else {
        zcopy_(&n, resid, &int1, &workd[ipj], &int1);
    }

LINE70:

     /*--------------------------------------------------*
     | Back from reverse communication if ORTH1 = .true. |
     | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
     *--------------------------------------------------*/

    V->aitr_orth1 = 0;

     /*-----------------------------*
     | Compute the B-norm of r_{j}. |
     *-----------------------------*/
    if (V->bmat)
    {
        *rnorm = cabsf(cdotc_(&n, resid, &int1, &workd[ipj], &int1));
    } else {
        *rnorm = dznrm2_(&n, resid, &int1);
    }

     /*----------------------------------------------------------*
     | STEP 5: Re-orthogonalization / Iterative refinement phase |
     | Maximum NITER_ITREF tries.                                |
     |                                                           |
     |          s      = V_{j}^T * B * r_{j}                     |
     |          r_{j}  = r_{j} - V_{j}*s                         |
     |          alphaj = alphaj + s_{j}                          |
     |                                                           |
     | The stopping criteria used for iterative refinement is    |
     | discussed in Parlett's book SEP, page 107 and in Gragg &  |
     | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
     | Determine if we need to correct the residual. The goal is |
     | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
     | The following test determines whether the sine of the     |
     | angle between  OP*x and the computed residual is less     |
     | than or equal to 0.7071.                                  |
     *----------------------------------------------------------*/

    if (*rnorm > sq2o2) { goto LINE100; }
    V->aitr_iter = 0;

     /*--------------------------------------------------*
     | Enter the Iterative refinement phase. If further  |
     | refinement is necessary, loop back here. The loop |
     | variable is ITER. Perform a step of Classical     |
     | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
     *--------------------------------------------------*/
LINE80:

     /*---------------------------------------------------*
     | Compute V_{j}^T * B * r_{j}.                       |
     | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
     *---------------------------------------------------*/

    TRANS = "T";
    zgemv_(TRANS, &n, &V->aitr_j, &cdbl1, v, &ldv, &workd[ipj], &int1, &cdbl0, &workd[irj], &int1);

     /*--------------------------------------------*
     | Compute the correction to the residual:     |
     | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1). |
     | The correction to H is v(:,1:J)*H(1:J,1:J)  |
     | + v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j.         |
     *--------------------------------------------*/
    TRANS = "N";
    tmp_dbl = -1.0;
    zgemv_(TRANS, &n, &V->aitr_j, &cdblm1, v, &ldv, &workd[irj], &int1, &cdbl1, resid, &int1);
    zaxpy_(&V->aitr_j, &cdbl1, &workd[irj], &int1, &h[ldh*V->aitr_j], &int1);
    V->aitr_orth2 = 1;
    if (V->bmat)
    {
        zcopy_(&n, resid, &int1, &workd[irj], &int1);
        ipntr[0] = irj;
        ipntr[1] = ipj;
        V->ido = 2;
         /*----------------------------------*
         | Exit in order to compute B*r_{j}. |
         | r_{j} is the corrected residual.  |
         *----------------------------------*/
        return;
    } else {
        zcopy_(&n, resid, &int1, &workd[ipj], &int1);
    }

LINE90:
     /*--------------------------------------------------*
     | Back from reverse communication if ORTH2 = .true. |
     *--------------------------------------------------*/

     /*----------------------------------------------------*
     | Compute the B-norm of the corrected residual r_{j}. |
     *----------------------------------------------------*/
    if (V->bmat)
    {
        V->aitr_rnorm1 = cabsf(cdotc_(&n, resid, &int1, &workd[ipj], &int1));
    } else {
        V->aitr_rnorm1 = dznrm2_(&n, resid, &int1);
    }

     /*----------------------------------------*
     | Determine if we need to perform another |
     | step of re-orthogonalization.           |
     *----------------------------------------*/
    if (V->aitr_rnorm1 > sq2o2)
    {
         /*--------------------------------------*
         | No need for further refinement.       |
         | The cosine of the angle between the   |
         | corrected residual vector and the old |
         | residual vector is greater than 0.717 |
         | In other words the corrected residual |
         | and the old residual vector share an  |
         | angle of less than arcCOS(0.717)      |
         *--------------------------------------*/
        *rnorm = V->aitr_rnorm1;

    } else {
         /*------------------------------------------*
         | Another step of iterative refinement step |
         | is required.                              |
         *------------------------------------------*/
        *rnorm = V->aitr_rnorm1;
        iter += 1;
        if (iter < 2) { goto LINE80; }

         /*------------------------------------------------*
         | Otherwise RESID is numerically in the span of V |
         *------------------------------------------------*/
        for (jj = 0; jj < n; jj++)
        {
            resid[jj] = ARPACK_cplxf(0.0, 0.0);
        }
        *rnorm = 0.0;
    }

LINE100:

    V->aitr_restart = 0;
    V->aitr_orth2 = 0;

     /*-----------------------------------*
     | STEP 6: Update  j = j+1;  Continue |
     *-----------------------------------*/
    V->aitr_j += 1;
    if (V->aitr_j >= V->nev + V->np)
    {
        V->ido = 99;
        int k = (V->nev > 1 ? V->nev : 1);
        for (i = k - 1; i < k - 2 + V->np; i++)
        {
             /*-------------------------------------------*
             | Check for splitting and deflation.         |
             | Use a standard test as in the QR algorithm |
             | REFERENCE: LAPACK subroutine dlahqr        |
             *-------------------------------------------*/
            tst1 = cabsf(h[i + ldh*i]) + cabsf(h[i+1 + ldh*(i+1)]);
            if (tst1 == 0.0)
            {
                tmp_int = k + V->np;
                // dlanhs(norm, n, a, lda, work) with "work" being float type
                // Recasting complex workspace to float
                tst1 = zlanhs_(NORM, &tmp_int, h, &ldh, (float*)&workd[n]);
            }
            if (cabsf(h[i+1 + ldh*i]) <= fmaxf(ulp*tst1, smlnum))
            {
                h[i+1 + ldh*i] = ARPACK_cplxf(0.0, 0.0);
            }
        }
        // 110
        return;
    }
    goto LINE1000;

}


void
cnapps(int n, int* kev, int np, ARPACK_CPLXF_TYPE* shift, ARPACK_CPLXF_TYPE* v,
       int ldv, ARPACK_CPLXF_TYPE* h, int ldh, ARPACK_CPLXF_TYPE* resid,
       ARPACK_CPLXF_TYPE* q, int ldq, ARPACK_CPLXF_TYPE* workl,
       ARPACK_CPLXF_TYPE* workd)
{
    int cconj;
    int i, j, jj, int1 = 1, istart, iend, tmp_int;
    int kplusp = *kev + np;

    float smlnum = unfl * ( n / ulp);
    float c, tst1;
    float tmp_dbl, dbl1 = 1.0, dbl0 = 0.0;
    // float u[3] = { 0.0 };
    ARPACK_CPLXF_TYPE f, g, h11, h21, sigma, s, r, t, tmp_cplx, tmp_cplx2;
    ARPACK_CPLXF_TYPE cdbl1 = ARPACK_cplxf(1.0, 0.0);
    ARPACK_CPLXF_TYPE cdbl0 = ARPACK_cplxf(0.0, 0.0);
    char *NORM = "1", *TRANS = "N";

     /*-------------------------------------------*
     | Initialize Q to the identity to accumulate |
     | the rotations and reflections              |
     *-------------------------------------------*/
    for (i = 0; i < kplusp; i++)
    {
        for (j = 0; j < kplusp; j++)
        {
            q[j + ldq*i] = ARPACK_cplxf(0.0, 0.0);
            if (i == j) { q[j + ldq*i] = ARPACK_cplxf(1.0, 0.0); }
        }
    }

     /*---------------------------------------------*
     | Quick return if there are no shifts to apply |
     *---------------------------------------------*/
    if (np == 0) { return; }

     /*---------------------------------------------*
     | Chase the bulge with the application of each |
     | implicit shift. Each shift is applied to the |
     | whole matrix including each block.           |
     *---------------------------------------------*/

    cconj = 0;
    for (jj = 0; jj < np; jj++)
    {
        sigma = shift[jj];
        istart = 0;

LINE20:

        for (i = istart; i < kplusp - 1; i++)
        {
             /*---------------------------------------*
             | Check for splitting and deflation. Use |
             | a standard test as in the QR algorithm |
             | REFERENCE: LAPACK subroutine zlahqr    |
             *---------------------------------------*/

            tst1 = fabs(crealf(h[i + ldh*i])) + fabs(cimagf(h[i + ldh*i])) +
                   fabs(crealf(h[i+1 + ldh*(i+1)])) + fabs(cimagf(h[i+1 + ldh*(i+1)]));
            if (tst1 == 0.0)
            {
                // clanhs(norm, n, a, lda, work)
                tmp_int = kplusp - jj + 1;
                zlanhs_(NORM, &tmp_int, h, &ldh, (float*)workl);
            }
            if (fabs(crealf(h[i+1 + ldh*i])) <= fmaxf(ulp*tst1, smlnum))
            {
                iend = i;
                h[i+1 + ldh*i] = ARPACK_cplxf(0.0, 0.0);
                break;
            }
        }
        // 30
        // Adjust iend if no break
        if (i == kplusp - 1) { iend = kplusp - 1; }
        // 40

         /*-----------------------------------------------*
         | No reason to apply a shift to block of order 1 |
         | or if the current block starts after the point |
         | of compression since we'll discard this stuff  |
         *-----------------------------------------------*/

        if ((istart == iend) || (istart > *kev))
        {
            // go to 100
            istart = iend + 1;
            if (iend < kplusp) { goto LINE20; }
            continue;
        }

        h11 = h[istart + ldh*istart];
        h21 = h[istart + 1 + ldh*istart];
        // Because MSVC is not C99 compliant hence does not support arithmetic ops.
        // f = h11 - sigma;
        f = ARPACK_cplxf(crealf(h11)-crealf(sigma), cimagf(h11)-cimagf(sigma));
        g = h21;

        for (i = istart; i <= iend-1; i++)
        {
             /*-----------------------------------------------------*
             | Construct the plane rotation G to zero out the bulge |
             *-----------------------------------------------------*/
            zlartg_(&f, &g, &c, &s, &r);
            if (i > istart)
            {
                h[i + ldh*(i-1)] = r;
                h[i + 1 + ldh*(i-1)] = ARPACK_cplxf(0.0, 0.0);
            }

             /*--------------------------------------------*
             | Apply rotation to the left of H;  H <- G'*H |
             *--------------------------------------------*/

            // More MSVC torture to get basic arithmetic ops
            // MSVC does not implement +,-,*,/ ops for complex numbers.
            for (j = i; j < kplusp; j++)
            {
                #if defined(_MSC_VER)
                    t = _FCmulcc(s, h[i + 1 + ldh*j]);
                    tmp_cplx = _FCmulcr(h[i + ldh*j], c);
                    t = ARPACK_cplxf(crealf(t) + crealf(tmp_cplx), cimagf(t) + cimagf(tmp_cplx));

                    tmp_cplx  = _FCmulcc(ARPACK_cplxf(-crealf(s),cimagf(s)), h[i + 1 + ldh*j]);
                    tmp_cplx2 = _FCmulcr(h[i + 1 + ldh*j], c);
                    h[i + 1 + ldh*j] = ARPACK_cplxf(crealf(tmp_cplx) + crealf(tmp_cplx2), cimagf(tmp_cplx)+cimagf(tmp_cplx2));

                    h[i + ldh*j]     = t;
                #else
                    t                =        c*h[i + ldh*j] + s*h[i + 1 + ldh*j];
                    h[i + 1 + ldh*j] = -conjf(s)*h[i + ldh*j] + c*h[i + 1 + ldh*j];
                    h[i + ldh*j]     = t;
                #endif
            }
            // 50

             /*--------------------------------------------*
             | Apply rotation to the right of H;  H <- H*G |
             *--------------------------------------------*/
            for (j = 0; j < (i+2 > iend ? iend : i+2); j++)
            {
                #if defined(_MSC_VER)
                    t = _FCmulcc(ARPACK_cplxf(crealf(s), -cimagf(s)), h[j + ldh*(i + 1)]);
                    tmp_cplx = _FCmulcr(h[j + ldh*i], c);
                    t = ARPACK_cplxf(crealf(t) + crealf(tmp_cplx), cimagf(t) + cimagf(tmp_cplx));

                    tmp_cplx  = _FCmulcc(ARPACK_cplxf(-crealf(s), -cimagf(s)), h[j + ldh*i]);
                    tmp_cplx2 = _FCmulcr(h[j + ldh*(i + 1)], c);
                    h[j + ldh*(i + 1)] = ARPACK_cplxf(crealf(tmp_cplx) + crealf(tmp_cplx2), cimagf(tmp_cplx) + cimagf(tmp_cplx2));

                    h[j + ldh*i]       = t;
                #else
                    t                  =  c*h[j + ldh*i] + conjf(s)*h[j + ldh*(i + 1)];
                    h[j + ldh*(i + 1)] = -s*h[j + ldh*i] +       c*h[j + ldh*(i + 1)];
                    h[j + ldh*i]       = t;
                #endif
            }
            // 60

             /*---------------------------------------------------*
             | Accumulate the rotation in the matrix Q;  Q <- Q*G |
             *---------------------------------------------------*/
            for (j = 0; j < (i+jj+1 > kplusp ? kplusp : i+jj+1); j++)
            {
                #if defined(_MSC_VER)
                    t = _FCmulcc(ARPACK_cplxf(crealf(s), -cimagf(s)), q[j + ldq*(i + 1)]);
                    tmp_cplx = _FCmulcr(q[j + ldq*i], c);
                    t = ARPACK_cplxf(crealf(t) + crealf(tmp_cplx), cimagf(t) + cimagf(tmp_cplx));

                    tmp_cplx  = _FCmulcc(ARPACK_cplxf(-crealf(s), -cimagf(s)), q[j + ldq*i]);
                    tmp_cplx2 = _FCmulcr(q[j + ldq*(i + 1)], c);
                    q[j + ldq*(i + 1)] = ARPACK_cplxf(crealf(tmp_cplx) + crealf(tmp_cplx2), cimagf(tmp_cplx) + cimagf(tmp_cplx2));

                    q[j + ldq*i]       = t;
                #else
                    t                  =  c*q[j + ldq*i] + conjf(s)*q[j + ldq*(i + 1)];
                    q[j + ldq*(i + 1)] = -s*q[j + ldq*i] +       c*q[j + ldq*(i + 1)];
                    q[j + ldq*i]       = t;
                #endif
            }
            // 70

             /*--------------------------*
             | Prepare for next rotation |
             *--------------------------*/
            if (i < iend-1)
            {
                f = h[i + 1 + ldh*i];
                g = h[i + 2 + ldh*i];
            }

        }
        //80


         /*--------------------------------------------------------*
         | Apply the same shift to the next block if there is any. |
         *--------------------------------------------------------*/
        istart = iend + 1;
        if (iend < kplusp) { goto LINE20; }
    }
    // 110

     /*-------------------------------------------------*
     | Perform a similarity transformation that makes   |
     | sure that H will have non negative sub diagonals |
     *-------------------------------------------------*/
    for (j = 0; j < *kev; j++)
    {
        if ((crealf(h[j+1 + ldh*j]) < 0.0) && (cimagf(h[j+1 + ldh*j]) != 0.0))
        {
            tmp_int = kplusp - j;
            tmp_dbl = cabsf(h[j+1 + ldh*j]);
            t = ARPACK_cplxf(crealf(h[j+1 + ldh*j]) / tmp_dbl,
                            cimagf(h[j+1 + ldh*j]) / tmp_dbl);
            tmp_cplx = conjf(t);

            zscal_(&tmp_int, &tmp_cplx, &h[j+1 + ldh*j], &ldh);
            tmp_int = (j+2 > kplusp ? kplusp : j+ 2) + 1;
            zscal_(&tmp_int, &t, &h[ldh*(j+1)], &int1);
            tmp_int = (j+np+1 > kplusp ? kplusp : j+np+1) + 1;
            zscal_(&tmp_int, &t, &q[ldq*(j+1)], &int1);
        }
    }
    // 120

    for (i = 0; i < *kev; i++)
    {
         /*-------------------------------------------*
         | Final check for splitting and deflation.   |
         | Use a standard test as in the QR algorithm |
         | REFERENCE: LAPACK subroutine zlahqr.       |
         | Note: Since the subdiagonals of the        |
         | compressed H are nonnegative real numbers, |
         | we take advantage of this.                 |
         *-------------------------------------------*/
        tst1 = fabs(crealf(h[i + ldh*i])) + fabs(crealf(h[i+1 + ldh*(i+1)])) +
               fabs(cimagf(h[i + ldh*i])) + fabs(cimagf(h[i+1 + ldh*(i+1)]));
        if (tst1 == 0.0)
        {
            tst1 = zlanhs_(NORM, kev, h, &ldh, (float*)workl);
        }
        if (crealf(h[i+1 + ldh*i]) <= fmaxf(ulp+tst1, smlnum))
        {
            h[i+1 + ldh*i] = ARPACK_cplxf(0.0, 0.0);
        }
    }
    // 130

     /*------------------------------------------------*
     | Compute the (kev+1)-st column of (V*Q) and      |
     | temporarily store the result in WORKD(N+1:2*N). |
     | This is needed in the residual update since we  |
     | cannot GUARANTEE that the corresponding entry   |
     | of H would be zero as in exact arithmetic.      |
     *------------------------------------------------*/
    if (crealf(h[*kev + ldh*(*kev-1)]) > 0.0)
    {
        zgemv_(TRANS, &n, &kplusp, &cdbl1, v, &ldv, &q[(*kev)*ldq], &int1, &cdbl0, &workd[n], &int1);
    }

     /*---------------------------------------------------------*
     | Compute column 1 to kev of (V*Q) in backward order       |
     | taking advantage of the upper Hessenberg structure of Q. |
     *---------------------------------------------------------*/
    for (i = 0; i < kev; i++)
    {
        tmp_int = kplusp - i;
        zgemv_(TRANS, &n, &tmp_int, &cdbl1, v, &ldv, &q[(*kev-i)*ldq], &int1, &cdbl0, workd, &int1);
        zcopy_(&n, workd, &int1, &v[(kplusp-i)*ldv], &int1);
    }

     /*------------------------------------------------*
     |  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev). |
     *------------------------------------------------*/

    // Use any letter other than "L", "U", hence TRANS here.
    zlacpy_(TRANS, &n, kev, &v[ldv*(kplusp - *kev)], &ldv, v, &ldv);

     /*-------------------------------------------------------------*
     | Copy the (kev+1)-st column of (V*Q) in the appropriate place |
     *-------------------------------------------------------------*/
    if (crealf(h[*kev + ldh*(*kev-1)]) > 0.0) {
        zcopy_(&n, &workd[n], &int1, &v[ldv*(*kev)], &int1);
    }

     /*------------------------------------*
     | Update the residual vector:         |
     |    r <- sigmak*r + betak*v(:,kev+1) |
     | where                               |
     |    sigmak = (e_{kplusp}'*Q)*e_{kev} |
     |    betak = e_{kev+1}'*H*e_{kev}     |
     *------------------------------------*/

    zscal_(&n, &q[kplusp-1 + ldq*(*kev-i)], resid, &int1);

    if (crealf(h[*kev + ldh*(*kev-1)]) > 0.0)
    {
        zaxpy_(&n, &h[*kev + ldh*(*kev-1)], &v[ldv*(*kev)], &int1, resid, &int1);
    }

    return;
}



void
cneigh(float* rnorm, int n, ARPACK_CPLXF_TYPE* h, int ldh, ARPACK_CPLXF_TYPE* ritz,
       ARPACK_CPLXF_TYPE* bounds, ARPACK_CPLXF_TYPE* q, int ldq, ARPACK_CPLXF_TYPE* workl,
       float* rwork, int* ierr)
{
    int select[1] = { 0 };
    int int1 = 1, j;
    float dbl1 = 1.0, temp;
    ARPACK_CPLXF_TYPE vl[1] = { 0.0 };
    ARPACK_CPLXF_TYPE c1 = ARPACK_cplxf(1.0, 0.0), c0 = ARPACK_cplxf(0.0, 0.0);
    char *UPLO = "A", *SIDE = "R", *HOWMNY = "B";


     /*---------------------------------------------------------*
     | 1. Compute the eigenvalues, the last components of the   |
     |    corresponding Schur vectors and the full Schur form T |
     |    of the current upper Hessenberg matrix H.             |
     |    zlahqr returns the full Schur form of H               |
     |    in WORKL(1:N**2), and the Schur vectors in q.         |
     *---------------------------------------------------------*/

    zlacpy_(UPLO, &n, &n, h, &ldh, workl, &n);
    zlaset_(UPLO, &n, &n, &c0, &c1, q, &ldq);
    for (j = 0; j < n-1; j++)
    {
        bounds[j] = ARPACK_cplxf(0.0, 0.0);
    }

    zlahqr_(&int1, &int1, &n, &int1, &n, workl, &ldh, ritz, &int1, &n, q, &ldq, ierr);
    if (*ierr != 0) { return; }

    zcopy_(&n, &q[n-2], &ldq, bounds, &int1);
     /*---------------------------------------------------------*
     | 2. Compute the eigenvectors of the full Schur form T and |
     |    apply the Schur vectors to get the corresponding      |
     |    eigenvectors.                                         |
     *---------------------------------------------------------*/
    ztrevc_(SIDE, HOWMNY, select, &n, workl, &n, vl, &n, q, &ldq, &n, &n, &workl[n*n], rwork, ierr);
    if (*ierr != 0) { return; }

     /*-----------------------------------------------*
     | Scale the returning eigenvectors so that their |
     | euclidean norms are all one. LAPACK subroutine |
     | ztrevc returns each eigenvector normalized so  |
     | that the element of largest magnitude has      |
     | magnitude 1; here the magnitude of a complex   |
     | number (x,y) is taken to be |x| + |y|.         |
     *-----------------------------------------------*/

    for (j = 0; j < n; j++)
    {
        temp = 1.0 / dznrm2_(&n, &q[j*ldq], &int1);
        zdscal_(&n, &temp, &q[j*ldq], &int1);
    }

     /*---------------------------*
     | Compute the Ritz estimates |
     *---------------------------*/

    zcopy_(&n, &q[n-1], &ldq, workl, &int1);
    zdscal_(&n, rnorm, bounds, &int1);

    return;
}


void
cngets(struct ARPACK_arnoldi_update_vars_s *V, int* kev, int* np,
       ARPACK_CPLXF_TYPE* ritz, ARPACK_CPLXF_TYPE* bounds)
{

    csortc(V->which, 1, *kev + *np, ritz, bounds);

    if (V->ishift == 1)
    {
         /*------------------------------------------------------*
         | Sort the unwanted Ritz values used as shifts so that  |
         | the ones with largest Ritz estimates are first        |
         | This will tend to minimize the effects of the         |
         | forward instability of the iteration when they shifts |
         | are applied in subroutine znapps.                     |
         | Be careful and use 'SM' since we want to sort BOUNDS! |
         *------------------------------------------------------*/
        csortc(which_SR, 1, *np, bounds, ritz);
    }

    return;
}


void
cgetv0(struct ARPACK_arnoldi_update_vars_s *V, int initv, int n, int j,
       ARPACK_CPLXF_TYPE* v, int ldv, ARPACK_CPLXF_TYPE* resid, float* rnorm,
       int* ipntr, ARPACK_CPLXF_TYPE* workd)
{
    int jj, int1 = 1, int0 = 0, intm1 = -1;
    char *TRANS = "C";
    const float sq2o2 = sqrtf(2.0) / 2.0;
    ARPACK_CPLXF_TYPE cnorm;
    ARPACK_CPLXF_TYPE c0 = ARPACK_cplxf(0.0, 0.0);
    ARPACK_CPLXF_TYPE c1 = ARPACK_cplxf(1.0, 0.0);
    ARPACK_CPLXF_TYPE cm1 = ARPACK_cplxf(-1.0, 0.0);

    if (V->ido == 0)
    {
        V->info = 0;
        V->getv0_iter = 0;
        V->getv0_first = 0;
        V->getv0_orth = 0;

         /*----------------------------------------------------*
         | Possibly generate a random starting vector in RESID |
         *----------------------------------------------------*/
        if (!(initv))
        {
            generate_random_vector_z(n, resid);
        }

         /*---------------------------------------------------------*
         | Force the starting vector into the range of OP to handle |
         | the generalized problem when B is possibly (singular).   |
         *---------------------------------------------------------*/

        if (V->getv0_itry == 1)
        {
            ipntr[0] = 0;
            ipntr[1] = n;
            zcopy_(&n, resid, &int1, workd, &int1);
            V->ido = -1;
            return;
        } else if ((V->getv0_itry > 1) && (V->bmat == 1))
        {
            zcopy_(&n, resid, &int1, &workd[n], &int1);
        }
    }

     /*----------------------------------------*
     | Back from computing OP*(initial-vector) |
     *----------------------------------------*/

    if (V->getv0_first) { goto LINE20; }

     /*-----------------------------------------------*
     | Back from computing OP*(orthogonalized-vector) |
     *-----------------------------------------------*/

    if (V->getv0_orth) { goto LINE40; }

     /*-----------------------------------------------------*
     | Starting vector is now in the range of OP; r = OP*r; |
     | Compute B-norm of starting vector.                   |
     *-----------------------------------------------------*/

    V->getv0_first = 1;
    if (V->getv0_itry == 1)
    {
        zcopy_(&n, &workd[n], &int1, resid, &int1);
    }
    if (V->bmat)
    {
        ipntr[0] = n;
        ipntr[1] = 0;
        V->ido = 2;
        return;
    }

LINE20:

    V->getv0_first = 0;
    if (V->bmat)
    {
        cnorm = cdotc_(&n, resid, &int1, workd, &int1);
        V->getv0_rnorm0 = sqrtf(cabsf(cnorm));
    } else {
        V->getv0_rnorm0 = dznrm2_(&n, resid, &int1);
    }
    *rnorm = V->getv0_rnorm0;

     /*--------------------------------------------%
     | Exit if this is the very first Arnoldi step |
     *--------------------------------------------*/

    if (j == 0)
    {
        V->ido = 99;
        return;
    }

     /*--------------------------------------------------------------*
     | Otherwise need to B-orthogonalize the starting vector against |
     | the current Arnoldi basis using Gram-Schmidt with iter. ref.  |
     | This is the case where an invariant subspace is encountered   |
     | in the middle of the Arnoldi factorization.                   |
     |                                                               |
     |       s = V^{T}*B*r;   r = r - V*s;                           |
     |                                                               |
     | Stopping criteria used for iter. ref. is discussed in         |
     | Parlett's book, page 107 and in Gragg & Reichel TOMS paper.   |
     *--------------------------------------------------------------*/

    V->getv0_orth = 1;

LINE30:
    zgemv_(TRANS, &n, &j, &c1, v, &ldv, workd, &int1, &c0, &workd[n], &int1);
    TRANS = "N";
    zgemv_(TRANS, &n, &j, &cm1, v, &ldv, workd, &int1, &c1, &workd[n], &int1);

     /*---------------------------------------------------------*
     | Compute the B-norm of the orthogonalized starting vector |
     *---------------------------------------------------------*/

    if (V->bmat)
    {
        zcopy_(&n, resid, &int1, &workd[n], &int1);
        ipntr[0] = n;
        ipntr[1] = 0;
        V->ido = 2;
        return;
    } else {
        zcopy_(&n, resid, &int1, workd, &int1);
    }

LINE40:
    if (V->bmat)
    {
        *rnorm = cabsf(cdotc_(&n, resid, &int1, workd, &int1));
    } else {
        *rnorm = dznrm2_(&n, resid, &int1);
    }

     /*-------------------------------------*
     | Check for further orthogonalization. |
     *-------------------------------------*/

    if (*rnorm > sq2o2*V->getv0_rnorm0)
    {
        V->ido = 99;
        return;
    }

    V->getv0_iter += 1;
    if (V->getv0_iter < 1)
    {
         /*----------------------------------*
         | Perform iterative refinement step |
         *----------------------------------*/
        V->getv0_rnorm0 = *rnorm;
        goto LINE30;
    } else {
         /*-----------------------------------*
         | Iterative refinement step "failed" |
         *-----------------------------------*/

        for (jj = 0; jj < n; jj++) { resid[jj] = ARPACK_cplxf(0.0, 0.0); }
        *rnorm = 0.0;
        V->info = -1;
    }

    V->ido = 99;

    return;
}


void
csortc(const enum ARPACK_which w, const int apply, const int n,  ARPACK_CPLXF_TYPE *x,  ARPACK_CPLXF_TYPE *y)
{
    int i, igap, j;
    ARPACK_CPLXF_TYPE temp;
    ARPACK_compare_cfunc *f;

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

    igap = n / 2;

    while (igap != 0)
    {
        j = 0;
        for (i = igap; i < n; i++)
        {
            while (f(x[j], y[j+igap]))
            {
                if (j < 0) { break; }
                temp = x[j];
                x[j] = x[j+igap];
                x[j+igap] = temp;

                if (apply)
                {
                    temp = y[j];
                    y[j] = y[j+igap];
                    y[j+igap] = temp;
                }
                j -= igap;
            }
            j = i - igap + 1;
        }
        igap = igap / 2;
    }
}


int sortc_LM(const ARPACK_CPLXF_TYPE x, const ARPACK_CPLXF_TYPE y) { return (cabsf(x) > cabsf(y)); }
int sortc_SM(const ARPACK_CPLXF_TYPE x, const ARPACK_CPLXF_TYPE y) { return (cabsf(x) < cabsf(y)); }
int sortc_LR(const ARPACK_CPLXF_TYPE x, const ARPACK_CPLXF_TYPE y) { return (crealf(x) > crealf(y)); }
int sortc_SR(const ARPACK_CPLXF_TYPE x, const ARPACK_CPLXF_TYPE y) { return (crealf(x) < crealf(y)); }
int sortc_LI(const ARPACK_CPLXF_TYPE x, const ARPACK_CPLXF_TYPE y) { return (cimagf(x) > cimagf(y)); }
int sortc_SI(const ARPACK_CPLXF_TYPE x, const ARPACK_CPLXF_TYPE y) { return (cimagf(x) < cimagf(y)); }
