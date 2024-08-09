#include "__lbfgsb.h"
#include "../_build_utils/src/npy_cblas.h"
#include "../_build_utils/src/fortran_defs.h"


enum Status {
    START        = 0,
    NEW_X        = 1,
    RESTART      = 2,      // RESTART FROM LINESEARCH
    FG           = 3,      // REQUEST FG FOR THE CURRENT X
    CONVERGENCE  = 4,
    STOP         = 5,
    WARNING      = 6,
    ERROR        = 7,
    ABNORMAL     = 8       // ABNORMAL TERMINATION IN LINE SEARCH (WHY NOT AN ERROR?)
};


enum StatusMsg {
    NO_MSG       = 0,

    FG_START     = 301,
    FG_LNSRCH    = 302,

    CONV_GRAD    = 401,    // CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL
    CONV_F       = 402,    // CONVERGENCE: PROGRESS <= FACTR*EPSMCH

    STOP_CPU     = 501,    // CPU EXCEEDING THE TIME LIMIT ?!
    STOP_ITER    = 502,    // TOTAL NO. OF F,G EVALUATIONS EXCEEDS LIMIT
    STOP_GRAD    = 503,    // PROJECTED GRADIENT IS SUFFICIENTLY SMALL

    WARN_ROUND   = 601,    // ROUNDING ERRORS PREVENT PROGRESS
    WARN_STPMAX  = 602,    // STP = STPMAX
    WARN_STPMIN  = 603,    // STP = STPMIN
    WARN_XTOL    = 604,    // XTOL TEST SATISFIED ?!

    ERROR_NOFEAS = 701,    // NO FEASIBLE SOLUTION
    ERROR_FACTR  = 702,    // FACTR < 0
    ERROR_FTOL   = 703,    // FTOL < 0
    ERROR_GTOL   = 704,    // GTOL < 0
    ERROR_XTOL   = 705,    // XTOL < 0
    ERROR_STP1   = 706,    // STP < STPMIN
    ERROR_STP2   = 707,    // STP > STPMAX
    ERROR_STPMIN = 708,    // STPMIN < 0
    ERROR_STPMAX = 709,    // STPMAX < STPMIN
    ERROR_INITG  = 710,    // INITIAL G >= 0
    ERROR_M      = 711,    // M <= 0
    ERROR_N      = 712,    // N <= 0
    ERROR_NBD    = 713     // INVALID NBD
};


enum StatusWord {
    SUBS_THREE_DASH,
    SUBS_CONVERGED,
    SUBS_BOUNDED,
    SUBS_TRUNCATED
};


// Internal functions

static void mainlb(
    int n, int m, double* x, double* l, double* u,
    int* nbd, double* f, double* g, double factr, double pgtol,
    double* ws, double* wy, double* sy, double* ss, double* wt, double* wn,
    double* snd, double* z, double* r, double* d, double* t, double* xp,
    double* wa, int* index, int* iwhere, int* indx2, int* task,
    int* task_msg, int* lsave, int* isave, double* dsave, int maxls
);
static void active(
    int n, double* l, double* u, int* nbd, double* x,
    int* iwhere, int* prjctd, int* cnstnd, int* boxed
);
static void bmv(
    int m, double* sy, double* wt, int col,
    double* v, double* p, int* info
);
static void cauchy(
    int n, double* x, double* l, double* u,
    int* nbd, double *g, int* iorder, int* iwhere, double* t,
    double* d, double* xcp, int m, double* wy, double* ws,
    double* sy, double* wt, double theta, int col,
    int head, double* p, double* c, double* wbp, double* v, int* nseg,
    double sbgnrm, int* info
);
static void cmprlb(
    int n, int m, double* x, double* g,
    double* ws, double* wy, double* sy, double* wt,
    double* z, double* r, double* wa, int* index,
    double theta, int col, int head, int nfree,
    int cnstnd, int* info
);
static void errclb(
    int n, int m, double factr, double* l,
    double* u, int* nbd, int* task, int* task_msg, int* info, int* k
);
static void formk(
    int n, int nsub, int* ind, int nenter,
    int ileave, int* indx2, int iupdat, int updatd,
    double* wn, double* wn1, int m, double* ws, double* wy,
    double* sy, double theta, int col, int head,
    int* info
);
static void formt(
    int m, double* wt, double* sy, double* ss, int col,
    double theta, int* info
);
static void freev(
    int n, int* nfree, int* idx, int* nenter, int* ileave, int* idx2,
    int* iwhere, int* wrk, int updatd, int cnstnd,
    int iter
);
static void hpsolb(int n, double* t, int* iorder, int iheap);
static void lnsrlb(
    int n, double* l, double* u, int* nbd, double* x,
    double f, double* fold, double *gd, double *gdold, double* g,
    double* d, double* r, double* t, double* z, double* stp, double* dnorm,
    double* dtd, double* xstep, double* stpmx, int iter, int* ifun,
    int* iback, int* nfgv, int* info, int* task, int* task_msg, int boxed,
    int cnstnd, int* isave, double* dsave
);
static void matupd(
    int n, int m, double* ws, double *wy, double* sy, double* ss,
    double* d, double* r, int* itail, int iupdat, int* col,
    int* head, double* theta, double rr, double dr, double stp,
    double dtd
);
static void projgr(
    int n, double* l, double* u, int* nbd,
    double* x, double* g, double* sbgnrm
);
static void subsm(
    int n, int m, int nsub, int* ind, double* l,
    double* u, int* nbd, double* x, double* d, double* xp,
    double* ws, double* wy, double theta, double* xx,
    double* gg, int col, int head, int* iword, double* wv,
    double* wn, int* info
);

// Line search functions
static void dcsrch(
    double f, double g, double* stp, double ftol,
    double gtol, double xtol, double stpmin,
    double stpmax, int* task, int* task_msg, int* isave, double* dsave
);
static void dcstep (
    double* stx, double* fx, double* dx, double* sty, double* fy, double* dy,
    double* stp, double fp, double dp, int* brackt,
    double stpmin, double stpmax
);

static double epsmach = 2.220446049250313e-016;  /* np.finfo(np.float64).eps  */


void
setulb(int n, int m, double* x, double* l, double* u, int* nbd,
       double* f, double* g, double factr, double pgtol,
       double* wa, int* iwa, int* task, int* task_msg, int iprint,
       int* lsave, int* isave, double* dsave, int maxls)
{
    int lws, lr, lz, lt, ld, lxp, lwa, lwy, lsy, lss, lwt, lwn, lsnd;

    if (*task == START)
    {
        // Pointer olympics
        isave[0]  = m*n;
        isave[1]  = m*m;
        isave[2]  = 4*m*m;
        isave[3]  = 1;                      // ws      m*n
        isave[4]  = isave[3]  + isave[0];   // wy      m*n
        isave[5]  = isave[4]  + isave[0];   // wsy     m**2
        isave[6]  = isave[5]  + isave[1];   // wss     m**2
        isave[7]  = isave[6]  + isave[1];   // wt      m**2
        isave[8]  = isave[7]  + isave[1];   // wn      4*m**2
        isave[9]  = isave[8]  + isave[2];   // wsnd    4*m**2
        isave[10] = isave[9]  + isave[3];   // wz      n
        isave[11] = isave[10] + n;          // wr      n
        isave[12] = isave[11] + n;          // wd      n
        isave[13] = isave[12] + n;          // wt      n
        isave[14] = isave[13] + n;          // wxp     n
        isave[15] = isave[14] + n;          // wa      8*m
    }

    lws  = isave[3];
    lwy  = isave[4];
    lsy  = isave[5];
    lss  = isave[6];
    lwt  = isave[7];
    lwn  = isave[8];
    lsnd = isave[9];
    lz   = isave[10];
    lr   = isave[11];
    ld   = isave[12];
    lt   = isave[13];
    lxp  = isave[14];
    lwa  = isave[15];

    mainlb(n, m, x, l, u, nbd, f, g, factr, pgtol, &wa[lws], &wa[lwy], &wa[lsy],
           &wa[lss], &wa[lwt], &wa[lwn], &wa[lsnd], &wa[lz], &wa[lr], &wa[ld],
           &wa[lt], &wa[lxp], &wa[lwa], iwa, &iwa[n], &iwa[2*n], task, task_msg,
           lsave, &isave[21], dsave, maxls);

    return;
}


void
mainlb(int n, int m, double* x, double* l, double* u,
       int* nbd, double* f, double* g, double factr, double pgtol,
       double* ws, double* wy, double* sy, double* ss, double* wt, double* wn,
       double* snd, double* z, double* r, double* d, double* t, double* xp,
       double* wa, int* index, int* iwhere, int* indx2, int* task,
       int* task_msg, int* lsave, int* isave, double* dsave, int maxls)
{

    int prjctd, cnstnd, boxed, updatd, wrk;
    enum StatusWord word;
    int i, k, nintol, itfile, iback, nskip, head, col, iter, itail, iupdat;
    int nseg, nfgv, info, ifun, iword, nfree, nact, ileave, nenter;
    double theta, fold, dr, rr, tol, xstep, sbgnrm, stpmx, ddum, dnorm, dtd;
    double gd, gdold, stp;

    int one_int = 1, mut_n;

    if (*task == START) {

        // Initialize counters and scalars
        col = 0;
        head = 1;
        theta = 1.0;
        iupdat = 0;
        updatd = 0;
        iback = 0;
        itail = 0;
        iword = 0;
        nact = 0;
        ileave = 0;
        nenter = 0;
        fold = 0.0;
        dnorm = 0.0;
        gd = 0.0;
        sbgnrm = 0.0;
        stp = 0.0;
        gdold = 0.0;
        dtd = 0.0;

        // For operation counts
        iter = 0;
        nfgv = 0;
        nseg = 0;
        nintol = 0;
        nskip = 0;
        nfree = n;
        ifun = 0;
        // For stopping tolerance
        tol = factr * epsmach;

        // 'word' records the status of subspace solutions.
        word = SUBS_THREE_DASH;

        // 'info' records the termination information.
        info = 0;

        // Check the input arguments for errors
        errclb(n, m, factr, l, u, nbd, task, task_msg, &info, &k);
        if (*task == ERROR) { return; }

        // Initialize iwhere & project x onto the feasible set
        active(n, l, u, nbd, x, iwhere, &prjctd, &cnstnd, &boxed);

        // End of the initialization
    } else {
        // Restore local variables
        prjctd = lsave[0];
        cnstnd = lsave[1];
        boxed = lsave[2];
        updatd = lsave[3];

        nintol = isave[0];
        itfile = isave[1];
        iback = isave[2];
        nskip = isave[3];
        head = isave[4];
        col = isave[5];
        iter = isave[6];
        itail = isave[7];
        iupdat = isave[8];
        nseg = isave[9];
        nfgv = isave[10];
        info = isave[11];
        ifun = isave[12];
        iword = isave[13];
        nfree = isave[14];
        nact = isave[15];
        ileave = isave[16];
        nenter = isave[17];

        theta = dsave[0];
        fold = dsave[1];
        tol = dsave[2];
        dnorm = dsave[3];
        // epsmach = dsave[4];
        // cpu1 = dsave[5];
        // cachyt = dsave[6];
        // sbtime = dsave[7];
        // lnscht = dsave[8];
        // time1 = dsave[9];
        gd = dsave[10];
        stpmx = dsave[11];
        sbgnrm = dsave[12];
        stp = dsave[13];
        gdold = dsave[14];
        dtd = dsave[15];

    // After returning from the driver go to the point where execution is to resume.
        if ((*task == FG) && (*task_msg == FG_LNSRCH)) { goto LINE666; }
        if (*task == NEW_X) { goto LINE777; }
        if ((*task == FG) && (*task_msg == FG_START)) { goto LINE111; }
        if (*task == STOP)
        {
            if (*task_msg == STOP_CPU)
            {
                mut_n = n;
                DCOPY(&mut_n, t, &one_int, x, &one_int);
                DCOPY(&mut_n, r, &one_int, g, &one_int);
                *f = fold;
            }
            goto LINE999;
        }
    }

    // Compute f0 and g0.
    *task = FG;
    *task_msg = FG_START;

    // Return to the driver to calculate f and g, then reenter at 111
    goto LINE1000;

LINE111:
    nfgv = 1;

    // Compute the infinity norm of the (-) projected gradient.
    projgr(n, l, u, nbd, x, g, &sbgnrm);

    if (sbgnrm <= pgtol)
    {
        *task = CONVERGENCE;
        *task_msg = CONV_GRAD;
        goto LINE999;
    }

    // ----------------- the beginning of the loop --------------------------
LINE222:
    iword = -1;

    if ((!cnstnd) && (col > 0))
    {
        mut_n = n;
        DCOPY(&mut_n, x, &one_int, z, &one_int);
        wrk = updatd;
        nseg = 0;
        goto LINE333;
    }

    /////////////////////////////////////////////////////
    //
    // Compute the Generalized Cauchy Point (GCP).
    //
    /////////////////////////////////////////////////////
    cauchy(n, x, l, u, nbd, g, indx2, iwhere, t, d, z, m, wy, ws, sy, wt, theta,
           col, head, wa, &wa[2*m], &wa[4*m], &wa[6*m], &nseg, sbgnrm, &info);
    if (info != 0)
    {
        // Singular triangular system detected; refresh the lbfgs memory.
        info = 0;
        col = 0;
        head = 0;
        theta = 1.0;
        iupdat = 0;
        updatd = 0;
        goto LINE222;
    }

    nintol += nseg;

    // Count the entering and leaving variables for iter > 0;
    // find the index set of free and active variables at the GCP.
    freev(n, &nfree, index, &nenter, &ileave, indx2, iwhere, &wrk, updatd,
          cnstnd, iter);
    nact = n - nfree;

LINE333:
    // If there are no free variables or B=theta*I, then skip the subspace
    // minimization.
    if ((nfree == 0) || (col == 0)) { goto LINE555; }

    /////////////////////////////////////////////////////
    //
    // Subspace minimization
    //
    /////////////////////////////////////////////////////

    // Form  the LEL^T factorization of the indefinite
    //   matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
    //                 [L_a -R_z           theta*S'AA'S ]
    //   where     E = [-I  0]
    //                 [ 0  I]
    if (wrk)
    {
        formk(n, nfree, index, nenter, ileave, indx2, iupdat, updatd, wn, snd,
              m, ws, wy, sy, theta, col, head, &info);

    }
    if (info != 0)
    {
        // nonpositive definiteness in Cholesky factorization;
        // refresh the lbfgs memory and restart the iteration.
        info = 0;
        col = 0;
        head = 0;
        theta = 1.0;
        iupdat = 0;
        updatd = 0;
        goto LINE222;
    }

    // compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x) from 'cauchy').
    cmprlb(n, m, x, g, ws, wy, sy, wt, z, r, wa, index, theta, col, head, nfree,
           cnstnd, &info);
    if (info != 0) { goto LINE444; }

    // Call the direct method.
    subsm(n, m, nfree, index, l, u, nbd, z, r, xp, ws, wy, theta, x, g, col,
          head, &iword, wa, wn, &info);

LINE444:
    if (info != 0)
    {
        // singular triangular system detected;
        // refresh the lbfgs memory and restart the iteration.
        info = 0;
        col = 0;
        head = 0;
        theta = 1.0;
        iupdat = 0;
        updatd = 0;
        goto LINE222;
    }

LINE555:

    /////////////////////////////////////////////////////
    //
    // Line search and optimality tests.
    //
    /////////////////////////////////////////////////////

    // Generate the search direction d := z - x.
    for (i = 0; i < n; i++)
    {
        d[i] = z[i] - x[i];
    }

LINE666:
    lnsrlb(n, l, u, nbd, x, *f, &fold, &gd, &gdold, g, d, r, t, z, &stp, &dnorm,
           &dtd, &xstep, &stpmx, iter, &ifun, &iback, &nfgv, &info, task,
           task_msg, boxed, cnstnd, isave, dsave);
    if ((info != 0) || (iback >= maxls))
    {
        mut_n = n;
        DCOPY(&mut_n, t, &one_int, x, &one_int);
        DCOPY(&mut_n, r, &one_int, g, &one_int);
        *f = fold;

        if (col == 0)
        {
            // Abnormal termination
            if (info == 0)
            {
                info = -9;
                // Restore the actual number of f and g evaluations etc.
                nfgv--;
                ifun--;
                iback--;
            }
            *task = ABNORMAL;
            *task_msg = NO_MSG;
            iter++;
            goto LINE999;
        } else {
            // Refresh the lbfgs memory and restart the iteration.
            if (info == 0) { nfgv--; }
            info = 0;
            col = 0;
            head = 0;
            theta = 1.0;
            iupdat = 0;
            updatd = 0;
            *task = RESTART;
            *task_msg = NO_MSG;
            goto LINE222;
        }
    } else if ((*task == FG) && (*task_msg == FG_LNSRCH))
    {
        // Return to the driver for calculating f and g; renter at 666.
        goto LINE1000;
    } else {
        // Calculate and print out the quantities related to the new X.
        iter++;

        // Compute the infinity norm of the projected (-)gradient.
        projgr(n, l, u, nbd, x, g, &sbgnrm);

        if (iword == 0)
        {
            word = SUBS_CONVERGED;
        } else if (iword == 1)
        {
            word = SUBS_BOUNDED;
        } else if (iword == 5)
        {
            word = SUBS_TRUNCATED;
        } else {
            word = SUBS_THREE_DASH;
        }

        goto LINE1000;
    }

LINE777:

    // Test for termination.
    if (sbgnrm <= pgtol)
    {
        // Terminate the algorithm.
        *task = CONVERGENCE;
        *task_msg = CONV_GRAD;
        goto LINE999;
    }
    ddum = fmax(fmax(fabs(fold), fabs(*f)), 1.0);
    if ((fold - *f) <= tol*ddum)
    {
        // Terminate the algorithm.
        *task = CONVERGENCE;
        *task_msg = CONV_F;
        if (iback >= 10) { info = -5; }
        // i.e. to issue a warning if iback > 10 in the line search.
        goto LINE999;
    }

    // Compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's.
    for (i = 0; i < n; i++)
    {
        r[i] = g[i] - r[i];
    }
    // 42
    mut_n = n;
    rr = pow(DNRM2(&mut_n, r, &one_int), 2.0);

    if (stp == 1.0)
    {
        dr = gd - gdold;
        ddum = -gdold;
    } else {
        dr = (gd - gdold)*stp;
        DSCAL(&mut_n, &stp, d, &one_int);
        ddum = -gdold*stp;
    }

    if (dr <= epsmach*ddum)
    {
        // Skip the L-BFGS update.
        nskip += 1;
        updatd = 0;
        goto LINE888;
    }

    /////////////////////////////////////////////////////
    //
    // Update the L-BFGS matrix.
    //
    /////////////////////////////////////////////////////
    updatd = 1;
    iupdat += 1;

    // Update matrices WS and WY and form the middle matrix in B.
    matupd(n, m, ws, wy, sy, ss, d, r, &itail, iupdat, &col, &head, &theta,
           rr, dr, stp, dtd);

    // Form the upper half of the pds T = theta*SS + L*D^(-1)*L';
    //    Store T in the upper triangular of the array wt;
    //    Cholesky factorize T to J*J' with
    //       J' stored in the upper triangular of wt.
    formt(m, wt, sy, ss, col, theta, &info);

    if (info != 0)
    {
        // Nonpositive definiteness in Cholesky factorization;
        // refresh the lbfgs memory and restart the iteration.
        info = 0;
        col = 0;
        head = 0;
        theta = 1.0;
        iupdat = 0;
        updatd = 0;
        goto LINE222;
    }

    // Now the inverse of the middle matrix in B is

    // [  D^(1/2)      O ] [ -D^(1/2)  D^(-1/2)*L' ]
    // [ -L*D^(-1/2)   J ] [  0        J'          ]
LINE888:
    // -------------------- the end of the loop -----------------------------
    goto LINE222;

LINE999:
    ;

LINE1000:
    // Save local variables
    lsave[0]  = prjctd;
    lsave[1]  = cnstnd;
    lsave[2]  = boxed;
    lsave[3]  = updatd;

    isave[0]  = nintol;
    isave[2]  = itfile;
    isave[3]  = iback;
    isave[4]  = nskip;
    isave[5]  = head;
    isave[6]  = col;
    isave[7]  = itail;
    isave[8]  = iter;
    isave[9] = iupdat;
    isave[11] = nseg;
    isave[12] = nfgv;
    isave[13] = info;
    isave[14] = ifun;
    isave[15] = iword;
    isave[16] = nfree;
    isave[17] = nact;
    isave[18] = ileave;
    isave[19] = nenter;

    dsave[0]  = theta;
    dsave[1]  = fold;
    dsave[2]  = tol;
    dsave[3]  = dnorm;
    // dsave[4]  = epsmch;
    // dsave[5]  = cpu1;
    // dsave[6]  = cachyt;
    // dsave[7]  = sbtime;
    // dsave[8]  = lnscht;
    // dsave[9] = time1;
    dsave[10] = gd;
    dsave[11] = stpmx;
    dsave[12] = sbgnrm;
    dsave[13] = stp;
    dsave[14] = gdold;
    dsave[15] = dtd;

    return;
}


void
active(int n, double* l, double* u, int* nbd, double* x,
       int* iwhere, int* prjctd, int* cnstnd, int* boxed)
{
    int nbdd, i;

    // Initialize nbdd, prjctd, cnstnd, and boxed.
    nbdd = 0;
    *prjctd = 0;
    *cnstnd = 0;
    *boxed = 1;

    // Project the initial x to the feasible set if necessary.

    for (i = 0; i < n; i++)
    {
        if (nbd[i] > 0)
        {
            if ((nbd[i] <= 2) && (x[i] <= l[i]))
            {
                if (x[i] < l[i])
                {
                    *prjctd = 1;
                    x[i] = l[i];
                }
                nbdd += 1;
            } else if ((nbd[i] >= 2) && (x[i] >= u[i])) {
                if (x[i] > u[i])
                {
                    *prjctd = 1;
                    x[i] = u[i];
                }
                nbdd += 1;
            }
        }
    }
    // 10

    // Initialize iwhere and assign values to cnstnd and boxed.
    for (i = 0; i < n; i++)
    {
        if (nbd[i] != 2) { boxed = 0; }
        if (nbd[i] == 0)
        {
            // This variable is always free
            iwhere[i] = -1;
        } else {
            // Otherwise set x(i)=mid(x(i), u(i), l(i)).
            *cnstnd = 1;
            if ((nbd[i] == 2) && (u[i] - l[i] <= 0.0))
            {
                // This variable is always fixed.
                iwhere[i] = 3;
            } else {
                iwhere[i] = 0;
            }
        }
    }
    // 20

    return;
}


void
bmv(int m, double* sy, double* wt, int col,
    double* v, double* p, int* info)
{
    int i, k, i2, one_int = 1;
    double ssum;
    char *uplo = "U", *trans = "T", *diag = "N";

    if (col == 0) { return; }

    // PART I: solve [  D^(1/2)      O ] [ p1 ] = [ v1 ]
    //               [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ].

    // solve Jp2=v2+LD^(-1)v1.
    p[col] = v[col];
    for (i = 1; i < col; i++)
    {
        i2 = col + i;
        ssum = 0.0;
        for (k = 0; k < i-1; k++)
        {
            ssum = ssum + sy[i + m*k] * v[k] / sy[k + m*k];
        }
        // 10
        p[i2] = v[i2] + ssum;
    }
    // 20

    // Solve the triangular system
    // dtrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info)
    DTRTRS(uplo, trans, diag, &col, &one_int, wt, &m, &p[col], &col, info);
    if (*info != 0) { return; }

    // Solve D^(1/2)p1=v1.
    for (i = 0; i < col; i++)
    {
        p[i] = v[i] / sqrt(sy[i + m*i]);
    }
    // 30

    // PART II: solve [ -D^(1/2)   D^(-1/2)*L'  ] [ p1 ] = [ p1 ]
    //                [  0         J'           ] [ p2 ]   [ p2 ].

    // solve J^Tp2=p2.
    // dtrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info)
    DTRTRS(uplo, trans, diag, &col, &one_int, wt, &m, &p[col], &col, info);
    if (*info != 0) { return; }

    // compute p1=-D^(-1/2)(p1-D^(-1/2)L'p2)
    //           =-D^(-1/2)p1+D^(-1)L'p2.

    for (i = 0; i < col; i++)
    {
        p[i] = -p[i] / sqrt(sy[i + m*i]);
    }
    // 40

    for (i = 0; i < col; i++)
    {
        ssum = 0.0;
        for (k = i+1; k < col; k++)
        {
            ssum = ssum + sy[k + m*i] * p[col + k] / sy[i + m*i];
        }
        // 50
        p[i] = p[i] + ssum;
    }
    // 60

    return;
}


void
cauchy(int n, double* x, double* l, double* u,
       int* nbd, double *g, int* iorder, int* iwhere, double* t,
       double* d, double* xcp, int m, double* wy, double* ws,
       double* sy, double* wt, double theta, int col,
       int head, double* p, double* c, double* wbp, double* v,
       int* nseg, double sbgnrm, int* info)
{
    int bnded, xlower, xupper, mut_n;
    int i, j, col2, nfree, nbreak, pointr, ibp, nleft, ibkmin, iter, one_int = 1;
    double f1, f2, dt, dtm, tsum, dibp, zibp, dibp2, bkmin, tu, tl, wmc, wmp;
    double wmw, tj, tj0, neggi, f2_org;

    if (sbgnrm <= 0.0)
    {
        mut_n = n;
        DCOPY(&mut_n, x, &one_int, xcp, &one_int);
        return;
    }
    bnded = 1;
    nfree = n;
    nbreak = -1;
    ibkmin = 0;
    bkmin = 0.0;
    col2 = 2*col;

    // We set p to zero and build it up as we determine d.
    for (i = 0; i < col2; i++)
    {
        p[i] = 0.0;
    }
    // 20

    // In the following loop we determine for each variable its boud status and
    // its breakpoint, and update p accordingly. Smallest breakpoint is identified.
    for (i = 0; i < n; i++)
    {
        neggi = -g[i];
        if ((iwhere[i] != 3) && (iwhere[i] != -1))
        {
            // If x[i] is not a constant and has bounds, compute the difference
            // between x[i] and its bounds.
            if (nbd[i] <= 2) { tl = x[i] - l[i]; }
            if (nbd[i] >= 2) { tu = u[i] - x[i]; }

            // If a varaible is close enough to a bound, we treat it as at bound.
            xlower = ((nbd[i] <= 2) && (tl <= 0.0));
            xupper = ((nbd[i] >= 2) && (tu <= 0.0));

            // Reset iwhere[i]
            iwhere[i] = 0;
            if (xlower) {
                if (neggi <= 0.0) { iwhere[i] = 1; }
            } else if (xupper) {
                if (neggi >= 0.0) { iwhere[i] = 2; }
            } else {
                if (fabs(neggi) <= 0.0) { iwhere[i] = -3; }
            }
        }
        pointr = head;
        if ((iwhere[i] != 0) && (iwhere[i] != -1))
        {
            d[i] = 0.0;
        } else {
            d[i] = neggi;
            f1 = f1 - neggi * neggi;
            // Calculate p := p - W'e_i* (g_i).
            for (j = 0; j < col; j++) {
                p[j] = p[j] + wy[i + n*pointr]*neggi;
                p[col + j] = p[col + j] + ws[i + n*pointr]*neggi;
                pointr = pointr % m;
            }
            // 40
            if ((nbd[i] <= 2) && (nbd[i] != 0) && (neggi < 0.0))
            {
                // x(i) + d(i) is bounded; compute t(i).
                nbreak += 1;
                iorder[nbreak] = i;
                t[nbreak] = tl / (-neggi);
                if ((nbreak == 1) || (t[nbreak] < bkmin))
                {
                    bkmin = t[nbreak];
                    ibkmin = nbreak;
                }
            } else if ((nbd[i] >= 2) && (neggi > 0.0)) {
                // x(i) + d(i) is bounded; compute t(i).
                nbreak += 1;
                iorder[nbreak] = i;
                t[nbreak] = tu / neggi;
                if ((nbreak == 1) || (t[nbreak] < bkmin))
                {
                    bkmin = t[nbreak];
                    ibkmin = nbreak;
                }
            } else {
                // x(i) + d(i) is not bounded.
                nfree -= 1;
                iorder[nfree] = i;
                if (fabs(neggi) > 0.0) { bnded = 0; }
            }
        }
    }
    // 50

    // The indices of the nonzero components of d are now stored
    // in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n).
    // The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin.
    if (theta != 1.0)
    {
        // Complete the initialization of p for theta not= one.
        DSCAL(&col, &theta, &p[col], &one_int);
    }
    // Initialize GCP xcp = x.
    mut_n = n;
    DCOPY(&mut_n, x, &one_int, xcp, &one_int);
    if ((nbreak == 0) && (nfree == n))
    {
        // is a zero vector, return with the initial xcp as GCP.
        return;
    }
    // Initialize c = W'(xcp - x) = 0.
    for (j = 0; j < col2; j++)
    {
        c[j] = 0.0;
    }
    // 60

    // Initialize derivative f2.
    f2 = -(theta) * f1;
    f2_org = f2;
    if (col > 0) {
        bmv(m, sy, wt, col, p, v, info);
        if (*info != 0) { return; }
        f2 = f2 - DDOT(&col2, v, &one_int, p, &one_int);
    }
    dtm = -f1 / f2;
    tsum = 0.0;
    *nseg = 1;

    // If there are no breakpoints, locate the GCP and return.
    if (nbreak == -1) { goto LINE888; }
    nleft = nbreak;
    iter = 0;
    tj = 0.0;

    // ------------------- the beginning of the loop -------------------------
    while (1)
    {
        // Find the next smallest breakpoint;
        // compute dt = t(nleft) - t(nleft + 1).
        tj0 = tj;
        if (iter == 0)
        {
            // Since we already have the smallest breakpoint we need not do
            // heapsort yet. Often only one breakpoint is used and the
            // cost of heapsort is avoided.
            tj = bkmin;
            ibp = iorder[ibkmin];
        } else {
            if (iter == 1)
            {
                // Replace the already used smallest breakpoint with the
                // breakpoint numbered nbreak > nlast, before heapsort call.
                if (ibkmin != nbreak)
                {
                    t[ibkmin] = t[nbreak];
                    iorder[ibkmin] = iorder[nbreak];
                }
            }
            // Update heap structure of breakpoints (if iter=2, initialize heap).
            hpsolb(nleft, t, iorder, iter - 1);
            tj = t[nleft];
            ibp = iorder[nleft];
        }

        dt = tj - tj0;

        // If a minimizer is within this interval, locate the GCP and return.
        if (dtm < dt) { goto LINE888; }

        // Otherwise fix one variable and reset the corresponding component
        // of d to zero.
        tsum = tsum + dt;
        nleft -= 1;
        iter += 1;
        dibp = d[ibp];
        d[ibp] = 0.0;
        if (dibp > 0.0) {
            zibp = u[ibp] - x[ibp];
            xcp[ibp] = u[ibp];
            iwhere[ibp] = 2;
        } else {
            zibp = l[ibp] - x[ibp];
            xcp[ibp] = l[ibp];
            iwhere[ibp] = 1;
        }
        if (nleft == -1 && nbreak == n) {
            // All n variables are fixed, return with xcp as GCP.
            dtm = dt;
            goto LINE999;
        }

        // Update the derivative information.
        nseg += 1;
        dibp2 = pow(dibp, 2.0);

        // Update f1 and f2

        // Temporarily set f1 and f2 for col=0.
        f1 = f1 + dt * f2 + dibp2 - theta*dibp*zibp;
        f2 = f2 - theta*dibp2;
        if (col > 0) {
            // Update c = c + dt*p.
            DAXPY(&col2, &dt, p, &one_int, c, &one_int);
            // Choose wbp, the row of W corresponding to the breakpoint encountered.
            pointr = head;

            for (j = 0; j < col; j++) {
                wbp[j] = wy[ibp + n*pointr];
                wbp[col + j] = theta * ws[ibp + n*pointr];
                pointr = pointr % m;
            }
            // 70

            // Compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'.
            bmv(m, sy, wt, col, wbp, v, info);
            if (*info != 0) { return; }
            wmc = DDOT(&col2, c, &one_int, v, &one_int);
            wmp = DDOT(&col2, p, &one_int, v, &one_int);
            wmw = DDOT(&col2, wbp, &one_int, v, &one_int);
            // Update p = p - dibp*wbp.
            dibp = -dibp;
            DAXPY(&col2, &dibp, wbp, &one_int, p, &one_int);
            dibp = -dibp;

            // Complete updating f1 and f2 while col > 0.
            f1 = f1 + dibp*wmc;
            f2 = f2 + dibp*2.0*wmp - dibp2*wmw;
        }

        f2 = fmax(epsmach*f2_org, f2);

        if (nleft >= 0) {
            // to repeat the loop for unsearched intervals.
            dtm = -f1 / f2;
        } else if (bnded) {
            f1 = 0.0;
            f2 = 0.0;
            dtm = 0.0;
            break;
        } else {
            dtm = -f1 / f2;
            break;
        }
    }

    // ------------------- the end of the loop -------------------------------

LINE888:
    if (dtm <= 0.0) { dtm = 0.0; }
    tsum = tsum + dtm;

    // Move free variables (i.e. the ones w/o breakpoints) and the variables
    // whose breakpoints haven't been reached.
    DAXPY(&n, &tsum, d, &one_int, xcp, &one_int);

LINE999:
    // Update c = c + dtm*p = W'(x^c - x)
    // which will be used in computing r = Z'(B(x^c - x) + g).

    if (col > 0) { DAXPY(&col2, &dtm, p, &one_int, c, &one_int); }

    return;
}


void
cmprlb(int n, int m, double* x, double* g,
       double* ws, double* wy, double* sy, double* wt,
       double* z, double* r, double* wa, int* index,
       double theta, int col, int head, int nfree,
       int cnstnd, int* info)
{
    int i, j, k, pointr;
    double a1, a2;

    if ((!cnstnd) && (col > 0))
    {
        for (i = 0; i < n; i++) { r[i] = -g[i]; }
        // 26
    } else {
        for (i = 0; i < nfree; i++)
        {
            k = index[i];
            r[i] = -theta*(z[k] - x[k]) - g[k];
        }
        // 30

        bmv(m, sy, wt, col, &wa[2*m], wa, info);
        if (*info != 0) { *info = -8; return; }
        pointr = head;
        for (j = 0; j < col; j++)
        {
            a1 = wa[j];
            a2 = theta*wa[col + j];
            for (i = 0; i < nfree; i++)
            {
                k = index[i];
                r[i] = r[i] + wy[k + n*pointr]*a1 + ws[k + n*pointr]*a2;
            }
            // 32
            pointr = pointr % m;
        }
        // 34
    }

    return;
}


void
errclb(int n, int m, double factr, double* l,
       double* u, int* nbd, int* task, int* task_msg, int* info,
       int* k)
{
    int i;

    // Check the input arguments for errors.
    if (n <= 0) { *task = ERROR; *task_msg = ERROR_N; }
    if (m <= 0) { *task = ERROR; *task_msg = ERROR_M; }
    if (factr < 0) { *task = ERROR; *task_msg = ERROR_FACTR; }

    // Check the validity of the arrays nbd(i), u(i), and l(i).
    for (i = 0; i < n; i++)
    {
        if ((nbd[i] < 0) || (nbd[i] > 3))
        {
            *task = ERROR;
            *task_msg = ERROR_NBD;
            *info = -6;
            *k = i;
        }
        if (nbd[i] == 2)
        {
            if (l[i] > u[i])
            {
                // Return
                *task = ERROR;
                *task_msg = ERROR_NOFEAS;
                *info = -7;
                *k = i;
            }
        }
    }

    return;
}


void
formk(int n, int nsub, int* ind, int nenter, int ileave,
      int* indx2, int iupdat, int updatd, double* wn,
      double* wn1, int m, double* ws, double* wy,
      double* sy, double theta, int col, int head, int* info)
{
    int m2, ipntr, jpntr, iy, is, jy, js, is1, js1, k1, i, k, col2, pbegin, pend;
    int dbegin, dend, upcl, one_int = 1, temp_int;
    double temp1, temp2, temp3, temp4;
    char *uplo = "U", *trans = "T", *diag = "N";

    // Form the lower triangular part of
    //           WN1 = [Y' ZZ'Y   L_a'+R_z']
    //                 [L_a+R_z   S'AA'S   ]
    //    where L_a is the strictly lower triangular part of S'AA'Y
    //          R_z is the upper triangular part of S'ZZ'Y.

    if (updatd)
    {
        if (iupdat > m)
        {
            //  Shift old part of WN1.
            for (jy = 0; jy < m - 1; jy++)
            {
                js = m + jy;
                temp_int = m - jy + 1;
                DCOPY(&temp_int, &wn1[(jy + 1) + 2*m*(jy + 1)], &one_int, &wn1[jy + 2*m*jy], &one_int);
                DCOPY(&temp_int, &wn1[(js + 1) + 2*m*(js + 1)], &one_int, &wn1[js + 2*m*js], &one_int);
                temp_int = m - 1;
                DCOPY(&temp_int, &wn1[(m + 1) + 2*m*(jy + 1)], &one_int, &wn1[m + 2*m*jy], &one_int);
            }
            // 10
        }

        // Put new rows in blocks (1,1), (2,1) and (2,2).
        pbegin = 0;     // Inclusive loop start
        pend = nsub;    // Exclusive loop end
        dbegin = nsub;  // Inclusive loop start
        dend = n;       // Exclusive loop end
        iy = col - 1;
        is = m + col - 1;
        ipntr = head + col - 1;
        if (ipntr >= m) { ipntr -= m; }
        jpntr = head;
        for (jy = 0; jy < col; jy++)
        {
            js = m + jy;
            temp1 = 0.0;
            temp2 = 0.0;
            temp3 = 0.0;

            // Compute element jy of row 'col' of Y'ZZ'Y.
            for (k = pbegin; k < pend; k++)
            {
                k1 = ind[k];
                temp1 = temp1 + wy[k1 + n*ipntr]*wy[k1 + n*jpntr];
            }
            // 15

            // Compute elements jy of row 'col' of L_a and S'AA'S.
            for (k = dbegin; k < dend; k++)
            {
                k1 = ind[k];
                temp2 = temp2 + ws[k1 + n*ipntr]*ws[k1 + n*jpntr];
                temp3 = temp3 + ws[k1 + n*ipntr]*wy[k1 + n*jpntr];
            }
            // 16
            wn1[iy + 2*m*jy] = temp1;
            wn1[is + 2*m*js] = temp2;
            wn1[is + 2*m*jy] = temp3;
            jpntr = (jpntr % m);
        }
        // 20

        // Put new column in block (2,1)
        jy = col - 1;
        jpntr = head + col - 1;
        if (jpntr >= m) { jpntr -= m; }
        ipntr = head;
        for (i = 0; i < col; i++)
        {
            is = m + i;
            temp3 = 0.0;
            // Compute element i of column 'col' of R_z.
            for (k = pbegin; k < pend; k++)
            {
                k1 = ind[k];
                temp3 = temp3 + ws[k1 + n*ipntr]*wy[k1 + n*jpntr];
            }
            // 25
            ipntr = (ipntr % m);
            wn1[is + 2*m*jy] = temp3;
        }
        // 30

        upcl = col - 1;
    } else {
        upcl = col;
    }

    // modify the old parts in blocks (1,1) and (2,2) due to changes
    // in the set of free variables.
    ipntr = head;
    for (iy = 0; iy < upcl; iy++)
    {
        is = m + iy;
        jpntr = head;
        for (jy = 0; jy < iy; jy++)
        {
            js = m + jy;
            temp1 = 0.0;
            temp2 = 0.0;
            temp3 = 0.0;
            temp4 = 0.0;

            for (k = 0; k < nenter; k++)
            {
                k1 = indx2[k];
                temp1 = temp1 + wy[k1 + n*ipntr]*wy[k1 + n*jpntr];
                temp2 = temp2 + ws[k1 + n*ipntr]*ws[k1 + n*jpntr];
            }
            // 35
            for (k = ileave; k < n; k++)  // CHECK THIS LINE MIGHT BE ileave - 1
            {
                k1 = indx2[k];
                temp3 = temp3 + wy[k1 + n*ipntr]*wy[k1 + n*jpntr];
                temp4 = temp4 + ws[k1 + n*ipntr]*ws[k1 + n*jpntr];
            }
            // 36

            wn1[iy + 2*m*jy] = wn1[iy + 2*m*jy] + temp1 - temp3;
            wn1[is + 2*m*js] = wn1[is + 2*m*js] - temp2 + temp4;
            jpntr = (jpntr % m);
        }
        // 40

        ipntr = (ipntr % m);
    }
    // 45

    // Modify the old parts in block (2,1)
    ipntr = head;
    for (is = m; is < m + upcl; is++)
    {
        jpntr = head;
        for (jy = 0; jy < upcl; jy++)
        {
            temp1 = 0.0;
            temp3 = 0.0;
            for (k = 0; k < nenter; k++)
            {
                k1 = indx2[k];
                temp1 = temp1 + ws[k1 + n*ipntr]*wy[k1 + n*jpntr];
            }
            // 50
            for (k = ileave; k < n; k++)
            {
                k1 = indx2[k];
                temp3 = temp3 + ws[k1 + n*ipntr]*wy[k1 + n*jpntr];
            }
            // 51
            if (is <= jy + m)
            {
                wn1[is + 2*m*jy] = wn1[is + 2*m*jy] + temp1 - temp3;
            } else {
                wn1[is + 2*m*jy] = wn1[is + 2*m*jy] - temp1 + temp3;
            }
            jpntr = (jpntr % m);
        }
        // 55
        ipntr = (ipntr % m);
    }
    // 60

    // Form the upper triangle of WN = [D+Y' ZZ'Y/theta   -L_a'+R_z' ]
    //                                 [-L_a +R_z        S'AA'S*theta]

    m2 = 2*m;
    for (iy = 0; iy < col; iy++)
    {
        is = col + iy;
        is1 = m + iy;
        for (jy = 0; jy < iy; jy++)
        {
            js = col + jy;
            js1 = m + jy;
            wn[jy + 2*m*iy] = wn1[iy + 2*m*jy]/theta;
            wn[js + 2*m*is] = wn1[is1 + 2*m*js1]*theta;
        }
        // 65
        for (jy = 0; jy < iy-1; jy++)
        {
            wn[jy + 2*m*is] = -wn1[is1 + 2*m*jy];
        }
        // 66
        for (jy = iy-1; jy < col; jy++)
        {
            wn[jy + 2*m*is] = wn1[is1 + 2*m*jy];
        }
        // 67
        wn[iy + 2*m*iy] = wn[iy + 2*m*iy] + sy[iy + m*iy];
    }
    // 70

    // Form the upper triangle of WN= [  LL'            L^-1(-L_a'+R_z')]
    //                                [(-L_a +R_z)L'^-1   S'AA'S*theta  ]
    //
    //    first Cholesky factor (1,1) block of wn to get LL'
    //                      with L' stored in the upper triangle of wn.

    // dpotrf(uplo, n, a, lda, info)
    DPOTRF(uplo, &col, wn, &m2, info);
    if (*info != 0) { *info = -1; return; }

    // Then form L^-1(-L_a'+R_z') in the (1,2) block.
    col2 = 2*col;

    // dtrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info)
    DTRTRS(uplo, trans, diag, &col, &col, wn, &m2, &wn[2*m*col], &col, info);

    // Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the upper
    //  triangle of (2,2) block of wn.
    for (is = col; is < col2; is++)
    {
        for (js = is; js < col2; js++)
        {
            wn[is + 2*m*js] = wn[is + 2*m*js] + DDOT(&col, &wn[2*m*is], &one_int, &wn[2*m*js], &one_int);
        }
        // 74
    }
    // 72

    // dpotrf(uplo, n, a, lda, info)
    DPOTRF(uplo, &col, wn[col+1 + 2*m*(col+1)], &m2, info);
    if (*info != 0) { *info = -2; }

    return;
}


void
formt(int m, double* wt, double* sy, double* ss, int col,
      double theta, int* info)
{
    int i, j, k, k1;
    double ddum;
    char *uplo = "U";
    // Form the upper half of  T = theta*SS + L*D^(-1)*L', store T in the upper
    // triangle of the array wt.
    for (j = 0; j < col; j++)
    {
        wt[m*j] = theta*ss[m*j];
    }
    // 52

    for (i = 1; i < col; i++)
    {
        for (j = i; j < col; j++)
        {
            k1 = (i > j ? j : i);  // CHECK THIS LINE, REMOVED -1
            ddum = 0.0;
            for (k = 0; k < k1; k++)
            {
                ddum = ddum + sy[i + m*k] * sy[j + m*k] / sy[k + m*k];
            }
            // 53
            wt[i + m*j] = ddum + theta*ss[i + m*j];
        }
        // 54
    }
    // 55

    // Cholesky factorize T to J*J' with J' stored in the upper triangle of wt.
    // dpotrf(uplo, n, a, lda, info)
    DPOTRF(uplo, &col, wt, &m, &info);
    if (info != 0) { *info = -3; }

    return;
}


void
freev(int n, int* nfree, int* idx, int* nenter, int* ileave, int* idx2,
      int* iwhere, int* wrk, int updatd, int cnstnd,
      int iter)
{
    int i, iact, k;

    *nenter = -1;
    *ileave = n;
    if ((iter > 0) && cnstnd)
    {
        // Count the entering and leaving variables.
        for (i = 0; i < *nfree; i++)
        {
            k = idx[i];

            if (iwhere[k] > 0)
            {
                *ileave -= 1;
                idx2[*ileave] = k;
            }
        }
        // 20

        for (i = *nfree; i < n; i++)
        {
            k = idx[i];
            if (iwhere[k] <= 0)
            {
                *nenter += 1;
                idx2[*nenter] = k;
            }
        }
        // 22
    }

    *wrk = ((*ileave < n) || (*nenter > -1) || updatd);

    // Find the index set of free and active variables at the GCP.
    *nfree = -1;
    iact = n;
    for (i = 0; i < n; i++)
    {
        if (iwhere[i] <= 0)
        {
            nfree += 1;
            idx[*nfree] = i;
        } else {
            iact -= 1;
            idx[iact] = i;
        }
    }
    // 24

    return;
}


void
hpsolb(int n, double* t, int* iorder, int iheap)
{
    int i, j, k, indxin, indxout;
    double ddum, out;

    if (iheap == 0)
    {
        // Rearrange the elements t(1) to t(n) to form a heap.
        for (k = 2; k <= n + 1; k++)
        {
            ddum = t[k-1];
            indxin = iorder[k-1];

            // Add ddum to the heap
            i = k;
            while (1)
            {
                if (i > 1)
                {
                    j = i / 2;
                    if (ddum < t[j-1])
                    {
                        t[i-1] = t[j-1];
                        iorder[i-1] = iorder[j-1];
                        i = j;
                    } else {
                        break;
                    }
                } else {
                    break;
                }
            }
            // 10

            t[i-1] = ddum;
            iorder[i-1] = indxin;
        }
        // 20
    }

    // Assign to 'out' the value of t(1), the least member of the heap, and
    // rearrange the remaining members to form a heap as elements 1 to n-1 of t.
    if (n > 0)
    {
        i = 1;
        out = t[0];
        indxout = iorder[0];
        ddum = t[n-1];
        indxin = iorder[n-1];

        // Restore the heap
        while (1)
        {
            j = i + i;
            if (j <= n-1)
            {
                if (t[j] < t[j-1]) { j++; }
                if (t[j-1] < ddum)
                {
                    t[i-1] = t[j-1];
                    iorder[i-1] = iorder[j-1];
                    i = j;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
        t[i-1] = ddum;
        iorder[i-1] = indxin;

        // Put the least member in t(n).
        t[n-1] = out;
        iorder[n-1] = indxout;
    }

    return;
}


void
lnsrlb(int n, double* l, double* u, int* nbd, double* x,
       double f, double* fold, double *gd, double *gdold, double* g,
       double* d, double* r, double* t, double* z, double* stp, double* dnorm,
       double* dtd, double* xstep, double* stpmx, int iter, int* ifun,
       int* iback, int* nfgv, int* info, int* task, int* task_msg,
       int boxed, int cnstnd, int* isave, double* dsave)
{
    int i, temp_task, temp_taskmsg, one_int = 1;
    double a1, a2, ftol = 1e-3, gtol = 0.9, xtol = 0.1;

    if (*task_msg == FG_LNSRCH) { goto LINE556; }

    *dnorm = DNRM2(&n, d, &one_int);
    *dtd = pow(*dnorm, 2.0);

    // Determine the maximum step length.
    *stpmx = 1.0e+10;
    if (cnstnd)
    {
        if (iter == 0)
        {
            *stpmx = 1.0;
        } else {
            for (i = 0; i < n; i++)
            {
                a1 = d[i];
                if (nbd[i] != 0)
                {
                    if ((a1 < 0.0) && (nbd[i] <= 2))
                    {
                        a2 = l[i] - x[i];
                        if (a2 >= 0.0)
                        {
                            *stpmx = 0.0;
                        } else if (a1*(*stpmx) < a2) {
                            *stpmx = a2 / a1;
                        }
                    } else if ((a1 > 0.0) && (nbd[i] >= 2)) {
                        a2 = u[i] - x[i];
                        if (a2 <= 0.0)
                        {
                            *stpmx = 0.0;
                        } else if (a1*(*stpmx) > a2) {
                            *stpmx = a2 / a1;
                        }
                    }
                }
            }
            // 43
        }
    }

    if ((iter == 0) && (!(boxed)))
    {
        *stp = fmin(1.0 / (*dnorm), *stpmx);
    } else {
        *stp = 1.0;
    }

    for (i = 0; i < n; i++)
    {
        t[i] = x[i];
        r[i] = g[i];
    }

    *fold = f;  // What's the point of this?
    ifun = 0;
    iback = 0;
    temp_task = START;
    temp_taskmsg = NO_MSG;
LINE556:
    *gd = DDOT(&n, g, &one_int, d, &one_int);
    if (ifun == 0)
    {
        *gdold = *gd;
        if (*gd >= 0.0)
        {
            // The directional derivative >= 0, line search is impossible.
            *info = -4;
            return;
        }
    }

    dcsrch(f, *gd, stp, ftol, gtol, xtol, 0.0, *stpmx, &temp_task, &temp_taskmsg,
           isave, dsave);

    *xstep = *stp * (*dnorm);
    if ((temp_task != CONVERGENCE) && (temp_task != WARNING))
    {
        *task = FG;
        *task_msg = FG_LNSRCH;
        *ifun += 1;
        *nfgv += 1;
        iback = ifun - 1;
        if (*stp == 1.0)
        {
            // for (i = 0; i < n; i++) { x[i] = z[i]; }
            DCOPY(&n, z, &one_int, x, &one_int);
        } else {
            for (i = 0; i < n; i++)
            {
                // Take step and prevent rounding error beyond bound.
                x[i] = (*stp) * d[i] + t[i];
                if ((nbd[i] == 1) || (nbd[i] == 2)) { x[i] = fmax(x[i], l[i]); }
                if ((nbd[i] == 2) || (nbd[i] == 3)) { x[i] = fmin(x[i], u[i]); }
            }
            // 41
        }
    } else {
        *task = NEW_X;
        *task_msg = NO_MSG;
    }

    return;
}


void
matupd(int n, int m, double* ws, double *wy, double* sy, double* ss,
       double* d, double* r, int* itail, int iupdat, int* col,
       int* head, double* theta, double rr, double dr, double stp,
       double dtd)
{
    int i, j, pointr, temp_int, one_int;

    // Set pointers for matrices WS and WY.
    if (iupdat <= m)
    {
        *col = iupdat;
        *itail = (*head + iupdat - 2) % m;  // CHECK THIS EXPRESSION!
    } else {
        *itail = (*itail % m);
        *head = (*head % m);
    }

    // Update matrices WS and WY.
    for (i = 0; i < n; i++)
    {
        ws[i + (*itail)*m] = d[i];
        wy[i + (*itail)*m] = r[i];
    }

    // Set theta=yy/ys.
    *theta = rr / dr;

    // Form the middle matrix in B.

    // update the upper triangle of SS (m, m), and the lower triangle of SY (m, m):
    if (iupdat > m)
    {
        for (j = 1; j < *col; j++)
        {
            DCOPY(&j, &ss[1 + m*j], &one_int, &ss[m*(j-1)], &one_int);
            temp_int = *col - j;
            DCOPY(&temp_int, &sy[j + m*j], &one_int, &sy[(j-1) + m*(j-1)], &one_int);
        }
        // 50
    }

    // add new information: the last row of SY and the last column of SS
    pointr = *head;
    for (j = 0; j < *col - 1; j++)
    {
        sy[*col - 1 + m*j] = DDOT(&n, &d[0], &one_int, &wy[pointr], &one_int);
        ss[j + m*(*col - 1)] = DDOT(&n, &ws[pointr], &one_int, &d[0], &one_int);
        pointr = (pointr % m);
    }
    // 51

    ss[(*col - 1) + m*(*col - 1)] = (stp == 1.0 ? dtd : stp*stp*dtd );
    sy[(*col - 1) + m*(*col - 1)] = dr;

    return;
}


void
projgr(int n, double* l, double* u, int* nbd,
       double* x, double* g, double* sbgnrm)
{
    int i;
    double gi;

    *sbgnrm = 0.0;
    for (i = 0; i < n; i++)
    {
        gi = g[i];
        if (gi != gi)
        {
            // NaN value in gradient: propagate it
            *sbgnrm = gi;
            return;
        }
        if (nbd[i] != 0)
        {
            if (gi < 0.0)
            {
                if (nbd[i] >= 2) { gi = fmax(x[i] - u[i], gi); }
            } else {
                if (nbd[i] <= 2) { gi = fmin(x[i] - l[i], gi); }
            }
        }
        *sbgnrm = fmax(*sbgnrm, fabs(gi));
    }
    // 15

    return;
}


void subsm(int n, int m, int nsub, int* ind,
           double* l, double* u, int* nbd, double* x,
           double* d, double* xp, double* ws, double* wy,
           double theta, double* xx, double* gg, int col,
           int head, int* iword, double* wv, double* wn, int* info)
{
    int pointr, m2, col2, ibd, jy, js, i, j, k, one_int = 1;
    double alpha, xk, dk, temp, temp1, temp2, dd_p;
    char *uplo = "U", *trans = "T", *diag = "N";

    // n, m, col, nsub are size variables.
    // ind, head are index variables/arrays.
    if (nsub <= 0) { return; }

    // Compute wv = W'Zd.
    pointr = head;
    for (i = 0; i < col; i++)
    {
        temp1 = 0.0;
        temp2 = 0.0;
        for (j = 0; j < nsub; j++)
        {
            k = ind[j];
            temp1 += wy[k + m*pointr]*d[j];
            temp2 += ws[k + m*pointr]*d[j];
        }
        wv[i] = temp1;
        wv[col + i] = theta*temp2;
        pointr = (pointr % m);
    }

    // Compute wv:=K^(-1)wv.
    m2 = 2*m;
    col2 = 2*col;

    // dtrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info)
    DTRTRS(uplo, trans, diag, &col2, &one_int, wn, &m2, wv, &col2, info);
    if (*info != 0) { return; }

    // Compute d = (1/theta)d + (1/theta**2)Z'W wv.
    pointr = head;
    for (jy = 0; jy < col; jy++)
    {
        js = col + jy;
        for (i = 0; i < nsub; i++)
        {
            k = ind[i];
            d[i] = d[i] + (wy[k + m*pointr] * wv[jy] / theta) +
                          (ws[k + m*pointr] * wv[js]);
        }
        // 30
        pointr = (pointr % m);
    }
    // 40

    temp = 1.0 / theta;
    // n, da, dx, incx
    DSCAL(&nsub, &temp, d, &one_int);

    // -----------------------------------------------------
    // Let us try the projection, d is the Newton direction.

    *iword = 0;
    // for (i = 0; i < n; i++) { xp[i] = x[i]; }
    DCOPY(&n, x, &one_int, xp, &one_int);

    for (i = 0; i < nsub; i++)
    {
        k = ind[i];
        dk = d[i];
        xk = x[k];
        if (nbd[k] != 0)
        {
            if (nbd[k] == 1)  // Lower bounds only
            {
                x[k] = fmax(l[k], xk + dk);
                if (x[k] == l[k]) { *iword = 1; }
            } else {
                if (nbd[k] == 2)  // Upper and lower bounds
                {
                    xk = fmax(l[k], xk + dk);
                    x[k] = fmin(u[k], xk);
                    if ((x[k] == l[k]) || (x[k] == u[k])) { *iword = 1; }
                } else {
                    if (nbd[k] == 3)  // Upper bounds only
                    {
                        x[k] = fmin(u[k], xk + dk);
                        if (x[k] == u[k]) { *iword = 1; }
                    }
                }
            }
        } else {  // Free variables
            x[k] = xk + dk;
        }
    }
    // 50

    if (*iword == 0) { return; }

    // Check sign of the directional derivative
    dd_p = 0.0;
    for (i = 0; i < n; i++) { dd_p = dd_p + (x[i] - xx[i])*gg[i]; }
    // 55

    if (dd_p > 0.0)
    {
        // for (i = 0; i < n; i++) { x[i] = xp[i]; }
        DCOPY(&n, xp, &one_int, x, &one_int);
    } else {
        return;
    }

    // -------------------------------------------------------------
    alpha = 1.0;
    temp1 = 1.0;
    ibd = 0;

    for (i = 0; i < nsub; i++)
    {
        k = ind[i];
        dk = d[i];
        if (nbd[k] != 0)
        {
            if ((dk < 0.0) && (nbd[k] <= 2))
            {
                temp2 = l[k] - x[k];
                if (temp2 >= 0.0)
                {
                    temp1 = 0.0;
                } else if (dk * alpha < temp2) {
                    temp1 = temp2 / dk;
                }
            } else if ((dk > 0.0) && (nbd[k] >= 2)) {
                temp2 = u[k] - x[k];
                if (temp2 <= 0.0) {
                    temp1 = 0.0;
                } else if (dk * alpha > temp2) {
                    temp1 = temp2 / dk;
                }
            }
            if (temp1 < alpha)
            {
                alpha = temp1;
                ibd = i;
            }
        }
    }
    // 60

    if (alpha < 1.0)
    {
        dk = d[ibd];
        k = ind[ibd];
        if (dk > 0.0)
        {
            x[k] = u[k];
            d[ibd] = 0.0;
        } else if (dk < 0.0) {
            x[k] = l[k];
            d[ibd] = 0.0;
        }
    }

    for (i = 0; i < nsub; i++)
    {
        k = ind[i];
        x[k] = x[k] + alpha*d[i];
    }
    // 70

    return;
}


void
dcsrch(double f, double g, double* stp, double ftol, double gtol,
       double xtol, double stpmin, double stpmax, int* task,
       int* task_msg, int* isave, double* dsave)
{

    int brackt, stage;
    double finit, ftest, fm, fx, fxm, fy, fym, ginit, gtest, gm, gx, gxm, gy;
    double gym, stx, sty, stmin, stmax, width, width1;
    double xtrapl = 1.1;
    double xtrapu = 4.0;


    // Initialization block.
    if (*task == START)
    {
        // Check the input arguments for errors.
        if (*stp < stpmin) { *task_msg = ERROR_STP1; }
        if (*stp > stpmax) { *task_msg = ERROR_STP2; }
        if (g >= 0.0) { *task_msg = ERROR_INITG; }
        if (ftol < 0.0) { *task_msg = ERROR_FTOL; }
        if (gtol < 0.0) { *task_msg = ERROR_GTOL; }
        if (xtol < 0.0) { *task_msg = ERROR_XTOL; }
        if (stpmin < 0.0) { *task_msg = ERROR_STPMIN; }
        if (stpmax < stpmin) { *task_msg = ERROR_STPMAX; }


        // Exit if there are errors on input.
        if (*task_msg > 700)
        {
            *task = ERROR;
            return;
        }

        // Initialize local variables.
        brackt = 0;
        stage = 1;
        finit = f;
        ginit = g;
        gtest = ftol * ginit;
        width = stpmax - stpmin;
        width1 = width / 0.5;

        // The variables stx, fx, gx contain the values of the step, function,
        // and derivative at the best step.
        // The variables sty, fy, gy contain the value of the step, function,
        // and derivative at sty.
        // The variables stp, f, g contain the values of the step, function,
        // and derivative at stp.
        stx = 0.0;
        fx = finit;
        gx = ginit;
        sty = 0.0;
        fy = finit;
        gy = ginit;
        stmin = 0.0;
        stmax = *stp + xtrapu*(*stp);

        *task = FG;
        *task_msg = NO_MSG;
        goto SAVEVARS;

    } else {

        // Restore local variables.
        if (isave[0] == 1) {
            brackt = 1;
        } else {
            brackt = 0;
        }
        stage = isave[1];
        ginit = dsave[0];
        gtest = dsave[1];
        gx = dsave[2];
        gy = dsave[3];
        finit = dsave[4];
        fx = dsave[5];
        fy = dsave[6];
        stx = dsave[7];
        sty = dsave[8];
        stmin = dsave[9];
        stmax = dsave[10];
        width = dsave[11];
        width1 = dsave[12];

    }

    // If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the algorithm
    // enters the second stage.
    ftest = finit + (*stp)*gtest;
    if ((stage == 1) && (f <= ftest) && (g >= 0.0))
    {
        stage = 2;
    }

    // Test for warnings.
    if ((brackt) && ((*stp <= stmin) || (*stp >= stmax)))
    {
        *task = WARNING;
        *task_msg = WARN_ROUND;
    }
    if ((brackt) && (stmax - stmin <= xtol * stmax))
    {
        *task = WARNING;
        *task_msg = WARN_XTOL;
    }
    if ((*stp == stpmax) && (f <= ftest) && (g <= gtest)) {
        *task = WARNING;
        *task_msg = WARN_STPMAX;
    }
    if ((*stp == stpmin) && ((f > ftest) || (g >= gtest))) {
        *task = WARNING;
        *task_msg = WARN_STPMIN;
    }
    // Test for convergence.
    if ((f <= ftest) && (fabs(g) <= gtol * (-ginit))) {
        *task = CONVERGENCE;
    }

    // Test for termination.
    if ((*task == WARNING) || (*task == CONVERGENCE))
    {
        goto SAVEVARS;
    }

    // A modified function is used to predict the step during the first stage if
    // a lower function value has been obtained but the decrease is not sufficient.

    if ((stage == 1) && (f <= fx) && (f > ftest))
    {
        // Define the modified function and derivative values.
        fm = f - (*stp)*gtest;
        fxm = fx - stx*gtest;
        fym = fy - sty*gtest;
        gm = g - gtest;
        gxm = gx - gtest;
        gym = gy - gtest;

        // Call dcstep to update stx, sty, and to compute the new step.
        dcstep(&stx, &fxm, &gxm, &sty, &fym, &gym, stp, fm, gm, &brackt, stmin, stmax);

        // Reset the function and derivative values for f.
        fx = fxm + stx*gtest;
        fy = fym + sty*gtest;
        gx = gxm + gtest;
        gy = gym + gtest;

    } else {
        // Call dcstep to update stx, sty, and to compute the new step.
        dcstep(&stx, &fx, &gx, &sty, &fy, &gy, stp, f, g, &brackt, stmin, stmax);
    }

    // Decide if a bisection step is needed.
    if (brackt)
    {
            if (fabs(sty - stx) >= 0.66*width1)
            {
                *stp = stx + 0.5*(sty - stx);
            }
        width1 = width;
        width = fabs(sty - stx);
    }

    // Set the minimum and maximum steps allowed for stp.
    if (brackt)
    {
        stmin = fmin(stx,sty);
        stmax = fmax(stx,sty);
    } else {
        stmin = *stp + xtrapl*(*stp - stx);
        stmax = *stp + xtrapu*(*stp - stx);
    }

    // Force the step to be within the bounds stpmax and stpmin.
    *stp = fmax(*stp, stpmin);
    *stp = fmin(*stp, stpmax);

    // If further progress is not possible, let stp be the best point obtained
    // during the search.

    if (((brackt) && (*stp <= stmin || *stp >= stmax)) ||
        ((brackt) && (stmax - stmin <= xtol * stmax)))
    {
        *stp = stx;
    }

    // Obtain another function and derivative.
    *task = FG;
    *task_msg = NO_MSG;

SAVEVARS:
    if (brackt)
    {
        isave[0] = 1;
    } else {
        isave[0] = 0;
    }

    isave[1] = stage;
    dsave[0] = ginit;
    dsave[1] = gtest;
    dsave[2] = gx;
    dsave[3] = gy;
    dsave[4] = finit;
    dsave[5] = fx;
    dsave[6] = fy;
    dsave[7] = stx;
    dsave[8] = sty;
    dsave[9] = stmin;
    dsave[10] = stmax;
    dsave[11] = width;
    dsave[12] = width1;

    return;
}


void
dcstep (double* stx, double* fx, double* dx, double* sty, double* fy, double* dy,
        double* stp, double fp, double dp, int* brackt,
        double stpmin, double stpmax)
{
    double gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta;
    sgnd = dp * (*dx / fabs(*dx));

    // First case: A higher function value. The minimum is bracketed.
    // If the cubic step is closer to stx than the quadratic step, the
    // cubic step is taken, otherwise the average of the cubic and
    // quadratic steps is taken.

    if (fp > *fx)
    {
        theta = 3.0 * (*fx - fp) / (*stp - *stx) + *dx + dp;
        s = fmax(fabs(theta), fmax(fabs(*dx), fabs(dp)));
        gamma = s * sqrt(pow(theta / s, 2.0) - (*dx / s)*(dp / s));
        if (*stp < *stx) { gamma = -gamma; }
        p = (gamma - *dx) + theta;
        q = ((gamma - *dx) + gamma) + dp;
        r = p / q;
        stpc = *stx + r * (*stp - *stx);
        stpq = *stx + ((*dx / ((*fx - fp) / (*stp - *stx) + *dx)) / 2.0) * (*stp - *stx);
        if (fabs(stpc - *stx) < fabs(stpq - *stx))
        {
            stpf = stpc;
        } else {
            stpf = stpc + (stpq - stpc) / 2.0;
        }
        *brackt = 1;

    } else if (sgnd < 0.0) {

        // Second case: A lower function value and derivatives of opposite
        // sign. The minimum is bracketed. If the cubic step is farther from
        // stp than the secant step, the cubic step is taken, otherwise the
        // secant step is taken.

        theta = 3.0 * (*fx - fp) / (*stp - *stx) + *dx + dp;
        s = fmax(fabs(theta), fmax(fabs(*dx), fabs(dp)));
        gamma = s * sqrt(pow(theta / s, 2.0) - (*dx / s)*(dp / s));
        if (*stp > *stx) { gamma = -gamma; }
        p = (gamma - dp) + theta;
        q = ((gamma - dp) + gamma) + *dx;
        r = p / q;
        stpc = *stp + r * (*stx - *stp);
        stpq = *stp + (dp / (dp - *dx)) * (stx - stp);
        if (fabs(stpc - *stp) > fabs(stpq - *stp))
        {
            stpf = stpc;
        } else {
            stpf = stpq;
        }
        *brackt = 1;

    } else if (fabs(dp) < fabs(*dx)) {

        // Third case: A lower function value, derivatives of the same sign,
        // and the magnitude of the derivative decreases.

        // The cubic step is computed only if the cubic tends to infinity
        // in the direction of the step or if the minimum of the cubic
        // is beyond stp. Otherwise the cubic step is defined to be the
        // secant step.

        theta = 3.0 * (*fx - fp) / (*stp - *stx) + *dx + dp;
        s = fmax(fabs(theta),fmax(fabs(*dx), fabs(dp)));

        // The case gamma = 0 only arises if the cubic does not tend
        // to infinity in the direction of the step.

        gamma = s*sqrt(fmax(0.0, pow(theta/s, 2.0) - (*dx / s) * (dp / s)));
        if (*stp > *stx) { gamma = -gamma; }
        p = (gamma - dp) + theta;
        q = (gamma + (*dx - dp)) + gamma;
        r = p / q;
        if ((r < 0.0) && (gamma != 0.0))
        {
            stpc = *stp + r*(*stx - *stp);
        } else if (*stp > *stx) {
            stpc = stpmax;
        } else {
            stpc = stpmin;
        }
        stpq = *stp + (dp / (dp - *dx)) * (*stx - *stp);

        if (*brackt)
        {
            // A minimizer has been bracketed. If the cubic step is
            // closer to stp than the secant step, the cubic step is
            // taken, otherwise the secant step is taken.

            if (fabs(stpc - *stp) < fabs(stpq - *stp))
            {
                stpf = stpc;
            } else {
                stpf = stpq;
            }
            if (*stp > *stx)
            {
                stpf = fmin(*stp + 0.66*(*sty - *stp), stpf);
            } else {
                stpf = fmax(*stp + 0.66*(*sty - *stp), stpf);
            }
        } else {
            // A minimizer has not been bracketed. If the cubic step is
            // farther from stp than the secant step, the cubic step is
            // taken, otherwise the secant step is taken.

            if (fabs(stpc - *stp) > fabs(stpq - *stp))
            {
                stpf = stpc;
            } else {
                stpf = stpq;
            }
            stpf = fmin(stpmax, stpf);
            stpf = fmax(stpmin, stpf);
        }

    } else {
        // Fourth case: A lower function value, derivatives of the same sign,
        // and the magnitude of the derivative does not decrease. If the
        // minimum is not bracketed, the step is either stpmin or stpmax,
        // otherwise the cubic step is taken.

        if (*brackt) {
            theta = 3.0 * (fp - *fy) / (*sty - *stp) + *dy + dp;
            s = fmax(fabs(theta), fmax(fabs(*dy), fabs(dp)));
            gamma = s*sqrt(pow(theta/s, 2.0) - (*dy / s) * (dp / s));
            if (*stp > *sty) { gamma = -gamma; }
            p = (gamma - dp) + theta;
            q = ((gamma - dp) + gamma) + *dy;
            r = p / q;
            stpc = *stp + r*(*sty - *stp);
            stpf = stpc;
        } else if (*stp > *stx) {
            stpf = stpmax;
        } else {
            stpf = stpmin;
        }
    }

    // Update the interval which contains a minimizer.

    if (fp > *fx) {
        *sty = *stp;
        *fy = fp;
        *dy = dp;
    } else {
        if (sgnd < 0.0) {
            *sty = *stx;
            *fy = *fx;
            *dy = *dx;
        }
         *stx = *stp;
         *fx = fp;
         *dx = dp;
    }

    // Compute the new step.

    *stp = stpf;

    return;
}
