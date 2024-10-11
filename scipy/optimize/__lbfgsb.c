/*
 * Copyright (C) 2024 SciPy developers
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * a. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * b. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * c. Names of the SciPy Developers may not be used to endorse or promote
 *    products derived from this software without specific prior written
 *    permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 *
 * This file is a C translation of the Fortran code written by Ciyou Zhu,
 * Richard Byrd, Jorge Nocedal and Jose Luis Morales with the original
 * description below.
 * Some of the functions have been modified and certain arguments are rendered
 * obsolete to be used in SciPy however original docstrings are kept at the top
 * of each function for reference.
 *
 */


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * ===========   L-BFGS-B (version 3.0.  April 25, 2011  ===================
 *
 *     This is a modified version of L-BFGS-B. Minor changes in the updated
 *     code appear preceded by a line comment as follows
 *
 *     c-jlm-jn
 *
 *     Major changes are described in the accompanying paper:
 *
 *         Jorge Nocedal and Jose Luis Morales, Remark on "Algorithm 778:
 *         L-BFGS-B: Fortran Subroutines for Large-Scale Bound Constrained
 *         Optimization"  (2011). To appear in  ACM Transactions on
 *         Mathematical Software,
 *
 *     The paper describes an improvement and a correction to Algorithm 778.
 *     It is shown that the performance of the algorithm can be improved
 *     significantly by making a relatively simple modication to the subspace
 *     minimization phase. The correction concerns an error caused by the use
 *     of routine dpmeps to estimate machine precision.
 *
 *     The total work space **wa** required by the new version is
 *
 *                  2*m*n + 11m*m + 5*n + 8*m
 *
 *     the old version required
 *
 *                  2*m*n + 12m*m + 4*n + 12*m
 *
 *
 *            J. Nocedal  Department of Electrical Engineering and
 *                        Computer Science.
 *                        Northwestern University. Evanston, IL. USA
 *
 *
 *           J.L Morales  Departamento de Matematicas,
 *                        Instituto Tecnologico Autonomo de Mexico
 *                        Mexico D.F. Mexico.
 *
 *                        March  2011
 *
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * The original BSD-3 licensed Fortran code can be found at
 *
 * https://users.iems.northwestern.edu/~nocedal/lbfgsb.html
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 * References:
 *
 * [1]: R. H. Byrd, P. Lu and J. Nocedal. A Limited Memory Algorithm for Bound
 *      Constrained Optimization, (1995), SIAM Journal on Scientific and
 *      Statistical Computing, 16, 5, pp. 1190-1208.
 * [2]: C. Zhu, R. H. Byrd and J. Nocedal. L-BFGS-B: Algorithm 778: L-BFGS-B,
 *      FORTRAN routines for large scale bound constrained optimization (1997),
 *      ACM Transactions on Mathematical Software, 23, 4, pp. 550 - 560.
 * [3]: J.L. Morales and J. Nocedal. L-BFGS-B: Remark on Algorithm 778:
 *      L-BFGS-B, FORTRAN routines for large scale bound constrained
 *      optimization (2011), ACM Transactions on Mathematical Software, 38, 1.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "__lbfgsb.h"


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
    STOP_ITERC   = 504,    // TOTAL NO. of ITERATIONS REACHED LIMIT
    STOP_CALLB   = 505,    // CALLBACK REQUESTED HALT

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


// Internal functions

static void mainlb(
    int n, int m, double* x, double* l, double* u,
    int* nbd, double* f, double* g, double factr, double pgtol,
    double* ws, double* wy, double* sy, double* ss, double* wt, double* wn,
    double* snd, double* z, double* r, double* d, double* t, double* xp,
    double* wa, int* index, int* iwhere, int* indx2, int* task, int* task_msg,
    int* lsave, int* isave, double* dsave, int maxls, int* ln_task, int* ln_taskmsg
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
    int cnstnd, int* isave, double* dsave, int* temp_task, int* temp_task_msg
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
setulb(int n, int m, double* x, double* l, double* u, int* nbd, double* f,
       double* g, double factr, double pgtol, double* wa, int* iwa, int* task,
       int* lsave, int* isave, double* dsave, int maxls, int* ln_task)
{
    //     ************
    //
    //     Subroutine setulb
    //
    //     This subroutine partitions the working arrays wa and iwa, and
    //       then uses the limited memory BFGS method to solve the bound
    //       constrained optimization problem by calling mainlb.
    //       (The direct method will be used in the subspace minimization.)
    //
    //     n is an integer variable.
    //       On entry n is the dimension of the problem.
    //       On exit n is unchanged.
    //
    //     m is an integer variable.
    //       On entry m is the maximum number of variable metric corrections
    //         used to define the limited memory matrix.
    //       On exit m is unchanged.
    //
    //     x is a double precision array of dimension n.
    //       On entry x is an approximation to the solution.
    //       On exit x is the current approximation.
    //
    //     l is a double precision array of dimension n.
    //       On entry l is the lower bound on x.
    //       On exit l is unchanged.
    //
    //     u is a double precision array of dimension n.
    //       On entry u is the upper bound on x.
    //       On exit u is unchanged.
    //
    //     nbd is an integer array of dimension n.
    //       On entry nbd represents the type of bounds imposed on the
    //         variables, and must be specified as follows:
    //         nbd(i)=0 if x(i) is unbounded,
    //                1 if x(i) has only a lower bound,
    //                2 if x(i) has both lower and upper bounds, and
    //                3 if x(i) has only an upper bound.
    //       On exit nbd is unchanged.
    //
    //     f is a double precision variable.
    //       On first entry f is unspecified.
    //       On final exit f is the value of the function at x.
    //
    //     g is a double precision array of dimension n.
    //       On first entry g is unspecified.
    //       On final exit g is the value of the gradient at x.
    //
    //     factr is a double precision variable.
    //       On entry factr >= 0 is specified by the user.  The iteration
    //         will stop when
    //
    //         (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
    //
    //         where epsmch is the machine precision, which is automatically
    //         generated by the code. Typical values for factr: 1.d+12 for
    //         low accuracy; 1.d+7 for moderate accuracy; 1.d+1 for extremely
    //         high accuracy.
    //       On exit factr is unchanged.
    //
    //     pgtol is a double precision variable.
    //       On entry pgtol >= 0 is specified by the user.  The iteration
    //         will stop when
    //
    //                 max{|proj g_i | i = 1, ..., n} <= pgtol
    //
    //         where pg_i is the ith component of the projected gradient.
    //       On exit pgtol is unchanged.
    //
    //     wa is a double precision working array of length
    //       (2mmax + 5)nmax + 12mmax^2 + 12mmax.
    //
    //     iwa is an integer working array of length 3nmax.
    //
    //     task is a working string of characters of length 60 indicating
    //       the current job when entering and quitting this subroutine.
    //
    //     iprint is an integer variable that must be set by the user.
    //       It controls the frequency and type of output generated:
    //        iprint<0    no output is generated;
    //        iprint=0    print only one line at the last iteration;
    //        0<iprint<99 print also f and |proj g| every iprint iterations;
    //        iprint=99   print details of every iteration except n-vectors;
    //        iprint=100  print also the changes of active set and final x;
    //        iprint>100  print details of every iteration including x and g;
    //
    //
    //     csave is a working string of characters of length 60.
    //
    //     lsave is a logical working array of dimension 4.
    //       On exit with 'task' = NEW_X, the following information is
    //                                                             available:
    //         If lsave(1) = .true.  then  the initial X has been replaced by
    //                                     its projection in the feasible set;
    //         If lsave(2) = .true.  then  the problem is constrained;
    //         If lsave(3) = .true.  then  each variable has upper and lower
    //                                     bounds;
    //
    //     isave is an integer working array of dimension 44.
    //       On exit with 'task' = NEW_X, the following information is
    //                                                             available:
    //         isave(22) = the total number of intervals explored in the
    //                         search of Cauchy points;
    //         isave(26) = the total number of skipped BFGS updates before
    //                         the current iteration;
    //         isave(30) = the number of current iteration;
    //         isave(31) = the total number of BFGS updates prior the current
    //                         iteration;
    //         isave(33) = the number of intervals explored in the search of
    //                         Cauchy point in the current iteration;
    //         isave(34) = the total number of function and gradient
    //                         evaluations;
    //         isave(36) = the number of function value or gradient
    //                                  evaluations in the current iteration;
    //         if isave(37) = 0  then the subspace argmin is within the box;
    //         if isave(37) = 1  then the subspace argmin is beyond the box;
    //         isave(38) = the number of free variables in the current
    //                         iteration;
    //         isave(39) = the number of active constraints in the current
    //                         iteration;
    //         n + 1 - isave(40) = the number of variables leaving the set of
    //                           active constraints in the current iteration;
    //         isave(41) = the number of variables entering the set of active
    //                         constraints in the current iteration.
    //
    //     dsave is a double precision working array of dimension 29.
    //       On exit with 'task' = NEW_X, the following information is
    //                                                             available:
    //         dsave(1) = current 'theta' in the BFGS matrix;
    //         dsave(2) = f(x) in the previous iteration;
    //         dsave(3) = factr*epsmch;
    //         dsave(4) = 2-norm of the line search direction vector;
    //         dsave(5) = the machine precision epsmch generated by the code;
    //         dsave(7) = the accumulated time spent on searching for
    //                                                         Cauchy points;
    //         dsave(8) = the accumulated time spent on
    //                                                 subspace minimization;
    //         dsave(9) = the accumulated time spent on line search;
    //         dsave(11) = the slope of the line search function at
    //                                  the current point of line search;
    //         dsave(12) = the maximum relative step length imposed in
    //                                                           line search;
    //         dsave(13) = the infinity norm of the projected gradient;
    //         dsave(14) = the relative step length in the line search;
    //         dsave(15) = the slope of the line search function at
    //                                 the starting point of the line search;
    //         dsave(16) = the square of the 2-norm of the line search
    //                                                      direction vector.
    //
    //     Subprograms called:
    //
    //       L-BFGS-B Library ... mainlb.
    //
    //
    //     References:
    //
    //       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
    //       memory algorithm for bound constrained optimization'',
    //       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
    //
    //       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
    //       limited memory FORTRAN code for solving bound constrained
    //       optimization problems'', Tech. Report, NAM-11, EECS Department,
    //       Northwestern University, 1994.
    //
    //       (Postscript files of these papers are available via anonymous
    //        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
    //
    //                           *  *  *
    //
    //     NEOS, November 1994. (Latest revision June 1996.)
    //     Optimization Technology Center.
    //     Argonne National Laboratory and Northwestern University.
    //     Written by
    //                        Ciyou Zhu
    //     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
    //
    //
    //     ************
    //
    int lws, lr, lz, lt, ld, lxp, lwa, lwy, lsy, lss, lwt, lwn, lsnd;

    if (task[0] == START)
    {
        // Pointer olympics
        isave[0]  = m*n;
        isave[1]  = m*m;
        isave[2]  = 4*m*m;
        isave[3]  = 0;                      // ws      m*n
        isave[4]  = isave[3]  + isave[0];   // wy      m*n
        isave[5]  = isave[4]  + isave[0];   // wsy     m**2
        isave[6]  = isave[5]  + isave[1];   // wss     m**2
        isave[7]  = isave[6]  + isave[1];   // wt      m**2
        isave[8]  = isave[7]  + isave[1];   // wn      4*m**2
        isave[9]  = isave[8]  + isave[2];   // wsnd    4*m**2
        isave[10] = isave[9]  + isave[2];   // wz      n
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
           &wa[lt], &wa[lxp], &wa[lwa], iwa, &iwa[n], &iwa[2*n], &task[0], &task[1],
           lsave, &isave[21], dsave, maxls, &ln_task[0], &ln_task[1]);

    return;
}


void
mainlb(int n, int m, double* x, double* l, double* u,
       int* nbd, double* f, double* g, double factr, double pgtol,
       double* ws, double* wy, double* sy, double* ss, double* wt, double* wn,
       double* snd, double* z, double* r, double* d, double* t, double* xp,
       double* wa, int* index, int* iwhere, int* indx2, int* task,
       int* task_msg, int* lsave, int* isave, double* dsave, int maxls,
       int* temp_task, int* temp_taskmsg)
{
    //     ************
    //
    //     Subroutine mainlb
    //
    //     This subroutine solves bound constrained optimization problems by
    //       using the compact formula of the limited memory BFGS updates.
    //
    //     n is an integer variable.
    //       On entry n is the number of variables.
    //       On exit n is unchanged.
    //
    //     m is an integer variable.
    //       On entry m is the maximum number of variable metric
    //          corrections allowed in the limited memory matrix.
    //       On exit m is unchanged.
    //
    //     x is a double precision array of dimension n.
    //       On entry x is an approximation to the solution.
    //       On exit x is the current approximation.
    //
    //     l is a double precision array of dimension n.
    //       On entry l is the lower bound of x.
    //       On exit l is unchanged.
    //
    //     u is a double precision array of dimension n.
    //       On entry u is the upper bound of x.
    //       On exit u is unchanged.
    //
    //     nbd is an integer array of dimension n.
    //       On entry nbd represents the type of bounds imposed on the
    //         variables, and must be specified as follows:
    //         nbd(i)=0 if x(i) is unbounded,
    //                1 if x(i) has only a lower bound,
    //                2 if x(i) has both lower and upper bounds,
    //                3 if x(i) has only an upper bound.
    //       On exit nbd is unchanged.
    //
    //     f is a double precision variable.
    //       On first entry f is unspecified.
    //       On final exit f is the value of the function at x.
    //
    //     g is a double precision array of dimension n.
    //       On first entry g is unspecified.
    //       On final exit g is the value of the gradient at x.
    //
    //     factr is a double precision variable.
    //       On entry factr >= 0 is specified by the user.  The iteration
    //         will stop when
    //
    //         (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
    //
    //         where epsmch is the machine precision, which is automatically
    //         generated by the code.
    //       On exit factr is unchanged.
    //
    //     pgtol is a double precision variable.
    //       On entry pgtol >= 0 is specified by the user.  The iteration
    //         will stop when
    //
    //                 max{|proj g_i | i = 1, ..., n} <= pgtol
    //
    //         where pg_i is the ith component of the projected gradient.
    //       On exit pgtol is unchanged.
    //
    //     ws, wy, sy, and wt are double precision working arrays used to
    //       store the following information defining the limited memory
    //          BFGS matrix:
    //          ws, of dimension n x m, stores S, the matrix of s-vectors;
    //          wy, of dimension n x m, stores Y, the matrix of y-vectors;
    //          sy, of dimension m x m, stores S'Y;
    //          ss, of dimension m x m, stores S'S;
    //          yy, of dimension m x m, stores Y'Y;
    //          wt, of dimension m x m, stores the Cholesky factorization
    //                                  of (theta*S'S+LD^(-1)L'); see eq.
    //                                  (2.26) in [3].
    //
    //     wn is a double precision working array of dimension 2m x 2m
    //       used to store the LEL^T factorization of the indefinite matrix
    //                 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
    //                     [L_a -R_z           theta*S'AA'S ]
    //
    //       where     E = [-I  0]
    //                     [ 0  I]
    //
    //     snd is a double precision working array of dimension 2m x 2m
    //       used to store the lower triangular part of
    //                 N = [Y' ZZ'Y   L_a'+R_z']
    //                     [L_a +R_z  S'AA'S   ]
    //
    //     z(n),r(n),d(n),t(n), xp(n),wa(8*m) are double precision working arrays.
    //       z  is used at different times to store the Cauchy point and
    //          the Newton point.
    //       xp is used to safeguard the projected Newton direction
    //
    //     sg(m),sgo(m),yg(m),ygo(m) are double precision working arrays.
    //
    //     index is an integer working array of dimension n.
    //       In subroutine freev, index is used to store the free and fixed
    //          variables at the Generalized Cauchy Point (GCP).
    //
    //     iwhere is an integer working array of dimension n used to record
    //       the status of the vector x for GCP computation.
    //       iwhere(i)=0 or -3 if x(i) is free and has bounds,
    //                 1       if x(i) is fixed at l(i), and l(i) .ne. u(i)
    //                 2       if x(i) is fixed at u(i), and u(i) .ne. l(i)
    //                 3       if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
    //                -1       if x(i) is always free, i.e., no bounds on it.
    //
    //     indx2 is an integer working array of dimension n.
    //       Within subroutine cauchy, indx2 corresponds to the array iorder.
    //       In subroutine freev, a list of variables entering and leaving
    //       the free set is stored in indx2, and it is passed on to
    //       subroutine formk with this information.
    //
    //     task is a working string of characters of length 60 indicating
    //       the current job when entering and leaving this subroutine.
    //
    //     iprint is an INTEGER variable that must be set by the user.
    //       It controls the frequency and type of output generated:
    //        iprint<0    no output is generated;
    //        iprint=0    print only one line at the last iteration;
    //        0<iprint<99 print also f and |proj g| every iprint iterations;
    //        iprint=99   print details of every iteration except n-vectors;
    //        iprint=100  print also the changes of active set and final x;
    //        iprint>100  print details of every iteration including x and g;
    //
    //
    //     csave is a working string of characters of length 60.
    //
    //     lsave is a logical working array of dimension 4.
    //
    //     isave is an integer working array of dimension 23.
    //
    //     dsave is a double precision working array of dimension 29.
    //
    //
    //     Subprograms called
    //
    //       L-BFGS-B Library ... cauchy, subsm, lnsrlb, formk,
    //
    //        errclb, prn1lb, prn2lb, prn3lb, active, projgr,
    //
    //        freev, cmprlb, matupd, formt.
    //
    //       Minpack2 Library ... timer
    //
    //       Linpack Library ... dcopy, ddot.
    //
    //
    //     References:
    //
    //       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
    //       memory algorithm for bound constrained optimization'',
    //       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
    //
    //       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
    //       Subroutines for Large Scale Bound Constrained Optimization''
    //       Tech. Report, NAM-11, EECS Department, Northwestern University,
    //       1994.
    //
    //       [3] R. Byrd, J. Nocedal and R. Schnabel "Representations of
    //       Quasi-Newton Matrices and their use in Limited Memory Methods'',
    //       Mathematical Programming 63 (1994), no. 4, pp. 129-156.
    //
    //       (Postscript files of these papers are available via anonymous
    //        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
    //
    //                           *  *  *
    //
    //     NEOS, November 1994. (Latest revision June 1996.)
    //     Optimization Technology Center.
    //     Argonne National Laboratory and Northwestern University.
    //     Written by
    //                        Ciyou Zhu
    //     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
    //
    //
    //     ************
    int prjctd, cnstnd, boxed, updatd, wrk;
    int i, k, nintol, iback, nskip, head, col, iter, itail, iupdat;
    int nseg, nfgv, info, ifun, iword, nfree, nact, ileave, nenter;
    double theta, fold, dr, rr, tol, xstep, sbgnrm, stpmx, ddum, dnorm, dtd;
    double gd, gdold, stp;
    int one_int = 1;

    stpmx = 0.0;

    if (*task == START) {

        // Initialize counters and scalars
        col = 0;      // size  : Actual number of variable metric corrections
        head = 0;     // index : Location of the first s- or y-vector in S or Y
        theta = 1.0;
        iupdat = 0;   // size  : Total number of BFGS updates made so far
        updatd = 0;   // bool  : Whether L-BFGS matrix is updated
        iback = 0;
        itail = 0;    // index
        iword = 0;
        nact = 0;
        ileave = 0;   // index  : idx2[ileave:n] are the variables leaving the free set
        nenter = 0;   // size   : No of variables entering the free set
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
        boxed  = lsave[2];
        updatd = lsave[3];

        nintol = isave[0];
        //     = isave[1];
        //     = isave[2];
        iback  = isave[3];
        nskip  = isave[4];
        head   = isave[5];
        col    = isave[6];
        itail  = isave[7];
        iter   = isave[8];
        iupdat = isave[9];
        //     = isave[10];
        nseg   = isave[11];
        nfgv   = isave[12];
        info   = isave[13];
        ifun   = isave[14];
        iword  = isave[15];
        nfree  = isave[16];
        nact   = isave[17];
        ileave = isave[18];
        nenter = isave[19];

        theta  = dsave[0];
        fold   = dsave[1];
        tol    = dsave[2];
        dnorm  = dsave[3];
        //     = dsave[4];
        //     = dsave[5];
        //     = dsave[6];
        //     = dsave[7];
        //     = dsave[8];
        //     = dsave[9];
        gd     = dsave[10];
        stpmx  = dsave[11];
        sbgnrm = dsave[12];
        stp    = dsave[13];
        gdold  = dsave[14];
        dtd    = dsave[15];

    // After returning from the driver go to the point where execution is to resume.
        if ((*task == FG) && (*task_msg == FG_LNSRCH)) { goto LINE666; }
        if (*task == NEW_X) { goto LINE777; }
        if ((*task == FG) && (*task_msg == FG_START)) { goto LINE111; }
        if (*task == STOP)
        {
            if (*task_msg == STOP_CPU)
            {
                dcopy_(&n, t, &one_int, x, &one_int);
                dcopy_(&n, r, &one_int, g, &one_int);
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
        dcopy_(&n, x, &one_int, z, &one_int);
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
           task_msg, boxed, cnstnd, &isave[21], &dsave[16], temp_task, temp_taskmsg);

    if ((info != 0) || (iback >= maxls))
    {
        dcopy_(&n, t, &one_int, x, &one_int);
        dcopy_(&n, r, &one_int, g, &one_int);
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
    } else if ((*task == FG) && (*task_msg == FG_LNSRCH)) {
        // Return to the driver for calculating f and g; renter at 666.
        goto LINE1000;
    } else {
        // Calculate and print out the quantities related to the new X.
        iter++;

        // Compute the infinity norm of the projected (-)gradient.
        projgr(n, l, u, nbd, x, g, &sbgnrm);

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
    rr = pow(dnrm2_(&n, r, &one_int), 2.0);

    if (stp == 1.0)
    {
        dr = gd - gdold;
        ddum = -gdold;
    } else {
        dr = (gd - gdold)*stp;
        dscal_(&n, &stp, d, &one_int);
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
    isave[3]  = iback;
    isave[4]  = nskip;
    isave[5]  = head;
    isave[6]  = col;
    isave[7]  = itail;
    isave[8]  = iter;
    isave[9]  = iupdat;
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
    // dsave[9]  = time1;
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
    //     ************
    //
    //     Subroutine active
    //
    //     This subroutine initializes iwhere and projects the initial x to
    //       the feasible set if necessary.
    //
    //     iwhere is an integer array of dimension n.
    //       On entry iwhere is unspecified.
    //       On exit iwhere(i)=-1  if x(i) has no bounds
    //                         3   if l(i)=u(i)
    //                         0   otherwise.
    //       In cauchy, iwhere is given finer gradations.
    //
    //
    //                           *  *  *
    //
    //     NEOS, November 1994. (Latest revision June 1996.)
    //     Optimization Technology Center.
    //     Argonne National Laboratory and Northwestern University.
    //     Written by
    //                        Ciyou Zhu
    //     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
    //
    //
    //     ************
    int i;

    // Initialize nbdd, prjctd, cnstnd, and boxed.
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
            } else if ((nbd[i] >= 2) && (x[i] >= u[i])) {
                if (x[i] > u[i])
                {
                    *prjctd = 1;
                    x[i] = u[i];
                }
            }
        }
    }
    // 10

    // Initialize iwhere and assign values to cnstnd and boxed.
    for (i = 0; i < n; i++)
    {
        if (nbd[i] != 2) { *boxed = 0; }
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
bmv(int m, double* sy, double* wt, int col, double* v, double* p, int* info)
{
    //     ************
    //
    //     Subroutine bmv
    //
    //     This subroutine computes the product of the 2m x 2m middle matrix
    //       in the compact L-BFGS formula of B and a 2m vector v;
    //       it returns the product in p.
    //
    //     m is an integer variable.
    //       On entry m is the maximum number of variable metric corrections
    //         used to define the limited memory matrix.
    //       On exit m is unchanged.
    //
    //     sy is a double precision array of dimension m x m.
    //       On entry sy specifies the matrix S'Y.
    //       On exit sy is unchanged.
    //
    //     wt is a double precision array of dimension m x m.
    //       On entry wt specifies the upper triangular matrix J' which is
    //         the Cholesky factor of (thetaS'S+LD^(-1)L').
    //       On exit wt is unchanged.
    //
    //     col is an integer variable.
    //       On entry col specifies the number of s-vectors (or y-vectors)
    //         stored in the compact L-BFGS formula.
    //       On exit col is unchanged.
    //
    //     v is a double precision array of dimension 2col.
    //       On entry v specifies vector v.
    //       On exit v is unchanged.
    //
    //     p is a double precision array of dimension 2col.
    //       On entry p is unspecified.
    //       On exit p is the product Mv.
    //
    //     info is an integer variable.
    //       On entry info is unspecified.
    //       On exit info = 0       for normal return,
    //                    = nonzero for abnormal return when the system
    //                                to be solved by dtrsl is singular.
    //
    //     Subprograms called:
    //
    //       Linpack ... dtrsl.
    //
    //
    //                           *  *  *
    //
    //     NEOS, November 1994. (Latest revision June 1996.)
    //     Optimization Technology Center.
    //     Argonne National Laboratory and Northwestern University.
    //     Written by
    //                        Ciyou Zhu
    //     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
    //
    //
    //     ************
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
        for (k = 0; k <= i-1; k++)
        {
            ssum = ssum + sy[i + m*k] * v[k] / sy[k + m*k];
        }
        // 10
        p[i2] = v[i2] + ssum;
    }
    // 20

    // Solve the triangular system
    // dtrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info)
    dtrtrs_(uplo, trans, diag, &col, &one_int, wt, &m, &p[col], &col, info);
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
    trans = "N";
    dtrtrs_(uplo, trans, diag, &col, &one_int, wt, &m, &p[col], &col, info);
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
    //     ************
    //
    //     Subroutine cauchy
    //
    //     For given x, l, u, g (with sbgnrm > 0), and a limited memory
    //       BFGS matrix B defined in terms of matrices WY, WS, WT, and
    //       scalars head, col, and theta, this subroutine computes the
    //       generalized Cauchy point (GCP), defined as the first local
    //       minimizer of the quadratic
    //
    //                  Q(x + s) = g's + 1/2 s'Bs
    //
    //       along the projected gradient direction P(x-tg,l,u).
    //       The routine returns the GCP in xcp.
    //
    //     n is an integer variable.
    //       On entry n is the dimension of the problem.
    //       On exit n is unchanged.
    //
    //     x is a double precision array of dimension n.
    //       On entry x is the starting point for the GCP computation.
    //       On exit x is unchanged.
    //
    //     l is a double precision array of dimension n.
    //       On entry l is the lower bound of x.
    //       On exit l is unchanged.
    //
    //     u is a double precision array of dimension n.
    //       On entry u is the upper bound of x.
    //       On exit u is unchanged.
    //
    //     nbd is an integer array of dimension n.
    //       On entry nbd represents the type of bounds imposed on the
    //         variables, and must be specified as follows:
    //         nbd(i)=0 if x(i) is unbounded,
    //                1 if x(i) has only a lower bound,
    //                2 if x(i) has both lower and upper bounds, and
    //                3 if x(i) has only an upper bound.
    //       On exit nbd is unchanged.
    //
    //     g is a double precision array of dimension n.
    //       On entry g is the gradient of f(x).  g must be a nonzero vector.
    //       On exit g is unchanged.
    //
    //     iorder is an integer working array of dimension n.
    //       iorder will be used to store the breakpoints in the piecewise
    //       linear path and free variables encountered. On exit,
    //         iorder(1),...,iorder(nleft) are indices of breakpoints
    //                                which have not been encountered;
    //         iorder(nleft+1),...,iorder(nbreak) are indices of
    //                                     encountered breakpoints; and
    //         iorder(nfree),...,iorder(n) are indices of variables which
    //                 have no bound constraits along the search direction.
    //
    //     iwhere is an integer array of dimension n.
    //       On entry iwhere indicates only the permanently fixed (iwhere=3)
    //       or free (iwhere= -1) components of x.
    //       On exit iwhere records the status of the current x variables.
    //       iwhere(i)=-3  if x(i) is free and has bounds, but is not moved
    //                 0   if x(i) is free and has bounds, and is moved
    //                 1   if x(i) is fixed at l(i), and l(i) .ne. u(i)
    //                 2   if x(i) is fixed at u(i), and u(i) .ne. l(i)
    //                 3   if x(i) is always fixed, i.e.,  u(i)=x(i)=l(i)
    //                 -1  if x(i) is always free, i.e., it has no bounds.
    //
    //     t is a double precision working array of dimension n.
    //       t will be used to store the break points.
    //
    //     d is a double precision array of dimension n used to store
    //       the Cauchy direction P(x-tg)-x.
    //
    //     xcp is a double precision array of dimension n used to return the
    //       GCP on exit.
    //
    //     m is an integer variable.
    //       On entry m is the maximum number of variable metric corrections
    //         used to define the limited memory matrix.
    //       On exit m is unchanged.
    //
    //     ws, wy, sy, and wt are double precision arrays.
    //       On entry they store information that defines the
    //                             limited memory BFGS matrix:
    //         ws(n,m) stores S, a set of s-vectors;
    //         wy(n,m) stores Y, a set of y-vectors;
    //         sy(m,m) stores S'Y;
    //         wt(m,m) stores the
    //                 Cholesky factorization of (theta*S'S+LD^(-1)L').
    //       On exit these arrays are unchanged.
    //
    //     theta is a double precision variable.
    //       On entry theta is the scaling factor specifying B_0 = theta I.
    //       On exit theta is unchanged.
    //
    //     col is an integer variable.
    //       On entry col is the actual number of variable metric
    //         corrections stored so far.
    //       On exit col is unchanged.
    //
    //     head is an integer variable.
    //       On entry head is the location of the first s-vector (or y-vector)
    //         in S (or Y).
    //       On exit col is unchanged.
    //
    //     p is a double precision working array of dimension 2m.
    //       p will be used to store the vector p = W^(T)d.
    //
    //     c is a double precision working array of dimension 2m.
    //       c will be used to store the vector c = W^(T)(xcp-x).
    //
    //     wbp is a double precision working array of dimension 2m.
    //       wbp will be used to store the row of W corresponding
    //         to a breakpoint.
    //
    //     v is a double precision working array of dimension 2m.
    //
    //     nseg is an integer variable.
    //       On exit nseg records the number of quadratic segments explored
    //         in searching for the GCP.
    //
    //     sg and yg are double precision arrays of dimension m.
    //       On entry sg  and yg store S'g and Y'g correspondingly.
    //       On exit they are unchanged.
    //
    //     iprint is an INTEGER variable that must be set by the user.
    //       It controls the frequency and type of output generated:
    //        iprint<0    no output is generated;
    //        iprint=0    print only one line at the last iteration;
    //        0<iprint<99 print also f and |proj g| every iprint iterations;
    //        iprint=99   print details of every iteration except n-vectors;
    //        iprint=100  print also the changes of active set and final x;
    //        iprint>100  print details of every iteration including x and g;
    //
    //
    //     sbgnrm is a double precision variable.
    //       On entry sbgnrm is the norm of the projected gradient at x.
    //       On exit sbgnrm is unchanged.
    //
    //     info is an integer variable.
    //       On entry info is 0.
    //       On exit info = 0       for normal return,
    //                    = nonzero for abnormal return when the system
    //                              used in routine bmv is singular.
    //
    //     Subprograms called:
    //
    //       L-BFGS-B Library ... hpsolb, bmv.
    //
    //       Linpack ... dscal dcopy, daxpy.
    //
    //
    //     References:
    //
    //       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
    //       memory algorithm for bound constrained optimization'',
    //       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
    //
    //       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: FORTRAN
    //       Subroutines for Large Scale Bound Constrained Optimization''
    //       Tech. Report, NAM-11, EECS Department, Northwestern University,
    //       1994.
    //
    //       (Postscript files of these papers are available via anonymous
    //        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
    //
    //                           *  *  *
    //
    //     NEOS, November 1994. (Latest revision June 1996.)
    //     Optimization Technology Center.
    //     Argonne National Laboratory and Northwestern University.
    //     Written by
    //                        Ciyou Zhu
    //     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
    //
    //
    //     ************
    int bnded, xlower, xupper;
    int i, j, col2, nfree, nbreak, pointr, ibp, nleft, ibkmin, iter, one_int = 1;
    double f1, f2, dt, dtm, tsum, dibp, zibp, dibp2, bkmin, tu, tl, wmc, wmp;
    double wmw, tj, tj0, neggi, f2_org;

    f1 = 0.0;
    tl = 0.0;
    tu = 0.0;

    if (sbgnrm <= 0.0)
    {
        dcopy_(&n, x, &one_int, xcp, &one_int);
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
                pointr = (pointr + 1) % m;
            }
            // 40
            if ((nbd[i] <= 2) && (nbd[i] != 0) && (neggi < 0.0))
            {
                // x(i) + d(i) is bounded; compute t(i).
                nbreak += 1;
                iorder[nbreak] = i;
                t[nbreak] = tl / (-neggi);
                if ((nbreak == 0) || (t[nbreak] < bkmin))
                {
                    bkmin = t[nbreak];
                    ibkmin = nbreak;
                }
            } else if ((nbd[i] >= 2) && (neggi > 0.0)) {
                // x(i) + d(i) is bounded; compute t(i).
                nbreak += 1;
                iorder[nbreak] = i;
                t[nbreak] = tu / neggi;
                if ((nbreak == 0) || (t[nbreak] < bkmin))
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
        dscal_(&col, &theta, &p[col], &one_int);
    }
    // Initialize GCP xcp = x.
    dcopy_(&n, x, &one_int, xcp, &one_int);
    if ((nbreak == -1) && (nfree == n))
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
        f2 = f2 - ddot_(&col2, v, &one_int, p, &one_int);
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
            daxpy_(&col2, &dt, p, &one_int, c, &one_int);
            // Choose wbp, the row of W corresponding to the breakpoint encountered.
            pointr = head;
            for (j = 0; j < col; j++) {
                wbp[j] = wy[ibp + n*pointr];
                wbp[col + j] = theta * ws[ibp + n*pointr];
                pointr = (pointr + 1) % m;
            }
            // 70

            // Compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'.
            bmv(m, sy, wt, col, wbp, v, info);
            if (*info != 0) { return; }
            wmc = ddot_(&col2, c, &one_int, v, &one_int);
            wmp = ddot_(&col2, p, &one_int, v, &one_int);
            wmw = ddot_(&col2, wbp, &one_int, v, &one_int);
            // Update p = p - dibp*wbp.
            dibp = -dibp;
            daxpy_(&col2, &dibp, wbp, &one_int, p, &one_int);
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
    daxpy_(&n, &tsum, d, &one_int, xcp, &one_int);

LINE999:
    // Update c = c + dtm*p = W'(x^c - x)
    // which will be used in computing r = Z'(B(x^c - x) + g).

    if (col > 0) { daxpy_(&col2, &dtm, p, &one_int, c, &one_int); }

    return;
}


void
cmprlb(int n, int m, double* x, double* g,
       double* ws, double* wy, double* sy, double* wt,
       double* z, double* r, double* wa, int* index,
       double theta, int col, int head, int nfree,
       int cnstnd, int* info)
{
    //     ************
    //
    //     Subroutine cmprlb
    //
    //       This subroutine computes r=-Z'B(xcp-xk)-Z'g by using
    //         wa(2m+1)=W'(xcp-x) from subroutine cauchy.
    //
    //     Subprograms called:
    //
    //       L-BFGS-B Library ... bmv.
    //
    //
    //                           *  *  *
    //
    //     NEOS, November 1994. (Latest revision June 1996.)
    //     Optimization Technology Center.
    //     Argonne National Laboratory and Northwestern University.
    //     Written by
    //                        Ciyou Zhu
    //     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
    //
    //
    //     ************
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
            pointr = (pointr + 1) % m;
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
    //     ************
    //
    //     Subroutine errclb
    //
    //     This subroutine checks the validity of the input data.
    //
    //
    //                           *  *  *
    //
    //     NEOS, November 1994. (Latest revision June 1996.)
    //     Optimization Technology Center.
    //     Argonne National Laboratory and Northwestern University.
    //     Written by
    //                        Ciyou Zhu
    //     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
    //
    //
    //     ************
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
    //     ************
    //
    //     Subroutine formk
    //
    //     This subroutine forms  the LEL^T factorization of the indefinite
    //
    //       matrix    K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
    //                     [L_a -R_z           theta*S'AA'S ]
    //                                                    where E = [-I  0]
    //                                                              [ 0  I]
    //     The matrix K can be shown to be equal to the matrix M^[-1]N
    //       occurring in section 5.1 of [1], as well as to the matrix
    //       Mbar^[-1] Nbar in section 5.3.
    //
    //     n is an integer variable.
    //       On entry n is the dimension of the problem.
    //       On exit n is unchanged.
    //
    //     nsub is an integer variable
    //       On entry nsub is the number of subspace variables in free set.
    //       On exit nsub is not changed.
    //
    //     ind is an integer array of dimension nsub.
    //       On entry ind specifies the indices of subspace variables.
    //       On exit ind is unchanged.
    //
    //     nenter is an integer variable.
    //       On entry nenter is the number of variables entering the
    //         free set.
    //       On exit nenter is unchanged.
    //
    //     ileave is an integer variable.
    //       On entry indx2(ileave),...,indx2(n) are the variables leaving
    //         the free set.
    //       On exit ileave is unchanged.
    //
    //     indx2 is an integer array of dimension n.
    //       On entry indx2(1),...,indx2(nenter) are the variables entering
    //         the free set, while indx2(ileave),...,indx2(n) are the
    //         variables leaving the free set.
    //       On exit indx2 is unchanged.
    //
    //     iupdat is an integer variable.
    //       On entry iupdat is the total number of BFGS updates made so far.
    //       On exit iupdat is unchanged.
    //
    //     updatd is a logical variable.
    //       On entry 'updatd' is true if the L-BFGS matrix is updatd.
    //       On exit 'updatd' is unchanged.
    //
    //     wn is a double precision array of dimension 2m x 2m.
    //       On entry wn is unspecified.
    //       On exit the upper triangle of wn stores the LEL^T factorization
    //         of the 2*col x 2*col indefinite matrix
    //                     [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
    //                     [L_a -R_z           theta*S'AA'S ]
    //
    //     wn1 is a double precision array of dimension 2m x 2m.
    //       On entry wn1 stores the lower triangular part of
    //                     [Y' ZZ'Y   L_a'+R_z']
    //                     [L_a+R_z   S'AA'S   ]
    //         in the previous iteration.
    //       On exit wn1 stores the corresponding updated matrices.
    //       The purpose of wn1 is just to store these inner products
    //       so they can be easily updated and inserted into wn.
    //
    //     m is an integer variable.
    //       On entry m is the maximum number of variable metric corrections
    //         used to define the limited memory matrix.
    //       On exit m is unchanged.
    //
    //     ws, wy, sy, and wtyy are double precision arrays;
    //     theta is a double precision variable;
    //     col is an integer variable;
    //     head is an integer variable.
    //       On entry they store the information defining the
    //                                          limited memory BFGS matrix:
    //         ws(n,m) stores S, a set of s-vectors;
    //         wy(n,m) stores Y, a set of y-vectors;
    //         sy(m,m) stores S'Y;
    //         wtyy(m,m) stores the Cholesky factorization
    //                                   of (theta*S'S+LD^(-1)L')
    //         theta is the scaling factor specifying B_0 = theta I;
    //         col is the number of variable metric corrections stored;
    //         head is the location of the 1st s- (or y-) vector in S (or Y).
    //       On exit they are unchanged.
    //
    //     info is an integer variable.
    //       On entry info is unspecified.
    //       On exit info =  0 for normal return;
    //                    = -1 when the 1st Cholesky factorization failed;
    //                    = -2 when the 2st Cholesky factorization failed.
    //
    //     Subprograms called:
    //
    //       Linpack ... dcopy, dpofa, dtrsl.
    //
    //
    //     References:
    //       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
    //       memory algorithm for bound constrained optimization'',
    //       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
    //
    //       [2] C. Zhu, R.H. Byrd, P. Lu, J. Nocedal, ``L-BFGS-B: a
    //       limited memory FORTRAN code for solving bound constrained
    //       optimization problems'', Tech. Report, NAM-11, EECS Department,
    //       Northwestern University, 1994.
    //
    //       (Postscript files of these papers are available via anonymous
    //        ftp to eecs.nwu.edu in the directory pub/lbfgs/lbfgs_bcm.)
    //
    //                           *  *  *
    //
    //     NEOS, November 1994. (Latest revision June 1996.)
    //     Optimization Technology Center.
    //     Argonne National Laboratory and Northwestern University.
    //     Written by
    //                        Ciyou Zhu
    //     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
    //
    //
    //     ************
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
                temp_int = m - (jy + 1);
                dcopy_(&temp_int, &wn1[(jy + 1) + 2*m*(jy + 1)], &one_int, &wn1[jy + 2*m*jy], &one_int);
                dcopy_(&temp_int, &wn1[(js + 1) + 2*m*(js + 1)], &one_int, &wn1[js + 2*m*js], &one_int);
                temp_int = m - 1;
                dcopy_(&temp_int, &wn1[(m + 1) + 2*m*(jy + 1)], &one_int, &wn1[m + 2*m*jy], &one_int);
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
            jpntr = (jpntr + 1) % m;
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
            ipntr = (ipntr + 1) % m;
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
        for (jy = 0; jy <= iy; jy++)
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
            for (k = ileave; k < n; k++)
            {
                k1 = indx2[k];
                temp3 = temp3 + wy[k1 + n*ipntr]*wy[k1 + n*jpntr];
                temp4 = temp4 + ws[k1 + n*ipntr]*ws[k1 + n*jpntr];
            }
            // 36

            wn1[iy + 2*m*jy] = wn1[iy + 2*m*jy] + temp1 - temp3;
            wn1[is + 2*m*js] = wn1[is + 2*m*js] - temp2 + temp4;
            jpntr = (jpntr + 1) % m;
        }
        // 40

        ipntr = (ipntr + 1) % m;
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
            jpntr = (jpntr + 1) % m;
        }
        // 55
        ipntr = (ipntr + 1) % m;
    }
    // 60

    // Form the upper triangle of WN = [D+Y' ZZ'Y/theta   -L_a'+R_z' ]
    //                                 [-L_a +R_z        S'AA'S*theta]

    m2 = 2*m;
    for (iy = 0; iy < col; iy++)
    {
        is = col + iy;
        is1 = m + iy;
        for (jy = 0; jy <= iy; jy++)
        {
            js = col + jy;
            js1 = m + jy;
            wn[jy + m2*iy] = wn1[iy + m2*jy] / theta;
            wn[js + m2*is] = wn1[is1 + m2*js1] * theta;
        }
        // 65
        for (jy = 0; jy < iy; jy++)
        {
            wn[jy + m2*is] = -wn1[is1 + m2*jy];
        }
        // 66
        for (jy = iy; jy < col; jy++)
        {
            wn[jy + m2*is] = wn1[is1 + m2*jy];
        }
        // 67
        wn[iy + m2*iy] = wn[iy + m2*iy] + sy[iy + m*iy];
    }
    // 70

    // Form the upper triangle of WN= [  LL'            L^-1(-L_a'+R_z')]
    //                                [(-L_a +R_z)L'^-1   S'AA'S*theta  ]
    //
    //    first Cholesky factor (1,1) block of wn to get LL'
    //                      with L' stored in the upper triangle of wn.

    // dpotrf(uplo, n, a, lda, info)
    dpotrf_(uplo, &col, wn, &m2, info);
    if (*info != 0) { *info = -1; return; }

    // Then form L^-1(-L_a'+R_z') in the (1,2) block.
    col2 = 2*col;

    // dtrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info)
    dtrtrs_(uplo, trans, diag, &col, &col, wn, &m2, &wn[m2*col], &m2, info);

    // Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the upper
    //  triangle of (2,2) block of wn.
    for (is = col; is < col2; is++)
    {
        for (js = is; js < col2; js++)
        {
            wn[is + m2*js] = wn[is + m2*js] + ddot_(&col, &wn[m2*is], &one_int, &wn[m2*js], &one_int);
        }
        // 74
    }
    // 72

    // dpotrf(uplo, n, a, lda, info)
    dpotrf_(uplo, &col, &wn[col + m2*col], &m2, info);
    if (*info != 0) { *info = -2; }

    return;
}


void
formt(int m, double* wt, double* sy, double* ss, int col,
      double theta, int* info)
{
    //     ************
    //
    //     Subroutine formt
    //
    //       This subroutine forms the upper half of the pos. def. and symm.
    //         T = theta*SS + L*D^(-1)*L', stores T in the upper triangle
    //         of the array wt, and performs the Cholesky factorization of T
    //         to produce J*J', with J' stored in the upper triangle of wt.
    //
    //     Subprograms called:
    //
    //       Linpack ... dpofa.
    //
    //
    //                           *  *  *
    //
    //     NEOS, November 1994. (Latest revision June 1996.)
    //     Optimization Technology Center.
    //     Argonne National Laboratory and Northwestern University.
    //     Written by
    //                        Ciyou Zhu
    //     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
    //
    //
    //     ************
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
            k1 = (i > j ? j : i);
            ddum = 0.0;
            for (k = 0; k < k1; k++)
            {
                ddum = ddum + ((sy[i + m*k] * sy[j + m*k]) / sy[k + m*k]);
            }
            // 53
            wt[i + m*j] = ddum + theta*ss[i + m*j];
        }
        // 54
    }
    // 55

    // Cholesky factorize T to J*J' with J' stored in the upper triangle of wt.
    // dpotrf(uplo, n, a, lda, info)
    dpotrf_(uplo, &col, wt, &m, info);
    if (*info != 0) { *info = -3; }

    return;
}


void
freev(int n, int* nfree, int* idx, int* nenter, int* ileave, int* idx2,
      int* iwhere, int* wrk, int updatd, int cnstnd,
      int iter)
{
    //     ************
    //
    //     Subroutine freev
    //
    //     This subroutine counts the entering and leaving variables when
    //       iter > 0, and finds the index set of free and active variables
    //       at the GCP.
    //
    //     cnstnd is a logical variable indicating whether bounds are present
    //
    //     index is an integer array of dimension n
    //       for i=1,...,nfree, index(i) are the indices of free variables
    //       for i=nfree+1,...,n, index(i) are the indices of bound variables
    //       On entry after the first iteration, index gives
    //         the free variables at the previous iteration.
    //       On exit it gives the free variables based on the determination
    //         in cauchy using the array iwhere.
    //
    //     indx2 is an integer array of dimension n
    //       On entry indx2 is unspecified.
    //       On exit with iter>0, indx2 indicates which variables
    //          have changed status since the previous iteration.
    //       For i= 1,...,nenter, indx2(i) have changed from bound to free.
    //       For i= ileave+1,...,n, indx2(i) have changed from free to bound.
    //
    //
    //                           *  *  *
    //
    //     NEOS, November 1994. (Latest revision June 1996.)
    //     Optimization Technology Center.
    //     Argonne National Laboratory and Northwestern University.
    //     Written by
    //                        Ciyou Zhu
    //     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
    //
    //
    //     ************
    int i, iact, k;

    // in formk
    // nenter will be used as exclusive upper bound
    // ileave will be used as the inclusive lower bound

    *nenter = 0;
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
                idx2[*nenter] = k;
                *nenter += 1;
            }
        }
        // 22
    }

    *wrk = ((*ileave < n) || (*nenter > 0) || updatd);

    // Find the index set of free and active variables at the GCP.
    *nfree = 0;
    iact = n;
    for (i = 0; i < n; i++)
    {
        if (iwhere[i] <= 0)
        {
            idx[*nfree] = i;
            *nfree += 1;
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
    //     ************
    //
    //     Subroutine hpsolb
    //
    //     This subroutine sorts out the least element of t, and puts the
    //       remaining elements of t in a heap.
    //
    //     n is an integer variable.
    //       On entry n is the dimension of the arrays t and iorder.
    //       On exit n is unchanged.
    //
    //     t is a double precision array of dimension n.
    //       On entry t stores the elements to be sorted,
    //       On exit t(n) stores the least elements of t, and t(1) to t(n-1)
    //         stores the remaining elements in the form of a heap.
    //
    //     iorder is an integer array of dimension n.
    //       On entry iorder(i) is the index of t(i).
    //       On exit iorder(i) is still the index of t(i), but iorder may be
    //         permuted in accordance with t.
    //
    //     iheap is an integer variable specifying the task.
    //       On entry iheap should be set as follows:
    //         iheap .eq. 0 if t(1) to t(n) is not in the form of a heap,
    //         iheap .ne. 0 if otherwise.
    //       On exit iheap is unchanged.
    //
    //
    //     References:
    //       Algorithm 232 of CACM (J. W. J. Williams): HEAPSORT.
    //
    //                           *  *  *
    //
    //     NEOS, November 1994. (Latest revision June 1996.)
    //     Optimization Technology Center.
    //     Argonne National Laboratory and Northwestern University.
    //     Written by
    //                        Ciyou Zhu
    //     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
    //
    //     ************
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
        ddum = t[n];
        indxin = iorder[n];

        // Restore the heap
        while (1)
        {
            j = i + i;
            if (j <= n)
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
        t[n] = out;
        iorder[n] = indxout;
    }

    return;
}


void
lnsrlb(int n, double* l, double* u, int* nbd, double* x,
       double f, double* fold, double *gd, double *gdold, double* g,
       double* d, double* r, double* t, double* z, double* stp, double* dnorm,
       double* dtd, double* xstep, double* stpmx, int iter, int* ifun,
       int* iback, int* nfgv, int* info, int* task, int* task_msg,
       int boxed, int cnstnd, int* isave, double* dsave, int* temp_task,
       int* temp_taskmsg)
{
    //     **********
    //
    //     Subroutine lnsrlb
    //
    //     This subroutine calls subroutine dcsrch from the Minpack2 library
    //       to perform the line search.  Subroutine dscrch is safeguarded so
    //       that all trial points lie within the feasible region.
    //
    //     Subprograms called:
    //
    //       Minpack2 Library ... dcsrch.
    //
    //       Linpack ... dtrsl, ddot.
    //
    //                           *  *  *
    //
    //     NEOS, November 1994. (Latest revision June 1996.)
    //     Optimization Technology Center.
    //     Argonne National Laboratory and Northwestern University.
    //     Written by
    //                        Ciyou Zhu
    //     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
    //
    //
    //     **********
    int i, one_int = 1;
    double a1, a2, ftol = 1e-3, gtol = 0.9, xtol = 0.1;
    if (*task_msg == FG_LNSRCH) { goto LINE556; }

    *dnorm = dnrm2_(&n, d, &one_int);
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

    *fold = f;  // Later used in mainlb, see control flow after returning from this function.
    *ifun = 0;
    *iback = 0;
    *temp_task = START;
    *temp_taskmsg = NO_MSG;
LINE556:
    *gd = ddot_(&n, g, &one_int, d, &one_int);
    if (*ifun == 0)
    {
        *gdold = *gd;
        if (*gd >= 0.0)
        {
            // The directional derivative >= 0, line search is impossible.
            *info = -4;
            return;
        }
    }
    dcsrch(f, *gd, stp, ftol, gtol, xtol, 0.0, *stpmx, temp_task, temp_taskmsg,
           isave, dsave);

    *xstep = (*stp) * (*dnorm);
    if ((*temp_task != CONVERGENCE) && (*temp_task != WARNING))
    {
        *task = FG;
        *task_msg = FG_LNSRCH;
        *ifun += 1;
        *nfgv += 1;
        *iback = *ifun - 1;
        if (*stp == 1.0)
        {
            // for (i = 0; i < n; i++) { x[i] = z[i]; }
            dcopy_(&n, z, &one_int, x, &one_int);
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
    //     ************
    //
    //     Subroutine matupd
    //
    //       This subroutine updates matrices WS and WY, and forms the
    //         middle matrix in B.
    //
    //     Subprograms called:
    //
    //       Linpack ... dcopy, ddot.
    //
    //
    //                           *  *  *
    //
    //     NEOS, November 1994. (Latest revision June 1996.)
    //     Optimization Technology Center.
    //     Argonne National Laboratory and Northwestern University.
    //     Written by
    //                        Ciyou Zhu
    //     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
    //
    //
    //     ************
    int j, pointr, temp_int, one_int = 1;

    // Set pointers for matrices WS and WY.
    if (iupdat <= m)
    {
        *col = iupdat;
        *itail = (*head + iupdat - 1) % m;
    } else {
        *itail = (*itail + 1) % m;
        *head = (*head + 1) % m;
    }

    // Update matrices WS and WY.
    dcopy_(&n, d, &one_int, &ws[(*itail) * n], &one_int);
    dcopy_(&n, r, &one_int, &wy[(*itail) * n], &one_int);

    // Set theta=yy/ys.
    *theta = rr / dr;

    // Form the middle matrix in B.
    // update the upper triangle of SS (m, m), and the lower triangle of SY (m, m):
    if (iupdat > m)
    {
        for (j = 1; j < *col; j++)
        {
            dcopy_(&j, &ss[1 + m*j], &one_int, &ss[m*(j-1)], &one_int);
            temp_int = *col - j;
            dcopy_(&temp_int, &sy[j + m*j], &one_int, &sy[(j-1) + m*(j-1)], &one_int);
        }
        // 50
    }

    // add new information: the last row of SY and the last column of SS
    pointr = *head;
    for (j = 0; j < *col - 1; j++)
    {
        sy[*col - 1 + m*j] = ddot_(&n, d, &one_int, &wy[pointr*n], &one_int);
        ss[j + m*(*col - 1)] = ddot_(&n, &ws[pointr*n], &one_int, &d[0], &one_int);
        pointr = (pointr + 1) % m;
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
    //     ************
    //
    //     Subroutine projgr
    //
    //     This subroutine computes the infinity norm of the projected
    //       gradient.
    //
    //
    //                           *  *  *
    //
    //     NEOS, November 1994. (Latest revision June 1996.)
    //     Optimization Technology Center.
    //     Argonne National Laboratory and Northwestern University.
    //     Written by
    //                        Ciyou Zhu
    //     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
    //
    //
    //     ************
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
    //     **********************************************************************
    //
    //     This routine contains the major changes in the updated version.
    //     The changes are described in the accompanying paper
    //
    //      Jose Luis Morales, Jorge Nocedal
    //      "Remark On Algorithm 788: L-BFGS-B: Fortran Subroutines for Large-Scale
    //       Bound Constrained Optimization". Decemmber 27, 2010.
    //
    //             J.L. Morales  Departamento de Matematicas,
    //                           Instituto Tecnologico Autonomo de Mexico
    //                           Mexico D.F.
    //
    //             J, Nocedal    Department of Electrical Engineering and
    //                           Computer Science.
    //                           Northwestern University. Evanston, IL. USA
    //
    //                           January 17, 2011
    //
    //      **********************************************************************
    //
    //
    //     Subroutine subsm
    //
    //     Given xcp, l, u, r, an index set that specifies
    //       the active set at xcp, and an l-BFGS matrix B
    //       (in terms of WY, WS, SY, WT, head, col, and theta),
    //       this subroutine computes an approximate solution
    //       of the subspace problem
    //
    //       (P)   min Q(x) = r'(x-xcp) + 1/2 (x-xcp)' B (x-xcp)
    //
    //             subject to l<=x<=u
    //                       x_i=xcp_i for all i in A(xcp)
    //
    //       along the subspace unconstrained Newton direction
    //
    //          d = -(Z'BZ)^(-1) r.
    //
    //       The formula for the Newton direction, given the L-BFGS matrix
    //       and the Sherman-Morrison formula, is
    //
    //          d = (1/theta)r + (1/theta*2) Z'WK^(-1)W'Z r.
    //
    //       where
    //                 K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
    //                     [L_a -R_z           theta*S'AA'S ]
    //
    //     Note that this procedure for computing d differs
    //     from that described in [1]. One can show that the matrix K is
    //     equal to the matrix M^[-1]N in that paper.
    //
    //     n is an integer variable.
    //       On entry n is the dimension of the problem.
    //       On exit n is unchanged.
    //
    //     m is an integer variable.
    //       On entry m is the maximum number of variable metric corrections
    //         used to define the limited memory matrix.
    //       On exit m is unchanged.
    //
    //     nsub is an integer variable.
    //       On entry nsub is the number of free variables.
    //       On exit nsub is unchanged.
    //
    //     ind is an integer array of dimension nsub.
    //       On entry ind specifies the coordinate indices of free variables.
    //       On exit ind is unchanged.
    //
    //     l is a double precision array of dimension n.
    //       On entry l is the lower bound of x.
    //       On exit l is unchanged.
    //
    //     u is a double precision array of dimension n.
    //       On entry u is the upper bound of x.
    //       On exit u is unchanged.
    //
    //     nbd is a integer array of dimension n.
    //       On entry nbd represents the type of bounds imposed on the
    //         variables, and must be specified as follows:
    //         nbd(i)=0 if x(i) is unbounded,
    //                1 if x(i) has only a lower bound,
    //                2 if x(i) has both lower and upper bounds, and
    //                3 if x(i) has only an upper bound.
    //       On exit nbd is unchanged.
    //
    //     x is a double precision array of dimension n.
    //       On entry x specifies the Cauchy point xcp.
    //       On exit x(i) is the minimizer of Q over the subspace of
    //                                                        free variables.
    //
    //     d is a double precision array of dimension n.
    //       On entry d is the reduced gradient of Q at xcp.
    //       On exit d is the Newton direction of Q.
    //
    //    xp is a double precision array of dimension n.
    //       used to safeguard the projected Newton direction
    //
    //    xx is a double precision array of dimension n
    //       On entry it holds the current iterate
    //       On output it is unchanged
    //    gg is a double precision array of dimension n
    //       On entry it holds the gradient at the current iterate
    //       On output it is unchanged
    //
    //     ws and wy are double precision arrays;
    //     theta is a double precision variable;
    //     col is an integer variable;
    //     head is an integer variable.
    //       On entry they store the information defining the
    //                                          limited memory BFGS matrix:
    //         ws(n,m) stores S, a set of s-vectors;
    //         wy(n,m) stores Y, a set of y-vectors;
    //         theta is the scaling factor specifying B_0 = theta I;
    //         col is the number of variable metric corrections stored;
    //         head is the location of the 1st s- (or y-) vector in S (or Y).
    //       On exit they are unchanged.
    //
    //     iword is an integer variable.
    //       On entry iword is unspecified.
    //       On exit iword specifies the status of the subspace solution.
    //         iword = 0 if the solution is in the box,
    //                 1 if some bound is encountered.
    //
    //     wv is a double precision working array of dimension 2m.
    //
    //     wn is a double precision array of dimension 2m x 2m.
    //       On entry the upper triangle of wn stores the LEL^T factorization
    //         of the indefinite matrix
    //
    //              K = [-D -Y'ZZ'Y/theta     L_a'-R_z'  ]
    //                  [L_a -R_z           theta*S'AA'S ]
    //                                                    where E = [-I  0]
    //                                                              [ 0  I]
    //       On exit wn is unchanged.
    //
    //     iprint is an INTEGER variable that must be set by the user.
    //       It controls the frequency and type of output generated:
    //        iprint<0    no output is generated;
    //        iprint=0    print only one line at the last iteration;
    //        0<iprint<99 print also f and |proj g| every iprint iterations;
    //        iprint=99   print details of every iteration except n-vectors;
    //        iprint=100  print also the changes of active set and final x;
    //        iprint>100  print details of every iteration including x and g;
    //
    //
    //     info is an integer variable.
    //       On entry info is unspecified.
    //       On exit info = 0       for normal return,
    //                    = nonzero for abnormal return
    //                                  when the matrix K is ill-conditioned.
    //
    //     Subprograms called:
    //
    //       Linpack dtrsl.
    //
    //
    //     References:
    //
    //       [1] R. H. Byrd, P. Lu, J. Nocedal and C. Zhu, ``A limited
    //       memory algorithm for bound constrained optimization'',
    //       SIAM J. Scientific Computing 16 (1995), no. 5, pp. 1190--1208.
    //
    //
    //
    //                           *  *  *
    //
    //     NEOS, November 1994. (Latest revision June 1996.)
    //     Optimization Technology Center.
    //     Argonne National Laboratory and Northwestern University.
    //     Written by
    //                        Ciyou Zhu
    //     in collaboration with R.H. Byrd, P. Lu-Chen and J. Nocedal.
    //
    //
    //     ************
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
            temp1 += wy[k + n*pointr]*d[j];
            temp2 += ws[k + n*pointr]*d[j];
        }
        wv[i] = temp1;
        wv[col + i] = theta*temp2;
        pointr = (pointr + 1) % m;
    }

    // Compute wv:=K^(-1)wv.
    m2 = 2*m;
    col2 = 2*col;

    // dtrtrs(uplo, trans, diag, n, nrhs, a, lda, b, ldb, info)
    dtrtrs_(uplo, trans, diag, &col2, &one_int, wn, &m2, wv, &m2, info);
    if (*info != 0) { return; }

    for (i = 0; i < col; i++) { wv[i] = -wv[i]; }

    trans = "N";
    dtrtrs_(uplo, trans, diag, &col2, &one_int, wn, &m2, wv, &m2, info);
    if (*info != 0) { return; }

    // Compute d = (1/theta)d + (1/theta**2)Z'W wv.
    pointr = head;
    for (jy = 0; jy < col; jy++)
    {
        js = col + jy;
        for (i = 0; i < nsub; i++)
        {
            k = ind[i];
            d[i] = d[i] + (wy[k + n*pointr] * wv[jy] / theta) +
                          (ws[k + n*pointr] * wv[js]);
        }
        // 30
        pointr = (pointr + 1) % m;
    }
    // 40

    temp = 1.0 / theta;
    // n, da, dx, incx
    dscal_(&nsub, &temp, d, &one_int);

    // -----------------------------------------------------
    // Let us try the projection, d is the Newton direction.

    *iword = 0;
    // for (i = 0; i < n; i++) { xp[i] = x[i]; }
    dcopy_(&n, x, &one_int, xp, &one_int);

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
        dcopy_(&n, xp, &one_int, x, &one_int);
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
    //     **********
    //
    //     Subroutine dcsrch
    //
    //     This subroutine finds a step that satisfies a sufficient
    //     decrease condition and a curvature condition.
    //
    //     Each call of the subroutine updates an interval with
    //     endpoints stx and sty. The interval is initially chosen
    //     so that it contains a minimizer of the modified function
    //
    //           psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).
    //
    //     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
    //     interval is chosen so that it contains a minimizer of f.
    //
    //     The algorithm is designed to find a step that satisfies
    //     the sufficient decrease condition
    //
    //           f(stp) <= f(0) + ftol*stp*f'(0),
    //
    //     and the curvature condition
    //
    //           abs(f'(stp)) <= gtol*abs(f'(0)).
    //
    //     If ftol is less than gtol and if, for example, the function
    //     is bounded below, then there is always a step which satisfies
    //     both conditions.
    //
    //     If no step can be found that satisfies both conditions, then
    //     the algorithm stops with a warning. In this case stp only
    //     satisfies the sufficient decrease condition.
    //
    //     A typical invocation of dcsrch has the following outline:
    //
    //     task = 'START'
    //  10 continue
    //        call dcsrch( ... )
    //        if (task .eq. 'FG') then
    //           Evaluate the function and the gradient at stp
    //           goto 10
    //           end if
    //
    //     NOTE: The user must no alter work arrays between calls.
    //
    //     The subroutine statement is
    //
    //        subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,
    //                          task,isave,dsave)
    //     where
    //
    //       f is a double precision variable.
    //         On initial entry f is the value of the function at 0.
    //            On subsequent entries f is the value of the
    //            function at stp.
    //         On exit f is the value of the function at stp.
    //
    //       g is a double precision variable.
    //         On initial entry g is the derivative of the function at 0.
    //            On subsequent entries g is the derivative of the
    //            function at stp.
    //         On exit g is the derivative of the function at stp.
    //
    //       stp is a double precision variable.
    //         On entry stp is the current estimate of a satisfactory
    //            step. On initial entry, a positive initial estimate
    //            must be provided.
    //         On exit stp is the current estimate of a satisfactory step
    //            if task = 'FG'. If task = 'CONV' then stp satisfies
    //            the sufficient decrease and curvature condition.
    //
    //       ftol is a double precision variable.
    //         On entry ftol specifies a nonnegative tolerance for the
    //            sufficient decrease condition.
    //         On exit ftol is unchanged.
    //
    //       gtol is a double precision variable.
    //         On entry gtol specifies a nonnegative tolerance for the
    //            curvature condition.
    //         On exit gtol is unchanged.
    //
    //       xtol is a double precision variable.
    //         On entry xtol specifies a nonnegative relative tolerance
    //            for an acceptable step. The subroutine exits with a
    //            warning if the relative difference between sty and stx
    //            is less than xtol.
    //         On exit xtol is unchanged.
    //
    //       stpmin is a double precision variable.
    //         On entry stpmin is a nonnegative lower bound for the step.
    //         On exit stpmin is unchanged.
    //
    //       stpmax is a double precision variable.
    //         On entry stpmax is a nonnegative upper bound for the step.
    //         On exit stpmax is unchanged.
    //
    //       task is a character variable of length at least 60.
    //         On initial entry task must be set to 'START'.
    //         On exit task indicates the required action:
    //
    //            If task(1:2) = 'FG' then evaluate the function and
    //            derivative at stp and call dcsrch again.
    //
    //            If task(1:4) = 'CONV' then the search is successful.
    //
    //            If task(1:4) = 'WARN' then the subroutine is not able
    //            to satisfy the convergence conditions. The exit value of
    //            stp contains the best point found during the search.
    //
    //            If task(1:5) = 'ERROR' then there is an error in the
    //            input arguments.
    //
    //         On exit with convergence, a warning or an error, the
    //            variable task contains additional information.
    //
    //       isave is an integer work array of dimension 2.
    //
    //       dsave is a double precision work array of dimension 13.
    //
    //     Subprograms called
    //
    //       MINPACK-2 ... dcstep
    //
    //     MINPACK-1 Project. June 1983.
    //     Argonne National Laboratory.
    //     Jorge J. More' and David J. Thuente.
    //
    //     MINPACK-2 Project. October 1993.
    //     Argonne National Laboratory and University of Minnesota.
    //     Brett M. Averick, Richard G. Carter, and Jorge J. More'.
    //
    //     **********
    int brackt, stage;
    double finit, ftest, fm, fx, fxm, fy, fym, ginit, gtest, gm, gx, gxm, gy;
    double gym, stx, sty, stmin, stmax, width, width1;
    double xtrapl = 1.1;
    double xtrapu = 4.0;

    // Initialization block.
    if (*task == START)
    {
        // Check the input arguments for errors.
        if (*stp < stpmin) { *task = ERROR; *task_msg = ERROR_STP1; }
        if (*stp > stpmax) { *task = ERROR; *task_msg = ERROR_STP2; }
        if (g >= 0.0) { *task = ERROR; *task_msg = ERROR_INITG; }
        if (ftol < 0.0) { *task = ERROR; *task_msg = ERROR_FTOL; }
        if (gtol < 0.0) { *task = ERROR; *task_msg = ERROR_GTOL; }
        if (xtol < 0.0) { *task = ERROR; *task_msg = ERROR_XTOL; }
        if (stpmin < 0.0) { *task = ERROR; *task_msg = ERROR_STPMIN; }
        if (stpmax < stpmin) { *task = ERROR; *task_msg = ERROR_STPMAX; }

        // Exit if there are errors on input.
        if (*task == ERROR) { return; }

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
        stage  = isave[1];
        ginit  = dsave[0];
        gtest  = dsave[1];
        gx     = dsave[2];
        gy     = dsave[3];
        finit  = dsave[4];
        fx     = dsave[5];
        fy     = dsave[6];
        stx    = dsave[7];
        sty    = dsave[8];
        stmin  = dsave[9];
        stmax  = dsave[10];
        width  = dsave[11];
        width1 = dsave[12];

    }

    // If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the algorithm
    // enters the second stage.
    ftest = finit + (*stp)*gtest;
    if ((stage == 1) && (f <= ftest) && (g >= 0.0)) { stage = 2; }

    // Test for warnings.
    if ((brackt) && ((*stp <= stmin) || (*stp >= stmax)))
    {
        *task = WARNING;
        *task_msg = WARN_ROUND;
    }
    if ((brackt) && ((stmax - stmin) <= xtol * stmax))
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

    isave[1]  = stage;
    dsave[0]  = ginit;
    dsave[1]  = gtest;
    dsave[2]  = gx;
    dsave[3]  = gy;
    dsave[4]  = finit;
    dsave[5]  = fx;
    dsave[6]  = fy;
    dsave[7]  = stx;
    dsave[8]  = sty;
    dsave[9]  = stmin;
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
    //     **********
    //
    //     Subroutine dcstep
    //
    //     This subroutine computes a safeguarded step for a search
    //     procedure and updates an interval that contains a step that
    //     satisfies a sufficient decrease and a curvature condition.
    //
    //     The parameter stx contains the step with the least function
    //     value. If brackt is set to .true. then a minimizer has
    //     been bracketed in an interval with endpoints stx and sty.
    //     The parameter stp contains the current step.
    //     The subroutine assumes that if brackt is set to .true. then
    //
    //           min(stx,sty) < stp < max(stx,sty),
    //
    //     and that the derivative at stx is negative in the direction
    //     of the step.
    //
    //     The subroutine statement is
    //
    //       subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,
    //                         stpmin,stpmax)
    //
    //     where
    //
    //       stx is a double precision variable.
    //         On entry stx is the best step obtained so far and is an
    //            endpoint of the interval that contains the minimizer.
    //         On exit stx is the updated best step.
    //
    //       fx is a double precision variable.
    //         On entry fx is the function at stx.
    //         On exit fx is the function at stx.
    //
    //       dx is a double precision variable.
    //         On entry dx is the derivative of the function at
    //            stx. The derivative must be negative in the direction of
    //            the step, that is, dx and stp - stx must have opposite
    //            signs.
    //         On exit dx is the derivative of the function at stx.
    //
    //       sty is a double precision variable.
    //         On entry sty is the second endpoint of the interval that
    //            contains the minimizer.
    //         On exit sty is the updated endpoint of the interval that
    //            contains the minimizer.
    //
    //       fy is a double precision variable.
    //         On entry fy is the function at sty.
    //         On exit fy is the function at sty.
    //
    //       dy is a double precision variable.
    //         On entry dy is the derivative of the function at sty.
    //         On exit dy is the derivative of the function at the exit sty.
    //
    //       stp is a double precision variable.
    //         On entry stp is the current step. If brackt is set to .true.
    //            then on input stp must be between stx and sty.
    //         On exit stp is a new trial step.
    //
    //       fp is a double precision variable.
    //         On entry fp is the function at stp
    //         On exit fp is unchanged.
    //
    //       dp is a double precision variable.
    //         On entry dp is the derivative of the function at stp.
    //         On exit dp is unchanged.
    //
    //       brackt is an logical variable.
    //         On entry brackt specifies if a minimizer has been bracketed.
    //            Initially brackt must be set to .false.
    //         On exit brackt specifies if a minimizer has been bracketed.
    //            When a minimizer is bracketed brackt is set to .true.
    //
    //       stpmin is a double precision variable.
    //         On entry stpmin is a lower bound for the step.
    //         On exit stpmin is unchanged.
    //
    //       stpmax is a double precision variable.
    //         On entry stpmax is an upper bound for the step.
    //         On exit stpmax is unchanged.
    //
    //     MINPACK-1 Project. June 1983
    //     Argonne National Laboratory.
    //     Jorge J. More' and David J. Thuente.
    //
    //     MINPACK-2 Project. October 1993.
    //     Argonne National Laboratory and University of Minnesota.
    //     Brett M. Averick and Jorge J. More'.
    //
    //     **********
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
        stpq = *stp + (dp / (dp - *dx)) * (*stx - *stp);
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
