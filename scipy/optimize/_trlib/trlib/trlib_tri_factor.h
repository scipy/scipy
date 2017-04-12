/* MIT License
 *
 * Copyright (c) 2016--2017 Felix Lenders
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#ifndef TRLIB_TRI_FACTOR_H
#define TRLIB_TRI_FACTOR_H

#define TRLIB_TTR_CONV_BOUND    (0)
#define TRLIB_TTR_CONV_INTERIOR (1)
#define TRLIB_TTR_HARD          (2)
#define TRLIB_TTR_NEWTON_BREAK  (3)
#define TRLIB_TTR_HARD_INIT_LAM (4)
#define TRLIB_TTR_ITMAX         (-1)
#define TRLIB_TTR_FAIL_FACTOR   (-2)
#define TRLIB_TTR_FAIL_LINSOLVE (-3)
#define TRLIB_TTR_FAIL_EIG      (-4)
#define TRLIB_TTR_FAIL_LM       (-5)
#define TRLIB_TTR_FAIL_HARD     (-10)

/** Solves tridiagonal trust region subproblem
 *
 *  Computes minimizer to
 * 
 *  :math:`\min \frac 12 \langle h, T h \rangle + \langle g, h \rangle`
 *  subject to the trust region constraint :math:`\Vert h \Vert \le r`,
 *  where :math:`T \in \mathbb R^{n \times n}` is symmetric tridiagonal.
 *
 *  Let :math:`T = \begin{pmatrix} T_1 & & \\ & \ddots & \\ & & T_\ell \end{pmatrix}`
 *  be composed into irreducible blocks :math:`T_i`.
 * 
 *  The minimizer is a global minimizer (modulo floating point).
 * 
 *  The algorithm is the Moré-Sorensen-Method as decribed as Algorithm 5.2 in [Gould1999]_.
 *
 *  **Convergence**
 *  
 *  Exit with success is reported in several cases:
 *
 *      - interior solution: the stationary point :math:`T h = - g` is suitable, iterative refinement is used in the solution of this system.
 *      - hard case: the smallest eigenvalue :math:`-\theta` was found to be degenerate and :math:`h = v + \alpha w` is returned with :math:`v` solution of :math:`(T + \theta I) v = -g`, :math:`w` eigenvector corresponding to :math:`- \theta` and :math:`\alpha` chosen to satisfy the trust region constraint and being a minimizer.
 *      - boundary and hard case: convergence in Newton iteration is reported if :math:`\Vert h(\lambda) \Vert - r \le \texttt{tol}\_\texttt{rel} \, r` with :math:`(T+\lambda I) h(\lambda) = -g`, Newton breakdown reported if :math:`\vert d\lambda \vert \le \texttt{macheps} \vert \lambda \vert`.
 *
 *  :param nirblk: number of irreducible blocks :math:`\ell`, ensure :math:`\ell > 0`
 *  :type nirblk: trlib_int_t, input
 *  :param irblk: pointer to indices of irreducible blocks, length :c:data:`nirblk+1`:
 *
 *                 - :c:data:`irblk[i]` is the start index of block :math:`i` in :c:data:`diag` and :c:data:`offdiag`
 *                 - :c:data:`irblk[i+1] - 1` is the stop index of block :math:`i`
 *                 - :c:data:`irblk[i+1] - irred[i]` the dimension :math:`n_\ell` of block :math:`i`
 *                 - :c:data:`irblk[nirred]` the dimension of :math:`T`
 *
 *  :type irblk: trlib_int_t, input
 *  :param diag: pointer to array holding diagonal of :math:`T`, length :c:data:`irblk[nirblk]`
 *  :type diag: trlib_flt_t, input
 *  :param offdiag: pointer to array holding offdiagonal of :math:`T`, length :c:data:`irblk[nirblk]`
 *  :type offdiag: trlib_flt_t, input
 *  :param neglin: pointer to array holding :math:`-g`, length :c:data:`n`
 *  :type neglin: trlib_flt_t, input
 *  :param radius: trust region constraint radius :math:`r`
 *  :type radius: trlib_flt_t, input
 *  :param itmax: maximum number of Newton iterations
 *  :type itmax: trlib_int_t, input
 *  :param tol_rel: relative stopping tolerance for residual in Newton iteration for :c:data:`lam`, good default may be :math:`\texttt{macheps}` (:c:macro:`TRLIB_EPS`)
 *  :type tol_rel: trlib_flt_t, input
 *  :param tol_newton_tiny: stopping tolerance for step in Newton iteration for :c:data:`lam`, good default may be :math:`10 \texttt{macheps}^{.75}`
 *  :type tol_newton_tiny: trlib_flt_t, input
 *  :param pos_def: set ``1`` if you know :math:`T` to be positive definite, otherwise ``0``
 *  :type pos_def: trlib_int_t, input
 *  :param equality: set ``1`` if you want to enfore trust region constraint as equality, otherwise ``0``
 *  :type equality: trlib_int_t, input
 *  :param warm0: set ``1`` if you provide a valid value in :c:data:`lam0`, otherwise ``0``
 *  :type warm0: trlib_int_t, input/output
 *  :param lam0: Lagrange multiplier such that :math:`T_0+ \mathtt{lam0} I` positive definite
 *                 and :math:`(T_0+ \mathtt{lam0} I) \mathtt{sol0} = \mathtt{neglin}`
 *                 - on entry: estimate suitable for warmstart
 *                 - on exit: computed multiplier
 *  :type lam0: trlib_flt_t, input/output
 *  :param warm: set ``1`` if you provide a valid value in :c:data:`lam`, otherwise ``0``
 *  :type warm: trlib_int_t, input/output
 *  :param lam: Lagrange multiplier such that :math:`T+ \mathtt{lam} I` positive definite
 *                 and :math:`(T+ \mathtt{lam} I) \mathtt{sol} = \mathtt{neglin}`
 *                 - on entry: estimate suitable for warmstart
 *                 - on exit: computed multiplier
 *  :type lam: trlib_flt_t, input/output
 *  :param warm_leftmost:
 *                 - on entry: set ``1`` if you provide a valid value in leftmost section in :c:data:`fwork`
 *                 - on exit: ``1`` if leftmost section in :c:data:`fwork` holds valid exit value, otherwise ``0``
 *  :type warm_leftmost: trlib_int_t, input/output
 *  :param ileftmost: index to block with smallest leftmost eigenvalue
 *  :type ileftmost: trlib_int_t, input/output
 *  :param leftmost: array holding leftmost eigenvalues of blocks
 *  :type leftmost: trlib_flt_t, input/output
 *  :param warm_fac0: set ``1`` if you provide a valid factoricization in :c:data:`diag_fac0`, :c:data:`offdiag_fac0`
 *  :type warm_fac0: trlib_int_t, input/output
 *  :param diag_fac0: pointer to array holding diagonal of Cholesky factorization of :math:`T_0 + \texttt{lam0} I`, length :c:data:`irblk[1]`
 *                 - on entry: factorization corresponding to provided estimation :c:data:`lam0` on entry
 *                 - on exit: factorization corresponding to computed multiplier :c:data:`lam0`
 *  :type diag_fac0: trlib_flt_t, input/output
 *  :param offdiag_fac0: pointer to array holding offdiagonal of Cholesky factorization of :math:`T_0 + \texttt{lam0} I`, length :c:data:`irblk[1]-1`
 *                 - on entry: factorization corresponding to provided estimation :c:data:`lam0` on entry
 *                 - on exit: factorization corresponding to computed multiplier :c:data:`lam0`
 *  :type offdiag_fac0: trlib_flt_t, input/output
 *  :param warm_fac: set ``1`` if you provide a valid factoricization in :c:data:`diag_fac`, :c:data:`offdiag_fac`
 *  :type warm_fac: trlib_int_t, input/output
 *  :param diag_fac: pointer to array of length :c:data:`n` that holds the following:
 *                 - let :math:`j = \texttt{ileftmost}` and :math:`\theta = - \texttt{leftmost[ileftmost]}`
 *                 - on position :math:`0, \ldots, \texttt{irblk[1]}`: diagonal of factorization of :math:`T_0 + \theta I`
 *                 - other positions have to be allocated
 *  :type diag_fac: trlib_flt_t, input/output
 *  :param offdiag_fac: pointer to array of length :c:data:`n-1` that holds the following:
 *                 - let :math:`j = \texttt{ileftmost}` and :math:`\theta = - \texttt{leftmost[ileftmost]}`
 *                 - on position :math:`0, \ldots, \texttt{irblk[1]}-1`: offdiagonal of factorization of :math:`T_0 + \theta I`
 *                 - other positions have to be allocated
 *  :type offdiag_fac: trlib_flt_t, input/output
 *  :param sol0: pointer to array holding solution, length :c:data:`irblk[1]`
 *  :type sol0: trlib_flt_t, input/output
 *  :param sol: pointer to array holding solution, length :c:data:`n`
 *  :type sol: trlib_flt_t, input/output
 *  :param ones: array with every value ``1.0``, length :c:data:`n`
 *  :type ones: trlib_flt_t, input
 *  :param fwork: floating point workspace, must be allocated memory on input of size :c:func:`trlib_tri_factor_memory_size` (call with argument :c:data:`n`) and can be discarded on output, memory layout:
 *
 *      ============== ============== ==========================================================
 *      start          end (excl)     description
 *      ============== ============== ==========================================================
 *                 0     :c:data:`n`  auxiliary vector
 *        :c:data:`n`  2 :c:data:`n`  holds diagonal of :math:`T + \lambda I`
 *      2 :c:data:`n`  4 :c:data:`n`  workspace for iterative refinement
 *      ============== ============== ==========================================================
 *
 *  :type fwork: trlib_flt_t, input/output
 *  :param refine: set to ``1`` if iterative refinement should be used on solving linear systems, otherwise to ``0``
 *  :type refine: trlib_int_t, input
 *  :param verbose: determines the verbosity level of output that is written to :c:data:`fout`
 *  :type verbose: trlib_int_t, input
 *  :param unicode: set to ``1`` if :c:data:`fout` can handle unicode, otherwise to ``0``
 *  :type unicode: trlib_int_t, input
 *  :param prefix: string that is printed before iteration output
 *  :type prefix: char, input
 *  :param fout: output stream
 *  :type fout: FILE, input
 *  :param timing: gives timing details, provide allocated zero initialized memory of length :c:func:`trlib_tri_timing_size`
 *
 *      ====== ================================
 *      block   description
 *      ====== ================================
 *      0       total duration
 *      1       timing of linear algebra calls
 *      2       :c:data:`timing` from :c:func:`trlib_leftmost_irreducible`
 *      3       :c:data:`timing` from :c:func:`trlib_eigen_inverse`
 *      ====== ================================
 *
 *  :type timing: trlib_int_t, input/output
 *  :param obj: objective function value at solution point
 *  :type obj: trlib_flt_t, output
 *  :param iter_newton: number of Newton iterations
 *  :type iter_newton: trlib_int_t, output
 *  :param sub_fail: status code of subroutine if failure occured in subroutines called
 *  :type sub_fail: trlib_int_t, output
 *
 *  :returns: status
 *
 *      - :c:macro:`TRLIB_TTR_CONV_BOUND`    success with solution on boundary
 *      - :c:macro:`TRLIB_TTR_CONV_INTERIOR` success with interior solution
 *      - :c:macro:`TRLIB_TTR_HARD`          success, but hard case encountered and solution may be approximate  
 *      - :c:macro:`TRLIB_TTR_NEWTON_BREAK`  most likely success with accurate result; premature end of Newton iteration due to tiny step
 *      - :c:macro:`TRLIB_TTR_HARD_INIT_LAM` hard case encountered without being able to find suitable initial :math:`\lambda` for Newton iteration, returned approximate stationary point that maybe suboptimal
 *      - :c:macro:`TRLIB_TTR_ITMAX`         iteration limit exceeded
 *      - :c:macro:`TRLIB_TTR_FAIL_FACTOR`   failure on matrix factorization
 *      - :c:macro:`TRLIB_TTR_FAIL_LINSOLVE` failure on backsolve
 *      - :c:macro:`TRLIB_TTR_FAIL_EIG`      failure on eigenvalue computation. status code of :c:func:`trlib_eigen_inverse` in :c:data:`sub_fail`
 *      - :c:macro:`TRLIB_TTR_FAIL_LM`       failure on leftmost eigenvalue computation. status code of :c:func:`trlib_leftmost_irreducible` in :c:data:`sub_fail`
 *  :rtype: trlib_int_t
 */

trlib_int_t trlib_tri_factor_min(
    trlib_int_t nirblk, trlib_int_t *irblk, trlib_flt_t *diag, trlib_flt_t *offdiag,
    trlib_flt_t *neglin, trlib_flt_t radius, 
    trlib_int_t itmax, trlib_flt_t tol_rel, trlib_flt_t tol_newton_tiny,
    trlib_int_t pos_def, trlib_int_t equality,
    trlib_int_t *warm0, trlib_flt_t *lam0, trlib_int_t *warm, trlib_flt_t *lam,
    trlib_int_t *warm_leftmost, trlib_int_t *ileftmost, trlib_flt_t *leftmost,
    trlib_int_t *warm_fac0, trlib_flt_t *diag_fac0, trlib_flt_t *offdiag_fac0,
    trlib_int_t *warm_fac, trlib_flt_t *diag_fac, trlib_flt_t *offdiag_fac,
    trlib_flt_t *sol0, trlib_flt_t *sol, trlib_flt_t *ones, trlib_flt_t *fwork,
    trlib_int_t refine,
    trlib_int_t verbose, trlib_int_t unicode, char *prefix, FILE *fout,
    trlib_int_t *timing, trlib_flt_t *obj, trlib_int_t *iter_newton, trlib_int_t *sub_fail);

/** Computes minimizer of regularized unconstrained problem
 *
 *  Computes minimizer of
 * 
 *  :math:`\min \frac 12 \langle h, (T + \lambda I) h \rangle + \langle g, h \rangle`,
 *  where :math:`T \in \mathbb R^{n \times n}` is symmetric tridiagonal and :math:`\lambda` such that :math:`T + \lambda I` is spd.
 *
 *  :param n: dimension, ensure :math:`n > 0`
 *  :type n: trlib_int_t, input
 *  :param diag: pointer to array holding diagonal of :math:`T`, length :c:data:`n`
 *  :type diag: trlib_flt_t, input
 *  :param offdiag: pointer to array holding offdiagonal of :math:`T`, length :c:data:`n-1`
 *  :type offdiag: trlib_flt_t, input
 *  :param neglin: pointer to array holding :math:`-g`, length :c:data:`n`
 *  :type neglin: trlib_flt_t, input
 *  :param lam: regularization parameter :math:`\lambda`
 *  :type lam: trlib_flt_t, input
 *  :param sol: pointer to array holding solution, length :c:data:`n`
 *  :type sol: trlib_flt_t, input/output
 *  :param ones: array with every value ``1.0``, length :c:data:`n`
 *  :type ones: trlib_flt_t, input
 *  :param fwork: floating point workspace, must be allocated memory on input of size :c:func:`trlib_tri_factor_memory_size` (call with argument :c:data:`n`) and can be discarded on output, memory layout:
 *
 *      ============== ============== ==========================================================
 *      start          end (excl)     description
 *      ============== ============== ==========================================================
 *                 0     :c:data:`n`  holds diagonal of :math:`T + \lambda I`
 *        :c:data:`n`  2 :c:data:`n`  holds diagonal of factorization :math:`T + \lambda I`
 *      2 :c:data:`n`  3 :c:data:`n`  holds offdiagonal of factorization :math:`T + \lambda I`
 *      3 :c:data:`n`  5 :c:data:`n`  workspace for iterative refinement
 *      ============== ============== ==========================================================
 *
 *  :type fwork: trlib_flt_t, input/output
 *  :param refine: set to ``1`` if iterative refinement should be used on solving linear systems, otherwise to ``0``
 *  :type refine: trlib_int_t, input
 *  :param verbose: determines the verbosity level of output that is written to :c:data:`fout`
 *  :type verbose: trlib_int_t, input
 *  :param unicode: set to ``1`` if :c:data:`fout` can handle unicode, otherwise to ``0``
 *  :type unicode: trlib_int_t, input
 *  :param prefix: string that is printed before iteration output
 *  :type prefix: char, input
 *  :param fout: output stream
 *  :type fout: FILE, input
 *  :param timing: gives timing details, provide allocated zero initialized memory of length :c:func:`trlib_tri_timing_size`
 *
 *      ====== ================================
 *      block   description
 *      ====== ================================
 *      0       total duration
 *      1       timing of linear algebra calls
 *      ====== ================================
 *
 *  :type timing: trlib_int_t, input/output
 *  :param norm_sol: norm of solution fector
 *  :type norm_sol: trlib_flt_t, output
 *  :param sub_fail: status code of subroutine if failure occured in subroutines called
 *  :type sub_fail: trlib_int_t, output
 *
 *  :returns: status
 *
 *      - :c:macro:`TRLIB_TTR_CONV_INTERIOR` success with interior solution
 *      - :c:macro:`TRLIB_TTR_FAIL_FACTOR`   failure on matrix factorization
 *      - :c:macro:`TRLIB_TTR_FAIL_LINSOLVE` failure on backsolve
 *
 *  :rtype: trlib_int_t
 */

trlib_int_t trlib_tri_factor_regularized_umin(
    trlib_int_t n, trlib_flt_t *diag, trlib_flt_t *offdiag,
    trlib_flt_t *neglin, trlib_flt_t lam,
    trlib_flt_t *sol,
    trlib_flt_t *ones, trlib_flt_t *fwork,
    trlib_int_t refine,
    trlib_int_t verbose, trlib_int_t unicode, char *prefix, FILE *fout,
    trlib_int_t *timing, trlib_flt_t *norm_sol, trlib_int_t *sub_fail);

/** Computes regularization parameter needed in trace
 *
 *  Let :math:`s(\lambda)` be solution of :math:`(T + \lambda I) s(\lambda) + g  = 0`,
 *  where :math:`T \in \mathbb R^{n \times n}` is symmetric tridiagonal and :math:`\lambda` such that :math:`T + \lambda I` is spd.
 *
 *  Then find :math:`\lambda` with :math:`\sigma_{\text l} \le \frac{\lambda}{s(\lambda)} \le \sigma_{\text u}`.
 *
 *  :param n: dimension, ensure :math:`n > 0`
 *  :type n: trlib_int_t, input
 *  :param diag: pointer to array holding diagonal of :math:`T`, length :c:data:`n`
 *  :type diag: trlib_flt_t, input
 *  :param offdiag: pointer to array holding offdiagonal of :math:`T`, length :c:data:`n-1`
 *  :type offdiag: trlib_flt_t, input
 *  :param neglin: pointer to array holding :math:`-g`, length :c:data:`n`
 *  :type neglin: trlib_flt_t, input
 *  :param lam: regularization parameter :math:`\lambda`, on input initial guess, on exit desired parameter
 *  :type lam: trlib_flt_t, input/output
 *  :param sigma: value that is used in root finding, e.g. :math:`\frac 12 (\sigma_{\text l} + \sigma_{\text u})`
 *  :type sigma: trlib_flt_t, input
 *  :param sigma_l: lower bound
 *  :type sigma_l: trlib_flt_t, input
 *  :param sigma_u: upper bound
 *  :type sigma_u: trlib_flt_t, input
 *  :param sol: pointer to array holding solution, length :c:data:`n`
 *  :type sol: trlib_flt_t, input/output
 *  :param ones: array with every value ``1.0``, length :c:data:`n`
 *  :type ones: trlib_flt_t, input
 *  :param fwork: floating point workspace, must be allocated memory on input of size :c:func:`trlib_tri_factor_memory_size` (call with argument :c:data:`n`) and can be discarded on output, memory layout:
 *
 *      ============== ============== ==========================================================
 *      start          end (excl)     description
 *      ============== ============== ==========================================================
 *                 0     :c:data:`n`  holds diagonal of :math:`T + \lambda I`
 *        :c:data:`n`  2 :c:data:`n`  holds diagonal of factorization :math:`T + \lambda I`
 *      2 :c:data:`n`  3 :c:data:`n`  holds offdiagonal of factorization :math:`T + \lambda I`
 *      3 :c:data:`n`  5 :c:data:`n`  workspace for iterative refinement
 *      5 :c:data:`n`  6 :c:data:`n`  auxiliary vector
 *      ============== ============== ==========================================================
 *
 *  :type fwork: trlib_flt_t, input/output
 *  :param refine: set to ``1`` if iterative refinement should be used on solving linear systems, otherwise to ``0``
 *  :type refine: trlib_int_t, input
 *  :param verbose: determines the verbosity level of output that is written to :c:data:`fout`
 *  :type verbose: trlib_int_t, input
 *  :param unicode: set to ``1`` if :c:data:`fout` can handle unicode, otherwise to ``0``
 *  :type unicode: trlib_int_t, input
 *  :param prefix: string that is printed before iteration output
 *  :type prefix: char, input
 *  :param fout: output stream
 *  :type fout: FILE, input
 *  :param timing: gives timing details, provide allocated zero initialized memory of length :c:func:`trlib_tri_timing_size`
 *
 *      ====== ================================
 *      block   description
 *      ====== ================================
 *      0       total duration
 *      1       timing of linear algebra calls
 *      ====== ================================
 *
 *  :type timing: trlib_int_t, input/output
 *  :param norm_sol: norm of solution fector
 *  :type norm_sol: trlib_flt_t, output
 *  :param sub_fail: status code of subroutine if failure occured in subroutines called
 *  :type sub_fail: trlib_int_t, output
 *
 *  :returns: status
 *
 *      - :c:macro:`TRLIB_TTR_CONV_INTERIOR` success with interior solution
 *      - :c:macro:`TRLIB_TTR_FAIL_FACTOR`   failure on matrix factorization
 *      - :c:macro:`TRLIB_TTR_FAIL_LINSOLVE` failure on backsolve
 *
 *  :rtype: trlib_int_t
 */

trlib_int_t trlib_tri_factor_get_regularization(
    trlib_int_t n, trlib_flt_t *diag, trlib_flt_t *offdiag,
    trlib_flt_t *neglin, trlib_flt_t *lam,
    trlib_flt_t sigma, trlib_flt_t sigma_l, trlib_flt_t sigma_u,
    trlib_flt_t *sol,
    trlib_flt_t *ones, trlib_flt_t *fwork,
    trlib_int_t refine,
    trlib_int_t verbose, trlib_int_t unicode, char *prefix, FILE *fout,
    trlib_int_t *timing, trlib_flt_t *norm_sol, trlib_int_t *sub_fail);

/** Compute diagonal regularization to make tridiagonal matrix positive definite
 *
 *  :param n: dimension, ensure :math:`n > 0`
 *  :type n: trlib_int_t, input
 *  :param diag: pointer to array holding diagonal of :math:`T`, length :c:data:`n`
 *  :type diag: trlib_flt_t, input
 *  :param offdiag: pointer to array holding offdiagonal of :math:`T`, length :c:data:`n-1`
 *  :type offdiag: trlib_flt_t, input
 *  :param tol_away: tolerance that diagonal entries in factorization should be away from zero, relative to previous entry. Good default :math:`10^{-12}`.
 *  :type tol_away: trlib_flt_t, input
 *  :param security_step: factor greater ``1.0`` that defines a margin to get away from zero in the step taken. Good default ``2.0``.
 *  :type security_step: trlib_flt_t, input
 *  :param regdiag: pointer to array holding regularization term, length :c:data:`n`
 *  :type regdiag: trlib_flt_t, input/output
 *
 *  :returns: ``0``
 *  :rtype: trlib_int_t
 */
trlib_int_t trlib_tri_factor_regularize_posdef(
    trlib_int_t n, trlib_flt_t *diag, trlib_flt_t *offdiag,
    trlib_flt_t tol_away, trlib_flt_t security_step, trlib_flt_t *regdiag);

/** Gives information on memory that has to be allocated for :c:func:`trlib_tri_factor_min`
 *  
 *  :param n: dimension, ensure :math:`n > 0`
 *  :type n: trlib_int_t, input
 *  :param fwork_size: size of floating point workspace fwork that has to be allocated for :c:func:`trlib_tri_factor_min`
 *  :type fwork_size: trlib_flt_t, output
 *  
 *  :returns: ``0``
 *  :rtype: trlib_int_t
 */

trlib_int_t trlib_tri_factor_memory_size(trlib_int_t n);

/** size that has to be allocated for :c:data:`timing` in :c:func:`trlib_tri_factor_min`
 */
trlib_int_t trlib_tri_timing_size(void);

#endif
