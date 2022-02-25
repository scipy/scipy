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

#ifndef TRLIB_EIGEN_INVERSE_H
#define TRLIB_EIGEN_INVERSE_H

#define TRLIB_EIR_CONV          (0)
#define TRLIB_EIR_ITMAX         (-1)
#define TRLIB_EIR_FAIL_FACTOR   (-2)
#define TRLIB_EIR_FAIL_LINSOLVE (-3)

#define TRLIB_EIR_N_STARTVEC    (5)

/** Computes eigenvector to provided eigenvalue of symmetric tridiagonal matrix :math:`T \in \mathbb R^{n\times n}`,
 *  using inverse iteration.
 *
 *  For a description of the method, see https://en.wikipedia.org/wiki/Inverse_iteration.
 *
 *  **Convergence**
 *
 *  Convergence is reported if :math:`\vert \frac{1}{\Vert w_{i+1} \Vert} - \texttt{pert} \vert \le \texttt{tol}\_\texttt{abs}`, where :math:`(T-\lambda I) w_{i+1} = v_i`, :math:`v_i` the current normalized iterate and :math:`\texttt{pert}` is the perturbation applied to the provided eigenvalue.
 *
 *  :param n: dimension, ensure :math:`n > 0`
 *  :type n: trlib_int_t, input
 *  :param diag: pointer to array holding diagonal of :math:`T`, length :c:data:`n`
 *  :type diag: trlib_flt_t, input
 *  :param offdiag: pointer to array holding offdiagonal of :math:`T`, length :c:data:`n-1`
 *  :type offdiag: trlib_flt_t, input
 *  :param lam_init: estimation of eigenvalue corresponding to eigenvector to compute
 *  :type lam_init: trlib_flt_t, input
 *  :param itmax: maximum number of iterations
 *  :type itmax: trlib_int_t, input
 *  :param tol_abs: absolute stopping tolerance in inverse iteration, good default may be :math:`\sqrt{\texttt{macheps}}` (:c:macro:`TRLIB_EPS_POW_5`)
 *  :type tol_abs: trlib_flt_t, input
 *  :param ones: array with every value ``1.0``, length :c:data:`n`
 *  :type ones: trlib_flt_t, input
 *  :param diag_fac: pointer to array holding diagonal of Cholesky factorization of :math:`T - \lambda I`, length :c:data:`n`
 *
 *      - on entry: allocated memory
 *      - on exit: factorization corresponding to computed eigenvalue @p lam
 *
 *  :type diag_fac: trlib_flt_t, input/output
 *  :param offdiag_fac: pointer to array holding offdiagonal of Cholesky factorization of :math:`T - \lambda I`, length :c:data:`n-1`
 *
 *      - on entry: allocated memory
 *      - on exit: factorization corresponding to computed eigenvalue @p lam
 *
 *  :type offidag_fac: trlib_flt_t, input/output
 *  :param eig: pointer to array holding eigenvector, length :c:data:`n`
 *  :type eig: trlib_flt_t, input/output
 *  :param verbose: determines the verbosity level of output that is written to :c:data:`fout`
 *  :type verbose: trlib_int_t, input
 *  :param unicode: set to ``1`` if :c:data:`fout` can handle unicode, otherwise to ``0``
 *  :type unicode: trlib_int_t, input
 *  :param prefix: string that is printed before iteration output
 *  :type prefix: char, input
 *  :param fout: output stream
 *  :type fout: FILE, input
 *  :param timing: gives timing details, provide allocated zero initialized memory of length :c:func:`trlib_eigen_timing_size`
 *
 *      ====== ================================
 *      block   description
 *      ====== ================================
 *      0       total duration
 *      1       timing of linear algebra calls
 *      ====== ================================
 *
 *  :type timing: trlib_int_t, input/output
 *  :param lam_pert: eigenvalue corresponding to eigenvector
 *  :type lam_pert: trlib_flt_t, output
 *  :param pert: perturbation applied to provided eigenvalue
 *  :type pert: trlib_flt_t, output
 *  :param iter_inv: number of inverse iterations
 *  :type iter_inv: trlib_int_t, output
 *
 *  :returns:      status
 *
 *                 - :c:macro:`TRLIB_EIR_CONV`          success
 *                 - :c:macro:`TRLIB_EIR_ITMAX`         iteration limit exceeded
 *                 - :c:macro:`TRLIB_EIR_FAIL_FACTOR`   failure on matrix factorization
 *                 - :c:macro:`TRLIB_EIR_FAIL_LINSOLVE` failure on backsolve
 *
 *  :rtype: trlib_int_t
 */

trlib_int_t trlib_eigen_inverse(
        trlib_int_t n, trlib_flt_t *diag, trlib_flt_t *offdiag,
        trlib_flt_t lam_init, trlib_int_t itmax, trlib_flt_t tol_abs,
        trlib_flt_t *ones, trlib_flt_t *diag_fac, trlib_flt_t *offdiag_fac,
        trlib_flt_t *eig, trlib_int_t verbose, trlib_int_t unicode, char *prefix, FILE *fout,
        trlib_int_t *timing, trlib_flt_t *lam_pert, trlib_flt_t *pert, trlib_int_t *iter_inv);

/** size that has to be allocated for :c:data:`timing` in :c:func:`trlib_eigen_inverse`
 */
trlib_int_t trlib_eigen_timing_size(void);

#endif
