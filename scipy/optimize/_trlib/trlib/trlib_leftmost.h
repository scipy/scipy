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

#ifndef TRLIB_LEFTMOST_H
#define TRLIB_LEFTMOST_H

#define TRLIB_LMR_CONV          (0)
#define TRLIB_LMR_ITMAX         (-1)
#define TRLIB_LMR_NEWTON_BREAK  (3)

/** Computes smallest eigenvalue of symmetric tridiagonal matrix
 *  :math:`T \in \mathbb R^{n\times n}`,
 *  using a iteration based on last-pivot function of Parlett and Reid.
 *
 *  Let :math:`T = \begin{pmatrix} T_1 & & \\ & \ddots & \\ & & T_\ell \end{pmatrix}`
 *  be composed into irreducible blocks :math:`T_i`.
 *  
 *  Calls :c:func:`trlib_leftmost_irreducible` on every irreducible block in case of coldstart,
 *  in case of warmstart just updates information on :math:`T_\ell`.
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
 *  :param warm: set :math:`\ge 1` if you want to update information about block :math:`\ell`, provide values in :c:data:`leftmost_minor`, :c:data:`ileftmost`, :c:data:`leftmost`; else ``0``
 *  :type warm: trlib_int_t, input
 *  :param leftmost_minor: smallest eigenvalue of principal :math:`(n_\ell-1)\times (n_\ell-1)` submatrix of :math:`T_\ell`
 *  :type leftmost_minor: trlib_flt_t, input
 *  :param itmax: maximum number of iterations
 *  :type itmax: trlib_int_t, input
 *  :param tol_abs: absolute stopping tolerance in Reid-Parlett zero finding, good default may be :math:`\sqrt{\texttt{macheps}}^{3/4}` (:c:macro:`TRLIB_EPS_POW_75`)
 *  :type tol_abs: trlib_flt_t, input
 *  :param verbose: determines the verbosity level of output that is written to :c:data:`fout`
 *  :type verbose: trlib_int_t, input
 *  :param unicode: set to ``1`` if :c:data:`fout` can handle unicode, otherwise to ``0``
 *  :type unicode: trlib_int_t, input
 *  :param prefix: string that is printed before iteration output
 *  :type prefix: char, input
 *  :param fout: output stream
 *  :type fout: FILE, input
 *  :param timing: gives timing details, provide allocated zero initialized memory of length :c:func:`trlib_leftmost_timing_size`
 *
 *      ====== ================================
 *      block   description
 *      ====== ================================
 *      0       total duration
 *      ====== ================================
 *
 *  :type timing: trlib_int_t, input/output
 *  :param ileftmost: index of block that corresponds to absolute smallest eigenvalue
 *  :type ileftmost: trlib_int_t, input/output
 *  :param leftmost: smallest eigenvalue of :math:`T`, length :math:`\ell`
 *
 *      - on entry: allocated memory
 *      - on exit: :c:data:`leftmost[i]` smallest eigenvalue of :math:`T_i`
 *
 *  :type leftmost: trlib_flt_t, input/output
 *  :returns:      status
 *
 *                 - :c:macro:`TRLIB_LMR_CONV`          success
 *                 - :c:macro:`TRLIB_LMR_ITMAX`         iteration limit exceeded
 *
 *  :rtype: trlib_int_t
 */

trlib_int_t trlib_leftmost(
        trlib_int_t nirblk, trlib_int_t *irblk, trlib_flt_t *diag, trlib_flt_t *offdiag,
        trlib_int_t warm, trlib_flt_t leftmost_minor, trlib_int_t itmax, trlib_flt_t tol_abs,
        trlib_int_t verbose, trlib_int_t unicode, char *prefix, FILE *fout,
        trlib_int_t *timing, trlib_int_t *ileftmost, trlib_flt_t *leftmost);

/** Computes smallest eigenvalue of irreducible symmetric tridiagonal matrix
 *  :math:`T \in \mathbb R^{n\times n}`,
 *  using a iteration based on last-pivot function of Parlett and Reid.
 *  
 *  Method is sketched on p. 516 in [Gould1999]_.
 *
 *  Note that this function most likely will fail in the case of a reducible matrix
 *  (:c:data:`offdiag` contains 0).
 *
 *  **Convergence**
 *  
 *  Convergence is reported if :math:`\texttt{up}-\texttt{low} \le \texttt{tol}\_\texttt{abs} * \max\{1, \vert \texttt{low} \vert, \vert \texttt{up} \vert \}` or :math:`\texttt{prlp} \le \texttt{tol}\_\texttt{abs}`, :math:`\texttt{low}` and :math:`\texttt{up}` denote bracket values enclosing the leftmost eigenvalue and :math:`\texttt{prlp}` denotes the last-pivot function value used in root finding.
 *
 *  :param n: dimension, ensure :math:`n > 0`
 *  :type n: trlib_int_t, input
 *  :param diag: pointer to array holding diagonal of :math:`T`, length :c:data:`n`
 *  :type diag: trlib_flt_t, input
 *  :param offdiag: pointer to array holding offdiagonal of :math:`T`, length :c:data:`n-1`
 *  :type offdiag: trlib_flt_t, input
 *  :param  warm: set :math:`\ge 1` if you provide a valid value in :c:data:`leftmost_minor`, else ``0``. Exact value determines which model will be used for zero-finding
 *
 *      ===== ======
 *      warm  model
 *      ===== ======
 *      0     linear model for rational function
 *      1     sensible heuristic model choice for lifted rational function
 *      2     asymptotic quadratic model :math:`\theta^2 + b \theta + c` for lifted rational function
 *      3     taylor quadratic model :math:`a \theta^2 + b \theta + c` for lifted rational function
 *      4     linear model :math:`b \theta + c` for lifted rational function
 *      ===== ======
 *
 *  :type warm: trlib_int_t, input
 *  :param leftmost_minor: smallest eigenvalue of principal :math:`(n-1)\times (n-1)` submatrix of :math:`T`
 *  :type leftmost_minor: trlib_flt_t, input
 *  :param itmax: maximum number of iterations
 *  :type itmax: trlib_int_t, input
 *  :param tol_abs: absolute stopping tolerance in Reid-Parlett zero finding, good default may be :math:`\sqrt{\texttt{macheps}}^{3/4}` (:c:macro:`TRLIB_EPS_POW_75`)
 *  :type tol_abs: trlib_flt_t, input
 *  :param verbose: determines the verbosity level of output that is written to :c:data:`fout`
 *  :type verbose: trlib_int_t, input
 *  :param unicode: set to ``1`` if :c:data:`fout` can handle unicode, otherwise to ``0``
 *  :type unicode: trlib_int_t, input
 *  :param prefix: string that is printed before iteration output
 *  :type prefix: char, input
 *  :param fout: output stream
 *  :type fout: FILE, input
 *  :param timing: gives timing details, provide allocated zero initialized memory of length :c:func:`trlib_leftmost_timing_size`
 *
 *      ====== ================================
 *      block   description
 *      ====== ================================
 *      0       total duration
 *      ====== ================================
 *
 *  :type timing: trlib_int_t, input/output
 *  :param leftmost: smallest eigenvalue of :math:`T`
 *  :type leftmost: trlib_flt_t, output
 *  :param iter_pr: number of Parlett-Reid iterations
 *  :type iter_pr: trlib_int_t, output
 *
 *  :returns:      status
 *
 *                 - :c:macro:`TRLIB_LMR_CONV`          success
 *                 - :c:macro:`TRLIB_LMR_ITMAX`         iteration limit exceeded
 *
 *  :rtype: trlib_int_t
 */

trlib_int_t trlib_leftmost_irreducible(
        trlib_int_t n, trlib_flt_t *diag, trlib_flt_t *offdiag,
        trlib_int_t warm, trlib_flt_t leftmost_minor, trlib_int_t itmax, trlib_flt_t tol_abs,
        trlib_int_t verbose, trlib_int_t unicode, char *prefix, FILE *fout,
        trlib_int_t *timing, trlib_flt_t *leftmost, trlib_int_t *iter_pr);

/** size that has to be allocated for :c:data:`timing` in :c:func:`trlib_leftmost_irreducible` and :c:func:`trlib_leftmost`
 */
trlib_int_t trlib_leftmost_timing_size(void);

#endif
