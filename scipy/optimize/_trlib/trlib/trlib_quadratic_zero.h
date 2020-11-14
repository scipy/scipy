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

#ifndef TRLIB_QUADRATIC_ZERO_H
#define TRLIB_QUADRATIC_ZERO_H

/** Computes real zeros of normalized quadratic polynomial.
 *
 *  :param c_abs: absolute coefficient
 *  :type c_abs: trlib_flt_t, input
 *  :param c_lin: coefficinet of linear term
 *  :type c_lin: trlib_flt_t, input
 *  :param tol: tolerance that indicates if ill-conditioning present, good default may be :math:`\texttt{macheps}^{3/4}` (:c:macro:`TRLIB_EPS_POW_75`)
 *  :type tol: trlib_flt_t, input
 *  :param verbose: determines the verbosity level of output that is written to :c:data:`fout`
 *  :type verbose: trlib_int_t, input
 *  :param unicode: set to ``1`` if :c:data:`fout` can handle unicode, otherwise to ``0``
 *  :type unicode: trlib_int_t, input
 *  :param prefix: string that is printed before iteration output
 *  :type prefix: trlib_int_t, input
 *  :param fout: output stream
 *  :type fout: FILE, input
 *  :param t1: first zero, :math:`\texttt{t1} \le \texttt{t2}`
 *  :type t2: trlib_flt_t, output
 *  :param t2: second zero, :math:`\texttt{t1} \le \texttt{t2}`
 *  :type t2: trlib_flt_t, output
 *
 *  :returns: number of zeros
 *  :rtype: trlib_int_t
 */

trlib_int_t trlib_quadratic_zero(trlib_flt_t c_abs, trlib_flt_t c_lin, trlib_flt_t tol,
        trlib_int_t verbose, trlib_int_t unicode, char *prefix, FILE *fout,
        trlib_flt_t *t1, trlib_flt_t *t2);

#endif
