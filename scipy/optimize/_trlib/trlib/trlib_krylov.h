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

#ifndef TRLIB_KRYLOV_H
#define TRLIB_KRYLOV_H

#define TRLIB_CLR_CONV_BOUND    (0)
#define TRLIB_CLR_CONV_INTERIOR (1)
#define TRLIB_CLR_APPROX_HARD   (2)
#define TRLIB_CLR_NEWTON_BREAK  (3)
#define TRLIB_CLR_HARD_INIT_LAM (4)
#define TRLIB_CLR_FAIL_HARD     (5)
#define TRLIB_CLR_UNBDBEL       (6)
#define TRLIB_CLR_UNLIKE_CONV   (7)
#define TRLIB_CLR_CONTINUE      (10)
#define TRLIB_CLR_ITMAX         (-1)
#define TRLIB_CLR_FAIL_FACTOR   (-3)
#define TRLIB_CLR_FAIL_LINSOLVE (-4)
#define TRLIB_CLR_FAIL_NUMERIC  (-5)
#define TRLIB_CLR_FAIL_TTR      (-7)
#define TRLIB_CLR_PCINDEF       (-8)
#define TRLIB_CLR_UNEXPECT_INT  (-9)

#define TRLIB_CLT_CG            (1)
#define TRLIB_CLT_L             (2)

#define TRLIB_CLA_TRIVIAL       (0)
#define TRLIB_CLA_INIT          (1)
#define TRLIB_CLA_RETRANSF      (2)
#define TRLIB_CLA_UPDATE_STATIO (3)
#define TRLIB_CLA_UPDATE_GRAD   (4)
#define TRLIB_CLA_UPDATE_DIR    (5)
#define TRLIB_CLA_NEW_KRYLOV    (6)
#define TRLIB_CLA_CONV_HARD     (7)
#define TRLIB_CLA_OBJVAL        (8)

#define TRLIB_CLS_INIT          (1)
#define TRLIB_CLS_HOTSTART      (2)
#define TRLIB_CLS_HOTSTART_G    (3)
#define TRLIB_CLS_HOTSTART_P    (4)
#define TRLIB_CLS_HOTSTART_R    (5)
#define TRLIB_CLS_HOTSTART_T    (6)
#define TRLIB_CLS_VEC_INIT      (7)
#define TRLIB_CLS_CG_NEW_ITER   (8)
#define TRLIB_CLS_CG_UPDATE_S   (9)
#define TRLIB_CLS_CG_UPDATE_GV  (10)
#define TRLIB_CLS_CG_UPDATE_P   (11)
#define TRLIB_CLS_LANCZOS_SWT   (12)
#define TRLIB_CLS_L_UPDATE_P    (13)
#define TRLIB_CLS_L_CMP_CONV    (14)
#define TRLIB_CLS_L_CMP_CONV_RT (15)
#define TRLIB_CLS_L_CHK_CONV    (16)
#define TRLIB_CLS_L_NEW_ITER    (17)
#define TRLIB_CLS_CG_IF_IRBLK_P (18)
#define TRLIB_CLS_CG_IF_IRBLK_C (19)
#define TRLIB_CLS_CG_IF_IRBLK_N (20)

#define TRLIB_CLC_NO_EXP_INV    (0)
#define TRLIB_CLC_EXP_INV_LOC   (1)
#define TRLIB_CLC_EXP_INV_GLO   (2)

#define TRLIB_CLT_CG_INT        (0)
#define TRLIB_CLT_CG_BOUND      (1)
#define TRLIB_CLT_LANCZOS       (2)
#define TRLIB_CLT_HOTSTART      (3)

/** 
 *  Solves trust region subproblem: Computes minimizer to
 * 
 *  :math:`\min \frac 12 \langle s, H s \rangle + \langle g_0, s \rangle`
 *  subject to the trust region constraint :math:`\Vert s \Vert_M \le r`,
 *  where 
 *
 *  - :math:`H` is available via matrix-vector products :math:`v \mapsto Hv`,
 *  - :math:`\Vert s \Vert_M = \sqrt{\langle s, Ms \rangle}` with :math:`M` positive definite,
 *  - :math:`M` is available via matrix-vector products of its inverse, :math:`v \mapsto M^{-1}v`.
 *
 *  The minimizer is a global minimizer (modulo floating point),
 *  as long as the hard case does not occur.
 *  The hard case is characterized by the fact that the eigenspace
 *  corresponding to the smallest eigenvalue is degenerate.
 *
 *  In that case a global minimizer in the Krylov subspaces sampled so far is returned,
 *  sampling the complete search space can be forced by setting :c:data:`ctl_invariant` to :c:macro:`TRLIB_CLC_EXP_INV_GLO`,
 *  but if this is desired factorization based method will most likely be much more efficient.
 *
 *  A preconditioned Conjugate Gradient / Lanczos method is used,
 *  that possibly employs a reduction to a tridiagnoal subproblem
 *  that is solved using Moré-Sorensens method by explicitly
 *  factorizing the tridiagonal matrix, calling :c:func:`trlib_tri_factor_min`.
 *
 *  The method builds upon algorithm 5.1 in [Gould1999]_, details of the implementation can be found [Lenders2016]_:
 *
 *  **Convergence**
 *
 *  The stopping criterion is based on the gradient of the Lagrangian. Convergence in iteration :math:`i` is reported as soon as
 *
 *  - interior case: :math:`\Vert g_{i+1} \Vert_{M^{-1}} \le \max\{ \texttt{tol}\_\texttt{abs}\_\texttt{i}, \eta_i \, \Vert g_0 \Vert_{M^{-1}} \}`.
 *  - boundary case: :math:`\Vert g_{i+1} \Vert_{M^{-1}} \le \max\{ \texttt{tol}\_\texttt{abs}\_\texttt{b}, \eta_b \, \Vert g_0 \Vert_{M^{-1}} \}`.
 *  - hard case: :c:data:`ctl_invariant` determines if this is the sole stopping criterion or if further invariant subspaces are going to be explored. See the documentation of this option.
 *
 *  Here 
 *  :math:`\eta_i = \begin{cases} \texttt{tol}\_\texttt{rel}\_\texttt{i} & \texttt{tol}\_\texttt{rel}\_\texttt{i} > 0 \\ \min\{ 0.5, \sqrt{\Vert g_{i+1} \Vert_{M^{-1}}} \} & \texttt{tol}\_\texttt{rel}\_\texttt{i} = -1 \\ \min\{ 0.5, \Vert g_{i+1} \Vert_{M^{-1}} \} & \texttt{tol}\_\texttt{rel}\_\texttt{i} = -2 \end{cases},`
 *  :math:`\eta_b = \begin{cases} \texttt{tol}\_\texttt{rel}\_\texttt{b} & \texttt{tol}\_\texttt{rel}\_\texttt{b} > 0 \\ \min\{ 0.5, \sqrt{\Vert g_{i+1} \Vert_{M^{-1}}} \} & \texttt{tol}\_\texttt{rel}\_\texttt{b} = -1 \\ \min\{ 0.5, \Vert g_{i+1} \Vert_{M^{-1}} \} & \texttt{tol}\_\texttt{rel}\_\texttt{b} = -2, \\ \max\{10^{-6}, \min\{ 0.5, \sqrt{\Vert g_{i+1} \Vert_{M^{-1}}} \}\} & \texttt{tol}\_\texttt{rel}\_\texttt{b} = -3 \\ \max\{10^{-6}, \min\{ 0.5, \Vert g_{i+1} \Vert_{M^{-1}}\} \} & \texttt{tol}\_\texttt{rel}\_\texttt{b} = -4 \end{cases}.`
 *
 *  Choice of :math:`\eta_i` and :math:`\eta_b` depend on the application, the author has good overall experience in unconstrained optimization with :math:`\texttt{tol}\_\texttt{rel}\_\texttt{i} = -2` and :math:`\texttt{tol}\_\texttt{rel}\_\texttt{b} = -3`.
 *  Some remarks to keep in mind when choosing the tolerances:
 *
 *  - Choosing fixed small values for :math:`\eta_i` and :math:`\eta_b`, for example :math:`\texttt{tol}\_\texttt{rel}\_\texttt{i} = 10^{-8}` and :math:`\texttt{tol}\_\texttt{rel}\_\texttt{b} = 10^{-6}` lead to quick convergence in a NLP algorithm with lowest number of function evaluations but at a possible excessive amount of matrix-vector products as this also solves problems accurately far away from the solution.
 *  - Choosing :math:`\eta \sim O(\sqrt{\Vert g \Vert_{M^{-1}}})` leads to superlinear convergence in a nonlinear programming algorithm.
 *  - Choosing :math:`\eta \sim O(\Vert g \Vert_{M^{-1}})` leads to quadratic convergence in a nonlinear programming algorithm.
 *  - It is questionable whether it makes sense to solve the boundary with very high accuracy, as the final solution of a nonlinear program should be interior.
 *
 *
 *  **Calling Scheme**
 *
 *  This function employs a reverse communication paradigm.
 *  The functions exits whenever there is an action to be performed
 *  by the user as indicated by the :c:data:`action`.
 *  The user should perform this action and continue with the other
 *  values unchanged as long as the return value is positive.
 *
 *  **User provided storage**
 *
 *  The user has to manage 5/6 vectors referred by the names
 *  :math:`g`, :math:`g_-`, :math:`v`, :math:`s`, :math:`p`, :math:`\textit{Hp}`
 *  and a matrix :math:`Q` with :c:data:`itmax` columns to store the
 *  Lanczos directions :math:`p_i`.
 *  The user must perform the actions on those as indicated by the return value.
 *
 *  In the case of trivial preconditioner, :math:`M = I`, it always
 *  holds :math:`v = g` and the vector :math:`v` is not necessary.
 *
 *  :math:`s` holds the solution candidate.
 *
 *  Note that the action :math:`s \leftarrow Q h` will only sometimes be requested
 *  in the very final iteration before convergence. If memory and not computational
 *  load is an issue, you may very well save :c:data:`iter`, :c:data:`ityp`, :c:data:`flt1`, :c:data:`flt2`
 *  :c:data:`flt2` instead of :math:`q_j` and when :math:`s \leftarrow Q s` is requested
 *  simultaniously recompute the directions :math:`q_j` and update the direction
 *  :math:`s` once :math:`q_j` has been computed as :math:`s \leftarrow h_j q_j`
 *  with initialization :math:`s \leftarrow 0`.
 *
 *  Furthermore the user has to provide a integer and floating point workspace, details
 *  are described in :c:data:`iwork` and :c:data:`fwork`.
 *
 *  **Resolves**
 *
 *  *Reentry with new radius*
 *  
 *  You can efficiently hotstart from old results if you have a new problem
 *  with *decreased* trust region radius. Just set `status` to :c:macro:`TRLIB_CLS_HOTSTART`.
 *  Furthermore hotstarting with increased trust region radius should be
 *  trivial as you should be able to just increase the radius and take off
 *  the computation from the previous stage. However as this is an untypical
 *  scenario, it has not been tested at all.
 *
 *  *Reentry to get minimizer of regularized problem*
 *
 *  You can reenter to compute the solution of a convexified trust region problem.
 *  You may be interested to this once you see the effect of bad conditioning or also as a cheap retry in a trust region algorithm
 *  before continuing to the next point after a discarded step.
 *
 *  In this case, call the function with `status` set to :c:macro:`TRLIB_CLS_HOTSTART_P`.
 *
 *  *Reentry to get unconstrained minimizer of constant regularized problem*
 *
 *  After a successful solve you can compute the norm of the unconstrained minimizer to the problem with regularized
 *  hessian :math:`H + \beta M` (:math:`\beta` large enough such that :math:`H + \beta M` is positive definite)
 *  in the previously expanded Krylov subspace. You will certainly be interested in this if you want to implement
 *  the TRACE algorithm described in [Curtis2016]_.
 *
 *  In this case, call the function with `status` set to :c:macro:`TRLIB_CLS_HOTSTART_R` and ``radius`` set to :math:`\beta`.
 *
 *  On exit, the ``obj`` field of :c:data:`fwork` yields the norm of the solution vector.
 *
 *  If you not only want the norm but also the unconstrained minimizer, you can get it by one step of :c:func:`TRLIB_CLA_RETRANSF`.
 *  However since this is usually not the case, you are not asked to do this in the reverse communication scheme.
 *
 *  *Reentry to find suitable TRACE regularization parameter*
 *
 *  For the aforementioned TRACE algorithm, there is also the need to compute :math:`\beta` such that :math:`\sigma_{\text l} \le \frac{\beta}{\Vert s(\beta) \Vert} \le \sigma_{\text u}`, where :math:`\Vert s(\beta) \Vert` denotes the norm of the regularized unconstrained minimizer.
 *
 *  To get this, set :c:data:`status` to :c:macro:`TRLIB_CLS_HOTSTART_T` and :c:data:`radius` to an initial guess for :math:`\beta`, :c:data:`tol_rel_i` to :math:`\sigma_{\text l}` and :c:data:`tol_rel_b` to :math:`\sigma_{\text u}`. The field ``obj`` of :c:data:`fwork` contains the desired solution norm on exit and :c:data:`flt1` is the desired regularization parameter :math:`\beta`.
 *
 *  If you not only want the norm but also the unconstrained minimizer, you can get it by one step of :c:macro:`TRLIB_CLA_RETRANSF`.
 *  However since this is usually not the case, you are not asked to do this in the reverse communication scheme.
 *
 *  *Reentry with new gradient*
 *  
 *  You can efficiently hotstart from old results if you have a new problem
 *  with changed :math:`g_0` where the previous sampled Krylov subspace is
 *  large enough to contain the solution to the new problem, convergence will
 *  *not* be checked:
 *
 *  - get pointer ``gt`` to negative gradient of tridiagonal problem in `fwork` with :c:func:`trlib_krylov_gt`
 *  - store ``gt[0]`` and overwrite :math:`\texttt{gt} \leftarrow - Q_i^T \texttt{grad}`
 *  - set :c:data:`status` to :c:macro:`TRLIB_CLS_HOTSTART_G` and start reverse communication process
 *  - to reset data for new problem make sure that you restore ``gt[0]`` and set ``gt[1:] = 0`` for those elements previously overwritten
 *
 *  **Hard case**
 *
 *  If you want to investigate the problem also if the hard case, you can either sample through the invariant subspaces (:c:data:`ctl_invariant` set to :c:macro:`TRLIB_CLC_EXP_INV_GLO`) or solve the problem with a gradient for which the hard does not occur and then hotstart with the actual gradient.

 *
 *  :param init: set to :c:macro:`TRLIB_CLS_INIT` on first call, to :c:macro:`TRLIB_CLS_HOTSTART` on hotstart with smaller radius and otherwise to ``0``.
 *  :type init: trlib_int_t, input
 *  :param radius: trust region radius
 *  :type radius: trlib_flt_t, input
 *  :param equality: set to ``1`` if trust region constraint should be enforced as equality
 *  :type equality: trlib_int_t, input
 *  :param itmax: maximum number of iterations
 *  :type itmax: trlib_int_t, input
 *  :param itmax_lanczos: maximum number of Lanczos type iterations.
 *                 You are strongly advised to set this to a small number (say 25) unless you know better.
 *                 Keep in mind that Lanczos type iteration are only performed when curvature
 *                 information is flat and Lanczos may amplify rounding errors without
 *                 reorthogonalization. If you allow a lot of Lanczos type iterations
 *                 consider reorthogonalizing the new direction against the vector storage.
 *
 *  :type itmax_lanczos: trlib_int_t, input
 *  :param tol_rel_i: relative stopping tolerance for interior solution
 *  :type tol_rel_i: trlib_flt_t, input
 *  :param tol_abs_i: absolute stopping tolerance for interior solution
 *  :type tol_abs_i: trlib_flt_t, input
 *  :param tol_rel_b: relative stopping tolerance for boundary solution
 *  :type tol_rel_b: trlib_flt_t, input
 *  :param tol_abs_b: absolute stopping tolerance for boundary solution
 *  :type tol_abs_b: trlib_flt_t, input
 *  :param zero: threshold that determines when :math:`\langle p, Hp \rangle` is considered to be zero and thus control eventual Lanczos switch
 *  :type zero: trlib_flt_t, input
 *  :param obj_lo: lower bound on objective, returns if a point is found with function value :math:`\le \texttt{obj}\_\texttt{lo}`. Set to a positive value to ignore this bound.
 *  :type obj_lo: trlib_flt_t, input
 *  :param ctl_invariant: 
 *
 *                 - set to :c:macro:`TRLIB_CLC_NO_EXP_INV` if you want to search only in the first invariant Krylov subspace
 *                 - set to :c:macro:`TRLIB_CLC_EXP_INV_LOC` if you want to continue to expand the Krylov subspaces but terminate if there is convergence indication in the subspaces sampled so far.
 *                 - set to :c:macro:`TRLIB_CLC_EXP_INV_GLO` if you want to continue to expand the Krylov subspaces even in the case of convergence to a local minimizer within in the subspaces sampled so far. Reverse communication continues as long as ctl_invariant is set to :c:macro:`TRLIB_CLC_EXP_INV_GLO`, so you should reset :c:data:`ctl_invariant` to either :c:macro:`TRLIB_CLC_EXP_INV_LOC` or :c:macro:`TRLIB_CLC_EXP_INV_LOC` if there is no reason to continue, for example because you cannot find a new nonzero random vector orthogonal to the sampled directions if :c:data:`action` is :c:macro:`TRLIB_CLA_NEW_KRYLOV`.
 *
 *  :type ctl_invariant: trlib_int_t, input
 *  :param convexify: set to `1` if you like to monitor if the tridiagonal solution and the backtransformed solution match and if not resolve with a convexified problem, else set to `0`
 *  :type convexify: trlib_int_t, input
 *  :param earlyterm: set to `1` if you like to terminate in the boundary case if it unlikely that much progress will be made fast but no convergence is reached, else set to `0`
 *  :type earlyterm: trlib_int_t, input
 *  :param g_dot_g: dot product :math:`\langle g, g \rangle`
 *  :type g_dot_g: trlib_flt_t, input
 *  :param v_dot_g: dot product :math:`\langle v, g \rangle`
 *  :type v_dot_g: trlib_flt_t, input
 *  :param p_dot_Hp: dot product :math:`\langle p, Hp \rangle`
 *  :type p_dot_Hp: trlib_flt_t, input
 *  :param iwork: integer workspace for problem, internal memory layout described in table below
 *
 *                 - on first call, provide allocated memory for :c:data:`iwork_size` entries provided by :c:func:`trlib_krylov_memory_size`
 *                 - leave untouched between iterations
 *
 *                 ====== ================= ==========================
 *                 start   end(exclusive)   description
 *                 ====== ================= ==========================
 *                 0       1                internal status flag
 *                 1       2                iteration counter
 *                 2       3                flag that indicates if :math:`H` is positive definite in sampled Krylov subspace
 *                 3       4                flag that indicates if interior solution is expected 
 *                 4       5                flag that indicates if warmstart information to :c:data:`leftmost` is available
 *                 5       6                index to block with smallest leftmost
 *                 6       7                flag that indicates if warmstart information to :math:`\lambda_0` is available
 *                 7       8                flag that indicates if warmstart information to :math:`\lambda` is available
 *                 8       9                iteration in which switched to Lanczos iteration, ``-1``: no switch occurred
 *                 9       10               return code from :c:func:`trlib_tri_factor_min`
 *                 10      11               :c:data:`sub_fail` exit code from subrotines called in :c:func:`trlib_tri_factor_min`
 *                 11      12               number of newton iterations needed in :c:func:`trlib_tri_factor_min`
 *                 12      13               last iteration in which a headline has been printed
 *                 13      14               kind of iteration headline that has been printed last
 *                 14      15               status variable if convexified version is going to be resolved
 *                 15      16               number of irreducible blocks
 *                 16      17 + ``itmax``   decomposition into irredecubile block, :c:data:`irblk` for :c:func:`trlib_tri_factor_min`
 *                 ====== ================= ==========================
 *
 *  :type iwork: trlib_int_t, input/output
 *  :param fwork: floating point workspace for problem, internal memory layout described in table below
 *
 *                 - on very first call, provide allocated memory for at least :c:data:`fwork_size` entries that has been initialized using :c:func:`trlib_prepare_memory`, :c:data:`fwork_size` can be obtained by calling :c:func:`trlib_krylov_memory_size`
 *                 - can be used for different problem instances with matching dimensions and itmax without reinitialization
 *                 - leave untouched between iterations
 *
 *                 ======================== ================================= =============================================
 *                 start                    end (exclusive)                   description
 *                 ======================== ================================= =============================================
 *                 0                        1                                 stopping tolerance in case of interior solution
 *                 1                        2                                 stopping tolerance in case of boundary solution
 *                 2                        3                                 dot product :math:`\langle v, g \rangle` in current iteration
 *                 3                        4                                 dot product :math:`\langle p, Hp \rangle` in current iteration
 *                 4                        5                                 ratio between projected CG gradient and Lanczos direction in current iteration
 *                 5                        6                                 ratio between  projected CG gradient and Lanczos direction in previous iteration
 *                 6                        7                                 Lagrange multiplier :math:`\lambda_0` for trust region constraint
 *                 7                        8                                 Lagrange multiplier :math:`\lambda` for trust region constraint
 *                 8                        9                                 objective function value in current iterate
 *                 9                        10                                :math:`\langle s_i, Mp_i \rangle`
 *                 10                       11                                :math:`\langle p_i, Mp_i \rangle`
 *                 11                       12                                :math:`\langle s_i, Ms_i \rangle`
 *                 12                       13                                :math:`\sigma_i`
 *                 13                       14                                max Rayleigh quotient, :math:`\max_i \frac{\langle p, Hp \rangle}{\langle p, M p \rangle}`
 *                 14                       15                                min Rayleigh quotient, :math:`\min_i \frac{\langle p, Hp \rangle}{\langle p, M p \rangle}`
 *                 15                       16 +    :c:data:`itmax`           :math:`\alpha_i, i \ge 0`, step length CG
 *                 16 +    :c:data:`itmax`  17 +  2 :c:data:`itmax`           :math:`\beta_i, i \ge 0`, step update factor CG
 *                 17 +  2 :c:data:`itmax`  18 +  3 :c:data:`itmax`           :c:data:`neglin` for :c:func:`trlib_tri_factor_min`, just given by :math:`- \gamma_0 e_1`
 *                 18 +  3 :c:data:`itmax`  19 +  4 :c:data:`itmax`           solution :math:`h_0` of tridiagonal subproblem provided as :c:data:`sol` by :c:func:`trlib_tri_factor_min`
 *                 19 +  4 :c:data:`itmax`  20 +  5 :c:data:`itmax`           solution :math:`h` of tridiagonal subproblem provided as :c:data:`sol` by :c:func:`trlib_tri_factor_min`
 *                 20 +  5 :c:data:`itmax`  21 +  6 :c:data:`itmax`           :math:`\delta_i, i \ge 0`, curvature in Lanczos, diagonal of :math:`T` in Lanczos tridiagonalization process
 *                 21 +  6 :c:data:`itmax`  22 +  7 :c:data:`itmax`           diagonal of Cholesky of :math:`T_0 + \lambda_0 I`
 *                 22 +  7 :c:data:`itmax`  23 +  8 :c:data:`itmax`           diagonal of Cholesky of :math:`T + \lambda I`
 *                 23 +  8 :c:data:`itmax`  23 +  9 :c:data:`itmax`           :math:`\gamma_i, i \ge 1`, norm of gradients in Lanczos; provides offdiagonal of :math:`T` in Lanczos tridiagonalization process
 *                 23 +  9 :c:data:`itmax`  23 + 10 :c:data:`itmax`           offdiagonal of Cholesky factorization of :math:`T_0 + \lambda_0 I`
 *                 23 + 10 :c:data:`itmax`  23 + 11 :c:data:`itmax`           offdiagonal of Cholesky factorization of :math:`T + \lambda I`
 *                 23 + 11 :c:data:`itmax`  24 + 12 :c:data:`itmax`           :c:data:`ones` for :c:func:`trlib_tri_factor_min` and :c:func:`trlib_eigen_inverse`
 *                 24 + 12 :c:data:`itmax`  25 + 13 :c:data:`itmax`           :c:data:`leftmost` for :c:func:`trlib_tri_factor_min`
 *                 25 + 13 :c:data:`itmax`  26 + 14 :c:data:`itmax`           :c:data:`regdiag` for :c:func:`trlib_tri_factor_regularize_posdef`
 *                 26 + 14 :c:data:`itmax`  27 + 15 :c:data:`itmax`           history of convergence criteria values
 *                 27 + 15 :c:data:`itmax`  27 + 15 :c:data:`itmax` + ``fmz`` :c:data:`fwork` for :c:func:`trlib_tri_factor_min`, ``fmz`` is given by :c:func:`trlib_tri_factor_memory_size` with argument ``itmax+1``
 *                 ======================== ================================= =============================================
 *
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
 *  :param timing: gives timing details, all values are multiples of nanoseconds, provide zero allocated memory of length :c:func:`trlib_krylov_timing_size`
 *                 
 *                 ====== ============
 *                 block   description
 *                 ====== ============
 *                 0       total duration
 *                 1       :c:data:`timing` of :c:func:`trlib_tri_factor_min`
 *                 ====== ============
 *
 *  :type timing: trlib_int_t, input/output
 *  :param action: The user should perform the following action depending on :c:data:`action` and :c:data:`ityp` on the vectors he manages, see the table below.
 *                 The table makes use of the notation explained in the *User provided storage* above and the following:
 *
 *                 - :math:`i`: :c:data:`iter`
 *                 - :math:`q_j`: :math:`j`-th column of :math:`Q`
 *                 - :math:`Q_i`: matrix consisting of the first :math:`i+1` columns of :math:`Q`, :math:`Q_i = (q_0, \ldots, q_i)`
 *                 - :math:`h_i`: vector of length :math:`i+1` stored in :c:data:`fwork` at start position :c:data:`h_pointer` provided by :c:func:`trlib_krylov_memory_size`
 *                 - :math:`p \leftarrow \perp_M Q_j`: optionally :math:`M`-reorthonormalize :math:`p` against :math:`Q_j`
 *                 - :math:`g \leftarrow \texttt{rand}\perp Q_j` find nonzero random vector that is orthogonal to :math:`Q_j`
 *                 - Note that :c:macro:`TRLIB_CLA_NEW_KRYLOV` is unlikely and only occurs on problems that employ the hard case and only if :c:data:`ctl_invariant` :math:`\neq` :c:macro:`TRLIB_CLC_NO_EXP_INV`.
 *                   If you want to safe yourself from the trouble implementing this and are confident that you don't need to expand several invariant Krylov subspaces, just ensure that :c:data:`ctl_invariant` is set to :c:macro:`TRLIB_CLC_NO_EXP_INV`.
 *
 *                 =================================== ================================================= =============================
 *                 action                              ityp                                              command
 *                 =================================== ================================================= =============================
 *                 :c:macro:`TRLIB_CLA_TRIVIAL`        :c:macro:`TRLIB_CLT_CG`, :c:macro:`TRLIB_CLT_L`   do nothing
 *                 :c:macro:`TRLIB_CLA_RETRANSF`       :c:macro:`TRLIB_CLT_CG`, :c:macro:`TRLIB_CLT_L`   :math:`s \leftarrow Q_i h_i`
 *                 :c:macro:`TRLIB_CLA_INIT`           :c:macro:`TRLIB_CLT_CG`, :c:macro:`TRLIB_CLT_L`   :math:`s \leftarrow 0`, :math:`g \leftarrow g_0`, :math:`g_- \leftarrow 0`, :math:`v \leftarrow M^{-1} g`, :math:`p \leftarrow -v`, :math:`\textit{Hp} \leftarrow Hp`, :math:`\texttt{g}\_\texttt{dot}\_\texttt{g} \leftarrow \langle g, g \rangle`, :math:`\texttt{v}\_\texttt{dot}\_\texttt{g} \leftarrow \langle v, g \rangle`, :math:`\texttt{p}\_\texttt{dot}\_\texttt{Hp} \leftarrow \langle p, \textit{Hp} \rangle`, :math:`q_0 \leftarrow \frac{1}{\sqrt{\texttt{v}\_\texttt{dot}\_\texttt{g}}} v`
 *                 :c:macro:`TRLIB_CLA_UPDATE_STATIO`  :c:macro:`TRLIB_CLT_CG`                           :math:`s \leftarrow s + \texttt{flt1} \, p`
 *                 :c:macro:`TRLIB_CLA_UPDATE_STATIO`  :c:macro:`TRLIB_CLT_L`                            do nothing
 *                 :c:macro:`TRLIB_CLA_UPDATE_GRAD`    :c:macro:`TRLIB_CLT_CG`                           :math:`q_i \leftarrow \texttt{flt2} \, v`, :math:`g_- \leftarrow g`, :math:`g \leftarrow g + \texttt{flt1} \, \textit{Hp}`, :math:`v \leftarrow M^{-1} g`, :math:`\texttt{g}\_\texttt{dot}\_\texttt{g} \leftarrow \langle g, g \rangle`, :math:`\texttt{v}\_\texttt{dot}\_\texttt{g} \leftarrow \langle v, g \rangle`
 *                 :c:macro:`TRLIB_CLA_UPDATE_GRAD`    :c:macro:`TRLIB_CLT_L`                            :math:`s \leftarrow \textit{Hp} + \texttt{flt1}\, g + \texttt{flt2}\, g_-`, :math:`g_- \leftarrow \texttt{flt3}\, g`, :math:`g \leftarrow s`, :math:`v \leftarrow M^{-1} g`, :math:`\texttt{v}\_\texttt{dot}\_\texttt{g} \leftarrow \langle v, g \rangle`
 *                 :c:macro:`TRLIB_CLA_UPDATE_DIR`     :c:macro:`TRLIB_CLT_CG`                           :math:`p \leftarrow \texttt{flt1} \, v + \texttt{flt2} \, p` with :math:`\texttt{flt1} = -1`, :math:`\textit{Hp} \leftarrow Hp`, :math:`\texttt{p}\_\texttt{dot}\_\texttt{Hp} \leftarrow \langle p, \textit{Hp} \rangle`
 *                 :c:macro:`TRLIB_CLA_UPDATE_DIR`     :c:macro:`TRLIB_CLT_L`                            :math:`p \leftarrow \texttt{flt1} \, v + \texttt{flt2} \, p` with :math:`\texttt{flt2} = 0`, :math:`p \leftarrow \perp_M Q_{i-1}`, :math:`\textit{Hp} \leftarrow Hp`, :math:`\texttt{p}\_\texttt{dot}\_\texttt{Hp} \leftarrow \langle p, \textit{Hp} \rangle`, :math:`q_i \leftarrow p`
 *                 :c:macro:`TRLIB_CLA_NEW_KRYLOV`     :c:macro:`TRLIB_CLT_CG`, :c:macro:`TRLIB_CLT_L`   :math:`g \leftarrow \texttt{rand}\perp Q_{i-1}`, :math:`g_- \leftarrow 0`, :math:`v \leftarrow M^{-1} g`, :math:`\texttt{v}\_\texttt{dot}\_\texttt{g} \leftarrow \langle v, g \rangle`, :math:`p \leftarrow \frac{1}{\sqrt{\texttt{v}\_\texttt{dot}\_\texttt{g}}} v`, :math:`\texttt{p}\_\texttt{dot}\_\texttt{Hp} \leftarrow \langle p, \textit{Hp} \rangle`, :math:`q_{i+1} \leftarrow p`
 *                 :c:macro:`TRLIB_CLA_CONV_HARD`      :c:macro:`TRLIB_CLT_CG`, :c:macro:`TRLIB_CLT_L`   :math:`\texttt{v}\_\texttt{dot}\_\texttt{g} \leftarrow \langle Hs+g_0+\texttt{flt1}\, Ms, M^{-1}(Hs+g_0) + \texttt{flt1} s \rangle`
 *                 :c:macro:`TRLIB_CLA_OBJVAL`         :c:macro:`TRLIB_CLT_CG`, :c:macro:`TRLIB_CLT_L`   :math:`\texttt{g}\_\texttt{dot}\_\texttt{g} \leftarrow \tfrac 12 \langle s, Hs \rangle + \langle g, s \rangle`
 *                 =================================== ================================================= =============================
 *
 *  :type action: trlib_int_t, output
 *  :param iter: iteration counter to tell user position in vector storage
 *  :type iter: trlib_int_t, output
 *  :param ityp: iteration type, see :c:data:`action`
 *  :type ityp: trlib_int_t, output
 *  :param flt1: floating point value that user needs for actions
 *  :type flt1: trlib_flt_t, output
 *  :param flt2: floating point value that user needs for actions
 *  :type flt2: trlib_flt_t, output
 *  :param flt3: floating point value that user needs for actions
 *  :type flt3: trlib_flt_t, output
 *
 *  :returns: status flag with following meaning:
 *
 *      - :c:macro:`TRLIB_CLR_CONTINUE` no convergence yet, continue in reverse communication
 *      - :c:macro:`TRLIB_CLR_CONV_BOUND` successful exit with converged solution on boundary, end reverse communication process
 *      - :c:macro:`TRLIB_CLR_CONV_INTERIOR` successful exit with converged interior solution, end reverse communication process
 *      - :c:macro:`TRLIB_CLR_APPROX_HARD` successful exit with approximate solution, hard case occurred, end reverse communication process
 *      - :c:macro:`TRLIB_CLR_NEWTON_BREAK` exit with breakdown in Newton iteration in :c:func:`trlib_tri_factor_min`, most likely converged to boundary solution
 *      - :c:macro:`TRLIB_TTR_HARD_INIT_LAM` hard case encountered without being able to find suitable initial :math:`\lambda` for Newton iteration, returned approximate stationary point that maybe suboptimal
 *      - :c:macro:`TRLIB_CLR_ITMAX` iteration limit exceeded, end reverse communication process
 *      - :c:macro:`TRLIB_CLR_UNBDBEL` problem seems to be unbounded from below, end reverse communication process
 *      - :c:macro:`TRLIB_CLR_FAIL_FACTOR` failure on factorization, end reverse communication process
 *      - :c:macro:`TRLIB_CLR_FAIL_LINSOLVE` failure on backsolve, end reverse communication process
 *      - :c:macro:`TRLIB_CLR_FAIL_NUMERIC` failure as one of the values :c:data:`v_dot_g` or :c:data:`p_dot_Hp` are not a floating point number
 *      - :c:macro:`TRLIB_CLR_UNLIKE_CONV` exit early as convergence seems to be unlikely
 *      - :c:macro:`TRLIB_CLR_PCINDEF` preconditioner apperas to be indefinite, end reverse communication process
 *      - :c:macro:`TRLIB_CLR_UNEXPECT_INT` unexpected interior solution found, expected boundary solution, end reverse communication process
 *      - :c:macro:`TRLIB_CLR_FAIL_TTR` failure occurred in :c:func:`trlib_tri_factor_min`, check :c:data:`iwork[7]` and :c:data:`iwork[8]` for details
 *      - :c:macro:`TRLIB_CLR_FAIL_HARD` failure due to occurrence of hard case: invariant subspace encountered without local convergence and request for early termination without exploring further invariant subspaces
 *
 *  :rtype: trlib_int_t
 *
 *
 */

trlib_int_t trlib_krylov_min(
    trlib_int_t init, trlib_flt_t radius, trlib_int_t equality, trlib_int_t itmax, trlib_int_t itmax_lanczos,
    trlib_flt_t tol_rel_i, trlib_flt_t tol_abs_i,
    trlib_flt_t tol_rel_b, trlib_flt_t tol_abs_b, trlib_flt_t zero, trlib_flt_t obj_lo,
    trlib_int_t ctl_invariant, trlib_int_t convexify, trlib_int_t earlyterm,
    trlib_flt_t g_dot_g, trlib_flt_t v_dot_g, trlib_flt_t p_dot_Hp,
    trlib_int_t *iwork, trlib_flt_t *fwork, trlib_int_t refine,
    trlib_int_t verbose, trlib_int_t unicode, char *prefix, FILE *fout, trlib_int_t *timing,
    trlib_int_t *action, trlib_int_t *iter, trlib_int_t *ityp,
    trlib_flt_t *flt1, trlib_flt_t *flt2, trlib_flt_t *flt3
);

/** Prepares floating point workspace for :c:func::`trlib_krylov_min`
 *  
 *  Initializes floating point workspace :c:data:`fwork` for :c:func:`trlib_krylov_min`
 *
 *  :param itmax: maximum number of iterations
 *  :type itmax: trlib_int_t, input
 *  :param fwork: floating point workspace to be used by :c:func:`trlib_krylov_min`
 *  :type fwork: trlib_flt_t, input/output
 *  
 *  :returns: ``0``
 *  :rtype: trlib_int_t
 */

trlib_int_t trlib_krylov_prepare_memory(trlib_int_t itmax, trlib_flt_t *fwork);

/** Gives information on memory that has to be allocated for :c:func::`trlib_krylov_min`
 *  
 *  :param itmax: maximum number of iterations
 *  :type itmax: trlib_int_t, input
 *  :param iwork_size: size of integer workspace iwork that has to be allocated for :c:func:`trlib_krylov_min`
 *  :type iwork_size: trlib_int_t, output
 *  :param fwork_size: size of floating point workspace fwork that has to be allocated for :c:func:`trlib_krylov_min`
 *  :type fwork_size: trlib_int_t, output
 *  :param h_pointer: start index of vector :math:`h` that has to be used in reverse communication in action :c:macro:`TRLIB_CLA_RETRANSF`
 *  :type h_pointer: trlib_int_t, output
 *  
 *  :returns: ``0``
 *  :rtype: trlib_int_t
 */

trlib_int_t trlib_krylov_memory_size(trlib_int_t itmax, trlib_int_t *iwork_size, trlib_int_t *fwork_size, trlib_int_t *h_pointer);

/** size that has to be allocated for :c:data:`timing` in :c:func:`trlib_krylov_min`
 *
 *  :returns: ``0``
 *  :rtype: trlib_int_t
 */
trlib_int_t trlib_krylov_timing_size(void);

/** Gives pointer to negative gradient of tridiagonal problem
 *  
 *  :param itmax: itmax maximum number of iterations
 *  :type itmax: trlib_int_t, input
 *  :param gt_pointer: pointer to negative gradient of tridiagonal subproblem
 *  :type gt_pointer: trlib_int_t, output
 *
 *  :returns: ``0``
 *  :rtype: trlib_int_t
 */

trlib_int_t trlib_krylov_gt(trlib_int_t itmax, trlib_int_t *gt_pointer);

#endif
