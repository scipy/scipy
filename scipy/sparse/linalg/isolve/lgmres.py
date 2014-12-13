# Copyright (C) 2009, Pauli Virtanen <pav@iki.fi>
# Distributed under the same license as Scipy.

from __future__ import division, print_function, absolute_import

import numpy as np
from scipy.lib.six import xrange
from scipy.linalg import get_blas_funcs
from .utils import make_system

__all__ = ['lgmres']


def norm2(q):
    q = np.asarray(q)
    nrm2 = get_blas_funcs('nrm2', dtype=q.dtype)
    return nrm2(q)


def lgmres(A, b, x0=None, tol=1e-5, maxiter=1000, M=None, callback=None,
           inner_m=30, outer_k=3, outer_v=None, store_outer_Av=True):
    """
    Solve a matrix equation using the LGMRES algorithm.

    The LGMRES algorithm [1]_ [2]_ is designed to avoid some problems
    in the convergence in restarted GMRES, and often converges in fewer
    iterations.

    Parameters
    ----------
    A : {sparse matrix, dense matrix, LinearOperator}
        The real or complex N-by-N matrix of the linear system.
    b : {array, matrix}
        Right hand side of the linear system. Has shape (N,) or (N,1).
    x0  : {array, matrix}
        Starting guess for the solution.
    tol : float
        Tolerance to achieve. The algorithm terminates when either the relative
        or the absolute residual is below `tol`.
    maxiter : int
        Maximum number of iterations.  Iteration will stop after maxiter
        steps even if the specified tolerance has not been achieved.
    M : {sparse matrix, dense matrix, LinearOperator}
        Preconditioner for A.  The preconditioner should approximate the
        inverse of A.  Effective preconditioning dramatically improves the
        rate of convergence, which implies that fewer iterations are needed
        to reach a given error tolerance.
    callback : function
        User-supplied function to call after each iteration.  It is called
        as callback(xk), where xk is the current solution vector.
    inner_m : int, optional
        Number of inner GMRES iterations per each outer iteration.
    outer_k : int, optional
        Number of vectors to carry between inner GMRES iterations.
        According to [1]_, good values are in the range of 1...3.
        However, note that if you want to use the additional vectors to
        accelerate solving multiple similar problems, larger values may
        be beneficial.
    outer_v : list of tuples, optional
        List containing tuples ``(v, Av)`` of vectors and corresponding
        matrix-vector products, used to augment the Krylov subspace, and
        carried between inner GMRES iterations. The element ``Av`` can
        be `None` if the matrix-vector product should be re-evaluated.
        This parameter is modified in-place by `lgmres`, and can be used
        to pass "guess" vectors in and out of the algorithm when solving
        similar problems.
    store_outer_Av : bool, optional
        Whether LGMRES should store also A*v in addition to vectors `v`
        in the `outer_v` list. Default is True.

    Returns
    -------
    x : array or matrix
        The converged solution.
    info : int
        Provides convergence information:

            - 0  : successful exit
            - >0 : convergence to tolerance not achieved, number of iterations
            - <0 : illegal input or breakdown

    Notes
    -----
    The LGMRES algorithm [1]_ [2]_ is designed to avoid the
    slowing of convergence in restarted GMRES, due to alternating
    residual vectors. Typically, it often outperforms GMRES(m) of
    comparable memory requirements by some measure, or at least is not
    much worse.

    Another advantage in this algorithm is that you can supply it with
    'guess' vectors in the `outer_v` argument that augment the Krylov
    subspace. If the solution lies close to the span of these vectors,
    the algorithm converges faster. This can be useful if several very
    similar matrices need to be inverted one after another, such as in
    Newton-Krylov iteration where the Jacobian matrix often changes
    little in the nonlinear steps.

    References
    ----------
    .. [1] A.H. Baker and E.R. Jessup and T. Manteuffel,
             SIAM J. Matrix Anal. Appl. 26, 962 (2005).
    .. [2] A.H. Baker, PhD thesis, University of Colorado (2003).
             http://amath.colorado.edu/activities/thesis/allisonb/Thesis.ps

    """
    from scipy.linalg.basic import lstsq
    A,M,x,b,postprocess = make_system(A,M,x0,b)

    if not np.isfinite(b).all():
        raise ValueError("RHS must contain only finite numbers")

    matvec = A.matvec
    psolve = M.matvec

    if outer_v is None:
        outer_v = []

    axpy, dot, scal = None, None, None

    b_norm = norm2(b)
    if b_norm == 0:
        b_norm = 1

    for k_outer in xrange(maxiter):
        r_outer = matvec(x) - b

        # -- callback
        if callback is not None:
            callback(x)

        # -- determine input type routines
        if axpy is None:
            if np.iscomplexobj(r_outer) and not np.iscomplexobj(x):
                x = x.astype(r_outer.dtype)
            axpy, dot, scal = get_blas_funcs(['axpy', 'dot', 'scal'],
                                              (x, r_outer))

        # -- check stopping condition
        r_norm = norm2(r_outer)
        if r_norm < tol * b_norm or r_norm < tol:
            break

        # -- inner LGMRES iteration
        vs0 = -psolve(r_outer)
        inner_res_0 = norm2(vs0)

        if inner_res_0 == 0:
            rnorm = norm2(r_outer)
            raise RuntimeError("Preconditioner returned a zero vector; "
                               "|v| ~ %.1g, |M v| = 0" % rnorm)

        vs0 = scal(1.0/inner_res_0, vs0)
        hs = []
        vs = [vs0]
        ws = []
        y = None

        for j in xrange(1, 1 + inner_m + len(outer_v)):
            # -- Arnoldi process:
            #
            #    Build an orthonormal basis V and matrices W and H such that
            #        A W = V H
            #    Columns of W, V, and H are stored in `ws`, `vs` and `hs`.
            #
            #    The first column of V is always the residual vector, `vs0`;
            #    V has *one more column* than the other of the three matrices.
            #
            #    The other columns in V are built by feeding in, one
            #    by one, some vectors `z` and orthonormalizing them
            #    against the basis so far. The trick here is to
            #    feed in first some augmentation vectors, before
            #    starting to construct the Krylov basis on `v0`.
            #
            #    It was shown in [BJM]_ that a good choice (the LGMRES choice)
            #    for these augmentation vectors are the `dx` vectors obtained
            #    from a couple of the previous restart cycles.
            #
            #    Note especially that while `vs0` is always the first
            #    column in V, there is no reason why it should also be
            #    the first column in W. (In fact, below `vs0` comes in
            #    W only after the augmentation vectors.)
            #
            #    The rest of the algorithm then goes as in GMRES, one
            #    solves a minimization problem in the smaller subspace
            #    spanned by W (range) and V (image).
            #
            #    XXX: Below, I'm lazy and use `lstsq` to solve the
            #    small least squares problem. Performance-wise, this
            #    is in practice acceptable, but it could be nice to do
            #    it on the fly with Givens etc.
            #

            #     ++ evaluate
            v_new = None
            if j < len(outer_v) + 1:
                z, v_new = outer_v[j-1]
            elif j == len(outer_v) + 1:
                z = vs0
            else:
                z = vs[-1]

            if v_new is None:
                v_new = psolve(matvec(z))
            else:
                # Note: v_new is modified in-place below. Must make a
                # copy to ensure that the outer_v vectors are not
                # clobbered.
                v_new = v_new.copy()

            #     ++ orthogonalize
            hcur = []
            for v in vs:
                alpha = dot(v, v_new)
                hcur.append(alpha)
                v_new = axpy(v, v_new, v.shape[0], -alpha)  # v_new -= alpha*v
            hcur.append(norm2(v_new))

            if hcur[-1] == 0:
                # Exact solution found; bail out.
                # Zero basis vector (v_new) in the least-squares problem
                # does no harm, so we can just use the same code as usually;
                # it will give zero (inner) residual as a result.
                bailout = True
            else:
                bailout = False
                v_new = scal(1.0/hcur[-1], v_new)

            vs.append(v_new)
            hs.append(hcur)
            ws.append(z)

            # XXX: Ugly: should implement the GMRES iteration properly,
            #      with Givens rotations and not using lstsq. Instead, we
            #      spare some work by solving the LSQ problem only every 5
            #      iterations.
            if not bailout and j % 5 != 1 and j < inner_m + len(outer_v) - 1:
                continue

            # -- GMRES optimization problem
            hess = np.zeros((j+1, j), x.dtype)
            e1 = np.zeros((j+1,), x.dtype)
            e1[0] = inner_res_0
            for q in xrange(j):
                hess[:(q+2),q] = hs[q]

            y, resids, rank, s = lstsq(hess, e1)
            inner_res = norm2(np.dot(hess, y) - e1)

            # -- check for termination
            if inner_res < tol * inner_res_0:
                break

        # -- GMRES terminated: eval solution
        dx = ws[0]*y[0]
        for w, yc in zip(ws[1:], y[1:]):
            dx = axpy(w, dx, dx.shape[0], yc)  # dx += w*yc

        # -- Store LGMRES augmentation vectors
        nx = norm2(dx)
        if store_outer_Av:
            q = np.dot(hess, y)
            ax = vs[0]*q[0]
            for v, qc in zip(vs[1:], q[1:]):
                ax = axpy(v, ax, ax.shape[0], qc)
            outer_v.append((dx/nx, ax/nx))
        else:
            outer_v.append((dx/nx, None))

        # -- Retain only a finite number of augmentation vectors
        while len(outer_v) > outer_k:
            del outer_v[0]

        # -- Apply step
        x += dx
    else:
        # didn't converge ...
        return postprocess(x), maxiter

    return postprocess(x), 0
