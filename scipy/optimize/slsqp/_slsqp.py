import numpy as np
import scipy.linalg as la
from scipy.optimize import nnls


def linmin(mode: int, ax: float, bx: float, f: float, tol: float):
    c = 0.3819660112501051  # (3-sqrt(5))/2
    e = d = u = 0
    a, b = ax, bx
    eps = 1.5e-8
    x = w = v = (c*(bx - ax) + ax)
    if mode == 0:
        return 1, x

    # if mode is 2, we need to start from the middle of the loop
    if mode == 2:
        jump_to_middle = True
    else:
        jump_to_middle = False

    fx = fv = fw = f
    while True:
        if not jump_to_middle:
            m = (a + b)*0.5
            tol1 = eps*np.abs(x) + tol
            tol2 = tol1 + tol1
            print(m, tol1, tol2)
            if np.abs(x - m) <= (tol2 - (b - a)*0.5):
                return 3, x

            p = q = r = 0

            # fit parabola
            if abs(e) > tol1:
                r = (x - w)*(fx - fv)
                q = (x - v)*(fx - fw)
                p = (x - v)*q - (x - w)*r
                q += (q - r)

                if q > 0:
                    p = -p
                if q < 0:
                    q = -q
                r = e
                e = d

            # is parabola acceptable
            if (np.abs(p) >= 0.5*np.abs(q*r)) or (p <= q*(a-x)) or (
                    p >= q*(b-x)):
                # golden section step
                e = b - x if (x < m) else a - x
                d = c*e
            else:
                # parabolic interpolation step
                d = p/q
                # f must not be evaluated too close to a or b
                if ((u - a) < tol2) or ((b - u) < tol2):
                    d = -tol1 if (m - x) < 0 else tol1
                    d = np.copysign(tol1, m - x)

            # f must not be evaluated too close to x
            if np.abs(d) < tol1:
                d = np.copysign(tol1, d)
            u = x + d
            return 2, u

        # in case we jumped here
        if jump_to_middle:
            jump_to_middle = False

        # 55
        fu = f

        # update a, b, v, w, and x
        if fu > fx:
            # 60
            if u < x:
                a = u
            else:
                b = u

            if (fu <= fw) or (w == x):
                # 70
                v = w
                fv = fw
                w = u
                fw = fu
            else:
                # 80
                if (fu <= fv) or (v == x) or (v == w):
                    v = u
                    fv = fu
        else:
            if u >= x:
                a = x
            else:
                b = x
            v = w
            fv = fw
            w = x
            fw = fx
            x = u
            fx = fu


def ldl_update(A, n, sigma, z):
    """In place update of A <- A + sigma * np.outer(z, z)"""
    eps = np.spacing(1.)
    if sigma == 0:
        return

    t = 1/sigma
    ij = 0

    if sigma < 0:
        w = z.copy()
        for i in range(n):  # DO - 170
            v = w[i]
            t += v*v/A[ij]
            for k in range(i+1, n):  # DO - 160
                ij += 1
                w[k] -= v*A[ij]
            ij += 1

        if t > 0:
            t = eps/sigma

        for i in range(n):
            j = n - i
            ij -= i
            u = w[j]
            w[j] = t
            t -= u*u/A[ij]

    for i in range(n):  # DO - 270
        v = z[i]
        delta = v / A[ij]
        if sigma < 0:
            tp = w[i]
        else:
            tp = t + delta*v
        alpha = tp / t
        A[ij] *= alpha
        if i == n - 1:
            return

        beta = delta / tp
        if alpha <= 4:
            for j in range(i+1, n):
                ij += 1
                z[j] -= v*A[ij]
        else:
            gamma = t/tp
            for j in range(i+1, n):
                ij += 1
                u = A[ij]
                A[ij] = gamma*u + beta*z[j]
                z[j] -= v*u

        ij += 1
        t = tp


def ldp(G, h):
    E = np.empty([G.shape[1]+1, G.shape[0]], dtype=G.dtype)
    E[:-1, :] = G.T
    E[-1, :] = h
    # E = np.vstack([G.T, h])
    f = np.zeros([E.shape[0]], dtype=E.dtype)
    f[-1] = 1
    u, unorm = nnls(E, f)
    r = E @ u - f
    if unorm == 0.:
        return np.zeros(len(r) - 1, dtype=r.dtype)

    return -r[:-1]/r[-1], 1


def lsi(A, b, G, h):
    """Assumes A full rank. Also b, h is assumed to be 1D"""
    eps = np.spacing(1.)
    Q, R = la.qr(A)
    m, n = A.shape
    # Only deal with full rank problems
    if (np.abs(np.diag(R)) < m*n*eps).any():
        return np.zeros(n, dtype=A.dtype), 5

    bb = Q.T @ b
    # Define x = inv(R) (z - ff)
    # Transform G and h:  X R = G -> R.T X.T = G.T
    GG = la.solve_triangular(R[:n, :].T, G.T,
                             lower=True, check_finite=False).T
    hh = h - GG @ bb[:n]

    z, mode = ldp(GG, hh)
    x = la.solve_triangular(R[:n, :], z + bb[:n], check_finite=False)
    return x, 1


def lsei(A, b, E, f, G, h):
    """
    A, E, G: 2D-ndarray
    b, f, h: 1D-ndarray

    We right triangulize E and apply Q.T to A and G, i.e.;

        [E]         [0  Et]   [  ]
        [A] @ Q.T = [A2 A1] = [AA]
        [G]         [G2 G1]   [GG]

    the indices are swapped due to scipy.linalg.rq returning Et on the right
    half as upper triangular.
    """
    eps = np.spacing(1.)
    Et, Q = la.rq(E)
    m1, n = E.shape
    nd = n - m1
    # Remove the zero part coming from RQ decomposition
    Et = Et[:, nd:]
    AA = A @ Q.T
    GG = G @ Q.T

    # Only deal with full rank problems
    if (np.abs(np.diag(Et)) < m1*n*eps).any():
        return np.zeros(E.shape[0], dtype=E.dtype), 6

    x1 = la.solve_triangular(Et, f)

    if G.size == 0:  # No inequality constraints
        x2, _, _, _ = la.lstsq(AA[:, :nd], b-AA[:, nd:]@x1)
    else:
        x2, mode = lsi(AA[:, :nd], b - AA[:, nd:]@x1,
                       GG[:, :nd], h - GG[:, nd:]@f[:m1])
        # If failed
        if mode != 1:
            return np.zeros(E.shape[0], dtype=E.dtype), mode

    # Note: it's [x2 x1] due to RQ swap
    x = Q.T @ np.hstack([x2, x1])

    return x, 1
