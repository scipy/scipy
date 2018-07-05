""" The module has the functions necessary to perform derivation on noisy discrete data.

It follows gaussian estimation to estimate the noise an reduce the noise in the derived data.

"""

import numpy as np
from scipy import linalg
from scipy import special


def derest(*args):
    r"""Estimates the smoothed signal of a noisy data along with it's first and second derivatives using gaussian
    estimation.

    :param args:
    t : list_like (double)
        The input 't' is a row or column list of the time steps of the measurements.
    y : list_like (double)
        The input 'y' is a row or column of the measurements taken to be derived.
    tt : array-like (double), optional
        The input 'tt' is an optional array of the desired time steps where the measurement derivatives are expected.
        Default : the input array 't'.
    nx : integer, optional
        Optional input parameter used to control the number of states in the state space model. Default : 3.

    :return:
    m : tuple_like (double)
        Returns the smoothed signal (m[0, :]), the first order differentiation of the smoothed signal (m[1, :]) and the
        second order differentiation of the smoothed signal (m[2, :]) at time steps 't'.
    mm : tuple_like (double)
        Returns the smoothed signal (mm[0, :]), the first order differentiation of the smoothed signal
        (mm[1, :]) and the second order differentiation of the smoothed signal (mm[2, :]) at time steps 'tt'.
    """
    pass

    # Pre define default parameters.
    nx = 3
    t = np.asarray(args[0])
    tt = t
    yy = np.asarray(args[1])

    # Assign input values from the user.
    if len(args) >= 3:
        tt = args[2]

    if len(args) >= 4:
        nx = args[3]

    n = np.size(t)

    # Check for validity of the input
    if n != np.size(yy):
        print('T and Y must be of the same length')
        quit()

    if any(np.diff(t)) <= 0:
        print('t values must be non descending and unique.')
        quit()

    if any(np.diff(tt)) <= 0:
        print('tt values must be non descending and unique.')
        quit()

    if n < 10:
        print('y needs atkeast 10 values')
        quit()

    # Prepare indices for dense output implementation.
    nk = len(t)
    dt = np.ediff1d(t)

    oldtt = np.unique(tt)
    tt = np.union1d(t, oldtt)

    dtt = np.ediff1d(tt)
    ind1 = __find_indices(tt, t)
    ind2 = __find_indices(tt, oldtt)

    if nx < 1:
        print('NX must be atleast 1')
        quit()

    h = np.c_[[1], np.zeros((1, (nx - 1)))]
    qbar = np.flip((np.flip((linalg.hilbert(nx)), 0)), 1)
    qbarsqrtin = np.linalg.cholesky(qbar)

    nfit = 10

    tfit = t[0:10] - t[0]
    yfit = yy[0:10]

    p, s, _, _, _ = np.polyfit(tfit, yfit, 1, full=True)

    if nx == 1:
        m = np.array(p[1])[np.newaxis]
    else:
        m = np.array(np.vstack(((np.vstack((p[1], p[0]))), (np.zeros(((nx - 2), 1))))))[np.newaxis]
    r = s/(nfit - 1)

    eps = np.finfo(float).eps
    p = eps*np.eye(nx)
    psqrt = np.sqrt(eps)*np.eye(nx)
    yspan = yy.max() - yy.min()
    tspan = t[-1] - t[0]
    qlist = (np.logspace(-3, 10, 21)*(yspan ** 2))/(tspan ** ((2 * nx) - 1))
    philist = np.zeros(np.size(qlist))

    for iq in range(0, np.size(qlist)):
        philist[iq], __, __, __, __, __ = __kalman(m, psqrt, dt, qlist[iq], qbarsqrtin, h, r, yy)

    iq_best = np.argmin(philist)
    q = qlist[iq_best]

    maxit_em = 100
    threshold = 1e-3
    x = np.zeros((nx, nk))
    mdense = np.zeros((nx, np.size(tt)))
    emstopflag = False

    xold = []

    for i_em in range(0, maxit_em):

        phi, m, psqrt, mkf, ppkfsqrt, pkfsqrt = __kalman(m, psqrt, dt, q, qbarsqrtin, h, r, yy)
        print("iteration ", i_em + 1, ", phi = ", phi)
        x[:, (nk - 1)] = m.ravel()
        r = 0
        hpsqrt = h.dot(psqrt)

        r += np.power((yy[nk - 1] - h.dot(m)), 2) + hpsqrt.dot(hpsqrt.T)

        dq = 0
        mdense[:, (np.size(tt) - 1)] = m.ravel()
        densek = np.size(dtt)

        for k in range((nk - 2), -1, -1):
            difft = dt[k]
            toeptemp = (np.power(difft, (range(0, nx))))/(special.gamma(range(1, nx + 1)))
            a = linalg.toeplitz(h, toeptemp)
            qtemp = np.power(difft, (np.array(range(((2 * nx) - 1), 0, -2)))/2)
            qbarsqrt = np.diag(qtemp).dot(qbarsqrtin)
            qsqrt = np.sqrt(q) * qbarsqrt

            ppkf = ppkfsqrt[:, :, (k + 1)].dot(ppkfsqrt[:, :, (k + 1)].T)
            pkf = pkfsqrt[:, :, k].dot(pkfsqrt[:, :, k].T)
            g = __mrdivide((pkf.dot(a.T)), ppkf)
            mkp1 = m
            psqrtkp1 = psqrt
            m = mkf[:, k][np.newaxis].T + g.dot((mkp1 - a.dot(mkf[:, k])[np.newaxis].T))
            x[:, k] = m.ravel()
            workmat = np.vstack([np.hstack([(pkfsqrt[:, :, k].T.dot(a.T)), pkfsqrt[:, :, k].T]),
                                    np.hstack([qsqrt.T, np.zeros((nx, nx))]),
                                    np.hstack([np.zeros((nx, nx)), psqrt.T.dot(g.T)])])
            __, temp = np.linalg.qr(workmat)
            psqrt = temp[nx:(2 * nx), nx:(2 * nx)].T
            dq1 = mkp1 - a.dot(m)
            dq2 = (a.dot(g)).dot(qsqrt)
            dq3 = (np.eye(nx) - a.dot(g)).dot((psqrtkp1 + a.dot(pkfsqrt[:, :, k])))
            dqhat = dq1.dot(dq1.T) + dq2.dot(dq2.T) + dq3.dot(dq3.T)
            ddq = np.trace(__mrdivide((np.linalg.solve(qbarsqrt, dqhat)), qbarsqrt.T))
            dq += ddq
            hpsqrt = h.dot(psqrt)

            r += np.power((yy[k] - h.dot(m)), 2) + hpsqrt.dot(hpsqrt.T)
            dttou = 0
            for j in range(ind1[k + 1], (ind1[k]), -1):
                dttou += dtt[(densek - 1)]
                toeptemp = (np.power((difft - dttou), (range(0, nx))))/(special.gamma(range(1, nx + 1)))
                atou = linalg.toeplitz(h, toeptemp)
                toeptemp = (np.power(dttou, (range(0, nx))))/(special.gamma(range(1, nx + 1)))
                a1tou = linalg.toeplitz(h, toeptemp)
                qtemp = np.power((difft - dttou), (np.array(range(((2 * nx) - 1), 0, -2))) / 2)
                qtoubarsqrt = np.diag(qtemp).dot(qbarsqrtin)
                qtousqrt = np.sqrt(q) * qtoubarsqrt
                qktou = qtousqrt.dot(qtousqrt.T)

                pktouk = atou.dot(pkf).dot(atou.T) + qktou
                mktou = atou.dot(mkf[:, k])
                gktou = __mrdivide(pktouk.dot(a1tou.T), ppkf)

                mkptou = mktou[np.newaxis].T + gktou.dot((mkp1 - a.dot(mkf[:, k][np.newaxis].T)))
                mdense[:, (densek - 1)] = mkptou.ravel()
                densek -= 1
            mdense[:, densek] = m.ravel()
        q = dq / (nk - 1) / nx
        r = r/n

        if i_em > 0:
            xdiff = np.linalg.norm((x[0, :] - xold))/np.linalg.norm(xold)
            emstopflag = (xdiff < threshold)

        xold = x[0, :]

        if emstopflag:
            print('number of EM iteration = ', i_em + 1)
            print(' q = ', q)
            print(' r = ', r)
            break
    mout = mdense[:, ind2]

    return x, mout


# Square root Kalman Filter implementation
def __kalman(m,
           psqrt,
           dt,
           q,
           qbarsqrtin,
           h,
           r,
           y):
    r"""Estimates the Kalman Filter for the signal using the square root Kalman Filter implementation.

    :param m: array-like (double)
        The 'mean' value. Size of [nx, 1]. - m(1, 0)
    :param psqrt: multidimensional-array-like (double)
        The square-root of the covariance matrix. Size of [nx, nx]. - p(1,0)^0.5
    :param dt: array-like (double)
        Array of the delta-time, or time difference between successive points in the time series.
    :param q: (double)
        Process noise variance.
    :param qbarsqrtin: multidimensional-array-like (double)
        Process covariance matrix. Size of [nx, nx].
    :param h: array-like (double)
        Measurement model matrix. Size of [1, nx]
    :param r: (double)
        Measurement noise.
    :param y: array-like (double)
        Measurement array.
    :return:
    phi: (double)
        Cost estimate.
    m: array-like (double)
        The 'mean' value. Size of [nx, 1]. - m(T, T)
    psqrt: multidimensional-array-like (double)
         The square-root of the covariance matrix. Size of [nx, nx]. - p(T, T)^0.5
    mkf: array-like (double)
        Kalman filter state.
    ppkfsqrt: multidimensional-array-like (double)
        Square root of filter prior state covariance.
    pkfsqrt: multidimensional-array-like (double)
        Square root of filter covariance.

    """
    nx = np.size(m)
    nk = np.size(y)

    # Pre allocate memory
    mkf = np.zeros((nx, nk))
    ppkfsqrt = np.zeros((nx, nx, nk))
    pkfsqrt = np.zeros((nx, nx, nk))

    phi = 0
    for kk in range(0, nk):
        ppkfsqrt[:, :, kk] = psqrt

        # Measurement update
        qt1 = np.c_[np.array(np.sqrt(r)), h.dot(psqrt)]
        qt2 = np.c_[np.zeros((nx, 1)), psqrt]
        qrmat = np.r_[qt1, qt2]
        __, temp = linalg.qr(qrmat.T)
        psqrt = temp[1:, 1:].T
        s = temp[0, 0] ** 2
        k = np.array(temp[0, 1:].T/temp[0, 0])[np.newaxis]
        k = k.T
        v = y[kk] - h.dot(m)
        m += k*v

        phi += (0.5 * (v ** 2))/s + 0.5 * np.log(2 * np.math.pi * s)
        mkf[:, kk] = m.ravel()
        pkfsqrt[:, :, kk] = psqrt

        # State update
        if kk < nk - 1:
            difft = dt[kk]
            a = linalg.toeplitz(h, ((np.power(difft, (range(0, nx))))/(special.gamma(range(1, nx + 1)))))
            rantemp = np.array(range(((2 * nx) - 1), 0, -2))
            diagels = np.power(difft, rantemp/2)
            qbarsqrt = np.diag(diagels).dot(qbarsqrtin)
            qsqrt = np.sqrt(q) * qbarsqrt
            am = a.dot(m)
            am = am.reshape(3, 1)
            m = am
            __, temp = linalg.qr(np.vstack(((psqrt.T.dot(a.T)), qsqrt.T)))
            psqrt = temp[0:nx, 0:nx].T

    return phi, m, psqrt, mkf, ppkfsqrt, pkfsqrt


# Find indices of elements of one array in another
def __find_indices(x,
                 y):
    r"""
    :param x: multidimensional-array-like (double)
        Array containing the elements where the indices are to be searched.
    :param y: multidimensional-array-like (double)
        Array containing the elements for those which indices are to be searched for.
    :return: multidimensional-array-like (integer)
        Indices of the elements of the second array in the first array.
    """
    xsorted = np.argsort(x)
    ypos = np.searchsorted(x[xsorted], y)
    inds = xsorted[ypos]

    return inds


# matrix right divide operator
def __mrdivide(a,
             b):
    r"""
    :param a: multidimensional-array-like (double)
        Matrix or column array with the same number of columns as the 'b' matrix or array.
    :param b: multidimensional-array-like (double)
        Matrix or column array with the same number of columns as the 'a' matrix or array.
    :return:
    x: multidimensional-array-lik (double)
        Returns the solution to a/b
    """
    # a/b
    xt = np.linalg.solve(b.T, a.T)
    x = xt.T
    return x
