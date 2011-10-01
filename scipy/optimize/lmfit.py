"""
====================================
Least-Squares Fitting (:mod:`lmfit`)
====================================

This module contains a class-based Levenberg-Marquardt fitter
`Fit` and a function `leastsq` which is a function-based
interface to the same functionality.  """
from __future__ import division
from numpy import (array, diagonal, dot, empty, eye, finfo, maximum, meshgrid,
    sqrt, zeros_like, diag as npdiag)
from itertools import count
from scipy.linalg import solve_triangular, qr, qr_multiply, norm
from _qrsolv import qrsolv

__all__ = ["Fit", "FitError", "InvalidParameter", "leastsq"]

class FitError(Exception):
    """ The error raised while fitting

    Parameters
    ----------
    exitcode : int
        The same as in the class `Fit`
    message : str
        A description of the error
    """
    def __init__(self, exitcode, message):
        Exception.__init__(self, message)
        self.exitcode = exitcode

class InvalidParameter(Exception):
    """ Indicate a parameter is invalid

    Parameters
    ----------
    number : int
        number is the number of the parameter, or 
        a list of numbers of parameters, which are erroneous.
        Set to None (the default) if no specific parameter is a problem.
    """
    def __init__(self, number=None):
        Exception.__init__(self)
        self.number = number

class Fit(object):
    """ Least-Squares Fitting using Levenberg-Marquardt

    This class contains methods to minimize the sum of the squares 
    of `M` nonlinear functions with `N` parameters by a modification of the
    Levenberg-Marquardt algorithm. The user must overwrite the method
    `func` that calculates the functions, and, optionally, the method
    `jacobian` to calculate its partial derivatives.

    Parameters
    ----------

    params : vector, length `N`
        contains the last successful set of parameters during fit.
        Useful if fit failed.
    running : boolean
        if set to False, fit will stop at next iteration. Use this to
        stop the fit from another thread.
    eps : float
        the stepsize used in the forward-difference approximation
    iterations : int
        number of iterations needed in last fit
    nfev : int
        the number of times the function has been evaluated
    fvec : vector, length `M`
        the function value at the last successful parameter optimization
    scale : vector, length `N`
        the scale of the changes of the parameters in each iteration 
    exitcode : int
        the reason the fitting has stopped. Possible values:
        (see also the parameter names of the method `fit`)

        :0:
            an unknown (typically user) error has occurred
        :1:
            both actual and predicted relative reductions
            in the sum of squares are at most `ftol`.

        :2:
            relative error between two consecutive iterates
            is at most `xtol`.

        :3:
            conditions for 1 and 2 both hold.

        :4:
            the cosine of the angle between function value and any
            column of the Jacobian is at most `gtol` in
            absolute value.

        :5:
            number of function evaluations has reached its limit

        :6:
            `ftol` is too small. no further reduction in
            the sum of squares is possible.

        :7:
            `xtol` is too small. No further improvement in
            the approximate solution `x` is possible.

        :8:
            `gtol` is too small. The function value is orthogonal to the
            columns of the jacobian to machine precision.

        :9:
            number of iterations has reached its limit
    """

    def fjac(self, j, x, f0):
        r""" Calculate the derivative of the function by one parameter.

        This method, which is used by `jacobian`, calculates the
        derivative

        .. math::
           \frac{\partial f}{\partial x_j}

        by a forward-difference approximation.

        The user should overwrite this method (or `jacobian`) if another
        way to calculate the derivative of the function is appropriate.

        Parameters
        ----------

        j : int
            The number of the function's parameter to derive by.
        x : vector, length `N`
            The position to evaluate the derivative at.
        f0 : vector, length `M`
            The value of the function at `x`.

        Returns
        -------
        out : vector, length `M`
            The derivative of the function at `x` by `x[j]`.
        """
        tmp = x[j]
        h = sqrt(self.eps) * abs(tmp)
        if h == 0:
            h = sqrt(self.eps)
        for c in [1, -1, 0.1, -0.1]:
            x[j] = tmp + c * h
            try:
                self.nfev += 1
                f1 = self.func(x)
            except InvalidParameter:
                continue
            x[j] = tmp
            return (f1 - f0) / (c * h)
        return zeros_like(f0)

    def jacobian(self, x, fvec, ret):
        r""" Calculate the Jacobian matrix

        The Jacobian matrix is a[i, j] with ::

            a[i, j] = \left.\frac{\partial f_i}{\partial x_j}\right|_x

        This default implementation uses `fjac` to calculate the
        derivatives. The user should overwrite this method (or
        `fjac`) if another method to calculate the Jacobian is
        more appropriate.

        Parameters
        ----------

        x : vector, length `N`
            is the value the Jacobian should be evaluated at
        fvec : vector, length `M`
            is the function value at `x`
        ret : matrix, `M` by `N`
            is a pre-allocated matrix into which the Jacobian should
            be written. It may be ignored.
        
        Returns
        -------
        ret : matrix, `M` by `N` 
            the Jacobian at `x` """
        for j in range(len(x)):
            ret[:, j] = self.fjac(j, x, fvec)
        return ret

    def func(self, x):
        """ The function to be fitted to

        This method is only a placeholder, it is to be overwritten by
        the user with the function to be fitted. Therefore, the following
        is a description of the function that the user has to supply,
        not the actual function which does nothing.

        Parameters
        ----------
        x : vector, length `N`
            the parameters where to calculate the function at

        Returns
        -------
        fvec : vector, length `M`
            the value of the function at `x`. 

        Raises
        ------
        InvalidParameter
            This method may raise the exception `InvalidParameter`, indicating
            that one or more of the values in `x` are not appropriate.
            The user may pass the number of the offending parameter as a 
            parameter to `InvalidParameter`, or a list of numbers in the
            case of multiple offending parameters.
        """
        raise NotImplemented("I don't know which function to fit!")

    def _qrsolv(self, r, ipvt, diag, qtb):
        r""" Solve the Levenberg-Marquart problem. Internal method.

        Given an `M` by `N` matrix `A`, an `N` by `N` diagonal matrix `D`,
        and an `M`-vector `b`, the problem is to determine an `x` which
        solves the system ::

            A x = b,    D x = 0 

        in the least squares sense.

        This subroutine completes the solution of the problem
        if it is provided with the necessary information from the
        QR factorization, with column pivoting, of `A`. That is, if
        `AP = QR`, where `P` is a permutation matrix, `Q` has orthogonal
        columns, and `R` is an upper triangular matrix with diagonal
        elements of nonincreasing magnitude, then `_qrsolv` expects
        the full upper triangle of `R`, the permutation matrix `P`,
        and the first n components of ``Q.T b``. The system
        ``Ax = b, Dx = 0``, is then equivalent to ::

              R z = Q.T b,  P.T D P z = 0,

        where ``x = Pz``. If this system does not have full rank,
        then a least squares solution is obtained. On output `_qrsolv`
        also provides an upper triangular matrix S such that ::

              P.T (A.T A + D D) P = S.T S.

        `S` is computed within `_qrsolv` and may be of separate interest.

        Parameters
        ----------
        r : ndarray, `N` by `N`
            On input the full upper triangle must contain the full upper
            triangle of the matrix `R`. 

        ipvt : integer array
            defines the permutation matrix `P` such that `A P = Q R`.
            Column `j` of `P` is column `ipvt[j]` of the identity matrix.

        diag : vector, length `N`
            contains the diagonal elements of the matrix `D`.

        qtb : vector, length `N`
            contains the first `N` elements of the vector ``Q.T b``.

        Returns
        -------
        x : vector, length `N`
            contains the least squares solution of the system
            ``A x = b, D x = 0``.

        s : matrix, `N` by `N`
            contains the aforementioned matrix `S`. """
        s = empty((r.shape[0], r.shape[1] + 1), dtype=r.dtype)
        s[:, :-1] = r
        s[:, -1] = qtb
        # Eliminate the diagonal matrix D using a givens rotation.
        s = qrsolv(s, diag[ipvt])
        condition = s.diagonal() == 0
        if condition.any():
            nsing = condition.nonzero()[0][0]
        else:
            nsing = len(ipvt)
        x = zeros_like(qtb)
        x[ipvt[:nsing]] = solve_triangular(s[:nsing, :nsing], s[:, -1])
        return x, s[:, :-1]
            
    def _lmpar(self, r, qtb, ipvt, diag, delta, par):
        r""" Calculate the Levenberg-Marquardt parameter. Internal method.

        Given an `M` by `N` matrix `A`, an `N` by `N` nonsingular diagonal
        matrix `D`, an `M`-vector `b`, and a positive number `delta`,
        the problem is to determine a value for the parameter
        `par` such that if `x` solves the system ::

              A x = b,    sqrt(par) D x = 0,

        in the least squares sense, and `dxnorm` is the Euclidean
        norm of `D x`, then either `par` is zero and ::

              10 (dxnorm - delta) <= delta,

        or `par` is positive and

              10 |dxnorm - delta| <= delta.

        This method completes the solution of the problem
        if it is provided with the necessary information from the
        QR factorization, with column pivoting, of `A`. That is, if
        `AP = QR`, where `P` is a permutation matrix, `Q` has orthogonal
        columns, and `R` is an upper triangular matrix with diagonal
        elements of nonincreasing magnitude, then `_lmpar` expects
        the full upper triangle of `R`, the permutation matrix `P`,
        and the first n components of ``Q.T b``. On output
        `_lmpar` also provides an upper triangular matrix `S` such that ::

              P.T (A.T A + par D D) P = S.T S.

        `S` is employed within `_lmpar` and may be of separate interest.

        Only a few iterations are generally needed for convergence
        of the algorithm. If, however, the limit of 10 iterations
        is reached, then the output par will contain the best
        value obtained so far.

        Parameters
        ----------
        r : matrix, `N` by `N`
            On input the full upper triangle must contain the full upper
            triangle of the matrix `R`.

        ipvt : integer array
            defines the permutation matrix `P` such that ``A P = Q R``.
            Column `j` of `P` is column ``ipvt[j]`` of the identity matrix.

        diag : vector, length `N`
            contains the diagonal elements of the matrix `D`.

        qtb : vector, length `N`
            contains the first `N` elements of the vector ``q.T b``.

        delta : float, > 0
            specifies an upper bound on the Euclidean norm of ``D x``.

        par : float, > 0
            contains an initial estimate of the Levenberg-Marquardt
            parameter.

        Returns
        -------
        par : float
            The optimal parameter found
        x : vector, length `N`
            contains the least squares solution of the system 
            ``A x = b,    sqrt(par) D x = 0``,
            for the output `par`. """
        dwarf = finfo(r.dtype).tiny
        # Compute and store in x the Gauss-Newton direction. If the
        # jacobian is rank-deficient, obtain a least squares solution.
        condition = diagonal(r) == 0
        if condition.any():
            nsing = condition.nonzero()[0][0]
        else:
            nsing = r.shape[1]
        x = zeros_like(qtb)
        x[ipvt[:nsing]] = solve_triangular(r[:nsing, :nsing], qtb[:nsing])
        # Evaluate the function at the origin, and test
        # for acceptance of the gauss-newton direction.
        dxnorm = norm(diag * x)
        fp = dxnorm - delta
        if fp < 0.1 * delta:
            return 0, x
        # If the jacobian is not rank deficient, the Newton
        # step provides a lower bound, parl, for the zero of
        # the function. Otherwise set this bound to zero.
        parl = 0
        if nsing == r.shape[0]:
            parl = fp / delta / (solve_triangular(r, 
                (diag ** 2 * x)[ipvt] / dxnorm, trans=1) ** 2).sum()
        # Calculate an upper bound, paru, for the zero of the function.
        gnorm = norm(dot(qtb, r) / diag[ipvt])
        paru = gnorm / delta
        if paru == 0:
            paru = dwarf / min(delta, 0.1)
        # If the input par lies outside of the interval (parl, paru),
        # set par to the closer endpoint.
        par = max(par, parl)
        par = min(par, paru)
        if par == 0:
            par = gnorm / dxnorm 
        for iter in count():
            # Evaluate the function at the current value of par.
            if par == 0:
                par = max(dwarf, 0.001 * paru)
            x, s = self._qrsolv(r, ipvt, sqrt(par) * diag, qtb)
            dxnorm = norm(diag * x)
            temp = fp
            fp = dxnorm - delta
            # If the function is small enough, accept the current value
            # of par. Also test for the exceptional cases where parl
            # is zero or the number of iterations has reached 10.
            if (abs(fp) <= 0.1 * delta or (parl == 0 and fp <= temp < 0)
                    or iter == 10):
                return par, x
            # Compute the Newton correction.
            parc = fp / delta / (
                solve_triangular(s, (diag ** 2 * x)[ipvt] / dxnorm, trans=1
                        ) ** 2).sum()
            if fp > 0:
                parl = max(parl, par)
            else:
                paru = min(paru, par)
            par = max(parl, par + parc)

    def plot(self, x, f, good):
        """ plot intermediate results

        This method, which does nothing in its default implemenation,
        is called by the fitting algorithm for every new parameters tested,
        and may be overwritten to plot the function while fitting.

        Parameters
        ----------
        x : vector, length `N`
            Current estimate of the fitting algorithm
        f : vector, length `M`
            The function at `x`
        good : boolean
            Whether this estimation is better (True) or worse (False) than
            the last estimate. """
        pass

    def fit(self, x, ftol=1.49012e-8, xtol=1.49012e-8, gtol=0.0, maxfev=None,
                epsfcn=None, factor=100, diagin=None, iterations=100):
        """ The Levenberg-Marquardt algorithm

        This method performs the actual Levenberg-Marquardt minimization
        this class is written for.

        Parameters
        ----------
        x : vector, length N
            an initial estimate of the solution vector

        ftol : float, optional
            Termination occurs when both the actual and predicted
            relative reductions in the sum of squares are at most
            `ftol`.  Therefore, `ftol` measures the relative error
            desired in the sum of squares.
     
        xtol : float, > 0, optional
            Termination occurs when the relative error between two
            consecutive iterates is at most `xtol`. Therefore, `xtol`
            measures the relative error desired in the approximate
            solution.
     
        gtol : float, > 0, optional
            Termination occurs when the cosine of the angle between
            the function and any column of the jacobian is at most
            `gtol` in absolute value. Therefore, `gtol` measures the
            orthogonality desired between the function vector and the
            columns of the Jacobian.
     
        iterations : int, optional
            is the maximum number of iterations.
     
        factor : float, > 0
            used in determining the initial step bound. this bound is
            set to the product of factor and the Euclidean norm of
            `diag * x` if nonzero, or else to factor itself. In most
            cases factor should lie in the interval 0.1 < x < 100.
            100 is a generally recommended value.

        diagin : vector, length N, optional
            contains positive entries that serve as multiplicative
            scale factors for the variables.  

        epsfcn : float, > 0, optional
            is an variable used in determining a suitable step length
            for the forward-difference approximation. This approximation
            assumes that the relative errors in the functions are of the
            order of epsfcn. """
        x = 1. * array(x)
        self.params = x
        self.nfev = 1
        self.fvec = self.func(x)
        fnorm = norm(self.fvec)
        par = 0
        epsmch = finfo(x.dtype).eps
        if epsfcn is None:
            epsfcn = epsmch
        self.eps = epsfcn
        self.running = True
        self.exitcode = 0
        fjac = empty((len(self.fvec), len(x)), dtype=self.fvec.dtype, order="F")

        for self.iterations in range(iterations):
            fjac = self.jacobian(x, self.fvec, fjac)
            
            fjacnorm = sqrt((fjac ** 2).sum(0))
            qtf, r, ipvt = qr_multiply(fjac, self.fvec, pivoting=True,
                overwrite_a=True, overwrite_c=True)

            if self.iterations == 0:
                if diagin is None:
                    self.scale = fjacnorm.copy()
                    self.scale[self.scale == 0] = 1
                else:
                    self.scale = diagin
                xnorm = norm(self.scale * x)
                delta = factor * xnorm
                if delta == 0:
                    delta = factor 
            fjacnorm[fjacnorm == 0] = -1
            gnorm = (abs(dot(r.T, qtf) / fnorm) / fjacnorm[ipvt]).max()
            if gnorm < gtol:
                self.exitcode = 4
                return x
            if diagin is None: 
                self.scale = maximum(self.scale, fjacnorm) 
            ratio = 0
            while ratio < 0.0001:
                par, p = self._lmpar(r, qtf, ipvt,
                                    self.scale, delta, par)
                testx = x - p
                pnorm = norm(self.scale * p)
                try:
                    self.nfev = 1
                    testf = self.func(testx)
                except InvalidParameter as error:
                    if error.number is not None:
                        self.scale[error.number] *= 2
                    else:
                        delta = 0.5 * min(delta, 10 * pnorm)
                        par *= 2
                    if self.running:
                        continue
                    else:
                        return x
                if self.iterations == 0:
                    delta = min(delta, pnorm)
                testfnorm = norm(testf)
                actual_reduction = -1
                if 0.1 * testfnorm < fnorm:
                    actual_reduction = 1 - (testfnorm / fnorm) ** 2
                temp1 = norm(dot(r, p[ipvt])) / fnorm
                temp2 = sqrt(par) * pnorm / fnorm
                predicted_reduction = temp1 ** 2 + 2 * temp2 ** 2
                directional_deriv = -(temp1 ** 2 + temp2 ** 2)
                ratio = 0
                if predicted_reduction != 0:
                    ratio = actual_reduction / predicted_reduction
                if ratio <= 0.25:
                    if actual_reduction > 0:
                        temp = 0.5
                    else:
                        temp = 0.5 * directional_deriv / (
                            directional_deriv + 0.5 * actual_reduction)
                    if 0.1 * testfnorm >= fnorm or temp < 0.1:
                        temp = 0.1
                    delta = temp * min(delta, 10 * pnorm)
                    par /= temp
                elif par == 0 or ratio > 0.75:
                    delta = 2 * pnorm
                    par = 0.5 * par
                if ratio >= 0.0001:
                    x = testx
                    self.params = x
                    self.fvec = testf
                    fnorm = testfnorm
                    xnorm = norm(self.scale * x)
                    self.plot(x, self.fvec, True)
                else:
                    self.plot(testx, testf, False)
                c1 = (abs(actual_reduction) <= ftol and
                      predicted_reduction <= ftol and ratio <= 2)
                c2 = delta <= (xtol * xnorm)
                self.exitcode = 2 * c2 + c1
                if c1 or c2 or not self.running:
                    return x
                if gnorm < epsmch:
                    self.exitcode = 8
                    raise FitError(8, "gtol=%f is too small, func(x) is " 
                        "orthogonal to the columns of the Jacobian to " 
                        "machine precision." % gtol)
                if delta < epsmch * xnorm:
                    self.exitcode = 7
                    raise FitError(7, "xtol=%f is too small, no further " 
                        "improvement in the approximate solution is possible."                          % xtol)
                if (abs(actual_reduction) < epsmch and
                        predicted_reduction < epsmch and 0.5 * ratio < 1):
                    self.exitcode = 6
                    raise FitError(6, "ftol=%f is too small, no further " 
                        "reduction in the sum of squares is possible." % ftol)
                if self.nfev > maxfev and maxfev is not None:
                    self.exitcode = 5
                    raise FitError(5,
                        "Number of function evaluations has reached %d." %
                        self.nfev)
        self.exitcode = 9
        raise FitError(9, "Number of iterations has reached %d." % iterations)

    def covar(self, x, eps=0):
        """ calculate the covariance matrix of the solution

        This method gives an estimation of the error of the fitting
        parameters.

        Parameters
        ----------
        x : vector, length N
            the value to calculate the errors at
        eps : float, optional
            the threshold below which a parameter is considered
            linear-dependent (ie, has no influence on a fitting result).
        f : vector, length M
            the function value at `p`, calculated internally if None """
        f = self.func(x)
        if not hasattr(self, "eps"):
            self.eps = finfo(f.dtype).eps
        fdjac = self.jacobian(x, f, empty((len(f), len(x)), dtype=f.dtype))
        r, jpvt = qr(fdjac, pivoting=True, mode="r")
        condition = abs(r.diagonal()) <= r[0, 0] * eps
        if condition.any():
            nsing = condition.nonzero()[0][0]
        else:
            nsing = len(jpvt)
        r = r[:nsing, :nsing]
        ri = eye(r.shape[0], dtype=r.dtype)
        ri = solve_triangular(r, ri.T, overwrite_b=True)
        cov = dot(ri.T, ri)
        ret = zeros_like(fdjac)
        j1, j2 = meshgrid(jpvt[:nsing], jpvt[:nsing])
        ret[j1, j2] = cov
        return ret

    def error(self, x=None, eps=0):
        """ calculate the error of the fit parameters

        Give an estimation on the error of the fit result.

        Parameters
        ----------
        x : vector, length `N`, optional
            the fit result to calculate the error at. Take the last fit
            result if not given.
        eps : float, optional
            threshold below which a parameter is considered linear-dependent
            (ie not relevant for the fit result)

        Returns
        -------
        err : vector, length `N`
            the estimated error for each parameter
        """
        if x is None:
            x = self.params
        return sqrt(diagonal(self.covar(x, eps)))

def leastsq(func, x0, args=(), Dfun=None, full_output=0, col_deriv=0, **kwargs):
    """ Function equivalent to the class `Fit`

    This function is supposed to be a drop-in replacement for 
    `scipy.optimize.leastsq`. See the documentation there.
    """
    class MyFit(Fit):
        def func(self, x):
            return func(x, *args)

        if Dfun is not None:
            def jacobian(self, x, fvec, ret):
                if col_deriv:
                    return Dfun(x)
                else:
                    return Dfun(x).T
    fit = MyFit()
    try:
        ret = fit.fit(x0, *kwargs)
    except FitError as e:
        mesg = e.message
    else:
        mesg = [
            "an unknown (typically user) error has occurred",
            "both actual and predicted relative reductions "
                "in the sum of squares are at most `ftol`.",
            "relative error between two consecutive "
                "iterates is at most `xtol`",
            "the actual and predicted relative reductions "
                "in the sum of squares are at most `ftol`, "
                "and relative error between two consecutive "
                "iterates is at most `xtol`",
            "the cosine of the angle between function value and any "
                "column of the Jacobian is at most `gtol` in absolute value."
            ][fit.exitcode]
    if full_output:
        return ret, fit.covar(ret), {"nfev": fit.nfev,
                "fvec": fit.fvec}, mesg, fit.exitcode
    else:
        return ret, fit.exitcode
