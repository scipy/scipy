from collections import namedtuple
from functools import partial

import numpy as np
import scipy.sparse as sps
from ._numdiff import approx_derivative, group_columns
from ._hessian_update_strategy import HessianUpdateStrategy
from scipy.sparse.linalg import LinearOperator
from scipy._lib._array_api import array_namespace
from scipy._lib import array_api_extra as xpx
from scipy._lib._util import _ScalarFunctionWrapper


FD_METHODS = ('2-point', '3-point', 'cs')


class _ScalarGradWrapper:
    """
    Wrapper class for gradient calculation
    """
    def __init__(
            self,
            grad,
            fun=None,
            args=None,
            finite_diff_options=None,
    ):
        self.fun = fun
        self.grad = grad
        self.args = [] if args is None else args
        self.finite_diff_options = finite_diff_options
        self.ngev = 0
        # number of function evaluations consumed by finite difference
        self.nfev = 0

    def __call__(self, x, f0=None, **kwds):
        # Send a copy because the user may overwrite it.
        # The user of this class might want `x` to remain unchanged.
        if callable(self.grad):
            g = np.atleast_1d(self.grad(np.copy(x), *self.args))
        elif self.grad in FD_METHODS:
            g, dct = approx_derivative(
                self.fun,
                x,
                f0=f0,
                **self.finite_diff_options,
            )
            self.nfev += dct['nfev']

        self.ngev += 1
        return g


class _ScalarHessWrapper:
    """
    Wrapper class for hess calculation via finite differences
    """
    def __init__(
            self,
            hess,
            x0=None,
            grad=None,
            args=None,
            finite_diff_options=None,
    ):
        self.hess = hess
        self.grad = grad
        self.args = [] if args is None else args
        self.finite_diff_options = finite_diff_options
        # keep track of any finite difference function evaluations for grad
        self.ngev = 0
        self.nhev = 0
        self.H = None
        self._hess_func = None

        if callable(hess):
            self.H = hess(np.copy(x0), *args)
            self.nhev += 1

            if sps.issparse(self.H):
                self._hess_func = self._sparse_callable
                self.H = sps.csr_array(self.H)
            elif isinstance(self.H, LinearOperator):
                self._hess_func = self._linearoperator_callable
            else:
                # dense
                self._hess_func = self._dense_callable
                self.H = np.atleast_2d(np.asarray(self.H))
        elif hess in FD_METHODS:
                self._hess_func = self._fd_hess

    def __call__(self, x, f0=None, **kwds):
        return self._hess_func(np.copy(x), f0=f0)

    def _fd_hess(self, x, f0=None, **kwds):
        self.H, dct = approx_derivative(
            self.grad, x, f0=f0, **self.finite_diff_options
        )
        self.ngev += dct["nfev"]
        return self.H

    def _sparse_callable(self, x, **kwds):
        self.nhev += 1
        self.H = sps.csr_array(self.hess(x, *self.args))
        return self.H

    def _dense_callable(self, x, **kwds):
        self.nhev += 1
        self.H = np.atleast_2d(
            np.asarray(self.hess(x, *self.args))
        )
        return self.H

    def _linearoperator_callable(self, x, **kwds):
        self.nhev += 1
        self.H = self.hess(x, *self.args)
        return self.H


class ScalarFunction:
    """Scalar function and its derivatives.

    This class defines a scalar function F: R^n->R and methods for
    computing or approximating its first and second derivatives.

    Parameters
    ----------
    fun : callable
        evaluates the scalar function. Must be of the form ``fun(x, *args)``,
        where ``x`` is the argument in the form of a 1-D array and ``args`` is
        a tuple of any additional fixed parameters needed to completely specify
        the function. Should return a scalar.
    x0 : array-like
        Provides an initial set of variables for evaluating fun. Array of real
        elements of size (n,), where 'n' is the number of independent
        variables.
    args : tuple, optional
        Any additional fixed parameters needed to completely specify the scalar
        function.
    grad : {callable, '2-point', '3-point', 'cs'}
        Method for computing the gradient vector.
        If it is a callable, it should be a function that returns the gradient
        vector:

            ``grad(x, *args) -> array_like, shape (n,)``

        where ``x`` is an array with shape (n,) and ``args`` is a tuple with
        the fixed parameters.
        Alternatively, the keywords  {'2-point', '3-point', 'cs'} can be used
        to select a finite difference scheme for numerical estimation of the
        gradient with a relative step size. These finite difference schemes
        obey any specified `bounds`.
    hess : {callable, '2-point', '3-point', 'cs', HessianUpdateStrategy}
        Method for computing the Hessian matrix. If it is callable, it should
        return the  Hessian matrix:

            ``hess(x, *args) -> {LinearOperator, spmatrix, array}, (n, n)``

        where x is a (n,) ndarray and `args` is a tuple with the fixed
        parameters. Alternatively, the keywords {'2-point', '3-point', 'cs'}
        select a finite difference scheme for numerical estimation. Or, objects
        implementing `HessianUpdateStrategy` interface can be used to
        approximate the Hessian.
        Whenever the gradient is estimated via finite-differences, the Hessian
        cannot be estimated with options {'2-point', '3-point', 'cs'} and needs
        to be estimated using one of the quasi-Newton strategies.
    finite_diff_rel_step : None or array_like
        Relative step size to use. The absolute step size is computed as
        ``h = finite_diff_rel_step * sign(x0) * max(1, abs(x0))``, possibly
        adjusted to fit into the bounds. For ``method='3-point'`` the sign
        of `h` is ignored. If None then finite_diff_rel_step is selected
        automatically,
    finite_diff_bounds : tuple of array_like
        Lower and upper bounds on independent variables. Defaults to no bounds,
        (-np.inf, np.inf). Each bound must match the size of `x0` or be a
        scalar, in the latter case the bound will be the same for all
        variables. Use it to limit the range of function evaluation.
    epsilon : None or array_like, optional
        Absolute step size to use, possibly adjusted to fit into the bounds.
        For ``method='3-point'`` the sign of `epsilon` is ignored. By default
        relative steps are used, only if ``epsilon is not None`` are absolute
        steps used.
    workers : map-like callable, optional
        A map-like callable, such as `multiprocessing.Pool.map` for evaluating
        any numerical differentiation in parallel.
        This evaluation is carried out as ``workers(fun, iterable)``, or
        ``workers(grad, iterable)``, depending on what is being numerically
        differentiated.
        Alternatively, if `workers` is an int the task is subdivided into `workers`
        sections and the function evaluated in parallel
        (uses `multiprocessing.Pool <multiprocessing>`).
        Supply -1 to use all available CPU cores.
        It is recommended that a map-like be used instead of int, as repeated
        calls to `approx_derivative` will incur large overhead from setting up
        new processes.

        .. versionadded:: 1.16.0

    Notes
    -----
    This class implements a memoization logic. There are methods `fun`,
    `grad`, hess` and corresponding attributes `f`, `g` and `H`. The following
    things should be considered:

        1. Use only public methods `fun`, `grad` and `hess`.
        2. After one of the methods is called, the corresponding attribute
           will be set. However, a subsequent call with a different argument
           of *any* of the methods may overwrite the attribute.
    """
    def __init__(self, fun, x0, args, grad, hess, finite_diff_rel_step=None,
                 finite_diff_bounds=(-np.inf, np.inf), epsilon=None, workers=None):

        if not callable(grad) and grad not in FD_METHODS:
            raise ValueError(
                f"`grad` must be either callable or one of {FD_METHODS}."
            )

        if not (callable(hess) or hess in FD_METHODS
                or isinstance(hess, HessianUpdateStrategy)):
            raise ValueError(
                f"`hess` must be either callable, HessianUpdateStrategy"
                f" or one of {FD_METHODS}."
            )

        if grad in FD_METHODS and hess in FD_METHODS:
            raise ValueError("Whenever the gradient is estimated via "
                             "finite-differences, we require the Hessian "
                             "to be estimated using one of the "
                             "quasi-Newton strategies.")
        self.xp = xp = array_namespace(x0)
        _x = xpx.atleast_nd(xp.asarray(x0), ndim=1, xp=xp)
        _dtype = xp.float64
        if xp.isdtype(_x.dtype, "real floating"):
            _dtype = _x.dtype

        # original arguments
        self._wrapped_fun = _ScalarFunctionWrapper(fun, args)
        self._orig_fun = fun
        self._orig_grad = grad
        self._orig_hess = hess
        self._args = args

        # promotes to floating
        self.x = xp.astype(_x, _dtype)
        self.x_dtype = _dtype
        self.n = self.x.size
        self.f_updated = False
        self.g_updated = False
        self.H_updated = False

        self._lowest_x = None
        self._lowest_f = np.inf

        # normalize workers
        workers = workers or map

        finite_diff_options = {}
        if grad in FD_METHODS:
            finite_diff_options["method"] = grad
            finite_diff_options["rel_step"] = finite_diff_rel_step
            finite_diff_options["abs_step"] = epsilon
            finite_diff_options["bounds"] = finite_diff_bounds
            finite_diff_options["workers"] = workers
            finite_diff_options["full_output"] = True
        if hess in FD_METHODS:
            finite_diff_options["method"] = hess
            finite_diff_options["rel_step"] = finite_diff_rel_step
            finite_diff_options["abs_step"] = epsilon
            finite_diff_options["as_linear_operator"] = True
            finite_diff_options["workers"] = workers
            finite_diff_options["full_output"] = True

        # Initial function evaluation
        self._nfev = 0
        self._update_fun()

        # Initial gradient evaluation
        self._wrapped_grad = _ScalarGradWrapper(
            grad,
            fun=self._wrapped_fun,
            args=args,
            finite_diff_options=finite_diff_options,
        )
        self._update_grad()

        # Hessian evaluation
        if isinstance(hess, HessianUpdateStrategy):
            self.H = hess
            self.H.initialize(self.n, 'hess')
            self.H_updated = True
            self.x_prev = None
            self.g_prev = None
            _FakeCounter = namedtuple('_FakeCounter', ['ngev', 'nhev'])
            self._wrapped_hess = _FakeCounter(ngev=0, nhev=0)
        else:
            if callable(hess):
                self._wrapped_hess = _ScalarHessWrapper(
                    hess,
                    x0=x0,
                    args=args,
                    finite_diff_options=finite_diff_options
                )
                self.H = self._wrapped_hess.H
                self.H_updated = True
            elif hess in FD_METHODS:
                self._wrapped_hess = _ScalarHessWrapper(
                    hess,
                    x0=x0,
                    args=args,
                    grad=self._wrapped_grad,
                    finite_diff_options=finite_diff_options
                )
                self._update_grad()
                self.H = self._wrapped_hess(self.x, f0=self.g)
                self.H_updated = True

    @property
    def nfev(self):
        return self._nfev + self._wrapped_grad.nfev

    @property
    def ngev(self):
        return self._wrapped_grad.ngev  #+ self._wrapped_hess.ngev

    @property
    def nhev(self):
        return self._wrapped_hess.nhev

    def _update_x(self, x):
        if isinstance(self._orig_hess, HessianUpdateStrategy):
            self._update_grad()
            self.x_prev = self.x
            self.g_prev = self.g
            # ensure that self.x is a copy of x. Don't store a reference
            # otherwise the memoization doesn't work properly.

            _x = xpx.atleast_nd(self.xp.asarray(x), ndim=1, xp=self.xp)
            self.x = self.xp.astype(_x, self.x_dtype)
            self.f_updated = False
            self.g_updated = False
            self.H_updated = False
            self._update_hess()
        else:
            # ensure that self.x is a copy of x. Don't store a reference
            # otherwise the memoization doesn't work properly.
            _x = xpx.atleast_nd(self.xp.asarray(x), ndim=1, xp=self.xp)
            self.x = self.xp.astype(_x, self.x_dtype)
            self.f_updated = False
            self.g_updated = False
            self.H_updated = False

    def _update_fun(self):
        if not self.f_updated:
            fx = self._wrapped_fun(self.x)
            self._nfev += 1
            if fx < self._lowest_f:
                self._lowest_x = self.x
                self._lowest_f = fx

            self.f = fx
            self.f_updated = True

    def _update_grad(self):
        if not self.g_updated:
            if self._orig_grad in FD_METHODS:
                self._update_fun()
            self.g = self._wrapped_grad(self.x, f0=self.f)
            self.g_updated = True

    def _update_hess(self):
        if not self.H_updated:
            if self._orig_hess in FD_METHODS:
                self._update_grad()
                self.H = self._wrapped_hess(self.x, f0=self.g)
            elif isinstance(self._orig_hess, HessianUpdateStrategy):
                self._update_grad()
                self.H.update(self.x - self.x_prev, self.g - self.g_prev)
            else:       # should be callable(hess)
                self.H = self._wrapped_hess(self.x)

            self.H_updated = True

    def fun(self, x):
        if not np.array_equal(x, self.x):
            self._update_x(x)
        self._update_fun()
        return self.f

    def grad(self, x):
        if not np.array_equal(x, self.x):
            self._update_x(x)
        self._update_grad()
        return self.g

    def hess(self, x):
        if not np.array_equal(x, self.x):
            self._update_x(x)
        self._update_hess()
        return self.H

    def fun_and_grad(self, x):
        if not np.array_equal(x, self.x):
            self._update_x(x)
        self._update_fun()
        self._update_grad()
        return self.f, self.g


def _VectorFunWrapper(fun, x):
    return np.atleast_1d(fun(x))


class VectorFunction:
    """Vector function and its derivatives.

    This class defines a vector function F: R^n->R^m and methods for
    computing or approximating its first and second derivatives.

    Notes
    -----
    This class implements a memoization logic. There are methods `fun`,
    `jac`, hess` and corresponding attributes `f`, `J` and `H`. The following
    things should be considered:

        1. Use only public methods `fun`, `jac` and `hess`.
        2. After one of the methods is called, the corresponding attribute
           will be set. However, a subsequent call with a different argument
           of *any* of the methods may overwrite the attribute.
    """
    def __init__(self, fun, x0, jac, hess,
                 finite_diff_rel_step=None, finite_diff_jac_sparsity=None,
                 finite_diff_bounds=None, sparse_jacobian=None, workers=None):
        if not callable(jac) and jac not in FD_METHODS:
            raise ValueError(f"`jac` must be either callable or one of {FD_METHODS}.")

        if not (callable(hess) or hess in FD_METHODS
                or isinstance(hess, HessianUpdateStrategy)):
            raise ValueError("`hess` must be either callable,"
                             f"HessianUpdateStrategy or one of {FD_METHODS}.")

        if jac in FD_METHODS and hess in FD_METHODS:
            raise ValueError("Whenever the Jacobian is estimated via "
                             "finite-differences, we require the Hessian to "
                             "be estimated using one of the quasi-Newton "
                             "strategies.")

        self.xp = xp = array_namespace(x0)
        _x = xpx.atleast_nd(xp.asarray(x0), ndim=1, xp=xp)
        _dtype = xp.float64
        if xp.isdtype(_x.dtype, "real floating"):
            _dtype = _x.dtype

        # promotes to floating
        self.x = xp.astype(_x, _dtype)
        self.x_dtype = _dtype

        self.n = self.x.size
        self.nfev = 0
        self.njev = 0
        self.nhev = 0
        self.f_updated = False
        self.J_updated = False
        self.H_updated = False

        # normalize workers
        workers = workers or map

        finite_diff_options = {}
        if jac in FD_METHODS:
            finite_diff_options["method"] = jac
            finite_diff_options["rel_step"] = finite_diff_rel_step
            if finite_diff_jac_sparsity is not None:
                sparsity_groups = group_columns(finite_diff_jac_sparsity)
                finite_diff_options["sparsity"] = (finite_diff_jac_sparsity,
                                                   sparsity_groups)
            finite_diff_options["bounds"] = finite_diff_bounds
            finite_diff_options["workers"] = workers
            finite_diff_options["full_output"] = True
            self.x_diff = np.copy(self.x)
        if hess in FD_METHODS:
            finite_diff_options["method"] = hess
            finite_diff_options["rel_step"] = finite_diff_rel_step
            finite_diff_options["as_linear_operator"] = True
            # workers is not useful for evaluation of the LinearOperator
            # produced by approx_derivative. Only two/three function
            # evaluations are used, and the LinearOperator may persist
            # outside the scope that workers is valid in.
            self.x_diff = np.copy(self.x)
        if jac in FD_METHODS and hess in FD_METHODS:
            raise ValueError("Whenever the Jacobian is estimated via "
                             "finite-differences, we require the Hessian to "
                             "be estimated using one of the quasi-Newton "
                             "strategies.")

        fun_wrapped = partial(_VectorFunWrapper, fun)

        def update_fun():
            self.nfev += 1
            self.f = fun_wrapped(self.x)

        self._update_fun_impl = update_fun
        update_fun()

        self.v = np.zeros_like(self.f)
        self.m = self.v.size

        # Jacobian Evaluation
        if callable(jac):
            self.J = jac(self.x)
            self.J_updated = True
            self.njev += 1

            if (sparse_jacobian or
                    sparse_jacobian is None and sps.issparse(self.J)):
                def jac_wrapped(x):
                    self.njev += 1
                    return sps.csr_array(jac(x))
                self.J = sps.csr_array(self.J)
                self.sparse_jacobian = True

            elif sps.issparse(self.J):
                def jac_wrapped(x):
                    self.njev += 1
                    return jac(x).toarray()
                self.J = self.J.toarray()
                self.sparse_jacobian = False

            else:
                def jac_wrapped(x):
                    self.njev += 1
                    return np.atleast_2d(jac(x))
                self.J = np.atleast_2d(self.J)
                self.sparse_jacobian = False

            def update_jac():
                self.J = jac_wrapped(self.x)

        elif jac in FD_METHODS:
            self.J, dct = approx_derivative(fun_wrapped, self.x, f0=self.f,
                                            **finite_diff_options)
            self.J_updated = True
            self.nfev += dct['nfev']

            if (sparse_jacobian or
                    sparse_jacobian is None and sps.issparse(self.J)):
                def update_jac():
                    self._update_fun()
                    self.J, dct = sps.csr_array(
                        approx_derivative(fun_wrapped, self.x, f0=self.f,
                                          **finite_diff_options))
                    self.nfev += dct['nfev']
                self.J = sps.csr_array(self.J)
                self.sparse_jacobian = True

            elif sps.issparse(self.J):
                def update_jac():
                    self._update_fun()
                    self.J, dct = approx_derivative(fun_wrapped, self.x, f0=self.f,
                                                    **finite_diff_options).toarray()
                    self.nfev += dct['nfev']
                self.J = self.J.toarray()
                self.sparse_jacobian = False

            else:
                def update_jac():
                    self._update_fun()
                    J, dct = approx_derivative(fun_wrapped, self.x, f0=self.f,
                                               **finite_diff_options)
                    self.J = np.atleast_2d(J)
                    self.nfev += dct['nfev']
                self.J = np.atleast_2d(self.J)
                self.sparse_jacobian = False

        self._update_jac_impl = update_jac

        # Define Hessian
        if callable(hess):
            self.H = hess(self.x, self.v)
            self.H_updated = True
            self.nhev += 1

            if sps.issparse(self.H):
                def hess_wrapped(x, v):
                    self.nhev += 1
                    return sps.csr_array(hess(x, v))
                self.H = sps.csr_array(self.H)

            elif isinstance(self.H, LinearOperator):
                def hess_wrapped(x, v):
                    self.nhev += 1
                    return hess(x, v)

            else:
                def hess_wrapped(x, v):
                    self.nhev += 1
                    return np.atleast_2d(np.asarray(hess(x, v)))
                self.H = np.atleast_2d(np.asarray(self.H))

            def update_hess():
                self.H = hess_wrapped(self.x, self.v)
        elif hess in FD_METHODS:
            def jac_dot_v(x, v):
                return jac_wrapped(x).T.dot(v)

            def update_hess():
                self._update_jac()
                self.H = approx_derivative(jac_dot_v, self.x,
                                           f0=self.J.T.dot(self.v),
                                           args=(self.v,),
                                           **finite_diff_options)

            update_hess()
            self.H_updated = True
        elif isinstance(hess, HessianUpdateStrategy):
            self.H = hess
            self.H.initialize(self.n, 'hess')
            self.H_updated = True
            self.x_prev = None
            self.J_prev = None

            def update_hess():
                self._update_jac()
                # When v is updated before x was updated, then x_prev and
                # J_prev are None and we need this check.
                if self.x_prev is not None and self.J_prev is not None:
                    delta_x = self.x - self.x_prev
                    delta_g = self.J.T.dot(self.v) - self.J_prev.T.dot(self.v)
                    self.H.update(delta_x, delta_g)

        self._update_hess_impl = update_hess

        if isinstance(hess, HessianUpdateStrategy):
            def update_x(x):
                self._update_jac()
                self.x_prev = self.x
                self.J_prev = self.J
                _x = xpx.atleast_nd(self.xp.asarray(x), ndim=1, xp=self.xp)
                self.x = self.xp.astype(_x, self.x_dtype)
                self.f_updated = False
                self.J_updated = False
                self.H_updated = False
                self._update_hess()
        else:
            def update_x(x):
                _x = xpx.atleast_nd(self.xp.asarray(x), ndim=1, xp=self.xp)
                self.x = self.xp.astype(_x, self.x_dtype)
                self.f_updated = False
                self.J_updated = False
                self.H_updated = False

        self._update_x_impl = update_x

    def _update_v(self, v):
        if not np.array_equal(v, self.v):
            self.v = v
            self.H_updated = False

    def _update_x(self, x):
        if not np.array_equal(x, self.x):
            self._update_x_impl(x)

    def _update_fun(self):
        if not self.f_updated:
            self._update_fun_impl()
            self.f_updated = True

    def _update_jac(self):
        if not self.J_updated:
            self._update_jac_impl()
            self.J_updated = True

    def _update_hess(self):
        if not self.H_updated:
            self._update_hess_impl()
            self.H_updated = True

    def fun(self, x):
        self._update_x(x)
        self._update_fun()
        return self.f

    def jac(self, x):
        self._update_x(x)
        self._update_jac()
        return self.J

    def hess(self, x, v):
        # v should be updated before x.
        self._update_v(v)
        self._update_x(x)
        self._update_hess()
        return self.H


class LinearVectorFunction:
    """Linear vector function and its derivatives.

    Defines a linear function F = A x, where x is N-D vector and
    A is m-by-n matrix. The Jacobian is constant and equals to A. The Hessian
    is identically zero and it is returned as a csr matrix.
    """
    def __init__(self, A, x0, sparse_jacobian):
        if sparse_jacobian or sparse_jacobian is None and sps.issparse(A):
            self.J = sps.csr_array(A)
            self.sparse_jacobian = True
        elif sps.issparse(A):
            self.J = A.toarray()
            self.sparse_jacobian = False
        else:
            # np.asarray makes sure A is ndarray and not matrix
            self.J = np.atleast_2d(np.asarray(A))
            self.sparse_jacobian = False

        self.m, self.n = self.J.shape

        self.xp = xp = array_namespace(x0)
        _x = xpx.atleast_nd(xp.asarray(x0), ndim=1, xp=xp)
        _dtype = xp.float64
        if xp.isdtype(_x.dtype, "real floating"):
            _dtype = _x.dtype

        # promotes to floating
        self.x = xp.astype(_x, _dtype)
        self.x_dtype = _dtype

        self.f = self.J.dot(self.x)
        self.f_updated = True

        self.v = np.zeros(self.m, dtype=float)
        self.H = sps.csr_array((self.n, self.n))

    def _update_x(self, x):
        if not np.array_equal(x, self.x):
            _x = xpx.atleast_nd(self.xp.asarray(x), ndim=1, xp=self.xp)
            self.x = self.xp.astype(_x, self.x_dtype)
            self.f_updated = False

    def fun(self, x):
        self._update_x(x)
        if not self.f_updated:
            self.f = self.J.dot(x)
            self.f_updated = True
        return self.f

    def jac(self, x):
        self._update_x(x)
        return self.J

    def hess(self, x, v):
        self._update_x(x)
        self.v = v
        return self.H


class IdentityVectorFunction(LinearVectorFunction):
    """Identity vector function and its derivatives.

    The Jacobian is the identity matrix, returned as a dense array when
    `sparse_jacobian=False` and as a csr matrix otherwise. The Hessian is
    identically zero and it is returned as a csr matrix.
    """
    def __init__(self, x0, sparse_jacobian):
        n = len(x0)
        if sparse_jacobian or sparse_jacobian is None:
            A = sps.eye_array(n, format='csr')
            sparse_jacobian = True
        else:
            A = np.eye(n)
            sparse_jacobian = False
        super().__init__(A, x0, sparse_jacobian)
