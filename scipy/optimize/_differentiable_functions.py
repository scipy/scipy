from __future__ import division, print_function, absolute_import
import numpy as np
import scipy.sparse as sps
from ._numdiff import approx_derivative, group_columns
from ._hessian_update_strategy import HessianUpdateStrategy
from scipy.sparse.linalg import LinearOperator
from copy import deepcopy


FD_METHODS = ('2-point', '3-point', 'cs')


class ScalarFunction(object):
    """Scalar function and its derivatives.

    This class defines a scalar function F: R^n->R and methods for
    computing or approximating its first and second derivatives.
    """

    def __init__(self, fun, x0, args, grad, hess, finite_diff_rel_step,
                 finite_diff_bounds):
        self.x = np.atleast_1d(x0).astype(float)
        self.n = self.x.size
        self.nfev, self.ngev, self.nhev = 0, 0, 0
        finite_diff_options = {}
        if grad in FD_METHODS:
            finite_diff_options["method"] = grad
            finite_diff_options["rel_step"] = finite_diff_rel_step
            finite_diff_options["bounds"] = finite_diff_bounds
            self.x_diff = np.copy(self.x)
        if hess in FD_METHODS:
            finite_diff_options["method"] = hess
            finite_diff_options["rel_step"] = finite_diff_rel_step
            finite_diff_options["as_linear_operator"] = True
            self.x_diff = np.copy(self.x)
        if grad in FD_METHODS and hess in FD_METHODS:
            raise ValueError("Whenever the gradient is estimated via "
                             "finite-differences, we require the Hessian "
                             "to be estimated using one of the "
                             "quasi-Newton strategies.")

        # Function evaluation
        def fun_wrapped(x):
            self.nfev += 1
            return fun(x, *args)
        self.f = fun(self.x, *args)
        self.f_updated = True
        self.nfev += 1

        def evaluate_fun():
            self.f = fun_wrapped(self.x)
            self.f_updated = True
            return self.f
        self.evaluate_fun = evaluate_fun

        # Gradient evaluation
        if callable(grad):
            def grad_wrapped(x):
                self.ngev += 1
                return np.atleast_1d(grad(x, *args))
            self.g = np.atleast_1d(grad(self.x, *args))
            self.g_updated = True
            self.ngev += 1

            def evaluate_grad():
                self.g = grad_wrapped(self.x)
                self.g_updated = True
                return self.g
        elif grad in FD_METHODS:
            self.g = approx_derivative(fun_wrapped, self.x, f0=self.f,
                                       **finite_diff_options)
            self.g_updated = True

            def evaluate_grad():
                if not self.f_updated:
                    self.evaluate_fun()
                self.g = approx_derivative(fun_wrapped, self.x, f0=self.f,
                                           **finite_diff_options)
                self.g_updated = True
                return self.g
        self.evaluate_grad = evaluate_grad

        # Hessian Evaluation
        if callable(hess):
            self.H = hess(x0, *args)
            self.H_updated = True
            self.nhev += 1

            if sps.issparse(self.H):
                def hess_wrapped(x):
                    self.nhev += 1
                    return sps.csr_matrix(hess(x, *args))
                self.H = sps.csr_matrix(self.H)

            elif isinstance(self.H, LinearOperator):
                def hess_wrapped(x):
                    self.nhev += 1
                    return hess(x, *args)

            else:
                def hess_wrapped(x):
                    self.nhev += 1
                    return np.atleast_2d(np.asarray(hess(x, *args)))
                self.H = np.atleast_2d(np.asarray(self.H))

            def evaluate_hess():
                self.H = hess_wrapped(self.x)
                self.H_updated = True
                return self.H

        elif hess in FD_METHODS:
            self.H = approx_derivative(grad_wrapped, self.x,
                                       f0=self.g, **finite_diff_options)
            self.H_updated = True

            def evaluate_hess():
                if not self.g_updated:
                    self.evaluate_grad()
                self.H = approx_derivative(grad_wrapped, self.x,
                                           f0=self.g, **finite_diff_options)
                self.H_updated = True
                return self.H
        elif isinstance(hess, HessianUpdateStrategy):
            self.H = hess
            self.H.initialize(len(x0), 'hess')
            self.H_updated = True

            def evaluate_hess():
                return self.H
        self.evaluate_hess = evaluate_hess

        # Update current point
        if isinstance(hess, HessianUpdateStrategy):
            def update_x(x):
                x_prev = deepcopy(self.x)
                g_prev = deepcopy(self.g)
                self.x = x
                self.f_updated = False
                self.evaluate_grad()
                self.g_updated = True
                delta_x = self.x - x_prev
                delta_grad = self.g - g_prev
                self.H.update(delta_x, delta_grad)
                self.H_updated = True

        else:
            def update_x(x):
                self.x = x
                self.f_updated = False
                self.g_updated = False
                self.H_updated = False
        self.update_x = update_x

    def fun(self, x):
        if not np.array_equal(x, self.x):
            self.update_x(x)
        if not self.f_updated:
            self.evaluate_fun()
        return self.f

    def grad(self, x):
        if not np.array_equal(x, self.x):
            self.update_x(x)
        if not self.g_updated:
            self.evaluate_grad()
        return self.g

    def hess(self, x):
        if not np.array_equal(x, self.x):
            self.update_x(x)
        if not self.H_updated:
            self.evaluate_hess()
        return self.H


class VectorFunction(object):
    """Vector function and its derivatives.

    This class defines a vector function F: R^n->R^m and methods for
    computing or approximating its first and second derivatives.
    """

    def __init__(self, fun, x0, jac, hess,
                 finite_diff_rel_step, finite_diff_jac_sparsity,
                 finite_diff_bounds, sparse_jacobian):
        self.x = np.atleast_1d(x0).astype(float)
        self.n = self.x.size
        self.nfev, self.njev, self.nhev = 0, 0, 0
        finite_diff_options = {}
        if jac in FD_METHODS:
            finite_diff_options["method"] = jac
            finite_diff_options["rel_step"] = finite_diff_rel_step
            if finite_diff_jac_sparsity is not None:
                sparsity_groups = group_columns(finite_diff_jac_sparsity)
                finite_diff_options["sparsity"] = (finite_diff_jac_sparsity,
                                                   sparsity_groups)
            finite_diff_options["bounds"] = finite_diff_bounds
            self.x_diff = np.copy(self.x)
        if hess in FD_METHODS:
            finite_diff_options["method"] = hess
            finite_diff_options["rel_step"] = finite_diff_rel_step
            finite_diff_options["as_linear_operator"] = True
            self.x_diff = np.copy(self.x)
        if jac in FD_METHODS and hess in FD_METHODS:
            raise ValueError("Whenever the Jacobian is estimated via "
                             "finite-differences, we require the Hessian to "
                             "be estimated using one of the quasi-Newton "
                             "strategies.")

        # Function evaluation
        def fun_wrapped(x):
            self.nfev += 1
            return np.atleast_1d(fun(x))
        self.f = np.atleast_1d(fun(self.x))
        self.f_updated = True
        self.nfev += 1
        self.v = np.zeros_like(self.f)
        self.m = self.v.size

        def evaluate_fun():
            self.f = fun_wrapped(self.x)
            self.f_updated = True
            return self.f
        self.evaluate_fun = evaluate_fun

        # Jacobian Evaluation
        if callable(jac):
            self.J = jac(self.x)
            self.J_updated = True
            self.njev += 1

            if sparse_jacobian or \
               (sparse_jacobian is None and sps.issparse(self.J)):
                def jac_wrapped(x):
                    self.njev += 1
                    return sps.csr_matrix(jac(x))
                self.J = sps.csr_matrix(self.J)
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

            def evaluate_jac():
                self.J = jac_wrapped(self.x)
                self.J_updated = True
                return self.J
        elif jac in FD_METHODS:
            self.J = approx_derivative(fun_wrapped, self.x, f0=self.f,
                                       **finite_diff_options)
            self.J_updated = True

            if sparse_jacobian or \
               (sparse_jacobian is None and sps.issparse(self.J)):
                def evaluate_jac():
                    if not self.f_updated:
                        self.evaluate_fun()
                    self.J = sps.csr_matrix(
                        approx_derivative(fun_wrapped, self.x, f0=self.f,
                                          **finite_diff_options))
                    self.J_updated = True
                    return self.J
                self.J = sps.csr_matrix(self.J)
                self.sparse_jacobian = True

            elif sps.issparse(self.J):
                def evaluate_jac():
                    if not self.f_updated:
                        self.evaluate_fun()
                    self.J = approx_derivative(fun_wrapped, self.x, f0=self.f,
                                               **finite_diff_options).toarray()
                    self.J_updated = True
                    return self.J
                self.J = self.J.toarray()
                self.sparse_jacobian = False

            else:
                def evaluate_jac():
                    if not self.f_updated:
                        self.evaluate_fun()
                    self.J = np.atleast_2d(
                        approx_derivative(fun_wrapped, self.x, f0=self.f,
                                          **finite_diff_options))
                    self.J_updated = True
                    return self.J
                self.J = np.atleast_2d(self.J)
                self.sparse_jacobian = False
        self.evaluate_jac = evaluate_jac

        # Define Hessian
        if callable(hess):
            self.H = hess(self.x, self.v)
            self.H_updated = True
            self.nhev += 1

            if sps.issparse(self.H):
                def hess_wrapped(x, v):
                    self.nhev += 1
                    return sps.csr_matrix(hess(x))
                self.H = sps.csr_matrix(self.H)

            elif isinstance(self.H, LinearOperator):
                def hess_wrapped(x, v):
                    self.nhev += 1
                    return hess(x, v)

            else:
                def hess_wrapped(x, v):
                    self.nhev += 1
                    return np.atleast_2d(np.asarray(hess(x, v)))
                self.H = np.atleast_2d(np.asarray(self.H))

            def evaluate_hess():
                self.H = hess_wrapped(self.x, self.v)
                self.H_updated = True
                return self.H

        elif hess in FD_METHODS:
            def jac_dot_v(x, v):
                return jac_wrapped(x).T.dot(v)
            self.H = approx_derivative(jac_dot_v, self.x,
                                       f0=self.J.T.dot(self.v),
                                       args=(self.v,),
                                       **finite_diff_options)
            self.H_updated = True

            def evaluate_hess():
                if not self.J_updated:
                    self.evaluate_jac()
                self.H = approx_derivative(jac_dot_v, self.x,
                                           f0=self.J.T.dot(self.v),
                                           args=(self.v,),
                                           **finite_diff_options)
                self.H_updated = True
                return self.H
        elif isinstance(hess, HessianUpdateStrategy):
            self.H = hess
            self.H.initialize(len(x0), 'hess')
            self.H_updated = True

            def evaluate_hess():
                return self.H
        self.evaluate_hess = evaluate_hess

        if isinstance(hess, HessianUpdateStrategy):
            def update_x(x):
                x_prev = deepcopy(self.x)
                J_prev = deepcopy(self.J)
                self.x = x
                self.f_updated = False
                self.evaluate_jac()
                self.J_updated = True
                delta_x = self.x - x_prev
                delta_grad = self.J.T.dot(self.v) - J_prev.T.dot(self.v)
                self.H.update(delta_x, delta_grad)
                self.H_updated = True

        else:
            def update_x(x):
                self.x = x
                self.f_updated = False
                self.J_updated = False
                self.H_updated = False

        def update_v(v):
            self.v = v
            self.H_updated = False
        self.update_x = update_x
        self.update_v = update_v

    def fun(self, x):
        if not np.array_equal(x, self.x):
            self.update_x(x)
        if not self.f_updated:
            self.evaluate_fun()
        return self.f

    def jac(self, x):
        if not np.array_equal(x, self.x):
            self.update_x(x)
        if not self.J_updated:
            self.evaluate_jac()
        return self.J

    def hess(self, x, v):
        # v should ALWAYS be updated before x
        # because of HessianUpdateStategy.
        if not np.array_equal(v, self.v):
            self.update_v(v)
        if not np.array_equal(x, self.x):
            self.update_x(x)
        if not self.H_updated:
            self.evaluate_hess()
        return self.H


class LinearVectorFunction(object):
    """Linear vector function and its derivatives.

    Defines a linear function F = A x, where x is n-dimensional vector and
    A is m-by-n matrix. The Jacobian is constant and equals to A. The Hessian
    is identically zero and it is returned as a csr matrix.
    """
    def __init__(self, A, x0, sparse_jacobian):
        if sparse_jacobian or sparse_jacobian is None and sps.issparse(A):
            self.J = sps.csr_matrix(A)
            self.sparse_jacobian = True
        elif sps.issparse(A):
            self.J = A.toarray()
            self.sparse_jacobian = False
        else:
            self.J = np.atleast_2d(A)
            self.sparse_jacobian = False

        self.m, self.n = self.J.shape

        self.x = np.atleast_1d(x0).astype(float)
        self.f = self.J.dot(self.x)
        self.f_updated = True

        self.v = np.zeros(self.m, dtype=float)
        self.H = sps.csr_matrix((self.n, self.n))

    def change_x(self, x):
        self.x = x
        self.f_updated = False

    def fun(self, x):
        if not np.array_equal(x, self.x):
            self.change_x(x)
        if not self.f_updated:
            self.f = self.J.dot(x)
            self.f_updated = True
        return self.f

    def jac(self, x):
        if not np.array_equal(x, self.x):
            self.change_x(x)
        return self.J

    def hess(self, x, v):
        if not np.array_equal(x, self.x):
            self.change_x(x)
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
            A = sps.eye(n, format='csr')
            sparse_jacobian = True
        else:
            A = np.eye(n)
            sparse_jacobian = False
        super(IdentityVectorFunction, self).__init__(A, x0, sparse_jacobian)
