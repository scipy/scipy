# Author: Matt Haberland

from __future__ import division, absolute_import, print_function

cimport cython
import numpy as np
cimport numpy as np
from scipy.linalg import (solve, lu_solve, lu_factor, solve_triangular,
                          LinAlgError)
from scipy.linalg.cython_blas cimport daxpy, dswap
try:
    from time import process_time as timer
except ImportError:
    from time import clock as timer

__all__ = ['LU', 'BGLU']

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void swap_rows(self, double[:, ::1] H, int i):
    """
    Swaps row i of H with next row; represents matrix product by PI_i
    matrix described after matrix 5.10
    """
    # Python
    # H[[i, i+1]] = H[[i+1, i]]
    # Cython, using BLAS
    cdef int n = H.shape[1]-i
    cdef int one = 1
    dswap(&n, &H[i, i], &one, &H[i+1, i], &one)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True) # not really important
cdef double row_subtract(self, double[:, ::1] H, int i):
    """
    Zeros first nonzero element of row i+1 of H by subtracting appropriate
    multiple of row i; represents matrix product by matrix 5.10. Returns
    factor g for storage.
    """
    cdef double g = H[i+1, i]/H[i,i]

    # Python
    # H[i+1, i:] -= g*H[i, i:]
    # Cython, using BLAS
    cdef int n = H.shape[1]-i
    cdef double ng = -g
    cdef int one = 1
    daxpy(&n, &ng, &H[i, i], &one, &H[i+1, i], &one)

    return g

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void hess_lu(self, double[:, ::1] H, int i, double[:,::1] ops):
    """
    Converts Hessenberg matrix H with first nonzero off-diagonal in
    column i to upper triangular, recording elementary row operations.
    That is, performs and records operations in Equation 5.9.
    """
    cdef int m = H.shape[1]
    cdef int j, k
    cdef double piv1, piv2
    cdef double g
    cdef double swap

    for k in range(i, m-1):
        j = k-i
        piv1, piv2 = abs(H[k, k]), abs(H[k+1, k]) # np.abs(H[k:k+2,k])
        swap = float(piv1 < piv2)
        # swap rows to ensure |g| <= 1
        if swap:
            swap_rows(self, H, k)
        g = row_subtract(self, H, k)

        ops[j, 0] = swap
        ops[j, 1] = g

@cython.boundscheck(False)
@cython.wraparound(False)
cdef void perform_ops(self, double[::1] y, double[:,::1] ops, bint rev = False):
    """
    Replays operations needed to convert Hessenberg matrix into upper
    triangular form on a vector y. Equivalent to matrix multlication by
    inverse of matrix 5.12.
    """
    cdef int i, j, k, m
    cdef double g, swap

    m = y.shape[0]
    i = m - ops.shape[0] - 1
    if not rev:
        for k in range(i, m-1):
            j = k - i
            swap = ops[j,0]
            g = ops[j, 1]
            if swap:
                y[k], y[k+1] = y[k+1], y[k]

            y[k+1] -= g*y[k]
    else:
        for k in range(m-2, i-1, -1):
            j = k -i
            swap = ops[j,0]
            g = ops[j, 1]
            y[k] -= g*y[k+1]
            if swap:
                y[k], y[k+1] = y[k+1], y[k]

def _consider_refactor(method):
    """
    This decorator records the time spent in the major BGLU
    routines - refactor, update, and solve - in order to
    calculate the average time required to solve a system.
    It also forces PLU factorization of the basis matrix from
    scratch to minimize the average solve time and to
    accumulation of roundoff error.

    Immediately after PLU factorization, the average solve time
    will be rather high because PLU factorization is slow. For
    some number of factor updates, the average solve time is
    expected to decrease because the updates and solves are fast.
    However, updates increase the compexity of the factorization,
    so solve times are expected to increase with each update.
    When the average solve time stops decreasing and begins
    increasing, we perform PLU factorization from scratch rather
    than updating. PLU factorization is also performed after the
    maximum permitted number of updates is reached to prevent
    further accumulation of roundoff error.
    """
    def f(self, *args, **kwargs):

        refactor_now = False
        out = None

        if method.__name__ == "update":
            # this will get zeroed if we refactor
            self.updates += 1

            # if average solve time is increasing, then it would
            # be faster to refactor from scratch
            slowing_down = (self.average_solve_times[1] >
                            self.average_solve_times[0])

            # if update limit is reached, we should refactor to
            # limit error buildup
            too_many_updates = self.updates >= self.max_updates

            if self.mast:
                refactor_now = (slowing_down or too_many_updates)
            else:
                refactor_now = too_many_updates

            if refactor_now:
                # update basis indices and factor from scratch
                self.update_basis(*args, **kwargs)
                out = self.refactor()  # time will be recorded

        # If refactor_now is True, then self.refactor() is called
        # We don't want to call method = self.update again here
        if not refactor_now:
            # record the time it took to call the method
            t0 = timer()
            out = method(self, *args, **kwargs)
            if isinstance(out, np.ndarray) and np.any(np.isnan(out)):
                raise LinAlgError("Nans in output")
            t1 = timer()
            self.bglu_time += (t1-t0)

        # calculate average solve time,
        # considering all significant method calls
        if method.__name__ == "solve":
            self.solves += 1
            avg = self.bglu_time/self.solves
            self.average_solve_times = [
                self.average_solve_times[1], avg]

        return out
    return f


cdef class LU(object):
    """
    Represents PLU factorization of a basis matrix with naive rank-one updates
    """

    cdef public np.ndarray A
    cdef public np.ndarray b
    cdef public np.ndarray B
    cdef public int m
    cdef public int n

    def __init__(self, A, b):
        """ Given matrix A and basis indices b, form basis matrix B """
        self.A = A
        self.b = b
        self.B = A[:, b]
        self.m, self.n = A.shape

    def update(self, i, j):
        """ Rank-one update to basis and basis matrix """
        self.b[i:self.m-1] = self.b[i+1:self.m]
        self.b[-1] = j
#        self.b[i] = j
        self.B = self.A[:, self.b]

    def solve(self, q, transposed = False):
        """
        Solve B @ v = q
        """
        v = solve(self.B, q, transposed=transposed)
        return v


cdef class BGLU(LU):
    """
    Represents PLU factorization with Golub rank-one updates from
    Bartels, Richard H. "A stabilization of the simplex method."
    Numerische Mathematik 16.5 (1971): 414-434.
    """

    cdef public tuple plu
    cdef public np.ndarray L
    cdef public np.ndarray U
    cdef public np.ndarray pi
    cdef public np.ndarray pit
    cdef public list ops_list
    cdef public double bglu_time
    cdef public int solves
    cdef public int updates
    cdef public int max_updates
    cdef public list average_solve_times
    cdef public bint mast

    def __init__(self, A, b, max_updates=10, mast=False):
        """
        Given matrix A and basis indices b, perform PLU factorization of
        basis matrix B
        """
        self.A = A
        self.b = b
        self.m, self.n = A.shape
        self.max_updates = max_updates  # maximum updates between refactor
        self.refactor()
        self.mast = mast

    @_consider_refactor
    def refactor(self):
        # Factor as in Equation 5.1
        self.B = self.A[:, self.b]  # get basis matrix
        self.plu = lu_factor(self.B)  # lu_factor tested faster than lu
        self.L = self.plu[0]  # L and U stored in the same matrix
        self.U = self.plu[0].copy()  # need to modify without changing L
        # indexing with self.pi is equivalent to PI matrix product
        self.pi = self.perform_perm(self.plu[1])  # permutation indices
        self.pit = np.zeros(self.m, dtype=int)  # permutation transposed
        self.pit[self.pi] = np.arange(self.m)
        self.ops_list = []  # elementary row operations in order

        self.bglu_time = 0  # cumulative time spent updating and solving
        self.solves = 0     # number of solves since refactoring
        self.updates = 0    # number of updates since refactoring
        self.average_solve_times = [np.inf, np.inf]  # current and last average solve time

    # ideally should time this, too, but I also want to call this
    # method in update below, which would double-count the time.
    def update_basis(self, i, j):
        self.b[i:self.m-1] = self.b[i+1:self.m]  # eliminate i from basis
        self.b[-1] = j  # add j to end of basis

    @_consider_refactor
    def update(self, i, j):
        """ Perform rank-one update to basis and factorization """
        self.update_basis(i, j)

        # calculate last column of Hessenberg matrix
        # TODO: share this calculation with simplex method
        pla = self.A[self.pi, j]
        um = solve_triangular(self.L, pla, lower = True,
                              check_finite=False, unit_diagonal=True)
        for ops in self.ops_list:
            perform_ops(self, um, ops) # modifies um in place

        # form Hessenberg matrix
        H = self.U
        H[:, i:self.m-1] = self.U[:, i+1:self.m]  # eliminate column i
        H[:, -1] = um  # add column corresponding with j

        # convert H to upper triangular, recording elementary row operations
        ops = np.zeros((self.m-1-i, 2))
        hess_lu(self, H, i, ops) # hess_lu modifies ops in place
        self.ops_list.append(ops)

        self.U = H

    @_consider_refactor
    def solve(self, q, transposed = False):
        """
        Solve B @ v = q efficiently using factorization
        """
        if not self.ops_list:
            # before any updates, solve according to Equation 5.2
            v = lu_solve(self.plu, q, trans=transposed)
        else:
            if not transposed:
                q = q[self.pi]  # paper skips this by making
                                # "inessential assumption" of no permutation

                # Equation 5.16
                t = solve_triangular(self.L, q, lower = True,
                                     check_finite=False, unit_diagonal=True)

                # Equation 5.17
                temp = t
                for ops in self.ops_list:
                    perform_ops(self, temp, ops) # modifies temp in place
                w = temp

                # Equation 5.18
                # Faster to use U.T and set trans=True due to array order
                v = solve_triangular(self.U.T, w, lower=True,
                                     trans=True, check_finite=False)

            else: # do everything transposed and in reverse order
                t = solve_triangular(self.U.T, q, lower=True,
                                     trans=False, check_finite=False)
                temp = t
                for ops in reversed(self.ops_list):
                    perform_ops(self, temp, ops, rev = True) # mod in place
                w = temp
                v = solve_triangular(self.L, w, lower = True, trans=True,
                                     check_finite=False, unit_diagonal=True)
                v = v[self.pit]

        return v

    def perform_perm(self, p):
        """
        Perform individual row swaps defined in p returned by factor_lu to
        generate final permutation indices pi
        """
        pi = np.arange(len(p))
        for i, row in enumerate(p):
            pi[i], pi[row] = pi[row], pi[i]
        return pi
