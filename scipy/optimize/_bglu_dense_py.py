#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  4 15:57:00 2018

@author: matthaberland
"""

from __future__ import absolute_import

import numpy as np
from scipy.linalg import (solve, lu_solve, lu_factor, solve_triangular,
                          LinAlgError)
import time

__all__ = ['LU', 'BGLU']


def consider_refactor(method):
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
            t0 = time.clock()
            out = method(self, *args, **kwargs)
            if isinstance(out, np.ndarray) and np.any(np.isnan(out)):
                raise LinAlgError("Nans in output")
            t1 = time.clock()
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


class LU(object):
    """
    Represents PLU factorization of a basis matrix with naive rank-one updates
    """

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

    def solve(self, q, transposed=False):
        """
        Solve B @ v = q
        """
        v = solve(self.B, q, transposed=transposed)
        return v


class BGLU(LU):
    """
    Represents PLU factorization with Golub rank-one updates from
    Bartels, Richard H. "A stabilization of the simplex method."
    Numerische Mathematik 16.5 (1971): 414-434.
    """

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

    @consider_refactor
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

    @consider_refactor
    def update(self, i, j):
        """ Perform rank-one update to basis and factorization """
        self.update_basis(i, j)

        # calculate last column of Hessenberg matrix
        # FIXME: share this calculation with simplex method
        pla = self.A[self.pi, j]
        um = solve_triangular(self.L, pla, lower=True,
                              check_finite=False, unit_diagonal=True)
        for ops in self.ops_list:
            um = self.perform_ops(um, ops)

        # form Hessenberg matrix
        H = self.U
        H[:, i:self.m-1] = self.U[:, i+1:self.m]  # eliminate column i
        H[:, -1] = um  # add column corresponding with j

        # convert H to upper triangular, recording elementary row operations
        self.ops_list.append(self.hess_lu(H, i))

        self.U = H

    @consider_refactor
    def solve(self, q, transposed=False):
        """
        Solve B @ v = q efficiently using factorization
        """
        if not self.ops_list:
            # before any updates, solve according to Equation 5.2
            v = lu_solve(self.plu, q, trans=transposed)
        else:
            if not transposed:
                # paper skips this by making "inessential assumption" of
                # no permutation
                q = q[self.pi]

                # Equation 5.16
                t = solve_triangular(self.L, q, lower=True,
                                     check_finite=False, unit_diagonal=True)

                # Equation 5.17
                temp = t
                for ops in self.ops_list:
                    temp = self.perform_ops(temp, ops)
                w = temp

                # Equation 5.18
                # For whatever reason, faster to use U.T and set trans=True
                v = solve_triangular(self.U.T, w, lower=True,
                                     trans=True, check_finite=False)

            else:  # do everything transposed and in reverse order
                t = solve_triangular(self.U.T, q, lower=True,
                                     trans=False, check_finite=False)
                temp = t
                for ops in reversed(self.ops_list):
                    temp = self.perform_ops(temp, ops, reverse=True)
                w = temp
                v = solve_triangular(self.L, w, lower=True, trans=True,
                                     check_finite=False, unit_diagonal=True)
                v = v[self.pit]

        return v

    def swap_rows(self, H, i):
        """
        Swaps row i of H with next row; represents matrix product by PI_i
        matrix described after matrix 5.10
        """
        if i == 0:
            H[i:i+2] = H[i+1::-1]
        else:
            H[i:i+2] = H[i+1:i-1:-1]

    def row_subtract(self, H, i):
        """
        Zeros first nonzero element of row i+1 of H by subtracting appropriate
        multiple of row i; represents matrix product by matrix 5.10. Returns
        factor g for storage.
        """
        g = H[i+1, i]/H[i, i]
        H[i+1, i:] -= g*H[i, i:]
        return g

    def hess_lu(self, H, i):
        """
        Converts Hessenberg matrix H with first nonzero off-diagonal in
        column i to upper triangular, recording elementary row operations.
        That is, performs and records operations in Equation 5.9.
        """
        m = H.shape[1]
        ops = []
        for k in range(i, m-1):
            piv1, piv2 = np.abs(H[k:k+2, k])
            swap = piv1 < piv2
            # swap rows to ensure |g| <= 1
            if swap:
                self.swap_rows(H, k)
            g = self.row_subtract(H, k)
            ops.append((swap, g))  # record elementary row operations
        return ops

    def perform_ops(self, y, ops, reverse=False):
        """
        Replays operations needed to convert Hessenberg matrix into upper
        triangular form on a vector y. Equivalent to matrix multlication by
        inverse of matrix 5.12.
        """
        m = len(y)
        i = m - len(ops) - 1
        if not reverse:
            for op, k in zip(ops, range(i, m-1)):
                swap, g = op
                if swap:
                    self.swap_rows(y, k)
                y[k+1] -= g*y[k]
        else:
            for op, k in zip(reversed(ops), reversed(range(i, m-1))):
                swap, g = op
                y[k] -= g*y[k+1]
                if swap:
                    self.swap_rows(y, k)
        return y

    def perform_perm(self, p):
        """
        Perform individual row swaps defined in p returned by factor_lu to
        generate final permutation indices pi
        """
        pi = np.arange(len(p))
        for i, row in enumerate(p):
            pi[i], pi[row] = pi[row], pi[i]
        return pi
