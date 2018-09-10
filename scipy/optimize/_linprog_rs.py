#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 18 22:07:15 2018

@author: matthaberland
"""

import numpy as np
from scipy.linalg import solve, norm
#from _bglu_dense import BGLU
from scipy.optimize._bglu_dense import BGLU
from scipy.optimize import linprog


def _phase_one(A, b, maxiter, tol, maxupdate):
    """
    The purpose of phase one is to find an initial basic feasible solution
    (BFS) to the original problem.

    Generates an auxiliary problem with a trivial BFS and an objective that
    minimizes infeasibility of the original problem. Solves the auxiliary
    problem using the main simplex routine (phase two). This either yields
    a BFS to the original problem or determines that the original problem is
    infeasible. If feasible, phase one detects redundant rows in the original
    constraint matrix and removes them, then chooses additional indices as
    necessary to complete a basis/BFS for the original problem.
    """

    m, n = A.shape
    status = 0

    # generate auxiliary problem to get initial BFS
    A, b, c, basis, x = _generate_auxiliary_problem(A, b)

    # solve auxiliary problem
    x, basis, status = _phase_two(c, A, x, basis, maxiter, tol, maxupdate)

    # check for infeasibility
    residual = c.dot(x)
    if status == 0 and residual > tol:
        status = 2

#    if status == 0:
    # detect redundancy
    B = A[:, basis]
    rank_revealer = solve(B, A[:, :n])
    z = _find_nonzero_rows(rank_revealer, tol)

    # eliminate redundancy
    A = A[z, :n]
    b = b[z]

    # form solution to original problem
    x = x[:n]
    m = A.shape[0]
    basis = basis[basis < n]

    # choose additional indices to complete basis
    if len(basis) < m:
        a = np.arange(m+n)
        bl = np.zeros(len(a), dtype=bool)
        bl[basis] = 1
        new_basis = a[~bl][:m-len(basis)]
        basis = np.concatenate((basis, new_basis))

    return x, basis, A, b, residual, status


def _generate_auxiliary_problem(A, b):
    """
    Modifies original problem to create an auxiliary problem with a trivial
    intial basic feasible solution and an objective that minimizes
    infeasibility in the original problem.

    Conceptually this is done by stacking an identity matrix on the right of
    the original constraint matrix, adding artificial variables to correspond
    with each of these new columns, and generating a cost vector that is all
    zeros except for ones corresponding with each of the new variables.

    A initial basic feasible solution is trivial: all variables are zero
    except for the artificial variables, which are set equal to the
    corresponding element of the right hand side `b`.

    Runnning the simplex method on this auxiliary problem drives all of the
    artificial variables - and thus the cost - to zero if the original problem
    is feasible. The original problem is declared infeasible otherwise.

    Much of the complexity below is to improve efficiency by using singleton
    columns in the original problem where possible and generating artificial
    variables only as necessary.
    """
    A = A.copy()  # FIXME: no need to copy once used with rest of linprog
    b = b.copy()
    m, n = A.shape

    A[b < 0] = -A[b < 0]  # express problem with RHS positive for trivial BFS
    b[b < 0] = -b[b < 0]  # to the auxiliary problem

    # chooses existing columns appropriate for inclusion in inital basis
    cols, rows = _select_singleton_columns(A, b)

    acols = np.arange(m-len(cols))          # indices of auxiliary columns
    arows = np.delete(np.arange(m), rows)   # indices of corresponding rows,
                                            #  that is, the row in each aux
                                            #  column with nonzero entry

    basis = np.concatenate((cols, n + acols))   # all initial basis columns
    basis_rows = np.concatenate((rows, arows))  # all intial basis rows

    # add auxiliary singleton columns
    A = np.hstack((A, np.zeros((m, m-len(cols)))))
    A[arows, n + acols] = 1

    # generate intial BFS
    x = np.zeros(m+n-len(cols))
    x[basis] = b[basis_rows]/A[basis_rows, basis]

    # generate costs to minimize infeasibility
    c = np.zeros(m+n-len(cols))
    c[basis] = 1

    return A, b, c, basis, x


def _select_singleton_columns(A, b):
    """
    Finds singleton columns for which the singleton entry is of the same sign
    as the right hand side; these columns are eligible for inclusion in an
    initial basis. Determines the rows in which the singleton entries are
    located. For each of these rows, returns the indices of the one singleton
    column and its corresponding row.
    """
    # find indices of all singleton columns and corresponding row indicies
    column_indices = np.nonzero(np.sum(np.abs(A) != 0, axis=0) == 1)[0]
    columns = A[:, column_indices]          # array of singleton columns
    row_indices = np.nonzero(columns)[0]    # corresponding row indicies

    # keep only singletons with entries that have same sign as RHS
    same_sign = A[row_indices, column_indices]*b[row_indices] >= 0
    column_indices = column_indices[same_sign]
    row_indices = row_indices[same_sign]
    # this is necessary because all elements of BFS must be non-negative

    # for each row, keep only one singleton column with an entry in that row
    unique_row_indices, first_columns = np.unique(row_indices,
                                                  return_index=True)
    return column_indices[first_columns], unique_row_indices


def _find_nonzero_rows(A, tol):
    """
    Returns logical array indicating the locations of rows with at least
    one nonzero element.
    """
    return np.any(np.abs(A) > tol, axis=1)


def _select_enter_pivot(c_hat, bl, a, rule="bland"):
    """
    Selects a pivot to enter the basis. Currently Bland's rule - the smallest
    index that has a negative reduced cost - is the default.
    """
    if rule.lower() == "mrc":  # index with minimum reduced cost
        return a[~bl][np.argmin(c_hat)]
    else:
        return a[~bl][c_hat < 0][0]  # smallest index w/ negative reduced cost


def _phase_two(c, A, x, b, maxiter, tol, maxupdate):
    """
    The heart of the simplex method. Beginning with a basic feasible solution,
    moves to adjacent basic feasible solutions successively lower reduced cost.
    Terminates when there are no basic feasible solutions with lower reduced
    cost or if the problem is determined to be unbounded.

    This implementation follows the revised simplex method based on LU
    decomposition. Rather than maintaining a tableau or an inverse of the
    basis matrix, we keep a factorization of the basis matrix that allows
    efficient solution of linear systems while avoiding stability issues
    associated with inverted matrices.
    """
    m, n = A.shape
    status = 0
    a = np.arange(n)                # indices of columns of A
    ab = np.arange(m)               # indices of columns of B
    B = BGLU(A, b, maxupdate)       # basis matrix factorization object
                                    # similar to B = A[:, b]

    for k in range(maxiter):
        bl = np.zeros(len(a), dtype=bool)
        bl[b] = 1

        N = A[:, ~bl]   # non-basis matrix
        xb = x[b]       # basic variables
        cb = c[b]       # basic costs

        v = B.solve(cb, transposed=True)    # similar to v = solve(B.T, cb)
        c_hat = c[~bl] - v.T.dot(N)         # reduced costs

        if np.all(c_hat >= -tol):  # all reduced costs positive -> terminate
            break

        j = _select_enter_pivot(c_hat, bl, a)
        u = B.solve(A[:, j])    # similar to u = solve(B, A[:, j])

        i = u > 0               # if there are none, it's unbounded
        if not np.any(i):
            status = 3
            break

        th = xb[i]/u[i]
        l = np.argmin(th)           # implicitly selects smallest subscript
        th_star = th[l]             # step size

        x[b] = x[b] - th_star*u     # take step
        x[j] = th_star
        B.update(ab[i][l], j)       # modify basis
        b = B.b                     # similar to b[ab[i][l]] = j
    else:
        status = 1

    return x, b, status


# FIXME: is maxiter for each phase?
def _linprog_rs(c, A, b, maxiter=1000, tol=1e-9, maxupdate = 20):
    """
    Performs the two phase simplex method to solve a linear programming
    problem in standard form.
    """
    messages = ["Optimization terminated successfully.",
                "Iteration limit reached.",
                "The problem appears infeasible, as the phase one auxiliary "
                "problem terminated successfully with a residual of {0:.1e}, "
                "greater than the tolerance {1} required for the solution to "
                "be considered feasible. Consider increasing the tolerance to "
                "be greater than {0:.1e}. If this tolerance is unnaceptably "
                "large, the problem is likely infeasible.",
                "The problem is unbounded, as the simplex algorithm found "
                "a basic feasible solution from which there is a direction "
                "with negative reduced cost in which all decision variables "
                "increase."]

    x, basis, A, b, residual, status = _phase_one(A, b, maxiter, tol, maxupdate)
    if status == 0:
        x, basis, status = _phase_two(c, A, x, basis, maxiter, tol, maxupdate)
    return x, status, messages[status].format(residual, tol)

#c = np.array([-10, -12, -12, 0, 0, 0])
#b = np.array([20, 20, 20])
#A = np.array([[1, 2, 2, 1, 0, 0], [2, 1, 2, 0, 1, 0], [2, 2, 1, 0, 0, 1]])
#x = np.array([0, 0, 0, 20, 20, 20])
#basis = np.array([3, 4, 5])

#c = np.array([-3/4, 20, -1/2, 6, 0, 0, 0])
#b = np.array([0, 0, 1])
#A = np.array([[1/4, -8, -1, 9, 1, 0, 0],
#              [1/2, -12, -1/2, 3, 0, 1, 0],
#              [0, 0, 1, 0, 0, 0, 1]], dtype=float)
#x = np.array([0, 0, 0, 0, 0, 0, 1], dtype=float)
#basis = np.array([4, 5, 6])

#c = np.array([1, 1, 1, 0])
#b = np.array([3, 2, 1])
#A = np.array([[1, 2, 3, 0],
#              [-1, 2, 6, 0],
#              [0, 0, 3, 1]])

#c = np.array([1, 1, 1, 0])
#b = np.array([3, 2, 5, 1])
#A = np.array([[1, 2, 3, 0],
#              [-1, 2, 6, 0],
#              [0, 4, 9, 0],
#              [0, 0, 3, 1]])

def myrand(shape, scale):
    return np.random.rand(*shape)*scale - scale/2
#
m = 10
n = round(2*m)
scale = 100
c = myrand((n,), scale)
A = myrand((m,n), scale)
b = myrand((m,), scale)
#A[5] = A[2] + A[3]
#A[6] = A[2] + A[4]

#A[:, 15] = 0
#A[5, 15] = 1
#A[:, 16] = 0
#A[5, 16] = 1
#A[:, 17] = 0
#A[6, 17] = 1

x, status, message = _linprog_rs(c, A, b, maxiter = 1000)

print(x)
print(message)

res = linprog(c, A_eq = A, b_eq = b, method = "simplex")
print(res.x)
print(res.message)

print(norm(x - res.x))