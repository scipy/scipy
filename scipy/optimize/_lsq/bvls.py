import numpy as np
from numpy.linalg import norm, lstsq
from scipy.optimize import OptimizeResult

from .common import print_header_linear, print_iteration_linear


def bvls(A, b, x_lsq, lb, ub, tol, max_iter, verbose):
    m, n = A.shape

    x = x_lsq.copy()
    on_bound = np.zeros(n)

    mask = x < lb
    x[mask] = lb[mask]
    on_bound[mask] = -1

    mask = x > ub
    x[mask] = ub[mask]
    on_bound[mask] = 1

    free_set = on_bound == 0
    active_set = ~free_set
    free_set, = np.where(free_set)

    r = A.dot(x) - b
    cost = 0.5 * np.dot(r, r)
    initial_cost = cost
    g = A.T.dot(r)

    cost_change = None
    step_norm = None
    iteration = 0

    if verbose == 2:
        print_header_linear()

    while free_set.size > 0:
        if verbose == 2:
            optimality = max(norm(g[free_set]), np.max(g * on_bound))
            print_iteration_linear(iteration, cost, cost_change, step_norm,
                                   optimality)

        iteration += 1
        x_free_old = x[free_set].copy()

        A_free = A[:, free_set]
        b_free = b - A.dot(x * active_set)
        z = lstsq(A_free, b_free)[0]

        lbv = z < lb[free_set]
        ubv = z > ub[free_set]
        v = lbv | ubv

        if np.any(lbv):
            ind = free_set[lbv]
            x[ind] = lb[ind]
            active_set[ind] = True
            on_bound[ind] = -1

        if np.any(ubv):
            ind = free_set[ubv]
            x[ind] = ub[ind]
            active_set[ind] = True
            on_bound[ind] = 1

        ind = free_set[~v]
        x[ind] = z[~v]

        r = A.dot(x) - b
        cost_new = 0.5 * np.dot(r, r)
        cost_change = cost - cost_new
        cost = cost_new
        g = A.T.dot(r)
        step_norm = norm(x[free_set] - x_free_old)

        if np.any(v):
            free_set = free_set[~v]
        else:
            break

    if max_iter is None:
        max_iter = n

    termination_status = None

    for iteration in range(iteration, max_iter):
        g_dir = g * on_bound
        move_to_free = np.argmax(g_dir)
        max_g_dir = g_dir[move_to_free]

        if max_g_dir <= 0:
            termination_status = 3

        if verbose == 2:
            g_norm = max(0, max_g_dir)
            if free_set.size > 0:
                g_norm = max(g_norm, norm(g[free_set], ord=np.inf))

            print_iteration_linear(iteration, cost, cost_change,
                                   step_norm, g_norm)

        if termination_status is not None:
            break

        on_bound[move_to_free] = 0
        free_set = on_bound == 0
        active_set = ~free_set
        free_set, = np.nonzero(free_set)

        x_free = x[free_set]
        x_free_old = x_free.copy()
        lb_free = lb[free_set]
        ub_free = ub[free_set]

        A_free = A[:, free_set]
        b_free = b - A.dot(x * active_set)

        z = lstsq(A_free, b_free)[0]

        lbv, = np.nonzero(z < lb_free)
        ubv, = np.nonzero(z > ub_free)
        v = np.hstack((lbv, ubv))

        if v.size > 0:
            alphas = np.hstack((
                lb_free[lbv] - x_free[lbv],
                ub_free[ubv] - x_free[ubv])) / (z[v] - x_free[v])

            i = np.argmin(alphas)
            i_free = v[i]
            alpha = alphas[i]

            x_free *= 1 - alpha
            x_free += alpha * z

            if i < lbv.size:
                x_free[i_free] = lb_free[i_free]
                on_bound[free_set[i_free]] = -1
            else:
                x_free[i_free] = ub_free[i_free]
                on_bound[free_set[i_free]] = 1
        else:
            x_free = z

        x[free_set] = x_free
        step_norm = norm(x_free - x_free_old)

        r = A.dot(x) - b
        cost_new = 0.5 * np.dot(r, r)
        cost_change = cost - cost_new

        if cost_change < tol * cost:
            termination_status = 2
        cost = cost_new

        g = A.T.dot(r)

    if termination_status is None:
        termination_status = 0

    g_norm = max(0, max_g_dir)
    if free_set.size > 0:
        g_norm = max(g_norm, norm(g[free_set], ord=np.inf))

    return OptimizeResult(
        x=x, fun=r, cost=cost, optimality=g_norm, active_mask=on_bound,
        nit=iteration + 1, status=termination_status,
        initial_cost=initial_cost)
