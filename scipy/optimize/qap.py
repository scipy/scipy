

import numpy as np
import math
from sklearn.utils import check_array


def quadratic_assignment(cost_matrix, dist_matrix, seed_cost = [], seed_dist = [], maximize = False):
    cost_matrix = check_array(cost_matrix, copy=True, ensure_2d=True)
    dist_matrix = check_array(dist_matrix, copy=True, ensure_2d=True)
    seed_cost = np.asarray(seed_cost)
    seed_dist = np.asarray(seed_dist)


    if cost_matrix.shape[0] != dist_matrix.shape[0]:
        msg = "Adjacency matrices must be of equal size"
        raise ValueError(msg)
    elif cost_matrix.shape[0] != cost_matrix.shape[1] or dist_matrix.shape[0] != dist_matrix.shape[1]:
        msg = "Adjacency matrix entries must be square"
        raise ValueError(msg)
    elif (cost_matrix >= 0).all() == False or (dist_matrix >= 0).all() == False:
        msg = "Adjacency matrix entries must be greater than or equal to zero"
        raise ValueError(msg)
    elif seed_cost.shape[0] != seed_dist.shape[0]:
        msg = "Seed arrays must be of equal size"
        raise ValueError(msg)
    elif seed_cost.shape[0] > cost_matrix.shape[0]:
        msg = "There cannot be more seeds than there are nodes"
        raise ValueError(msg)
    elif (seed_cost >= 0).all() == False or (seed_dist >= 0).all() == False:
        msg = "Seed array entries must be greater than or equal to zero"
        raise ValueError(msg)
    elif (seed_cost <= (cost_matrix.shape[0] - 1)).all() == False or (
            seed_dist <= (cost_matrix.shape[0] - 1)
    ).all() == False:
        msg = "Seed array entries must be less than or equal to n-1"
        raise ValueError(msg)

    n = cost_matrix.shape[0]  # number of vertices in graphs
    n_seeds = seed_cost.shape[0]  # number of seeds
    n_unseed = n - n_seeds

    score = math.inf
    perm_inds = np.zeros(n)

    obj_func_scalar = 1
    if maximize:
        obj_func_scalar = -1
        score = 0

    if self.shuffle_input:
        W2_c = np.random.permutation(np.array([x for x in range(n) if x not in seed_dist]))
        # shuffle_input to avoid results from inputs that were already matched
    else:
        W2_c = np.array([x for x in range(n) if x not in seed_dist])

    W1_c = np.array([x for x in range(n) if x not in seed_cost])
    p_A = np.concatenate([seed_cost, W1_c], axis=None).astype(int)
    p_B = np.concatenate([seed_dist, W2_c], axis=None).astype(int)
    A = cost_matrix[np.ix_(p_A, p_A)]
    B = dist_matrix[np.ix_(p_B, p_B)]

    A11 = A[:n_seeds, :n_seeds]
    A12 = A[:n_seeds, n_seeds:]
    A21 = A[n_seeds:, :n_seeds]
    A22 = A[n_seeds:, n_seeds:]
    B11 = B[:n_seeds, :n_seeds]
    B12 = B[:n_seeds, n_seeds:]
    B21 = B[n_seeds:, :n_seeds]
    B22 = B[n_seeds:, n_seeds:]
    A11T = np.transpose(A11)
    A12T = np.transpose(A12)
    A22T = np.transpose(A22)
    B21T = np.transpose(B21)
    B22T = np.transpose(B22)

    for i in range(self.n_init):

        # setting initialization matrix
        if self.init_method == "rand":
            sk = SinkhornKnopp()
            K = np.random.rand(
                n_unseed, n_unseed
            )  # generate a nxn matrix where each entry is a random integer [0,1]
            for i in range(10):  # perform 10 iterations of Sinkhorn balancing
                K = sk.fit(K)
            J = np.ones((n_unseed, n_unseed)) / float(
                n_unseed
            )  # initialize J, a doubly stochastic barycenter
            P = (K + J) / 2
        elif self.init_method == "barycenter":
            P = np.ones((n_unseed, n_unseed)) / float(n_unseed)

        const_sum = A21 @ np.transpose(B21) + np.transpose(A12) @ B12
        grad_P = math.inf  # gradient of P
        n_iter = 0  # number of FW iterations

        # OPTIMIZATION WHILE LOOP BEGINS
        while grad_P > self.eps and n_iter < self.max_iter:
            delta_f = (
                    const_sum + A22 @ P @ B22T + A22T @ P @ B22
            )  # computing the gradient of f(P) = -tr(APB^tP^t)
            rows, cols = linear_sum_assignment(
                obj_func_scalar * delta_f
            )  # run hungarian algorithm on gradient(f(P))
            Q = np.zeros((n_unseed, n_unseed))
            Q[rows, cols] = 1  # initialize search direction matrix Q

            def f(x):  # computing the original optimization function
                return obj_func_scalar * (
                        np.trace(A11T @ B11)
                        + np.trace(np.transpose(x * P + (1 - x) * Q) @ A21 @ B21T)
                        + np.trace(np.transpose(x * P + (1 - x) * Q) @ A12T @ B12)
                        + np.trace(
                    A22T
                    @ (x * P + (1 - x) * Q)
                    @ B22
                    @ np.transpose(x * P + (1 - x) * Q)
                )
                )

            alpha = minimize_scalar(
                f, bounds=(0, 1), method="bounded"
            ).x  # computing the step size
            P_i1 = alpha * P + (1 - alpha) * Q  # Update P
            grad_P = np.linalg.norm(P - P_i1)
            P = P_i1
            n_iter += 1
        # end of FW optimization loop

        row, col = linear_sum_assignment(
            -P
        )  # Project onto the set of permutation matrices
        perm_inds_new = np.concatenate(
            (np.arange(n_seeds), np.array([x + n_seeds for x in col]))
        )

        score_new = np.trace(
            np.transpose(A) @ B[np.ix_(perm_inds_new, perm_inds_new)]
        )  # computing objective function value

        if obj_func_scalar * score_new < obj_func_scalar * score:  # minimizing
            score = score_new
            perm_inds = np.array([0] * n)
            perm_inds[p_A] = p_B[perm_inds_new]

    perm_inds = perm_inds.astype(int)
    p_A_unshuffle = np.array(range(n))
    p_A_unshuffle[p_A] = np.array(range(n))
    A = A[np.ix_(p_A_unshuffle, p_A_unshuffle)]
    p_B_unshuffle = np.array(range(n))
    p_B_unshuffle[p_B] = np.array(range(n))
    B = B[np.ix_(p_B_unshuffle, p_B_unshuffle)]
    score = np.trace(np.transpose(A) @ B[np.ix_(perm_inds, perm_inds)])

    self.perm_inds_ = perm_inds  # permutation indices
    self.score_ = score  # objective function value
    return self