

import numpy as np
import math
from sklearn.utils import check_array


def quadratic_assignment(cost_matrix, dist_matrix, seed_cost = [], seed_dist = [], maximize = False):
    """Solve the quadratic assignment problem

        This class solves the Quadratic Assignment Problem and the Graph Matching Problem
        (QAP) through an implementation of the Fast Approximate QAP Algorithm (FAQ) (these
        two problems are the same up to a sign change) [1].

        This algorithm can be thought of as finding an alignment of the vertices of two
        graphs which minimizes the number of induced edge disagreements, or, in the case
        of weighted graphs, the sum of squared differences of edge weight disagreements.
        The option to add seeds (known vertex correspondence between some nodes) is also
        available [2].


        Parameters
        ----------
        cost_matrix : 2d-array, square, positive
            A square adjacency matrix

        dist_matrix : 2d-array, square, positive
            A square adjacency matrix

        seed_cost : 1d-array, shape (m , 1) where m <= number of nodes (default = [])
            An array where each entry is an index of a node in `cost_matrix`.

        seeds_dist : 1d-array, shape (m , 1) where m <= number of nodes (default = [])
            An array where each entry is an index of a node in `dist_matrix` The elements of
            `seed_cost` and `seed_dist` are vertices which are known to be matched, that is,
            `seed_cost[i]` is matched to vertex `seed_dist[i]`.

        maximize : bool (default = False)
            Gives users the option to solve the Graph Matching Problem (GMP) rather than QAP.
            This is accomplished through trivial negation of the objective function.

        n_init : int, positive (default = 1)
            Number of random initializations of the starting permutation matrix that
            the FAQ algorithm will undergo. n_init automatically set to 1 if
            init_method = 'barycenter'

        init_method : string (default = 'barycenter')
            The initial position chosen

            "barycenter" : the non-informative “flat doubly stochastic matrix,”
            :math:`J=1*1^T /n` , i.e the barycenter of the feasible region

            "rand" : some random point near :math:`J, (J+K)/2`, where K is some random doubly
            stochastic matrix

        max_iter : int, positive (default = 30)
            Integer specifying the max number of Franke-Wolfe iterations.
            FAQ typically converges with modest number of iterations.

        shuffle_input : bool (default = True)
            Gives users the option to shuffle the nodes of A matrix to avoid results
            from inputs that were already matched.

        eps : float (default = 0.1)
            A positive, threshold stopping criteria such that FW continues to iterate
            while Frobenius norm of :math:`(P_{i}-P_{i+1}) > eps`

        gmp : bool (default = True)
            Gives users the option to solve QAP rather than the Graph Matching Problem
            (GMP). This is accomplished through trivial negation of the objective function.

        Attributes
        ----------

        perm_inds_ : array, size (n,) where n is the number of vertices in the fitted graphs.
            The indices of the optimal permutation (with the fixed seeds given) on the nodes of B,
            to best minimize the objective function :math:`f(P) = trace(A^T PBP^T )`.


        score_ : float
            The objective function value of for the optimal permutation found.


        References
        ----------
        .. [1] J.T. Vogelstein, J.M. Conroy, V. Lyzinski, L.J. Podrazik, S.G. Kratzer,
            E.T. Harley, D.E. Fishkind, R.J. Vogelstein, and C.E. Priebe, “Fast
            approximate quadratic programming for graph matching,” PLOS one, vol. 10,
            no. 4, p. e0121002, 2015.

        .. [2] D. Fishkind, S. Adali, H. Patsolic, L. Meng, D. Singh, V. Lyzinski, C. Priebe,
            Seeded graph matching, Pattern Recognit. 87 (2019) 203–215



        """

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

    seed_dist_c = np.setdiff1d(range(n), seed_dist)
    if self.shuffle_input:
        seed_dist_c = np.random.permutation(seed_dist_c)
        # shuffle_input to avoid results from inputs that were already matched

    seed_cost_c = np.setdiff1d(range(n), seed_cost)
    permutation_cost = np.concatenate([seed_cost, seed_cost_c], axis=None).astype(int)
    permutation_dist = np.concatenate([seed_dist, seed_dist_c], axis=None).astype(int)
    cost_matrix = cost_matrix[np.ix_(permutation_cost, permutation_cost)]
    dist_matrix = dist_matrix[np.ix_(permutation_dist, permutation_dist)]

    # definitions according to Seeded Graph Matching [2].
    A11 = cost_matrix[:n_seeds, :n_seeds]
    A12 = cost_matrix[:n_seeds, n_seeds:]
    A21 = cost_matrix[n_seeds:, :n_seeds]
    A22 = cost_matrix[n_seeds:, n_seeds:]
    B11 = dist_matrix[:n_seeds, :n_seeds]
    B12 = dist_matrix[:n_seeds, n_seeds:]
    B21 = dist_matrix[n_seeds:, :n_seeds]
    B22 = dist_matrix[n_seeds:, n_seeds:]
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
            np.transpose(cost_matrix) @ dist_matrix[np.ix_(perm_inds_new, perm_inds_new)]
        )  # computing objective function value

        if obj_func_scalar * score_new < obj_func_scalar * score:  # minimizing
            score = score_new
            perm_inds = np.zeros(n, dtype=int)
            perm_inds[permutation_cost] = permutation_dist[perm_inds_new]

    permutation_cost_unshuffle = _unshuffle(permutation_cost, n)
    cost_matrix = cost_matrix[np.ix_(permutation_cost_unshuffle, permutation_cost_unshuffle)]
    permutation_dist_unshuffle = _unshuffle(permutation_dist, n)
    dist_matrix = dist_matrix[np.ix_(permutation_dist_unshuffle, permutation_dist_unshuffle)]
    score = np.trace(np.transpose(cost_matrix) @ dist_matrix[np.ix_(perm_inds, perm_inds)])

    self.perm_inds_ = perm_inds  # permutation indices
    self.score_ = score  # objective function value
    return self

def _unshuffle(array, n):
    unshuffle = np.array(range(n))
    unshuffle[array] = np.array(range(n))
    return unshuffle
