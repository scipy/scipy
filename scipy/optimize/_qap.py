import numpy as np
import operator
from . import linear_sum_assignment, minimize_scalar, OptimizeResult


def quadratic_assignment(
    cost_matrix,
    dist_matrix,
    maximize=False,
    options=None
):
    r"""
    Solve the quadratic assignment problem.

    This function solves the Quadratic Assignment Problem (QAP) and the
    Graph Matching Problem through an implementation of the Fast
    Approximate QAP Algorithm (FAQ) (these two problems are the same up
    to a sign change) [1]_.

    Quadratic Assignment solves problems of the following form:

    .. math::

        \min_P & \ {\ \text{trace}(APB^T P^T)}\\
        \mbox{s.t. } & {P \ \epsilon \ \mathcal{P}}\\

    where :math:`\mathcal{P}` is the set of all permutation matrices,
    and :math:`A` and :math:`B` are adjacency matrices.

    This algorithm can be thought of as finding an alignment of the
    vertices of two graphs which minimizes the number of induced edge
    disagreements, or, in the case of weighted graphs, the sum of squared
    differences of edge weight disagreements. The option to add seeds
    (known vertex correspondence between some nodes) is also available
    [2]_.

    Note that the quadratic assignment problem is NP-hard, is not
    known to be solvable in polynomial time, and is computationally
    intractable. Therefore, the results given are approximations,
    not guaranteed to be exact solutions.


    Parameters
    ----------
    cost_matrix : 2d-array, square, non-negative
        A square adjacency matrix. In this implementation, :math:`A` =
        `cost-matrix` in the objective function above.

    dist_matrix : 2d-array, square, non-negative
        A square adjacency matrix.  In this implementation, :math:`B` =
        `dist-matrix` in the objective function above.

    seed : 2d-array, optional, (default = None)
        Allows the user apply a seed, fixing part of the matching between
        the two adjacency matrices.
        For column 1, each entry is an index of a node in `cost_matrix`.
        For column 2, each entry is an index of a node in `dist_matrix`.
        The elements of ``seed[:, 0]`` and ``seed[:, 1]`` are vertices
        which are known to be matched, that is, ``seed[i, 0]`` is matched to
        vertex ``seed[i, 1]``. Array shape ``(m , 2)`` where ``m <= number of
        nodes``.

    maximize : bool (default = False)
        Gives users the option to solve the Graph Matching Problem (GMP)
        rather than QAP. This is accomplished through trivial negation
        of the objective function.


    options : dict, optional
        A dictionary of solver options. All methods accept the following
        options:

            init : 2d-array (default = 'barycenter')
                The algorithm may be sensitive to the initial permutation
                matrix (or search position) chosen due to the possibility
                of several local minima within the feasible region.
                With only 1 initialization, a barycenter init will
                likely return a more accurate permutation.

                Choosing several random initializations (through
                `init_weight` and `init_n`) as opposed to the non-informative
                barycenter will likely result in a more accurate result at
                the cost of higher runtime.

                The initial position chosen:

                "barycenter" : the non-informative "flat doubly stochastic
                matrix," :math:`J=1*1^T /n` , i.e the barycenter of the
                feasible region (where :math:`n` is the number of nodes and
                :math:`1` is a ``(n, 1)`` array of ones).

                If an ndarray is passed, it should have the same shape as
                `cost_matrix` and `dist_matrix`, and its rows and columns
                must sum to 1 (doubly stochastic).
            init_weight : float, positive in range [0,1] (Default = None)
                Allows the user to specify the amount of random perturbation
                from the starting position of `init`
                At each initialization, the initial permutation matrix
                is some random point near :math:`J`, defined as
                :math:`(\alpha J + (1- \alpha) K`, where :math: `\alpha`
                is `init_weight` and :math:`K` is some random doubly
                stochastic matrix.
            init_n : int, positive (default = 1)
                Number of random initializations of the starting
                permutation matrix that the FAQ algorithm will undergo.
            maxiter : int, positive (default = 30)
                Integer specifying the max number of Franke-Wolfe iterations.
                FAQ typically converges with modest number of iterations.
            shuffle_input : bool (default = True)
                To avoid artificially high or low matching due to inherent
                sorting of input adjacency matrices, gives users the option
                to shuffle the nodes of `cost_matrix`. Results are then
                unshuffled so that returned `col_ind` matches the node order
                of inputs.
            eps : float (default = 0.05)
                A positive, threshold stopping criteria such that Franke-
                Wolfe continues to iterate while Frobenius norm of
                :math:`(P_{i}-P_{i+1}) > eps`, where :math:`i` is the
                iteration number.

    Returns
    -------
    res : OptimizeResult
        A :class:`scipy.optimize.OptimizeResult` consisting of the fields:

            col_ind : 1-D array
                An array of column indices corresponding to the optimal
                permutation (with the fixed seeds given) of the
                nodes of `dist_matrix`, to best minimize the objective
                function.
            score : float
                The optimal value of the objective function.
            nit : int
                The total number of Franke-Wolfe iterations performed during
                optimization.

    References
    ----------
    .. [1] J.T. Vogelstein, J.M. Conroy, V. Lyzinski, L.J. Podrazik,
           S.G. Kratzer, E.T. Harley, D.E. Fishkind, R.J. Vogelstein, and
           C.E. Priebe, "Fast approximate quadratic programming for graph
           matching," PLOS one, vol. 10, no. 4, p. e0121002, 2015.

    .. [2] D. Fishkind, S. Adali, H. Patsolic, L. Meng, D. Singh, V. Lyzinski,
           C. Priebe, "Seeded graph matching", Pattern Recognit. 87 (2019):
           203-215.

    Examples
    --------

    >>> cost = np.array([[0, 80, 150, 170], [80, 0, 130, 100],
    ...         [150, 130, 0, 120], [170, 100, 120, 0]])
    >>> dist = np.array([[0, 5, 2, 7], [0, 0, 3, 8],
    ...         [0, 0, 0, 3], [0, 0, 0, 0]])
    >>> from scipy.optimize import quadratic_assignment
    >>> res = quadratic_assignment(cost, dist)
    >>> print(res)
     col_ind: array([0, 3, 2, 1])
         nit: 9
       score: 3260

    To demonstrate explicitly how the `score` value
    :math:`f(P) = trace(A^T PBP^T )` is calculated, one may construct the
    permutation matrix, and perform the necessary algebra.

    >>> n = cost.shape[0]
    >>> P = np.zeros((n, n))
    >>> P[np.arange(n), res['col_ind']] = 1
    >>> score = int(np.trace(cost.T @ P @ dist @ P.T))
    >>> print(score)
    3260

    As you can see, the value here matches res['score'] reported above.
    Alternatively, to avoid constructing the permutation matrix, one can also
    perform the following calculation.

    >>> score = np.trace(cost.T @ dist[np.ix_(res['col_ind'], res['col_ind'])])
    >>> print(score)
    3260

    Here, we are simply permuting the distance matrix.

    """

    if options is None:
        options = {}

    return _quadratic_assignment_faq(cost_matrix, dist_matrix, maximize,
                                     **options)


def _quadratic_assignment_faq(
        cost_matrix,
        dist_matrix,
        maximize=False,
        partial_match=None,
        init="barycenter",
        init_weight=None,
        init_n=1,
        maxiter=30,
        shuffle_input=True,
        eps=0.05
):

    cost_matrix = np.asarray(cost_matrix)
    dist_matrix = np.asarray(dist_matrix)

    if partial_match is None:
        partial_match = np.array([[], []]).T
    partial_match = np.asarray(partial_match)
    init_n = operator.index(init_n)
    maxiter = operator.index(maxiter)

    # ValueError check
    msg = None
    if cost_matrix.shape[0] != cost_matrix.shape[1]:
        msg = "'cost_matrix' must be square"
    elif dist_matrix.shape[0] != dist_matrix.shape[1]:
        msg = "'dist_matrix' must be square"
    elif cost_matrix.shape != dist_matrix.shape:
        msg = "Adjacency matrices must be of equal size"
    elif (cost_matrix < 0).any() or (dist_matrix < 0).any():
        msg = "Adjacency matrix contains negative entries"
    elif partial_match.shape[0] > cost_matrix.shape[0]:
        msg = "There cannot be more seeds than there are nodes"
    elif partial_match.shape[1] != 2:
        msg = "`partial_match` must have two columns"
    elif (partial_match < 0).any():
        msg = "`partial_match` contains negative entries"
    elif (partial_match >= len(cost_matrix)).any():
        msg = "`partial_match` entries must be less than the number of nodes"
    elif not len(set(partial_match[:, 0])) == len(partial_match[:, 0]) or not \
            len(set(partial_match[:, 1])) == len(partial_match[:, 1]):
        msg = "`partial_match` column entries must be unique"
    elif isinstance(init, str) and init not in {'barycenter'}:
        msg = "Invalid 'init' parameter string"
    elif init_weight is not None and (init_weight < 0 or init_weight > 1):
        msg = "'init_weight' must be in range [0, 1]"
    elif init_n <= 0:
        msg = "'n_init' must be a positive integer"
    elif maxiter <= 0:
        msg = "'maxiter' must be a positive integer"
    if msg is not None:
        raise ValueError(msg)

    # TypeError check
    if type(shuffle_input) is not bool:
        msg = "'shuffle_input' must be a boolean"
    elif eps <= 0 or type(eps) is not float:
        msg = "'eps' must be a positive float"
    elif type(maximize) is not bool:
        msg = "'maximize' must be a boolean"
    if msg is not None:
        raise TypeError(msg)

    rng = np.random.RandomState()
    n = cost_matrix.shape[0]  # number of vertices in graphs
    n_seeds = partial_match.shape[0]  # number of seeds
    n_unseed = n - n_seeds

    perm_inds = np.zeros(n)

    obj_func_scalar = 1
    if maximize:
        obj_func_scalar = -1
    score = obj_func_scalar * np.inf

    seed_dist_c = np.setdiff1d(range(n), partial_match[:, 1])
    if shuffle_input:
        seed_dist_c = rng.permutation(seed_dist_c)
        # shuffle_input to avoid results from inputs that were already matched

    seed_cost_c = np.setdiff1d(range(n), partial_match[:, 0])
    permutation_cost = np.concatenate([partial_match[:, 0],
                                       seed_cost_c], axis=None).astype(int)
    permutation_dist = np.concatenate([partial_match[:, 1],
                                       seed_dist_c], axis=None).astype(int)
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

    for i in range(init_n):
        # setting initialization matrix
        if isinstance(init, str) and init == 'barycenter':
            J = np.ones((n_unseed, n_unseed)) / float(n_unseed)
        else:
            _check_init_input(init, n_unseed)
            J = init
        if init_weight is not None:
            # generate a nxn matrix where each entry is a random number [0, 1]
            K = rng.rand(n_unseed, n_unseed)
            # Sinkhorn balancing
            K = _doubly_stochastic(K)
            # initialize J, a doubly stochastic barycenter
            P = J * init_weight + (1 - init_weight) * K
        else:
            P = J
        const_sum = A21 @ B21.T + A12.T @ B12
        grad_P = np.inf  # gradient of P
        n_iter = 0  # number of FW iterations

        # OPTIMIZATION WHILE LOOP BEGINS
        while grad_P > eps and n_iter < maxiter:
            # computing the gradient of f(P) = -tr(APB^tP^t)
            delta_f = (const_sum + A22 @ P @ B22.T + A22.T @ P @ B22)
            # run hungarian algorithm on gradient(f(P))
            rows, cols = linear_sum_assignment(obj_func_scalar * delta_f)
            Q = np.zeros((n_unseed, n_unseed))
            Q[rows, cols] = 1  # initialize search direction matrix Q

            def f(x):  # computing the original optimization function
                return obj_func_scalar * (
                    np.trace(A11.T @ B11)
                    + np.trace(np.transpose(x * P + (1 - x) * Q) @ A21 @ B21.T)
                    + np.trace(np.transpose(x * P + (1 - x) * Q) @ A12.T @ B12)
                    + np.trace(
                        A22.T
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
        ).astype(int)

        score_new = np.trace(
            np.transpose(cost_matrix)
            @ dist_matrix[np.ix_(perm_inds_new, perm_inds_new)]
        )  # computing objective function value

        if obj_func_scalar * score_new < obj_func_scalar * score:  # minimizing
            score = score_new
            perm_inds = np.zeros(n, dtype=int)
            perm_inds[permutation_cost] = permutation_dist[perm_inds_new]

    permutation_cost_inv = np.argsort(permutation_cost)
    cost_matrix = cost_matrix[
        np.ix_(permutation_cost_inv, permutation_cost_inv)
    ]
    permutation_dist_inv = np.argsort(permutation_dist)
    dist_matrix = dist_matrix[
        np.ix_(permutation_dist_inv, permutation_dist_inv)
    ]

    score = np.trace(
        np.transpose(cost_matrix) @ dist_matrix[np.ix_(perm_inds, perm_inds)]
    )

    res = {"col_ind": perm_inds, "score": score, "nit": n_iter}

    return OptimizeResult(res)


def _check_init_input(init, n):
    row_sum = np.round(np.sum(init, axis=0), decimals=2)
    col_sum = np.round(np.sum(init, axis=1), decimals=2)
    msg = None
    if init.shape != (n, n):
        msg = "`init` matrix must have same shape as A and B"
    elif (row_sum != 1.).any() or (col_sum != 1.).any() or (init < 0).any():
        msg = "`init` matrix must be doubly stochastic"
    if msg is not None:
        raise ValueError(msg)


def _doubly_stochastic(P, eps=1e-3):
    # cleaner implementation of btaba/sinkhorn_knopp
    # Title: sinkhorn_knopp Source Code
    # Author: Tabanpour, B
    # Date: 2018
    # Code version:  0.2
    # Availability: https://pypi.org/project/sinkhorn_knopp/
    #

    max_iter = 1000
    c = 1 / P.sum(axis=0)
    r = 1 / (P @ c)
    P_eps = P

    for it in range(max_iter):
        if ((np.abs(P_eps.sum(axis=1) - 1) < eps).all() and
                (np.abs(P_eps.sum(axis=0) - 1) < eps).all()):
            # All column/row sums ~= 1 within threshold
            break

        c = 1 / (r @ P)
        r = 1 / (P @ c)
        P_eps = r[:, None] * P * c

    return P_eps
