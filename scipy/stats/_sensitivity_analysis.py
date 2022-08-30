import numpy as np
from scipy.stats import qmc


def f_ishigami(x):
    """Ishigami function.

    :param array_like x: [n, [x1, x2, x3]].
    :return: Function evaluation.
    :rtype: array_like (n, 1)
    """
    x = np.atleast_2d(x)
    return np.array([np.sin(xi[0]) + 7 * np.sin(xi[1])**2 + \
                     0.1 * (xi[2]**4) * np.sin(xi[0]) for xi in x]).reshape(-1, 1)


def sample_A_B_AB(d, n, seed):
    rng = np.random.default_rng(seed)

    # A and B
    A_B = qmc.Sobol(d=2*d, seed=rng, bits=64).random(n)
    A, B = A_B[:, :d], A_B[:, d:]

    # AB: columns of B into A
    AB = np.empty((int(d*n), d))
    for i in range(d):
        AB[i*n:(i+1)*n, :] = np.column_stack((A[:, :i], B[:, i], A[:, i+1:]))

    # # BA: columns of A into B
    # BA = np.empty((int(d*n), d))
    # for i in range(d):
    #     BA[i*n:(i+1)*n, :] = np.column_stack((B[:, :i], A[:, i], B[:, i+1:]))

    return A, B, AB


def sobol_saltelli(f_A, f_B, f_AB):
    var = np.var(np.vstack([f_A, f_B]), axis=0)
    # var = np.clip(var, a_min=0, a_max=None)

    # Total effect of pairs of factors: generalization of Saltenis
    # st2 = np.zeros((dim, dim))
    # for j, i in itertools.combinations(range(0, dim), 2):
    #     st2[i, j] = 1 / (2 * n_sample) * np.sum((f_AB[i] - f_AB[j]) ** 2, axis=0) / var

    s = np.mean(f_B * (f_AB - f_A.flatten()).T, axis=0) / var
    st = 1 / 2 * np.mean((f_A.flatten() - f_AB).T ** 2, axis=0) / var

    return s, st


def sobol_indices(func, *, n, d, l_bounds, u_bounds=None, seed=None):
    """Sobol' indices.

    The total number of function call is N(p+2).
    Three matrices are required for the computation of
    the indices: A, B and a permutation matrix AB based
    on both A and B.

    """
    A, B, AB = sample_A_B_AB(d, n, seed)

    if (l_bounds is not None) or (u_bounds is not None):
        A = qmc.scale(A, l_bounds=l_bounds, u_bounds=u_bounds)
        B = qmc.scale(B, l_bounds=l_bounds, u_bounds=u_bounds)
        AB = qmc.scale(AB, l_bounds=l_bounds, u_bounds=u_bounds)

    f_A = func(A)
    f_B = func(B)

    f_AB = func(AB)
    f_AB = f_AB.reshape((d, n))

    s, st= sobol_saltelli(f_A, f_B, f_AB)

    return s, st


# print(f_ishigami([[0, 1, 3], [2, 4, 2]]))

indices = sobol_indices(
    f_ishigami, n=1024, d=3,
    l_bounds=[-np.pi, -np.pi, -np.pi], u_bounds=[np.pi, np.pi, np.pi]
)
print(indices)

pass
