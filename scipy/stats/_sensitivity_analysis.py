import numpy as np
from scipy.stats import qmc, bootstrap


def f_ishigami(x):
    """Ishigami function.

    :param array_like x: [n, [x1, x2, x3]].
    :return: Function evaluation.
    :rtype: array_like (n, 1)
    """
    x = np.atleast_2d(x)
    return np.array([np.sin(xi[0]) + 7 * np.sin(xi[1])**2 +
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

    return A, B, AB


def sobol_saltelli(f_A, f_B, f_AB):
    f_AB = f_AB.reshape(-1, f_A.shape[0])

    var = np.var(np.vstack([f_A, f_B]))
    # var = np.clip(var, a_min=0, a_max=None)

    s = np.mean(f_B * (f_AB - f_A.flatten()).T, axis=0) / var
    st = 0.5 * np.mean((f_A.flatten() - f_AB).T ** 2, axis=0) / var

    return s, st


def sobol_indices(func, *, n, d, l_bounds, u_bounds=None, seed=None):
    """Sobol' indices.

    The total number of function call is ``d(d+2)``.
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

    s, st = sobol_saltelli(f_A, f_B, f_AB)

    def sobol_saltelli_(idx):
        f_A_ = f_A[idx]
        f_B_ = f_B[idx]
        f_AB_ = f_AB[idx]
        return sobol_saltelli(f_A_, f_B_, f_AB_)

    ci = bootstrap(
        [np.arange(n)], sobol_saltelli_, method="BCa",
        batch=int(n*0.7), n_resamples=99
    )

    return s, st, ci


indices = sobol_indices(
    f_ishigami, n=1024, d=3,
    l_bounds=[-np.pi, -np.pi, -np.pi], u_bounds=[np.pi, np.pi, np.pi]
)
print(indices)
