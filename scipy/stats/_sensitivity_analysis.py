import numpy as np
from scipy.stats import qmc, bootstrap
from scipy.stats.sampling import NumericalInversePolynomial
import scipy.stats as stats


def f_ishigami(x):
    """Ishigami function.

    :param array_like x: [n, [x1, x2, x3]].
    :return: Function evaluation.
    :rtype: array_like (n, 1)
    """
    x = np.atleast_2d(x)
    return np.array([np.sin(xi[0]) + 7 * np.sin(xi[1])**2 +
                     0.1 * (xi[2]**4) * np.sin(xi[0]) for xi in x]).reshape(-1, 1)


def sample_A_B(d, n, *, dists=None, seed=None):
    rng = np.random.default_rng(seed)
    A_B = qmc.Sobol(d=2*d, seed=rng, bits=64).random(n)

    A, B = A_B[:, :d], A_B[:, d:]

    if dists is not None:
        for d_, dist in enumerate(dists):
            dist_rng = NumericalInversePolynomial(dist, random_state=rng)
            A[:, d_] = dist_rng.ppf(A[:, d_])
            B[:, d_] = dist_rng.ppf(B[:, d_])

    return A, B


def sample_A_B_AB(d, n, *, dists=None, seed=None):
    # A and B
    A, B = sample_A_B(d=d, n=n, dists=dists, seed=seed)

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


def sobol_indices(*, func, n, d, dists=None, l_bounds=None, u_bounds=None, seed=None):
    """Sobol' indices.

    The total number of function call is ``d(d+2)``.
    Three matrices are required for the computation of
    the indices: A, B and a permutation matrix AB based
    on both A and B.

    """
    A, B, AB = sample_A_B_AB(d=d, n=n, dists=dists, seed=seed)

    if (l_bounds is not None) or (u_bounds is not None):
        A = qmc.scale(A, l_bounds=l_bounds, u_bounds=u_bounds)
        B = qmc.scale(B, l_bounds=l_bounds, u_bounds=u_bounds)
        AB = qmc.scale(AB, l_bounds=l_bounds, u_bounds=u_bounds)

    f_A = func(A)
    f_B = func(B)
    f_AB = func(AB)

    # Y = (Y - Y.mean()) / Y.std()

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


# indices = sobol_indices(
#     func=f_ishigami, n=1024, d=3,
#     l_bounds=[-np.pi, -np.pi, -np.pi], u_bounds=[np.pi, np.pi, np.pi]
# )
# print(indices[:2])

# aggregated indices
# sum(Vik)/sum(var(Vk))
# Vik = Sik*Vk


def discrepancy_indices(*, func, n, d, dists=None, l_bounds=None, u_bounds=None, seed=None):
    rng = np.random.default_rng(seed)
    sample = qmc.Sobol(d=d, seed=rng, bits=64).random(n)

    if dists is not None:
        for d_, dist in enumerate(dists):
            dist_rng = NumericalInversePolynomial(dist, random_state=rng)
            sample[:, d_] = dist_rng.ppf(sample[:, d_])

    if (l_bounds is not None) or (u_bounds is not None):
        sample_scaled = qmc.scale(sample, l_bounds=l_bounds, u_bounds=u_bounds)
    else:
        sample_scaled = sample

    f_sample = func(sample_scaled)

    f_sample = qmc.scale(
        f_sample, l_bounds=np.min(f_sample), u_bounds=np.max(f_sample), reverse=True
    )

    indices = [
        qmc.discrepancy(
            np.concatenate([sample[:, i].reshape(-1, 1), f_sample], axis=1)
        )
        for i in range(d)
    ]

    # l_discrepancy = qmc.discrepancy(sample[:, :2])
    # u_discepancy = qmc.discrepancy(np.concatenate([sample[:, 0].reshape(-1, 1), np.ones((1024, 1))*1], axis=1))
    # indices_scaled = qmc.scale(np.array(indices).reshape(-1, 1), l_bounds=l_discrepancy, u_bounds=u_discepancy, reverse=True)

    indices = indices / np.sum(indices)

    return indices


# indices = discrepancy_indices(
#     func=f_ishigami, n=1024, d=3,
#     l_bounds=[-np.pi, -np.pi, -np.pi], u_bounds=[np.pi, np.pi, np.pi]
# )
# print(indices)

def discontinuous(x):
    res = np.zeros_like(x)
    res[np.where(x > 0.5)] = 1
    return res


functions = [
    ("Linear", lambda x: x),
    ("Quadratic", lambda x: x**2),
    ("Cubic", lambda x: x**3),
    ("Exponential", lambda x: (np.exp(x) - 1) / (np.exp(1) - 1)),
    ("Periodic", lambda x: np.sin(2 * np.pi * x) / 2 + 0.5),
    ("Discontinuous", discontinuous),
    ("No effect", lambda x: x * 0),
    ("Non monotonic", lambda x: 4 * (x - 0.5) ** 2),
    ("Inverse", lambda x: (10 - 1 / 1.1) ** (-1) * (x + 0.1) ** (- 1) - 0.1),
    ("Trigonometric", lambda x: np.cos(x)),
    ("Piecewise large", lambda x: ((-1) ** (4 * x).astype(int) * (0.125 - np.mod(x, 0.25)) + 0.125)),
    ("Piecewise small", lambda x: ((-1) ** (32 * x).astype(int) * (0.03125 - 2 * np.mod(x, 0.03125)) + 0.03125) / 2),
    ("Oscillation", lambda x: x ** 2 - 0.2 * np.cos(7 * np.pi * x)),
]

distributions = [
    ("uniform", stats.uniform()),
    ("normal", stats.norm(loc=0.5, scale=0.15)),
    ("beta", stats.beta(a=8, b=2)),
    ("beta2", stats.beta(a=2, b=8)),
    ("beta3", stats.beta(a=2, b=0.8)),
    ("beta4", stats.beta(a=0.8, b=2)),
    ("logitnormal", stats.logistic(loc=0, scale=3.16)),
]


def mixture_phi(n, rng):
    dists = [
        stats.norm(loc=0, scale=np.sqrt(0.5)),
        stats.norm(loc=0, scale=np.sqrt(5)),
    ]
    coefficients = np.array([0.7, 0.3])

    n_dists = len(dists)
    data = np.zeros((n, n_dists))
    for i, dist in enumerate(dists):
        data[:, i] = dist.rvs(n, random_state=rng)
    random_idx = rng.choice(np.arange(n_dists), size=(n,), p=coefficients)
    sample = data[np.arange(n), random_idx]

    return sample





class MetaModel:

    def __init__(self, d, seed=None):

        rng = np.random.default_rng(seed)

        n_funcs = len(functions)

        self.d = d

        # first order
        u = rng.integers(low=0, high=n_funcs, size=d)
        self.u = np.array(functions)[u]
        # phi[0] * u * u

        d2 = int(0.5 * d)
        self.V = rng.integers(low=0, high=d, size=(d2, 2))
        # phi[1] * u[V[0]] * u[V[1]]

        d3 = int(0.2 * d)
        self.W = rng.integers(low=0, high=d, size=(d3, 3))
        # phi[2] * u[W[0]] * u[W[1]] * u[W[2]]

        phi = mixture_phi(d + d2 + d3, rng)
        self.phi = [phi[:d], phi[d:d+d2], phi[d+d2:]]

    def __call__(self, x):

        first_order = np.zeros((len(x),))
        for d_ in range(self.d):
            first_order += self.phi[0][d_] * self.u[d_][1](x[:, d_])

        second_order = np.zeros((len(x),))
        for i, v_ in enumerate(self.V):
            second_order += self.phi[1][i] \
                            * self.u[v_[0]][1](x[:, v_[0]]) \
                            * self.u[v_[1]][1](x[:, v_[1]])

        third_order = np.zeros((len(x),))
        for i, w_ in enumerate(self.W):
            third_order += self.phi[2][i] \
                           * self.u[w_[0]][1](x[:, w_[0]]) \
                           * self.u[w_[1]][1](x[:, w_[1]]) \
                           * self.u[w_[2]][1](x[:, w_[2]])

        return (first_order + second_order + third_order).reshape(-1, 1)


rng = np.random.default_rng()
metamodel = MetaModel(d=4, seed=rng)
sample = rng.random((100, 4))
metamodel(sample)

pass
