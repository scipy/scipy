import os

import numpy as np
from numpy.testing import assert_allclose

from .common import Benchmark, is_xslow, safe_import

with safe_import():
    from scipy.optimize import milp


# MIPLIB 2017 benchmarks included with permission of the authors
# The MIPLIB benchmark problem set was downloaded from https://miplib.zib.de/.
# An MPS converter (scikit-glpk) was used to load the data into Python. The
# arrays were arranged to the format required by `milp` and saved to `npz`
# format using `np.savez`. The reduced case below is derived from piperout-27
# so the default suite retains a sizable MILP benchmark without timing out.
milp_problems = ["piperout-27-reduced", "piperout-27"]
milp_problem_data = {
    "piperout-27-reduced": ("piperout-27", 128),
    "piperout-27": ("piperout-27", None),
}


def magic_square(n):
    """
    Generate a linear program for which integer solutions represent an
    n x n magic square.
    """
    rng = np.random.default_rng(92350948245690509234)
    M = n * (n**2 + 1) / 2

    numbers = np.arange(n**4) // n**2 + 1
    numbers = numbers.reshape(n**2, n, n)

    zeros = np.zeros((n**2, n, n))

    A_list = []
    b_list = []

    # Rule 1: use every number exactly once
    for i in range(n**2):
        A_row = zeros.copy()
        A_row[i, :, :] = 1
        A_list.append(A_row.flatten())
        b_list.append(1)

    # Rule 2: only one number per square
    for i in range(n):
        for j in range(n):
            A_row = zeros.copy()
            A_row[:, i, j] = 1
            A_list.append(A_row.flatten())
            b_list.append(1)

    # Rule 3: sum of rows is M
    for i in range(n):
        A_row = zeros.copy()
        A_row[:, i, :] = numbers[:, i, :]
        A_list.append(A_row.flatten())
        b_list.append(M)

    # Rule 4: sum of columns is M
    for i in range(n):
        A_row = zeros.copy()
        A_row[:, :, i] = numbers[:, :, i]
        A_list.append(A_row.flatten())
        b_list.append(M)

    # Rule 5: sum of diagonals is M
    A_row = zeros.copy()
    A_row[:, range(n), range(n)] = numbers[:, range(n), range(n)]
    A_list.append(A_row.flatten())
    b_list.append(M)
    A_row = zeros.copy()
    A_row[:, range(n), range(-1, -n - 1, -1)] = \
        numbers[:, range(n), range(-1, -n - 1, -1)]
    A_list.append(A_row.flatten())
    b_list.append(M)

    A = np.array(np.vstack(A_list), dtype=float)
    b = np.array(b_list, dtype=float)
    c = rng.random(A.shape[1])

    return A, b, c, numbers, M


class MilpMiplibBenchmarks(Benchmark):
    params = [milp_problems]
    param_names = ['problem']

    def setup(self, prob):
        if prob == "piperout-27" and not is_xslow():
            raise NotImplementedError("skipped")

        if not hasattr(self, 'data'):
            dir_path = os.path.dirname(os.path.realpath(__file__))
            datafile = os.path.join(dir_path, "linprog_benchmark_files",
                                    "milp_benchmarks.npz")
            self.data = np.load(datafile, allow_pickle=True)

        data_name, max_ub_constraints = milp_problem_data[prob]
        c, A_ub, b_ub, A_eq, b_eq, bounds, integrality = self.data[data_name]

        lb = [li for li, ui in bounds]
        ub = [ui for li, ui in bounds]

        cons = []
        if A_ub is not None:
            if max_ub_constraints is not None:
                # The full MIPLIB instance has timed out in regular ASV runs.
                # Use a deterministic subset of inequality constraints to
                # keep this benchmark suitable for the default suite.
                A_ub = A_ub[:max_ub_constraints, :]
                b_ub = b_ub[:max_ub_constraints]
            cons.append((A_ub, -np.inf, b_ub))
        if A_eq is not None:
            cons.append((A_eq, b_eq, b_eq))

        self.c = c
        self.constraints = cons
        self.bounds = (lb, ub)
        self.integrality = integrality

    def time_milp(self, prob):
        res = milp(c=self.c, constraints=self.constraints, bounds=self.bounds,
                   integrality=self.integrality)
        assert res.success


class MilpMagicSquare(Benchmark):

    params = [[3, 4, 5, 6]]
    param_names = ['size']

    def setup(self, n):
        if n > 4 and not is_xslow():
            raise NotImplementedError("skipped")

        A_eq, b_eq, self.c, self.numbers, self.M = magic_square(n)
        self.constraints = (A_eq, b_eq, b_eq)

    def time_magic_square(self, n):
        res = milp(c=self.c*0, constraints=self.constraints,
                   bounds=(0, 1), integrality=True)
        assert res.status == 0
        x = np.round(res.x)
        s = (self.numbers.flatten() * x).reshape(n**2, n, n)
        square = np.sum(s, axis=0)
        assert_allclose(square.sum(axis=0), self.M)
        assert_allclose(square.sum(axis=1), self.M)
        assert_allclose(np.diag(square).sum(), self.M)
        assert_allclose(np.diag(square[:, ::-1]).sum(), self.M)
