import os

import numpy as np
from numpy.testing import assert_allclose

from .common import Benchmark, is_xslow, safe_import

with safe_import():
    from scipy.optimize import Bounds, milp

with safe_import():
    from scipy.optimize.tests.test_linprog import magic_square


# MIPLIB 2017 benchmarks included with permission of the authors
# The MIPLIB benchmark problem set was downloaded from https://miplib.zib.de/.
# An MPS converter (scikit-glpk) was used to load the data into Python. The
# arrays were arranged to the format required by `milp` and saved to `npz`
# format using `np.savez`. The reduced case below keeps most of the
# piperout-27 inequality constraints so the default suite retains a sizable
# MILP benchmark without timing out.
milp_problems = ["piperout-27-reduced", "piperout-27"]
milp_problem_data = {
    "piperout-27-reduced": ("piperout-27", 14336),
    "piperout-27": ("piperout-27", None),
}


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
    """Benchmark MILP feasibility with magic squares.

    A magic square arranges the integers from 1 through n**2 so that every
    row, column, and the two main diagonals have the same sum.
    See https://en.wikipedia.org/wiki/Magic_square
    """

    params = [[3, 4, 5, 6]]
    param_names = ['size']

    @staticmethod
    def _bounds(n, fix_top_left=False):
        """Return binary bounds, optionally fixing the top-left cell."""
        if not fix_top_left:
            return (0, 1)

        lb = np.zeros(n**4)
        ub = np.ones(n**4)
        row = col = 0

        # A magic square uses each integer from 1 through n**2, and
        # `magic_square` encodes its placement as x[value - 1, row, col].
        # Fixing the top-left cell to n**2 - 1, verified feasible for sizes 3-6,
        # reduces symmetric solutions and makes the problem easier to solve.
        value = n**2 - 1
        lb[(value - 1) * n**2 + row * n + col] = 1

        return Bounds(lb, ub)

    def setup(self, n):
        if n > 4 and not is_xslow():
            raise NotImplementedError("skipped")

        A_eq, b_eq, self.c, self.numbers, self.M = magic_square(n)
        self.constraints = (A_eq, b_eq, b_eq)

        # Break symmetry in the feasibility problem so the xslow sizes remain
        # under ASV's timeout.
        self.bounds = self._bounds(n, fix_top_left=True)

    def time_magic_square(self, n):
        res = milp(c=self.c*0, constraints=self.constraints,
                   bounds=self.bounds, integrality=True)
        assert res.status == 0
        x = np.round(res.x)
        s = (self.numbers.flatten() * x).reshape(n**2, n, n)
        square = np.sum(s, axis=0)
        assert_allclose(square.sum(axis=0), self.M)
        assert_allclose(square.sum(axis=1), self.M)
        assert_allclose(np.diag(square).sum(), self.M)
        assert_allclose(np.diag(square[:, ::-1]).sum(), self.M)
