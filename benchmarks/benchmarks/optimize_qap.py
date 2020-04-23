import numpy as np
from .common import Benchmark
import math
import os

try:
    from scipy.optimize import quadratic_assignment
except ImportError:
    pass


class QuadraticAssignment(Benchmark):
    params = ["bur26a", "bur26b", "bur26c", "bur26d", "bur26e", "bur26f", "bur26g", "bur26h",
              "chr12a", "chr12b", "chr12c", "chr15a", "chr15b", "chr15c", "chr18a", "chr18b",
              "chr20a", "chr20b", "chr20c", "chr22a", "chr22b", "chr25a", "els19",
              "esc16a", "esc16b", "esc16c", "esc16d", "esc16e", "esc16g", "esc16h", "esc16i",
              "esc16j", "esc32e", "esc32g", "esc128",
              "had12", "had14", "had16", "had18", "had20", "kra30a", "kra30b", "kra32",
              "lipa20a", "lipa20b", "lipa30a", "lipa30b", "lipa40a", "lipa40b", "lipa50a",
              "lipa50b", "lipa60a", "lipa60b", "lipa70a", "lipa70b", "lipa80a", "lipa90a",
              "lipa90b",
              "nug12", "nug14", "nug16a", "nug16b", "nug17", "nug18", "nug20", "nug21", "nug22",
              "nug24", "nug25", "nug27", "nug28", "nug30",
              "rou12", "rou15", "rou20",
              "scr12", "scr15", "scr20",
              "sko42", "sko49", "sko56", "sko64", "sko72", "sko81", "sko90", "sko100a", "sko100b",
              "sko100c", "sko100d", "sko100e", "sko100f", "ste36b", "ste36c",
              "tai12a", "tai12b", "tai15a", "tai15b", "tai17a", "tai20a", "tai20b", "tai25a",
              "tai25b", "tai30a", "tai30b", "tai35a", "tai40a", "tai40b", "tai50a", "tai50b",
              "tai60a", "tai60b", "tai64c", "tai80a", "tai100a", "tai100b", "tai150b", "tai256c",
              "tho30", "tho40", "tho150", "wil50", "wil100"]
    param_names = ['QAP Problem']

    def setup(self, qap_prob):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        datafile = os.path.join(dir_path, "qapdata", qap_prob + ".dat")
        slnfile = os.path.join(dir_path, "qapdata/qapsoln", qap_prob + ".sln")
        with open(datafile) as f:
            f = [int(elem) for elem in f.read().split()]

            # adjusting
            f = np.array(f[1:])
            n = int(math.sqrt(len(f) / 2))
            f = f.reshape(2 * n, n)
            self.cost_matrix = f[:n, :]
            self.dist_matrix = f[n:, :]
        with open(slnfile) as f:
            f = [int(elem) for elem in f.read().split()]
            self.opt_solution = f[1]

    def time_evaluation(self, qap_prob):
        quadratic_assignment(self.cost_matrix, self.dist_matrix)

    def track_score(self, qap_prob):
        row, col = quadratic_assignment(self.cost_matrix, self.dist_matrix)
        score = int(np.trace(np.transpose(self.cost_matrix) @ self.dist_matrix[np.ix_(col, col)]))
        percent_diff = (score - self.opt_solution) / self.opt_solution
        return percent_diff
