from __future__ import division, absolute_import, print_function

import numpy as np

from .common import run_monitored, set_mem_rlimit, Benchmark

try:
    from scipy.stats import spearmanr
except ImportError:
    pass

try:
    from scipy.interpolate import PPoly
except ImportError:
    pass


class Leaks(Benchmark):
    unit = "relative increase with repeats"

    def track_leaks(self):
        set_mem_rlimit()

        # Setup temp file, make it fit in memory
        repeats = [2, 5, 10, 50, 200]
        peak_mems = []

        for repeat in repeats:
            code = """
            import numpy as np
            from scipy.interpolate import griddata

            def func(x, y):
                return x*(1-x)*np.cos(4*np.pi*x) * np.sin(4*np.pi*y**2)**2

            grid_x, grid_y = np.mgrid[0:1:100j, 0:1:200j]
            points = np.random.rand(1000, 2)
            values = func(points[:,0], points[:,1])

            for t in range(%(repeat)d):
                for method in ['nearest', 'linear', 'cubic']:
                    griddata(points, values, (grid_x, grid_y), method=method)

            """ % dict(repeat=repeat)

            _, peak_mem = run_monitored(code)
            peak_mems.append(peak_mem)

        corr, p = spearmanr(repeats, peak_mems)
        if p < 0.05:
            print("*"*79)
            print("PROBABLE MEMORY LEAK")
            print("*"*79)
        else:
            print("PROBABLY NO MEMORY LEAK")

        return max(peak_mems) / min(peak_mems)


class BenchPPoly(Benchmark):

    def setup(self):
        np.random.seed(1234)
        m, k = 55, 3
        x = np.sort(np.random.random(m+1))
        c = np.random.random((3, m))
        self.pp = PPoly(c, x)

        npts = 100
        self.xp = np.linspace(0, 1, npts)

    def time_evaluation(self):
        self.pp(self.xp)
