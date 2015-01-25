# Posix-only benchmark
from __future__ import division, absolute_import, print_function

import os
import sys
import re
import subprocess
import time
import textwrap

from numpy.testing import dec
from scipy.stats import spearmanr

import numpy as np
from numpy.testing import Tester


@dec.skipif(not sys.platform.startswith('linux'), "Memory benchmark works only on Linux")
def bench_leaks():
    mem_info = get_mem_info()
    set_mem_rlimit(int(mem_info['memtotal'] * 0.7))

    # Setup temp file, make it fit in memory
    print_table_row(['repeats', 'peak memory (MB)'])

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

        print_table_row(["%d" % repeat, "%.1f" % (peak_mem/1e6,)])

    print("")
    corr, p = spearmanr(repeats, peak_mems)
    if p < 0.05:
        print("*"*79)
        print("PROBABLE MEMORY LEAK")
        print("*"*79)
        raise AssertionError("Probable memory leak")
    else:
        print("PROBABLY NO MEMORY LEAK")


def print_table_row(columns):
    print(" | ".join("%-20s" % x for x in columns))


def run_monitored(code):
    """
    Run code in a new Python process, and monitor peak memory usage.

    Returns
    -------
    duration : float
        Duration in seconds (including Python startup time)
    peak_memusage : float
        Peak memory usage (rough estimate only) in bytes

    """
    code = textwrap.dedent(code)
    process = subprocess.Popen([sys.executable, '-c', code],
                               cwd=os.path.dirname(__file__))

    peak_memusage = -1

    start = time.time()
    while True:
        ret = process.poll()
        if ret is not None:
            break

        with open('/proc/%d/status' % process.pid, 'r') as f:
            procdata = f.read()

        m = re.search('VmRSS:\s*(\d+)\s*kB', procdata, re.S | re.I)
        if m is not None:
            memusage = float(m.group(1)) * 1e3
            peak_memusage = max(memusage, peak_memusage)

        time.sleep(0.01)

    process.wait()

    duration = time.time() - start

    if process.returncode != 0:
        raise AssertionError("Running failed:\n%s" % code)

    return duration, peak_memusage


def get_mem_info():
    """Get information about available memory"""
    info = {}
    with open('/proc/meminfo', 'r') as f:
        for line in f:
            p = line.split()
            info[p[0].strip(':').lower()] = float(p[1]) * 1e3
    return info


def set_mem_rlimit(max_mem):
    """
    Set rlimit to 80% of total system memory, to avoid grinding halt
    because of swapping.
    """
    import resource
    cur_limit = resource.getrlimit(resource.RLIMIT_AS)
    if cur_limit[0] > 0:
        max_mem = min(max_mem, cur_limit[0])

    resource.setrlimit(resource.RLIMIT_AS, (max_mem, cur_limit[1]))


if __name__ == "__main__":
    Tester().bench()
