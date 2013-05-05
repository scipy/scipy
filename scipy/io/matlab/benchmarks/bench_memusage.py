# Posix-only benchmark
from __future__ import division, absolute_import, print_function

import os
import sys
import re
import subprocess
import time
import textwrap
import tempfile
import warnings

from numpy.testing import dec

import numpy as np
from scipy.io import savemat, loadmat


@dec.skipif(not sys.platform.startswith('linux'), "Memory benchmark works only on Linux")
def bench_run():
    mem_info = get_mem_info()
    set_mem_rlimit(int(mem_info['memtotal'] * 0.7))

    # Setup temp file, make it fit in memory
    f = tempfile.NamedTemporaryFile(suffix='.mat')
    os.unlink(f.name)

    max_size = int(mem_info['memtotal'] * 0.7)//4
    sizes = [1e6, 10e6, 100e6, 300e6, 500e6, 1000e6]

    print_table_row(['** loadmat benchmark'])
    print_table_row(['size (MB)', 'compression', 'time (s)',
                     'peak memory (MB)', 'mem factor'])

    for size in sizes:
        for compressed in (False, True):
            if size > max_size:
                print_table_row(["%.1f" % (size/1e6,), compressed, "SKIP"])
                continue

            try:
                x = np.random.rand(size//8).view(dtype=np.uint8)
                savemat(f.name, dict(x=x), do_compression=compressed, oned_as='row')
                del x
            except MemoryError:
                x = None
                print_table_row(["%.1f" % (size/1e6,), compressed, "FAIL"])
                continue

            code = """
            from scipy.io import loadmat
            loadmat('%s')
            """ % (f.name,)
            time, peak_mem = run_monitored(code)

            print_table_row(["%.1f" % (size/1e6,), compressed, time,
                             "%.1f" % (peak_mem/1e6,),
                             "%.2f x" % (peak_mem/size,)])

    print_table_row(['** savemat memory benchmark'])
    print_table_row(['size (MB)', 'compression', 'time (s)',
                     'peak memory (MB)', 'mem factor'])

    for size in sizes:
        for compressed in (False, True):
            if size > max_size:
                print_table_row(["%.1f" % (size/1e6,), compressed, "SKIP"])
                continue

            code = """
            import numpy as np
            from scipy.io import savemat
            x = np.random.rand(%d//8).view(dtype=np.uint8)
            savemat('%s', dict(x=x), do_compression=%r, oned_as='row')
            """ % (size, f.name, compressed)
            try:
                time, peak_mem = run_monitored(code)
            except AssertionError:
                print_table_row(["%.1f" % (size/1e6,), compressed, "FAIL"])
                continue

            print_table_row(["%.1f" % (size/1e6,), compressed, time,
                             "%.1f" % (peak_mem/1e6,),
                             "%.2f x" % (peak_mem/size,)])


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
    process = subprocess.Popen([sys.executable, '-c', code])

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
    bench_run()
