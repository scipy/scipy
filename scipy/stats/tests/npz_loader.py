"""
This is a helper script. Its purpose is to populate an single npz file from
several other files (e.g. for nist test cases)
"""
import os
import sys
import numpy as np
try:
    from scipy._lib.six import u
except ImportError:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                    os.pardir, 'scipy', 'lib'))
    from six import u

# These are the nist ANOVA files. They can be found at:
# http://www.itl.nist.gov/div898/strd/anova/anova.html
filenames = [u('SiRstv.dat'), u('SmLs01.dat'), u('SmLs02.dat'), u('SmLs03.dat'),
             u('AtmWtAg.dat'), u('SmLs04.dat'), u('SmLs05.dat'),
             u('SmLs06.dat'), u('SmLs07.dat'), u('SmLs08.dat'), u('SmLs09.dat')]


def getnist(filename):
    fname = os.path.abspath(os.path.join(os.path.dirname(__file__), filename))
    content = file(fname, 'r').read().split('\n')
    certified = [line.split() for line in content[40:48] if line != '\r']
    dataf = np.loadtxt(fname, skiprows=60)
    y, x = dataf.T
    y = y.astype(int)
    caty = np.unique(y)
    f = float(certified[0][-1])
    R2 = float(certified[2][-1])
    resstd = float(certified[4][-1])
    dfbn = int(certified[0][-4])
    dfwn = int(certified[1][-3])
    return [y, x, f, dfbn, dfwn, R2, resstd, certified, caty]

if __name__ == "__main__":
    res = {}
    for fn in filenames:
        res[fn] = (getnist(fn))

    np.savez(os.path.abspath(os.path.join(os.path.dirname(__file__),
                                          'nist_anova.npz')),
             **res)
