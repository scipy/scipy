'''Use cython wrapper to solve problem described by MPS file.'''

import pathlib

from pyHiGHS import linprog_mps

if __name__ == '__main__':

    #mpsfile = str(pathlib.Path(__file__).parent / '25fv47.mps')
    #linprog_mps(mpsfile)

    mpsfile = str(pathlib.Path(__file__).parent.parent.parent / 'test.mps')
    linprog_mps(mpsfile)
