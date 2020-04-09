'''Solve LP given numpy matrices.'''

from pyHiGHS import highs_wrapper

import numpy as np

if __name__ == '__main__':

    c = np.array([5, 4, 3]).astype('double')
    A = np.array([
        [2, 3, 1],
        [4, 1, 2],
        [3, 4, 2],
    ]).astype('double')
    b = np.array([5, 11, 8]).astype('double')
    options = {
        'presolve': True,
        'sense': -1,
        'solver': 'simplex',
        'parallel': True,
        'time_limit': 1,
        'message_level': 0,
        'write_solution_to_file': False,
        'solution_file': 'test.sol',
        'write_solution_pretty': True,
    }
    res = highs_wrapper(c, A, b, options=options)
    print(res)
