from pyprima.common.linalg import matprod, inprod
import numpy as np
import timeit
import sys

for numrows in [2, 3, 5, 10, 20, 50, 100]:
    np.random.seed(0)
    sys.modules['pyprima'].common.linalg.COMPARING = False
    np_result = timeit.timeit(lambda: matprod(np.random.random((numrows, numrows)), np.random.random((numrows, numrows))), number=1000)
    np.random.seed(0)
    sys.modules['pyprima'].common.linalg.COMPARING = True
    naive_result = timeit.timeit(lambda: matprod(np.random.random((numrows, numrows)), np.random.random((numrows, numrows))), number=1000)
    speedup = naive_result / np_result
    print(f"For {numrows} rows, numpy implementation takes {np_result} seconds, naive implementation takes {naive_result} seconds. Numpy is {speedup:.2f}x faster.")


# for numrows in [2, 3, 5, 10, 20, 50, 100]:
#     np.random.seed(0)
#     sys.modules['pyprima'].common.linalg.COMPARING = False
#     np_result = timeit.timeit(lambda: inprod(np.random.random(numrows), np.random.random(numrows)), number=1000)
#     np.random.seed(0)
#     sys.modules['pyprima'].common.linalg.COMPARING = True
#     naive_result = timeit.timeit(lambda: inprod(np.random.random(numrows), np.random.random(numrows)), number=1000)
#     speedup = naive_result / np_result
#     print(f"For {numrows} rows, numpy implementation takes {np_result} seconds, naive implementation takes {naive_result} seconds. Numpy is {speedup:.2f}x faster.")