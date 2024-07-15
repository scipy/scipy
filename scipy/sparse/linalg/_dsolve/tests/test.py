import scipy
from scipy.sparse.linalg import splu
import time
import numpy as np
rng = np.random.default_rng()

n = 1000
a = scipy.sparse.random_array((n, n), density=0.01, random_state=rng)
m = (a.T @ a) + scipy.sparse.identity(n)
print("sparsity: ", float(m.nnz)/n**2)
st = time.time()
B = splu(m,diag_pivot_thresh = 0.01,
         options = dict(SymmetricMode = True),
         permc_spec = "METIS_AT_PLUS_A")
print(time.time()-st)
