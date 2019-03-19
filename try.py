from scipy.sparse import csr_matrix,coo_matrix
import numpy as np

row  = np.array([0, 3, 1, 0])
col  = np.array([0, 3, 1, 2])
data1 = np.array([0, 5, 7, 9])
data2 = np.array([0, 4, 6, 8])
m1 = coo_matrix((data1, (row, col)), shape=(4, 4))
m2 = coo_matrix((data2, (row, col)), shape=(4, 4))
mm = m1 - m2
assert(mm.nnz == 4)