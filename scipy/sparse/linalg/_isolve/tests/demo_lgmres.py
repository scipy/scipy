import scipy.sparse.linalg as la
import scipy.io as io
import numpy as np
import sys

#problem = "SPARSKIT/drivcav/e05r0100"
problem = "SPARSKIT/drivcav/e05r0200"
#problem = "Harwell-Boeing/sherman/sherman1"
#problem = "misc/hamm/add32"

mm = np.lib._datasource.Repository('https://math.nist.gov/pub/MatrixMarket2/')
f = mm.open(f'{problem}.mtx.gz')
Am = io.mmread(f).tocsr()
f.close()

f = mm.open(f'{problem}_rhs1.mtx.gz')
b = np.array(io.mmread(f)).ravel()
f.close()

bnorm = np.linalg.norm(b)
count = [0]


def matvec(v):
    count[0] += 1
    sys.stderr.write(f"{count[0]}\r")
    return Am@v


A = la.LinearOperator(matvec=matvec, shape=Am.shape, dtype=Am.dtype)

M = 100

print("MatrixMarket problem %s" % problem)
print(f"Invert {Am.shape[0]} x {Am.shape[1]} matrix; nnz = {Am.nnz}")

count[0] = 0
x0, info = la.gmres(A, b, restrt=M, tol=1e-14)
count_0 = count[0]
err0 = np.linalg.norm(Am@x0 - b) / bnorm
print(f"GMRES({M}): {count_0} matvecs, relative residual: {err0}")
if info != 0:
    print("Didn't converge")

count[0] = 0
x1, info = la.lgmres(A, b, inner_m=M-6*2, outer_k=6, tol=1e-14)
count_1 = count[0]
err1 = np.linalg.norm(Am@x1 - b) / bnorm
print(f"LGMRES({M - 2*6}, 6) [same memory req.]: {count_1} "
      f"matvecs, relative residual: {err1}")
if info != 0:
    print("Didn't converge")

count[0] = 0
x2, info = la.lgmres(A, b, inner_m=M-6, outer_k=6, tol=1e-14)
count_2 = count[0]
err2 = np.linalg.norm(Am@x2 - b) / bnorm
print(f"LGMRES({M - 6}, 6) [same subspace size]: {count_2} "
      f"matvecs, relative residual: {err2}")
if info != 0:
    print("Didn't converge")
