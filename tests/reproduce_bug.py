import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

print("scipy", __import__("scipy").__version__, "numpy", np.__version__)

#fixed seed for reproducible results
A = sp.random_array((5, 5), density=0.5, format="csr",
                    rng=np.random.default_rng(0))

dtypes = [np.float32, np.float64, np.complex64, np.complex128]

print("\nnumpy dense, ord=1")
for dtype in dtypes:
    print(f"in: {dtype.__name__}, out: "
          f"{np.linalg.norm(A.astype(dtype).toarray(), ord=1).dtype}")

print("\nsparse norm, default ord")
for dtype in dtypes:
    print(f"in: {dtype.__name__}, out: {spla.norm(A.astype(dtype)).dtype}")

print("\nsparse norm, ord=1  (bug: float32 -> float64)")
for dtype in dtypes:
    print(f"in: {dtype.__name__}, out: {spla.norm(A.astype(dtype), 1).dtype}")