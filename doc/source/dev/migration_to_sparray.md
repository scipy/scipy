## Migration from spmatrix to sparray

This document provides guidance for converting code from sparse matrices
to sparse arrays in SciPy Sparse.

The change from sparse matrices to sparse arrays mirrors conversion from
`np.matrix` to `np.array`. Essentially we must move from an all-2D matrix
multiplication centric `matrix` object to a 1D or 2D "array" object that
supports operand matrix multiplication and elementwise computation.

Notation: For this guide we denote the sparse array classes generally
as `sparray` and the sparse matrix classes `spmatrix`. Dense numpy
arrays are denoted `np.array` and dense matrix classes are `np.matrix`.
Supported sparse formats are denoted BSR, COO, CSC, CSR, DIA, DOK, LIL
and all formats are supported by both sparray and spmatrix.
The term `sparse` refers to either `sparray` or `spmatrix`, while
`dense` refers to either `np.array` or `np.matrix`.

### Overview and big picture:
  - the constructor names `*_matrix`, e.g. `csr_matrix`, are changed to `*_array`.
  - spmatrix `M` is always 2D (rows x columns) even e.g. `M.min(axis=0)`.
    sparray `A` can be 1D or 2D. Numpy scalars are used as 0D values where needed.
  - Operators change:
    - scalar operations with `*` are not implemented for some cases (see below).
    - matrices use operators `*` and @ for matrix multiplication, while arrays
      use `*` for elementwise multiplication and `@` for matrix multiplication.
    - scalar exponents for matrices are matrix power, arrays use elementwise
      power. To get matrix power for arrays use `sp.sparse.linalg.matrix_power(A, n)`
        - operators that change:  `*, @, *=, @=, **`
    -
  - Checking the sparse type and format:
    - `issparse(A)` returns `True` for any sparse array/matrix.
    - `isspmatrix(M)` returns `True` for any sparse matrix.
    - `issparse(A) and A.format == 'csr'` looks for a CSR sparse array/matrix

### Recommended steps for migration:
  - First Pass (leaving spmatrix in the code)
    - In your spmatrix code, change `*` to `@` for matrix multiplication.
      Note that scalar multiplication with sparse should use `*`.
    - scalar powers, e.g. `M**3`, should be converted to
      `sp.sparse.linalg.matrix_power(M, 3)` or `M @ M @ M` is the power is known.
    - convert construction functions like `hstack` and `triu`
      (see below: `Construction Functions:`)
    - implement alternatives to unsupported functions/methods
      like `A.getnnz()` -> `A.nnz` (see below: `Remove Methods`)
    - change any logic regarding `issparse()` and `isspmatrix()` as needed.
    - run all your tests on the resulting code.
      You are still using spmatrix, not sparray.
      But your code and tests are prepared for the change.

  - Second Pass (switching to sparray)
    - Check all functions/methods for which migration causes 1D return values.
      These are mostly indexing and the reduction functions.
      (see below: `Shape changes and reductions`)
    - Find and change places where your code makes use of `np.matrix` features.
      Convert those to np.array features.
    - Rename all `*_matrix` constructor calls to `*_array`.
    - Test your code. And **read** your code.
      You have migrated to sparray.

### Details:

#### Construction Functions:

##### New functions
```
def block(blocks, format=None, dtype=None):
def diags_array(diagonals, /, *, offsets=0, shape=None, format=None, dtype=None):
def eye_array(m, n=None, *, k=0, dtype=float, format=None):
def random_array(m, n, density=0.01, format='coo', dtype=None, random_state=None, data_random_state=None):
```

##### Existing functions you should be careful while migrating
These functions return sparray or spmatrix dependent on their input.
Use of these should be examined and inputs adjusted to ensure return
values are sparrays. And in turn the outputs should be treated as sparrays.

```
def kron(A, B, format=None):
def kronsum(A, B, format=None):
def hstack(blocks, format=None, dtype=None):
def vstack(blocks, format=None, dtype=None):
def block_diag(mats, format=None, dtype=None):
def tril(A, k=0, format=None):
def triu(A, k=0, format=None):
```

##### Functions that you have to convert while migrating

        Function | New function
             --- | ---
        eye      | eye_array
        identity | eye_array
        diags    | diags_array
        spdiags  | diags_array
        bmat     | block
        rand     | random_array
        random   | random_array

#### Shape changes and reductions:
  - Construction using 1d-list of values:
      - `csr_matrix([1, 2, 3]).shape == (1, 3)` creates 2D matrix
      - `csr_array([1, 2, 3]).shape == (3,)` creates 1D array
  - Indexing:
    - Indexing of sparray allows 1D objects which can be made 2D using
      `np.newaxis` and/or `None`. E.g. `A[3, None, :]` gives a 2D row.
      Indexing of 2D sparray with implicit (not given) column index
      gives a 1D result, e.g. `A[3]`.  If you need a 2D result, use
      `np.newaxis`, or `None` in your index, or wrap the integer index
      as a list for which fancy indexing gives 2D, e.g. `A[[3], :]`
    - Iteration over sparse object:
        `next(M)` -> sparse 2D row matrix
        `next(A)` -> sparse 1D array
  - Reduction operations along an axis reduce the shape:
      - `M.sum(axis=1)` makes a 2D row matrix by summing along axis 1.
      - `A.sum(axis=1)` makes a 1D `coo_array` summing along axis 1.
    Some reductions return dense array/matrices instead of sparse:

        Method | result
        ---  | ---
        sum(axis)    | dense
        mean(axis)   | dense
        argmin(axis) | dense
        argmax(axis) | dense
        min(axis)    | sparse
        max(axis)    | sparse
        nanmin(axis) | sparse
        nanmax(axis) | sparse

    Generally, 2D `sparray` inputs lead to 1D results. 2D `spmatrix` inputs lead to 2D.

  - Some reductions return a scalar. Those should behave as they did before
    and shouldn't need to be considered during migration. E.g. `A.sum()`

#### Removed methods:
  - getrow, getcol, asfptype, getnnz, getH. Attributes M.A and M.H.
    It is recommended that you replace these functions with alternatives
    before starting the shift to sparray.

    Function  | Alternative
          --- | ---
    `M.get_shape()`  | `M.shape`
    `M.getformat()`  | `M.format`
    `M.asfptype(...)`| `M.astype(...)`
    `M.getmaxprint()`| `M.maxprint`
    `M.getnnz()`     | `M.nnz`
    `M.getnnz(axis)` | `M.count_nonzero(axis)`
    `M.getH()`       | `M.conj().T`
    `M.getrow(i)`    | `M[i, :]`
    `M.getcol(j)`    | `M[:, j]`
    `M.A`            | `M.toarray()`
    `M.H`            | `M.conj().T`

  - Shape assignment (`M.shape = (2, 6)`) is not permitted for sparray.
    Instead you should use `A.reshape`.

#### 
  - If you provide `(data, indices, indptr)` to a constructor, sparrays
    set the index dtype (`idxdtype`) without checking the content of
    the indices. See gh-18509
  - Binary operations with sparse and dense operands:
    Let `A` denote a `dense` object and `B` denote a `sparse` object.
    With a `sparse` and `dense` operand, some binary operations return
    `sparse` and others return `dense`. The logic is that operations
    that are likely to retain many zeros in the result return `sparse`.
    Switching `A <op> B` to `B <op> A` produces the same result type
    except for the operator `/` where sparray cannot be in the denominator.

      operator | result
         ---   |  ---
      A + B    | dense
      A - B    | dense
      A @ B    | dense
      A > B    | dense
      A != B   | dense
      A * B    | sparse
      A / B    | Not supported

   - Binary operations with two sparse operands:

      operator | result
         ---   |  ---
      B + B    | sparse
      B - B    | sparse
      B @ B    | sparse
      B > B    | sparse
      B != B   | sparse
      B * B    | sparse
      B / B    | dense
