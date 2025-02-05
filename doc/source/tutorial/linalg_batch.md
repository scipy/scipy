---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.4
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
orphan: true
---

+++ {"tags": ["jupyterlite_sphinx_strip"]}

```{eval-rst}
.. notebooklite:: linalg_batch.md
   :new_tab: True
```

(linalg_batch)=

+++

# Batched Linear Operations

Almost all of SciPy's linear algebra functions now support N-dimensional array input. These operations have not been mathematically generalized to higher-order tensors; rather, the indicated operation is performed on a *batch* (or "stack") of input scalars, vectors, and/or matrices.

Consider the `linalg.det` function, which maps a matrix to a scalar.

```{code-cell} ipython3
import numpy as np
from scipy import linalg
A = np.eye(3)
linalg.det(A)
```

Sometimes we need the determinant of a batch of matrices of the same dimensionality.

```{code-cell} ipython3
batch = [i*np.eye(3) for i in range(1, 4)]
batch
```

We could perform the operation for each element of the batch in a loop or list comprehension:

```{code-cell} ipython3
[linalg.det(A) for A in batch]
```

However, just as we might use NumPy broadcasting and vectorization rules to create the batch of matrices in the first place:

```{code-cell} ipython3
i = np.arange(1, 4).reshape(-1, 1, 1)
batch = i * np.eye(3)
batch
```

we might also wish to perform the determinant operation on all of the matrices in one function call.

```{code-cell} ipython3
linalg.det(batch)
```

In SciPy, we prefer the term "batch" instead of "stack" because the idea is generalized to N-dimensional batches. Suppose the input is a 2 x 4 batch of 3 x 3 matrices.

```{code-cell} ipython3
batch_shape = (2, 4)
i = np.arange(np.prod(batch_shape)).reshape(*batch_shape, 1, 1)
input = i * np.eye(3)
```

In this case, we say that the *batch shape* is `(2, 4)`, and the *core shape* of the input is `(3, 3)`. The net shape of the input is the sum (concatenation) of the batch shape and core shape.

```{code-cell} ipython3
input.shape
```

Since each 3 x 3 matrix is converted to a zero-dimensional scalar, we say that the core shape of the outuput is `()`. The shape of the output is the sum of the batch shape and core shape, so the result is a 2 x 4 array.

```{code-cell} ipython3
output = linalg.det(input)
output
```

```{code-cell} ipython3
output.shape
```

Not all linear algebra functions map to scalars. For instance, the {func}`scipy.linalg.expm` function maps from a matrix to a matrix with the same shape.

```{code-cell} ipython3
A = np.eye(3)
linalg.expm(A)
```

In this case, the core shape of the output is `(3, 3)`, so with a batch shape of `(2, 4)`, we expect an output of shape `(2, 4, 3, 3)`.

```{code-cell} ipython3
output = linalg.expm(input)
output.shape
```

Generalization of these rules to functions with multiple inputs and outputs is straightforward. For instance, the {func}`scipy.linalg.eig` function produces two outputs by default, a vector and a matrix.

```{code-cell} ipython3
evals, evecs = linalg.eig(A)
evals.shape, evecs.shape
```

In this case, the core shape of the output vector is `(3,)` and the core shape of the output matrix is `(3, 3)`. The shape of each output is the batch shape plus the core shape as before.

```{code-cell} ipython3
evals, evecs = linalg.eig(input)
evals.shape, evecs.shape
```

When there is more than one input, there is no complication if the input shapes are identical.

```{code-cell} ipython3
evals, evecs = linalg.eig(input, b=input)
evals.shape, evecs.shape
```

The rules when the shapes are not identical follow logically. Each input can have its own batch shape as long as the shapes are broadcastable according to [NumPy's broadcasting rules](#array-broadcasting-in-numpy). The net batch shape is the broadcasted shape of the individual batch shapes, and the shape of each output is the net batch shape plus its core shape.

```{code-cell} ipython3
rng = np.random.default_rng(2859239482)

# Define input core shapes
m = 3
core_shape_a = (m, m)
core_shape_b = (m, m)

# Define broadcastable batch shapes
batch_shape_a = (2, 4)
batch_shape_b = (5, 1, 4)

# Define output core shapes
core_shape_evals = (m,)
core_shape_evecs = (m, m)

# Predict shapes of outputs: broadcast batch shapes,
# and append output core shapes
net_batch_shape = np.broadcast_shapes(batch_shape_a, batch_shape_b)
output_shape_evals = net_batch_shape + core_shape_evals
output_shape_evecs = net_batch_shape + core_shape_evecs
output_shape_evals, output_shape_evecs
```

```{code-cell} ipython3
# Check predictions
input_a = rng.random(batch_shape_a + core_shape_a)
input_b = rng.random(batch_shape_b + core_shape_b)
evals, evecs = linalg.eig(input_a, b=input_b)
evals.shape, evecs.shape
```

There are a few functions for which the core dimensionality (i.e., the length of the core shape) of an argument or output can be either 1 or 2. In these cases, the core dimensionality is taken to be 1 if the array has only one dimension and 2 if the array has two or more dimensions. For instance, consider the following calls to {func}`scipy.linalg.solve`. The simplest case is a single square matrix `A` and a single vector `b`:

```{code-cell} ipython3
A = np.eye(5)
b = np.arange(5)
linalg.solve(A, b)
```

In this case, the core dimensionality of `A` is 2 (shape `(5, 5)`), the core dimensionality of `b` is 1  (shape `(5,)`), and the core dimensionality of the output is 1  (shape `(5,)`).

However, `b` can also be a two-dimensional array in which the *columns* are taken to be one-dimensional vectors.

```{code-cell} ipython3
b = np.empty((5, 2))
b[:, 0] = np.arange(5)
b[:, 1] = np.arange(5, 10)
linalg.solve(A, b)
```

```{code-cell} ipython3
b.shape
```

At first glance, it might seem that the core shape of `b` is still `(5,)`, and we have simply performed the operation with a batch shape of `(2,)`. However, if this were the case, the batch shape of `b` would be *prepended* to the core shape, resulting in `b` and the output having shape `(2, 5)`. Thinking more carefully, it is correct to consider the core dimensionality of both inputs and the output to be 2; the batch shape is `()`.

Likewise, whenever `b` has more than two dimensions, the core dimensionality of `b` and the output is considered to be 2. For example, to solve a batch of three entirely separate linear systems, each with only one right hand side, `b` must be provided as a three-dimensional array: one dimensions for the batch shape (`(3,)`) and two for the core shape (`(5, 1)`).

```{code-cell} ipython3
A = rng.random((3, 5, 5))
b = rng.random((3, 5, 1))  # batch shape (3,), core shape (5, 1)
linalg.solve(A, b).shape
```
