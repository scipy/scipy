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

```{eval-rst}
.. jupyterlite:: ../_contents/linalg_batch.ipynb
   :new_tab: True
```

(linalg_batch)=
# Batched Linear Operations

Some of SciPy's linear algebra functions support N-dimensional array input. These operations have not been mathematically generalized to higher-order tensors; rather, the indicated operation is performed on a *batch* (or "stack") of input scalars, vectors, and/or matrices.

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
