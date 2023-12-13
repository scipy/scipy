.. _sparsetutorial:

Sparse Arrays (`scipy.sparse`)
========================================================

.. sectionauthor:: Levi John Wolf <levi.john.wolf@gmail.com>

.. currentmodule:: scipy.sparse

Introduction
-------------

``scipy.sparse`` and its submodules provide tools for working with *sparse arrays*. Sparse arrays are arrays where only a few locations in the array have any data, most of the locations are considered as "empty". Sparse arrays are useful because they allow for simpler, faster, and/or less memory-intensive algorithms for linear algebra (`scipy.sparse.linalg`) or graph-based computations (`scipy.sparse.csgraph`), but they are generally less flexible for operations like slicing, reshaping, or assignment. This guide will introduce the basics of sparse arrays in `scipy.sparse`, explain the unique aspects of sparse data structures, and refer onward for other sections of the user guide explaining `sparse linear algebra <https://docs.scipy.org/doc/scipy/tutorial/arpack.html>`_ and `graph methods <https://docs.scipy.org/doc/scipy/tutorial/csgraph.html>`_. 

Getting started with sparse arrays
-----------------------------------

Sparse arrays are a special kind of array where only a few locations in the array have data. This allows for *compressed* representations of the data to be used, where only the locations where data exists are recorded. There are many different sparse array formats, each of which makes a different tradeoff between compression and functionality. To start, let's build a very simple sparse array, the Coordinate (COO) array (:func:`coo_array`) and compare it to a dense array:

   >>> import scipy as sp
   >>> import numpy
   >>> dense = numpy.array([[1, 0, 0, 2], [0, 4, 1, 0], [0, 0, 5, 0]])
   >>> sparse = sp.sparse.coo_array(dense)
   >>> dense
   array([[1, 0, 0, 2],
       [0, 4, 1, 0],
       [0, 0, 5, 0]])
   >>> sparse
   <3x4 sparse array of type '<class 'numpy.int64'>'
         with 5 stored elements in COOrdinate format>

Note that in our dense array, we have five nonzero values. For example, ``2`` is at location ``0,3``, and ``4`` is at location ``1,1``. All of the other values are zero. The sparse array records these five values *explicitly* (see the ``5 stored elements in COOrdinate format``), and then represents all of the remaining zeros as *implicit* values. 

Most sparse array methods work in a similar fashion to dense array methods: 

   >>> sparse.max()
   5
   >>> dense.max()
   5
   >>> sparse.argmax()
   10
   >>> dense.argmax()
   10
   >>> sparse.mean()
   1.0833333333333333
   >>> dense.mean()
   1.0833333333333333

A few "extra" properties, such as ``.nnz`` which returns the number of stored values, are present on sparse arrays as well: 

   >>> sparse.nnz
   5 

Most of the reduction operations, such as ``.mean()``, ``.sum()``, or ``.max()`` will return a numpy array when applied over an axis of the sparse array: 

   >>> sparse.mean(axis=1)
   array([0.75, 1.25, 1.25])

This is because reductions over sparse arrays are often dense.  

Understanding sparse array formats
-----------------------------------

Different kinds of sparse arrays have different capabilities. For example, COO arrays cannot be subscripted or sliced:  
   
   >>> dense[2, 2]
   5
   >>> sparse[2, 2]
   Traceback (most recent call last):
     File "<stdin>", line 1, in <module>
   TypeError: 'coo_array' object is not subscriptable

But, other formats, such as the Compressed Sparse Row (CSR) :func:`csr_array()` support slicing and element indexing:
   
   >>> sparse.tocsr()[2, 2]
   5

Sometimes, `scipy.sparse` will return a different sparse matrix format than the input sparse matrix format. For example, the dot product of two sparse arrays in COO format will be a CSR format array: 

   >>> sparse @ sparse.T
   <3x3 sparse array of type '<class 'numpy.int64'>'
        with 5 stored elements in Compressed Sparse Row format>

This change occurs because `scipy.sparse` will change the format of input sparse arrays in order to use the most efficient computational method. 

The `scipy.sparse` module contains the following formats, each with their own distinct advantages and disadvantages: 

- Block Sparse Row (BSR) arrays :func:`scipy.sparse.bsr_array()`, which are most appropriate when the parts of the array with data occur in contiguous blocks.  
- Coordinate (COO) arrays :func:`scipy.sparse.coo_array()`, which provide a simple way to construct sparse arrays and modify them in place. COO can also be quickly converted into other formats, such CSR, CSC, or BSR. 
- Compressed Sparse Row (CSR) arrays :func:`scipy.sparse.csr_array()`, which are most useful for fast arithmetic, vector products, and slicing by row. 
- Compressed Sparse Column (CSC) arrays :func:`scipy.sparse.csc_array()`, which are most useful for fast arithmetic, vector products, and slicing by column. 
- Diagonal (DIA) arrays :func:`scipy.sparse.dia_array()`, which are useful for efficient storage and fast arithmetic so long as the data primarily occurs along diagonals of the array.
- Dictionary of Keys (DOK) arrays :func:`scipy.sparse.dok_array()`, which are useful for fast construction and single-element access.
- List of Lists (LIL) arrays :func:`scipy.sparse.lil_array()`, which are useful for fast construction and modification of sparse arrays. 

More information on the strengths and weaknesses of each of the sparse array formats can be found in `their documentation <https://docs.scipy.org/doc/scipy/reference/sparse.html#sparse-array-classes>`_.

All formats of `scipy.sparse` arrays can be constructed directly from a `numpy.ndarray`. However, some sparse formats can be constructed in different ways, too. Each sparse array format has different strengths, and these strengths are documented in each class. For example, one of the most common methods for constructing sparse arrays is to build a sparse array from the individual ``row``, ``column``, and ``data`` values. For our array from before: 

   >>> dense
   array([[1, 0, 0, 2],
       [0, 4, 1, 0],
       [0, 0, 5, 0]])

The ``row``, ``column``, and ``data`` arrays describe the rows, columns, and values where our sparse array has entries: 

   >>> row = [0,0,1,1,2]
   >>> col = [0,3,1,2,2]
   >>> data = [1,2,4,1,5]

Using these, we can now define a sparse array without building a dense array first: 

   >>> csr = sp.sparse.csr_array((data, (row, col)))
   >>> csr
   <3x4 sparse array of type '<class 'numpy.int64'>'
        with 5 stored elements in Compressed Sparse Row format>

Different classes have different constructors, but the :func:`scipy.sparse.csr_array`, :func:`scipy.sparse.csc_array`, and :func:`scipy.sparse.coo_array` allow for this style of construction. 

Sparse arrays, implicit zeros, and duplicates
----------------------------------------------

Sparse arrays are useful because they represent much of their values *implicitly*, without storing an actual placeholder value. In `scipy.sparse`, the value used to represent "no data" is an *implicit zero*. This can be confusing when *explicit zeros* are required. For example, in `graph methods <https://docs.scipy.org/doc/scipy/tutorial/csgraph.html>`_ from `scipy.sparse.csgraph`, we often need to be able to distinguish between (A) a link connecting nodes ``i`` and ``j`` with zero weight and (B) no link between ``i`` and ``j``. Sparse matrices can do this, so long as we keep the *explicit* and *implicit* zeros in mind. 

For example, in our previous ``csr`` array, we could include an explicit zero by including it in the ``data`` list. Let's treat the final entry in the array at the bottom row and last column as an *explicit zero*: 

   >>> row = [0,0,1,1,2,2]
   >>> col = [0,3,1,2,2,3]
   >>> data = [1,2,4,1,5,0]

Then, our sparse array will have *six* stored elements, not five: 

   >>> csr = sp.sparse.csr_array((data, (row, col)))
   >>> csr
   <3x4 sparse array of type '<class 'numpy.int64'>'
        with 6 stored elements in Compressed Sparse Row format>

The "extra" element is our *explicit zero*. The two are still identical when converted back into a dense array, because dense arrays represent *everything* explicitly: 

   >>> csr.todense()
   array([[1, 0, 0, 2],
       [0, 4, 1, 0],
       [0, 0, 5, 0]])
   >>> dense
   array([[1, 0, 0, 2],
       [0, 4, 1, 0],
       [0, 0, 5, 0]])

But, for sparse arithmetic, linear algebra, and graph methods, the value at ``2,3`` will be considered an *explicit zero*. To remove this explicit zero, we can use the ``csr.eliminate_zeros()`` method. This operates on the sparse array *in place*, and removes any zero-value stored elements: 

   >>> csr
   <3x4 sparse array of type '<class 'numpy.int64'>'
        with 6 stored elements in Compressed Sparse Row format>
   >>> csr.eliminate_zeros()
   >>> csr
   <3x4 sparse array of type '<class 'numpy.int64'>'
        with 5 stored elements in Compressed Sparse Row format>

Before ``csr.eliminate_zeros()``, there were six stored elements. After, there are only five stored elements. 

Another point of complication arises from how *duplicates* are processed when constructing a sparse array. A *duplicate* can occur when we have two or more entries at ``row,col`` when constructing a sparse array. This often occurs when building sparse arrays using the ``data``, ``row``, and ``col`` vectors. For example, we might represent our previous array with a duplicate value at ``1,1``:

   >>> row = [0,0,1,1,1,2]
   >>> col = [0,3,1,1,2,2]
   >>> data = [1,2,1,3,1,5]

In this case, we can see that there are *two* ``data`` values that correspond to the ``1,1`` location in our final array. `scipy.sparse` will store these values separately: 

   >>> dupes = sp.sparse.coo_array((data, (row, col)))
   >>> dupes
   <3x4 sparse array of type '<class 'numpy.int64'>'
        with 6 stored elements in COOrdinate format>

Note that there are six stored elements in this sparse array, despite only having five unique locations where data occurs. When these arrays are converted back to dense arrays, the duplicate values are summed. So, at location ``1,1``, the dense array will contain the sum of duplicate stored entries, ``1 + 3``: 

   >>> dupes.todense()
   array([[1, 0, 0, 2],
         [0, 4, 1, 0],
         [0, 0, 5, 0]])

To remove duplicate values within the sparse array itself and thus reduce the number of stored elements, we can use the ``.sum_duplicates()`` method:

   >>> dupes.sum_duplicates()
   >>> dupes
   <3x4 sparse array of type '<class 'numpy.int64'>'
        with 5 stored elements in COOrdinate format>

Now there are only five stored elements in our sparse array, and it is identical to the array we have been working with throughout this guide: 
   
   >>> dupes.todense()
   array([[1, 0, 0, 2],
          [0, 4, 1, 0],
          [0, 0, 5, 0]])

Canonical formats
-----------------

Several sparse array formats have "canonical formats" to allow for more efficient operations.
Generally these consist of added restrictions like:

- No duplicate entries for any value
- Sorted indices

Classes with a canonical form include: :func:`coo_array`, :func:`csr_array`, :func:`csc_array`, and :func:`bsr_array`.
See the docstrings of these classes for details on each canonical representation.

To check if an instance of these classes is in canonical form, use the ``.has_canonical_format`` attribute:

   >>> coo = sp.sparse.coo_array(([1, 1, 1], ([0, 2, 1], [0, 1, 2])))
   >>> coo.has_canonical_format
   False

To convert an instance to canonical form, use the ``.sum_duplicates()`` method:

   >>> coo.sum_duplicates()
   >>> coo.has_canonical_format
   True

Next steps with sparse arrays 
-------------------------------

Sparse array types are most helpful when working with large, nearly empty arrays. Specifically, `sparse linear algebra <https://docs.scipy.org/doc/scipy/tutorial/arpack.html>`_ and `sparse graph methods <https://docs.scipy.org/doc/scipy/tutorial/csgraph.html>`_ see the largest improvements in efficiency in these circumstances. 

