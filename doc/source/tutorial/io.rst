File IO (:mod:`scipy.io`)
=========================

.. sectionauthor:: Matthew Brett

.. currentmodule:: scipy.io

.. seealso:: `NumPy IO routines <https://www.numpy.org/devdocs/reference/routines.io.html>`__

MATLAB files
------------

.. autosummary::

   loadmat
   savemat
   whosmat

The basic functions
```````````````````

We'll start by importing :mod:`scipy.io` and calling it ``sio`` for
convenience:

   >>> import scipy.io as sio

If you are using IPython, try tab-completing on ``sio``. Among the many
options, you will find::

   sio.loadmat
   sio.savemat
   sio.whosmat

These are the high-level functions you will most likely use when working
with MATLAB files. You'll also find::

   sio.matlab

This is the package from which ``loadmat``, ``savemat``, and ``whosmat``
are imported. Within ``sio.matlab``, you will find the ``mio`` module
This module contains the machinery that ``loadmat`` and ``savemat`` use.
From time to time you may find yourself re-using this machinery.

How do I start?
```````````````

You may have a ``.mat`` file that you want to read into SciPy. Or, you
want to pass some variables from SciPy / NumPy into MATLAB.

To save us using a MATLAB license, let's start in Octave_. Octave has
MATLAB-compatible save and load functions. Start Octave (``octave`` at
the command line for me):

.. sourcecode:: octave

  octave:1> a = 1:12
  a =

     1   2   3   4   5   6   7   8   9  10  11  12

  octave:2> a = reshape(a, [1 3 4])
  a =

  ans(:,:,1) =

     1   2   3

  ans(:,:,2) =

     4   5   6

  ans(:,:,3) =

     7   8   9

  ans(:,:,4) =

     10   11   12

  octave:3> save -6 octave_a.mat a % MATLAB 6 compatible
  octave:4> ls octave_a.mat
  octave_a.mat

Now, to Python:

  >>> mat_contents = sio.loadmat('octave_a.mat')
  >>> mat_contents
  {'a': array([[[  1.,   4.,   7.,  10.],
          [  2.,   5.,   8.,  11.],
          [  3.,   6.,   9.,  12.]]]),
   '__version__': '1.0',
   '__header__': 'MATLAB 5.0 MAT-file, written by
   Octave 3.6.3, 2013-02-17 21:02:11 UTC',
   '__globals__': []}
  >>> oct_a = mat_contents['a']
  >>> oct_a
  array([[[  1.,   4.,   7.,  10.],
          [  2.,   5.,   8.,  11.],
          [  3.,   6.,   9.,  12.]]])
  >>> oct_a.shape
  (1, 3, 4)

Now let's try the other way round:

  >>> import numpy as np
  >>> vect = np.arange(10)
  >>> vect.shape
  (10,)
  >>> sio.savemat('np_vector.mat', {'vect':vect})

Then back to Octave:

.. sourcecode:: octave

   octave:8> load np_vector.mat
   octave:9> vect
   vect =

     0  1  2  3  4  5  6  7  8  9

   octave:10> size(vect)
   ans =

       1   10

If you want to inspect the contents of a MATLAB file without reading the
data into memory, use the ``whosmat`` command:

   >>> sio.whosmat('octave_a.mat')
   [('a', (1, 3, 4), 'double')]

``whosmat`` returns a list of tuples, one for each array (or other object)
in the file. Each tuple contains the name, shape and data type of the
array.

MATLAB structs
``````````````

MATLAB structs are a little bit like Python dicts, except the field
names must be strings. Any MATLAB object can be a value of a field. As
for all objects in MATLAB, structs are, in fact, arrays of structs, where
a single struct is an array of shape (1, 1).

.. sourcecode:: octave

   octave:11> my_struct = struct('field1', 1, 'field2', 2)
   my_struct =
   {
     field1 =  1
     field2 =  2
   }

   octave:12> save -6 octave_struct.mat my_struct

We can load this in Python:

   >>> mat_contents = sio.loadmat('octave_struct.mat')
   >>> mat_contents
   {'my_struct': array([[([[1.0]], [[2.0]])]],
         dtype=[('field1', 'O'), ('field2', 'O')]), '__version__': '1.0', '__header__': 'MATLAB 5.0 MAT-file, written by Octave 3.6.3, 2013-02-17 21:23:14 UTC', '__globals__': []}
   >>> oct_struct = mat_contents['my_struct']
   >>> oct_struct.shape
   (1, 1)
   >>> val = oct_struct[0,0]
   >>> val
   ([[1.0]], [[2.0]])
   >>> val['field1']
   array([[ 1.]])
   >>> val['field2']
   array([[ 2.]])
   >>> val.dtype
   dtype([('field1', 'O'), ('field2', 'O')])

In the SciPy versions from 0.12.0, MATLAB structs come back as NumPy
structured arrays, with fields named for the struct fields. You can see
the field names in the ``dtype`` output above. Note also:

   >>> val = oct_struct[0,0]

and:

.. sourcecode:: octave

  octave:13> size(my_struct)
  ans =

     1   1

So, in MATLAB, the struct array must be at least 2-D, and we replicate
that when we read into SciPy. If you want all length 1 dimensions
squeezed out, try this:

   >>> mat_contents = sio.loadmat('octave_struct.mat', squeeze_me=True)
   >>> oct_struct = mat_contents['my_struct']
   >>> oct_struct.shape
   ()

Sometimes, it's more convenient to load the MATLAB structs as Python
objects rather than NumPy structured arrays - it can make the access
syntax in Python a bit more similar to that in MATLAB.  In order to do
this, use the ``struct_as_record=False`` parameter setting to ``loadmat``.

   >>> mat_contents = sio.loadmat('octave_struct.mat', struct_as_record=False)
   >>> oct_struct = mat_contents['my_struct']
   >>> oct_struct[0,0].field1
   array([[ 1.]])

``struct_as_record=False`` works nicely with ``squeeze_me``:

   >>> mat_contents = sio.loadmat('octave_struct.mat', struct_as_record=False, squeeze_me=True)
   >>> oct_struct = mat_contents['my_struct']
   >>> oct_struct.shape # but no - it's a scalar
   Traceback (most recent call last):
     File "<stdin>", line 1, in <module>
   AttributeError: 'mat_struct' object has no attribute 'shape'
   >>> type(oct_struct)
   <class 'scipy.io.matlab.mio5_params.mat_struct'>
   >>> oct_struct.field1
   1.0

Saving struct arrays can be done in various ways. One simple method is
to use dicts:

   >>> a_dict = {'field1': 0.5, 'field2': 'a string'}
   >>> sio.savemat('saved_struct.mat', {'a_dict': a_dict})

loaded as:

.. sourcecode:: octave

   octave:21> load saved_struct
   octave:22> a_dict
   a_dict =

     scalar structure containing the fields:

       field2 = a string
       field1 =  0.50000

You can also save structs back again to MATLAB (or Octave in our case)
like this:

   >>> dt = [('f1', 'f8'), ('f2', 'S10')]
   >>> arr = np.zeros((2,), dtype=dt)
   >>> arr
   array([(0.0, ''), (0.0, '')],
         dtype=[('f1', '<f8'), ('f2', 'S10')])
   >>> arr[0]['f1'] = 0.5
   >>> arr[0]['f2'] = 'python'
   >>> arr[1]['f1'] = 99
   >>> arr[1]['f2'] = 'not perl'
   >>> sio.savemat('np_struct_arr.mat', {'arr': arr})

MATLAB cell arrays
``````````````````

Cell arrays in MATLAB are rather like Python lists, in the sense that
the elements in the arrays can contain any type of MATLAB object. In
fact, they are most similar to NumPy object arrays, and that is how we
load them into NumPy.

.. sourcecode:: octave

   octave:14> my_cells = {1, [2, 3]}
   my_cells =
   {
     [1,1] =  1
     [1,2] =

        2   3

   }

   octave:15> save -6 octave_cells.mat my_cells

Back to Python:

   >>> mat_contents = sio.loadmat('octave_cells.mat')
   >>> oct_cells = mat_contents['my_cells']
   >>> print(oct_cells.dtype)
   object
   >>> val = oct_cells[0,0]
   >>> val
   array([[ 1.]])
   >>> print(val.dtype)
   float64

Saving to a MATLAB cell array just involves making a NumPy object array:

   >>> obj_arr = np.zeros((2,), dtype=np.object)
   >>> obj_arr[0] = 1
   >>> obj_arr[1] = 'a string'
   >>> obj_arr
   array([1, 'a string'], dtype=object)
   >>> sio.savemat('np_cells.mat', {'obj_arr':obj_arr})

.. sourcecode:: octave

   octave:16> load np_cells.mat
   octave:17> obj_arr
   obj_arr =
   {
     [1,1] = 1
     [2,1] = a string
   }

IDL files
---------

.. autosummary::

   readsav

Matrix Market files
-------------------

.. autosummary::

   mminfo
   mmread
   mmwrite

Wav sound files (:mod:`scipy.io.wavfile`)
-----------------------------------------

.. currentmodule:: scipy.io.wavfile

.. autosummary::

   read
   write

Arff files (:mod:`scipy.io.arff`)
---------------------------------

.. currentmodule:: scipy.io.arff

.. autosummary::

   loadarff

Netcdf
------

.. currentmodule:: scipy.io

.. autosummary::

   netcdf_file

Allows reading of  NetCDF files (version of pupynere_ package)

.. _pupynere: https://pypi.org/project/pupynere/
.. _octave: https://www.gnu.org/software/octave
.. _matlab: https://www.mathworks.com/
