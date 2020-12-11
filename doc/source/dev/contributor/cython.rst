:orphan:

.. _adding-cython:

Adding Cython to SciPy
======================

As written on the `Cython website`_:

 Cython is an optimising static
 compiler for both the Python programming language and the extended
 Cython programming language (based on Pyrex). It makes writing C
 extensions for Python as easy as Python itself.

If your code currently performs a lot of loops in Python, it might
benefit from compilation with Cython. This document is intended to be a
very brief introduction: just enough to see how to use Cython with
SciPy. Once you have your code compiling, you can learn more about how
to optimize it by reviewing the `Cython documentation`_.

There are only two things you need to do in order for SciPy compile your
code with Cython:

#. Include your code in a file with a ``.pyx``
   extension rather than a ``.py`` extension. All files with a ``.pyx``
   extension are automatically converted by Cython to ``.c`` files when
   SciPy is built.

#. Add an extension from this ``.c`` file to the
   configuration of the subpackage in which your code lives. Typically,
   this is very easy: add a single, formulaic line to the subpackage’s
   ``setup.py`` file. Once added as an extension, the ``.c`` code will be
   compiled by your C compiler to machine code when SciPy is built.

Example
-------

|linprog-rs|_ contains the implementation of the
revised simplex method for ``scipy.optimize.linprog``. The revised
simplex method performs many elementary row operations on matrices, and
so it was a natural candidate to be Cythonized.

Note that ``scipy/optimize/_linprog_rs.py`` imports the ``BGLU`` and
``LU`` classes from ``._bglu_dense`` exactly as if they were regular
Python classes. But they’re not. ``BGLU`` and ``LU`` are Cython classes
defined in |bglu-dense|_. There is nothing
about the way they are imported or used that suggests that they are
written in Cython; the only way so far that we can tell they are Cython
classes is that they are defined in a file with a ``.pyx`` extension.

Even in ``/scipy/optimize/_bglu_dense.pyx``, most of the code resembles
Python. The most notable differences are the presence of |cimport|_,
|cdef|_, and `Cython decorators`_. None of these are strictly
necessary. Without them, the pure Python code can still be compiled by
Cython. The Cython language extensions are \*just\* tweaks to improve
performance. This ``.pyx`` file is automatically converted to a ``.c``
file by Cython when SciPy is built.

The only thing left is to add an extension from this ``.c`` file using
|distutils|_. This takes just a single line in |optimize-setup|_:
``config.add_extension('_bglu_dense', sources=['_bglu_dense.c'])``,
``_bglu_dense.c`` is the source and ``_bglu_dense`` is the name of the
extension (for consistency). When SciPy is built, ``_bglu_dense.c`` will
be compiled to machine code, and we will be able to import the ``LU``
and ``BGLU`` classes from the extension ``_bglu_dense``.

Exercise
--------

*See a video run-through of this exercise:* \ `Cythonizing SciPy Code`_ \

#. Update Cython and create a new branch
   (e.g., ``git checkout -b cython_test``) in which to make some
   experimental changes to SciPy

#. Add some simple Python code in a ``.py`` file in the
   ``/scipy/optimize`` directory, say ``/scipy/optimize/mypython.py``.
   For example:

   ::

      def myfun():
          i = 1
          while i < 10000000:
              i += 1
          return i

#. Let’s see how long this pure-Python loop takes so we can compare the
   performance of Cython. For example, in an IPython console in Spyder:

   ::

      from scipy.optimize.mypython import myfun
      %timeit myfun()

   I get something like:

   ::

      715 ms ± 10.7 ms per loop

#. Save your ``.py`` file to a ``.pyx`` file, e.g. \ ``mycython.pyx``.

#. Build SciPy. Note that a ``.c`` file has been added to the
   ``/scipy/optimize`` directory.

#. Somewhere near similar lines, add an extension from your ``.c`` file
   to ``/scipy/optimize/setup.py``. e.g.:

   ::

      config.add_extension('_group_columns', sources=['_group_columns.c'],)  # was already here
      config.add_extension('mycython', sources=['mycython.c'],)  # this was new
      config.add_extension('_bglu_dense', sources=['_bglu_dense.c'])  # was already there

#. Rebuild SciPy. Note that a ``.so`` file has been added to the
   ``/scipy/optimize`` directory.

#. Time it:

   ::

      from scipy.optimize.mycython import myfun
      %timeit myfun()

   I get something like:

   ::

      359 ms ± 6.98 ms per loop

   Cython sped up the pure Python code by a factor of ~2.

#.  That’s not much of an improvement in the scheme of things. To see
    why, it helps to have Cython create an “annotated” version of our
    code to show bottlenecks. In a terminal window, call Cython on your
    ``.pyx`` file with the ``-a`` flag:

    ::

       cython -a scipy/optimize/mycython.pyx

    Note that this creates a new ``.html`` file in the
    ``/scipy/optimize`` directory. Open the ``.html`` file in any
    browser.

#.  The yellow-highlighted lines in the file indicate potential
    interaction between the compiled code and Python, which slows things
    down considerably. The intensity of the highlighting indicates the
    estimated severity of the interaction. In this case, much of the
    interaction can be avoided if we define the variable ``i`` as an
    integer so that Cython doesn’t have to consider the possibility of
    it being a general Python object:

    ::

       def myfun():
           cdef int i = 1  # our first line of Cython code
           while i < 10000000:
               i += 1
           return i

    Recreating the annotated ``.html`` file shows that most of the
    Python interaction has disappeared.

#. Rebuild SciPy, open an fresh IPython console, and ``%timeit``:

::

   from scipy.optimize.mycython import myfun
   %timeit myfun()

I get something like: ``68.6 ns ± 1.95 ns per loop``. The Cython code ran
about 10 million times faster than the original Python code.

In this case, the compiler probably optimized away the loop, simply
returning the final result. This sort of speedup is not typical for real
code, but this exercise certainly illustrates the power of Cython when
the alternative is many low-level operations in Python.

.. _Cython website: https://cython.org/
.. _Cython documentation: http://docs.cython.org/en/latest/

.. |cimport| replace:: ``cimport``
.. _cimport: https://cython.readthedocs.io/en/latest/src/userguide/sharing_declarations.html

.. |cdef| replace:: ``cdef``
.. _cdef: https://github.com/scipy/scipy/blob/master/scipy/optimize/setup.py

.. _Cython decorators: https://cython.readthedocs.io/en/latest/src/userguide/numpy_tutorial.html

.. |linprog-rs| replace:: ``scipy.optimize._linprog_rs.py``
.. _linprog-rs: https://github.com/scipy/scipy/blob/master/scipy/optimize/_linprog_rs.py

.. |bglu-dense| replace:: ``/scipy/optimize/_bglu_dense.pyx``
.. _bglu-dense: https://github.com/scipy/scipy/blob/master/scipy/optimize/_bglu_dense.pyx

.. |distutils| replace:: ``numpy.distutils``
.. _distutils: https://docs.scipy.org/doc/numpy/reference/distutils.html

.. |optimize-setup| replace:: ``scipy/optimize/setup.py``
.. _optimize-setup: https://github.com/scipy/scipy/blob/master/scipy/optimize/setup.py

.. _Cythonizing SciPy Code: https://youtu.be/K9bF7cjUJ7c
