:orphan:

.. _other-languages:

=============
Beyond Python
=============

This is a small collection of thoughts related to the inclusion of code written
in languages other than Python. Currently, the only option for languages other
than Python that is explicitly documented is :ref:`Cython<adding-cython>`.

*Can I use a programming language other than Python to speed up my code?*

Yes.  The languages used in SciPy are Python, Cython, C, C++ and Fortran.  All
of these have their pros and cons.  If Python really doesn't offer enough
performance, one of those languages can be used.  Important concerns when
using compiled languages are maintainability and portability.  For
maintainability, Cython is clearly preferred over C/C++/Fortran.  Cython and C
are more portable than C++/Fortran.  A lot of the existing C and Fortran code
in SciPy is older, battle-tested code that was only wrapped in (but not
specifically written for) Python/SciPy.  Therefore the basic advice is: use
Cython.  If there's specific reasons why C/C++/Fortran should be preferred,
please discuss those reasons first.

*Can I use Numba or Pythran?*

Not yet, but we're considering it for the future.


*How do I debug code written in C/C++/Fortran inside SciPy?*

The easiest way to do this is to first write a Python script that
invokes the C code whose execution you want to debug. For instance
``mytest.py``::

    from scipy.special import hyp2f1
    print(hyp2f1(5.0, 1.0, -1.8, 0.95))

Now, you can run::

    gdb --args python runtests.py -g --python mytest.py

If you didn't compile with debug symbols enabled before, remove the
``build`` directory first. While in the debugger::

    (gdb) break cephes_hyp2f1
    (gdb) run

The execution will now stop at the corresponding C function and you
can step through it as usual. Instead of plain ``gdb`` you can of
course use your favorite alternative debugger; run it on the
``python`` binary with arguments ``runtests.py -g --python mytest.py``.
