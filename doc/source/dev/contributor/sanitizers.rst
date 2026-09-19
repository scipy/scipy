.. _scipy-sanitizers:

==========
Sanitizers
==========

Sanitizers are compiler features that allow us to add run-time checks to detect
multiple types of bugs. They usually consist of an instrumentation module, which
adds error detection mechanisms to the generated code, and a supporting run-time
library.

There are multiple kinds of sanitizers, with support varying between compilers,
but the most common is the address sanitizer — often abbreviated to "ASAN" —
which checks for memory bugs, such as buffer overflows.


Side-effects
------------

Enabling sanitizers will almost always result in very significant slowdowns
and size increases to the generated objects, due to the added code
instrumentation. It can also make the build workflow substantially more complex.

Moreover, the resulting artifacts often require a carefully setup execution
environment, as well as a substantially more complex build workflow, especially
when multiple components are used together (eg. the Python interpreter and
extension modules). Even in the correct execution environment, the run-time
operation of the sanitizers has some sharp edges, as they're intended to be a
developer tool.

Consequently, sanitizers are generally only suitable to be used during
development. Even so, in larger projects like SciPy, it's usually best not to
use them as a part of your general development workflow. Rather, only to
validate code when necessary.


:ref:`CI <continuous-integration>` Coverage
-------------------------------------------

.. list-table::

    * - Type
      - Environment
      - Toolchain
      - Job

    * - Address Sanitizer
      - macOS
      - LLVM
      - .. image:: https://img.shields.io/badge/macOS.yml-clang__ASAN-dodgerblue?logo=github
            :alt: Badge showing the "macOS | clang_ASAN" github actions job, and workflow status
            :target: https://github.com/search?q=clang_ASAN+repo:scipy/scipy+path:.github/workflows/macos.yml


.. _scipy-asan:

Address Sanitizer (abbv. "ASAN")
================================

As mentioned above, the address sanitizer helps detect memory bugs.

To setup the LLVM address sanitizer in SciPy, the following workflow is
recommended:

#. Set ``ASAN_OPTIONS=detect_leaks=0:symbolize=1:strict_init_order=true:allocator_may_return_null=1:use_sigaltstack=0``
#. Compile CPython_ with ``--with-address-sanitizer`` and use it to setup a test environment
#. Build NumPy_ with the ``setup-args`` config setting set to ``-Db_sanitize=address``
    * Eg. ``pip install numpy --no-cache-dir --no-binary numpy -Csetup-args="-Db_sanitize=address" -v``
#. In the project root directory, set the compiler flags as follows:
    * ``CFLAGS="-fsanitize=address -fno-omit-frame-pointer -fsanitize-ignorelist=$(pwd)/tools/asan-ignore.txt"``
    * ``CXXFLAGS="-fsanitize=address -fno-omit-frame-pointer -fsanitize-ignorelist=$(pwd)/tools/asan-ignore.txt"``
#. Build SciPy using ``spin build --clean``
#. Run tests using ``-s`` (``--capture=no``), otherwise the ASAN failures won't be shown

Alternatively, the cpython_sanity_ project, provides container images with
CPython_ and NumPy_ pre-built with the address sanitizer (eg.
``ghcr.io/nascheme/numpy-asan:3.14t``).

To learn more about the address sanitizer usage, refer to the
`relevant LLVM documentation <https://clang.llvm.org/docs/AddressSanitizer.html>`_.


Known Errors
------------

The ``tools/asan-ignore.txt`` file specifies address sanitizer suppressions for
known issues.

You can learn more at the
`LLVM issue suppression documentation <https://clang.llvm.org/docs/SanitizerSpecialCaseList.html>`_


Common Issues
-------------


Meson not finding NumPy
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block::

    (...)
    numpy-config found: NO
    Run-time dependency numpy found: NO (tried pkgconfig and config-tool)
    (...)

This usually means that ``detect_leaks`` is not disabled in ``ASAN_OPTIONS``.

You can verify this by running ``numpy-config``. If it fails with the address
sanitizer reporting a memory leak, like below, set the ``ASAN_OPTIONS``
environment as specified in the workflow instructions.

.. code-block:: sh

    $ numpy-config
    (...)
    SUMMARY: AddressSanitizer: 3027 byte(s) leaked in 98 allocation(s).


Meson reports a nonfunctional C/C++ compiler
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-bloack::

    meson.build:1:0: ERROR: Compiler clang-20 cannot compile programs.

This may occur if the ``asan-ignore.txt`` file passed to
``-fsanitize-ignorelist`` is not specified using a absolute path.

You can verify this in the Meson log.


THe ASAN runtime failing to intercept libc++
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some cases, such as when loading an extension module that links against
libc++, the ASAN runtime has some trouble intercepting symbols like
``__cxa_throw`` to hook into them. This then leads to the interceptor handling
code to fail.

.. code-block::

    AddressSanitizer: CHECK failed: asan_interceptors.cpp:463 "((__interception::real___cxa_throw)) != (0)" (0x0, 0x0) (tid=77864)

The best solution would be linking the extension modules against an instrumented
build of libc++, however this can make the execution environment extremely hard
to setup.

However, a commonly viable work-around is to have the linker load libc++ before
anything else.

.. code-block::

    $ export LD_PRELOAD=/usr/lib/libstdc++.so


Incompatible Runtimes
~~~~~~~~~~~~~~~~~~~~~

.. code-block::

    ==82818==Your application is linked against incompatible ASan runtimes.

This most likely means that SciPy, NumPy_, and/or CPython_ were not compiled
with same toolchain.

During SciPy development, this usually occurs due to having an outdated build
directory. Try rebuilding the project with ``--clean`` (``spin build --clean``)
— if ``--clean`` is not specified, the old build configuration will be used.


Missing ASAN Runtime
~~~~~~~~~~~~~~~~~~~~

.. code-block::

    ==270868==ASan runtime does not come first in initial library list; you should either link runtime to your application or manually preload it with LD_PRELOAD.

This usually means that your CPython_ installation was not compiled with
``--with-address-sanitizer``.


.. _CPython: https://github.com/python/cpython
.. _NumPy: https://github.com/numpy/numpy
.. _cpython_sanity: https://github.com/nascheme/cpython_sanity
