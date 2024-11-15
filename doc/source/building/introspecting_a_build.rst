.. _meson-introspection:

Introspecting build steps
=========================

When you have an issue with a particular Python extension module or other build
target, there are a number of ways to figure out what the build system is doing
exactly. Beyond looking at the ``meson.build`` content for the target of
interest, these include:

1. Reading the generated ``build.ninja`` file in the build directory,
2. Using ``meson introspect`` to learn more about build options, dependencies
   and flags used for the target,
3. Reading ``<build-dir>/meson-info/*.json`` for details on discovered
   dependencies, where Meson plans to install files to, etc.

These things are all available after the configure stage of the build (i.e.,
``meson setup``) has run. It is typically more effective to look at this
information, rather than running the build and reading the full build log.


The ``ninja.build`` file
------------------------

As an example, let's say we are interested in ``scipy.linalg._decomp_update``.
From ``scipy/linalg/meson.build`` we learn that this extension is written in
templated Cython code, and there are no special compilation flags used nor
include directories beyond the ``numpy`` one. So the next step is to look at
``build.ninja``. Open that file in an editor and search for ``_decomp_update``.
You will find this set of generic and target-specific rules that apply (note,
comments in this code block are not present in ``build.ninja`` but only added
in this doc section to explain what is happening):

.. note that Pygments doesn't support Ninja syntax, so using Bash as an
   approximation here.

.. code-block:: bash

    # These rules are usually not needed to understand the problem, but can be looked up at the top of the file:
    rule c_COMPILER
     command = /home/username/anaconda3/envs/scipy-dev/bin/x86_64-conda-linux-gnu-cc $ARGS -MD -MQ $out -MF $DEPFILE -o $out -c $in
     deps = gcc
     depfile = $DEPFILE_UNQUOTED
     description = Compiling C object $out

    rule c_LINKER
     command = /home/username/anaconda3/envs/scipy-dev/bin/x86_64-conda-linux-gnu-cc $ARGS -o $out $in $LINK_ARGS
     description = Linking target $out

    # step 1: `.pyx.in` to `.pyx` code generation with Tempita
    build scipy/linalg/_decomp_update.pyx: CUSTOM_COMMAND ../scipy/linalg/_decomp_update.pyx.in | ../scipy/_build_utils/tempita.py /home/username/anaconda3/envs/scipy-dev/bin/python3.10
     COMMAND = /home/username/anaconda3/envs/scipy-dev/bin/python3.10 ../scipy/_build_utils/tempita.py ../scipy/linalg/_decomp_update.pyx.in -o scipy/linalg
     description = Generating$ scipy/linalg/_decomp_update$ with$ a$ custom$ command

    # step 2: `.pyx` to `.c` compilation with Cython
    build scipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so.p/_decomp_update.c: CUSTOM_COMMAND scipy/linalg/_decomp_update.pyx | /home/username/code/scipy/scipy/_build_utils/cythoner.py scipy/__init__.py scipy/linalg/__init__.py scipy/linalg/cython_blas.pyx
     DESC = Generating$ 'scipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so.p/_decomp_update.c'.
     COMMAND = /home/username/anaconda3/envs/scipy-dev/bin/python3.10 /home/username/code/scipy/scipy/_build_utils/cythoner.py scipy/linalg/_decomp_update.pyx scipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so.p/_decomp_update.c

    # step 3: use C compiler to go from `.c` to object file (`.o`)
    build scipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so.p/meson-generated__decomp_update.c.o: c_COMPILER scipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so.p/_decomp_update.c
     DEPFILE = scipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so.p/meson-generated__decomp_update.c.o.d
     DEPFILE_UNQUOTED = scipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so.p/meson-generated__decomp_update.c.o.d
     ARGS = -Iscipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so.p -Iscipy/linalg -I../scipy/linalg -I/home/username/anaconda3/envs/scipy-dev/lib/python3.10/site-packages/numpy/core/include -I/home/username/anaconda3/envs/scipy-dev/include/python3.10 -fvisibility=hidden -fdiagnostics-color=always -D_FILE_OFFSET_BITS=64 -Wall -Winvalid-pch -std=c99 -O2 -g -Wno-unused-but-set-variable -Wno-unused-function -Wno-conversion -Wno-misleading-indentation -fPIC -Wno-cpp

    # step 4: generate a symbol file (uses `meson --internal symbolextractor`); you can safely ignore this step
    build scipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so.p/_decomp_update.cpython-310-x86_64-linux-gnu.so.symbols: SHSYM scipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so
     IMPLIB = scipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so

    # step 5: link the `.o` file to obtain the file extension module (`.so`)
    build scipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so: c_LINKER scipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so.p/meson-generated__decomp_update.c.o | /home/username/anaconda3/envs/scipy-dev/x86_64-conda-linux-gnu/sysroot/lib64/libm-2.12.so /home/username/anaconda3/envs/scipy-dev/x86_64-conda-linux-gnu/sysroot/usr/lib64/libm.a
     LINK_ARGS = -L/home/username/anaconda3/envs/scipy-dev/lib -Wl,--as-needed -Wl,--allow-shlib-undefined -shared -fPIC -Wl,--start-group -lm -Wl,--end-group -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--disable-new-dtags -Wl,--gc-sections -Wl,--allow-shlib-undefined -Wl,-rpath,/home/username/anaconda3/envs/scipy-dev/lib -Wl,-rpath-link,/home/username/anaconda3/envs/scipy-dev/lib

Using ``meson introspect``
--------------------------

If we want to look at ``_decomp_update`` from another perspective, we can use
(for example) ``meson introspect --targets -i <build-dir> > targets.json`` to
generate readable JSON. Searching that generated file for our target of
interest shows:

.. code-block:: json

    {
        "name": "_decomp_update",
        "id": "b4ac6f0@@_decomp_update@cus",
        "type": "custom",
        "defined_in": "/home/username/code/scipy/scipy/linalg/meson.build",
        "filename": [
            "/home/username/code/scipy/build/scipy/linalg/_decomp_update.pyx"
        ],
        "build_by_default": false,
        "target_sources": [
            {
                "language": "unknown",
                "compiler": [
                    "/home/username/anaconda3/envs/scipy-dev/bin/python3.10",
                    "/home/username/code/scipy/scipy/_build_utils/tempita.py",
                    "@INPUT@",
                    "-o",
                    "@OUTDIR@"
                ],
                "parameters": [],
                "sources": [
                    "/home/username/code/scipy/scipy/linalg/_decomp_update.pyx.in"
                ],
                "generated_sources": []
            }
        ],
        "extra_files": [],
        "subproject": null,
        "installed": false
    },
    {
        "name": "_decomp_update.cpython-310-x86_64-linux-gnu",
        "id": "b4ac6f0@@_decomp_update.cpython-310-x86_64-linux-gnu@sha",
        "type": "shared module",
        "defined_in": "/home/username/code/scipy/scipy/linalg/meson.build",
        "filename": [
            "/home/username/code/scipy/build/scipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so"
        ],
        "build_by_default": true,
        "target_sources": [
            {
                "language": "c",
                "compiler": [
                    "/home/username/anaconda3/envs/scipy-dev/bin/x86_64-conda-linux-gnu-cc"
                ],
                "parameters": [
                    "-I/home/username/code/scipy/build/scipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so.p",
                    "-I/home/username/code/scipy/build/scipy/linalg",
                    "-I/home/username/code/scipy/scipy/linalg",
                    "-I/home/username/anaconda3/envs/scipy-dev/lib/python3.10/site-packages/numpy/core/include",
                    "-I/home/username/anaconda3/envs/scipy-dev/include/python3.10",
                    "-fvisibility=hidden",
                    "-fdiagnostics-color=always",
                    "-D_FILE_OFFSET_BITS=64",
                    "-Wall",
                    "-Winvalid-pch",
                    "-std=c99",
                    "-O2",
                    "-g",
                    "-Wno-unused-but-set-variable",
                    "-Wno-unused-function",
                    "-Wno-conversion",
                    "-Wno-misleading-indentation",
                    "-fPIC",
                    "-Wno-cpp"
                ],
                "sources": [],
                "generated_sources": [
                    "/home/username/code/scipy/build/scipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so.p/_decomp_update.c"
                ]
            }
        ],
        "extra_files": [],
        "subproject": null,
        "installed": true,
        "install_filename": [
            "/home/username/code/scipy/build-install/lib/python3.10/site-packages/scipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so"
        ]
    },

This tells us a lot of things, like which include directories will be used,
where the Cython-generated C code can be found, and what compile flags are
used. ``meson introspect --help`` has good documentation on the full range of
capabilities and how to use them.

``meson-info`` JSON files
-------------------------

There are a number of different JSON files in ``<build-dir>/meson-info/``.
These have descriptive names, hinting at their content. For example, where the
final ``_decomp_update`` extension gets installed to is described in
``intro-install_plan.json`` (note, these files aren't pretty-printed, running
them through a JSON formatter helps):

.. code-block:: json

      "/home/username/code/scipy/build/scipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so":{
         "destination":"{py_platlib}/scipy/linalg/_decomp_update.cpython-310-x86_64-linux-gnu.so",
         "tag":"runtime"
      },

We may also be interested in knowing what dependencies were detected by the
configure stage of the build. So we look in ``intro-dependencies.json``:

.. code-block:: json

    [
       {
          "name":"python",
          "version":"3.10",
          "compile_args":[
             "-I/home/username/anaconda3/envs/scipy-dev/include/python3.10"
          ],
          "link_args":[

          ]
       },
       {
          "name":"openblas",
          "version":"0.3.20",
          "compile_args":[
             "-I/home/username/anaconda3/envs/scipy-dev/include"
          ],
          "link_args":[
             "/home/username/anaconda3/envs/scipy-dev/lib/libopenblas.so"
          ]
       },
       {
          "name":"threads",
          "version":"unknown",
          "compile_args":[
             "-pthread"
          ],
          "link_args":[
             "-pthread"
          ]
       }
    ]

This tells us that we have three dependencies that were found. Note: ``numpy``
and a few other build-time dependencies are missing here because we do not
(yet) search for those with the builtin ``dependency()`` Meson command.
