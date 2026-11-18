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
templated Cython code, and that beyond the ``numpy`` dependency there are no
special compilation flags used nor extra include directories. So the next step
is to look at ``build.ninja``. Open that file in an editor and search for
``_decomp_update``.
You will find this set of generic and target-specific rules that apply (note,
comments in this code block are not present in ``build.ninja`` but only added
in this doc section to explain what is happening):

.. note that Pygments doesn't support Ninja syntax, so using Bash as an
   approximation here.

.. code-block:: bash

    # These rules are usually not needed to understand the problem, but can be looked up at the top of the file:
    rule c_COMPILER
     command = /home/username/mambaforge/envs/scipy-dev/bin/x86_64-conda-linux-gnu-cc $ARGS -MD -MQ $out -MF $DEPFILE -o $out -c $in
     deps = gcc
     depfile = $DEPFILE_UNQUOTED
     description = Compiling C object $out

    rule c_LINKER
     command = /home/username/mambaforge/envs/scipy-dev/bin/x86_64-conda-linux-gnu-cc $ARGS -o $out $in $LINK_ARGS
     description = Linking target $out

    # step 1: `.pyx.in` to `.pyx` code generation with Tempita
    build scipy/linalg/_decomp_update.pyx: CUSTOM_COMMAND ../scipy/linalg/_decomp_update.pyx.in | /home/username/code/scipy/tools/building/tempita.py
     COMMAND = /home/username/mambaforge/envs/scipy-dev/bin/python3.14 /home/username/code/scipy/tools/building/tempita.py ../scipy/linalg/_decomp_update.pyx.in -o scipy/linalg
     description = Generating$ scipy/linalg/_decomp_update$ with$ a$ custom$ command

    # step 2: `.pyx` to `.c` compilation with Cython
    build scipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so.p/_decomp_update.c: CUSTOM_COMMAND scipy/linalg/_decomp_update.pyx | /home/username/mambaforge/envs/scipy-dev/bin/cython scipy/__init__.py scipy/linalg/__init__.pxd scipy/linalg/__init__.py scipy/linalg/cython_blas.pxd
     DESC = Generating$ 'scipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so.p/_decomp_update.c'
     COMMAND = /home/username/mambaforge/envs/scipy-dev/bin/cython -3 --fast-fail --output-file scipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so.p/_decomp_update.c --include-dir . scipy/linalg/_decomp_update.pyx -Xfreethreading_compatible=True --shared=scipy._cyutility

    # step 3: use C compiler to go from `.c` to object file (`.o`)
    build scipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so.p/meson-generated__decomp_update.c.o: c_COMPILER scipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so.p/_decomp_update.c
     DEPFILE = scipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so.p/meson-generated__decomp_update.c.o.d
     DEPFILE_UNQUOTED = scipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so.p/meson-generated__decomp_update.c.o.d
     ARGS = -Iscipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so.p -Iscipy/linalg -I../scipy/linalg -I/home/username/mambaforge/envs/scipy-dev/lib/python3.14/site-packages/numpy/_core/include -I/home/username/mambaforge/envs/scipy-dev/include/python3.14 -fvisibility=hidden -fdiagnostics-color=always -D_FILE_OFFSET_BITS=64 -Wall -Winvalid-pch -std=c17 -O2 -g -Wno-unused-but-set-variable -Wno-unused-function -Wno-conversion -Wno-misleading-indentation -fPIC -DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION -DCYTHON_CCOMPLEX=0

    # step 4: generate a symbol file (uses `meson --internal symbolextractor`); you can safely ignore this step
    build scipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so.p/_decomp_update.cpython-314-x86_64-linux-gnu.so.symbols: SHSYM scipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so
     IMPLIB = scipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so

    # step 5: link the `.o` file to obtain the file extension module (`.so`)
    build scipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so: c_LINKER scipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so.p/meson-generated__decomp_update.c.o
     LINK_ARGS = -L/home/username/mambaforge/envs/scipy-dev/lib -Wl,--as-needed -Wl,--allow-shlib-undefined -shared -fPIC -lm -Wl,-O2 -Wl,--sort-common -Wl,--as-needed -Wl,-z,relro -Wl,-z,now -Wl,--allow-shlib-undefined -Wl,-rpath,/home/username/mambaforge/envs/scipy-dev/lib -Wl,-rpath-link,/home/username/mambaforge/envs/scipy-dev/lib -Wl,--version-script=/home/username/code/scipy/tools/building/link-version-pyinit.map

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
                    "/home/username/mambaforge/envs/scipy-dev/bin/python3.14",
                    "/home/username/code/scipy/tools/building/tempita.py",
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
        "dependencies": [],
        "depends": [],
        "installed": false
    },
    {
        "name": "_decomp_update.cpython-314-x86_64-linux-gnu",
        "id": "b4ac6f0@@_decomp_update.cpython-314-x86_64-linux-gnu.so@sha",
        "type": "shared module",
        "defined_in": "/home/username/code/scipy/scipy/linalg/meson.build",
        "filename": [
            "/home/username/code/scipy/build/scipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so"
        ],
        "build_by_default": true,
        "target_sources": [
            {
                "language": "c",
                "machine": "host",
                "compiler": [
                    "/home/username/mambaforge/envs/scipy-dev/bin/x86_64-conda-linux-gnu-cc"
                ],
                "parameters": [
                    "-I/home/username/code/scipy/build/scipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so.p",
                    "-I/home/username/code/scipy/build/scipy/linalg",
                    "-I/home/username/code/scipy/scipy/linalg",
                    "-I/home/username/mambaforge/envs/scipy-dev/lib/python3.14/site-packages/numpy/_core/include",
                    "-I/home/username/mambaforge/envs/scipy-dev/include/python3.14",
                    "-fvisibility=hidden",
                    "-fdiagnostics-color=always",
                    "-D_FILE_OFFSET_BITS=64",
                    "-Wall",
                    "-Winvalid-pch",
                    "-std=c17",
                    "-O2",
                    "-g",
                    "-Wno-unused-but-set-variable",
                    "-Wno-unused-function",
                    "-Wno-conversion",
                    "-Wno-misleading-indentation",
                    "-fPIC",
                    "-DNPY_NO_DEPRECATED_API=NPY_1_9_API_VERSION",
                    "-DCYTHON_CCOMPLEX=0"
                ],
                "sources": [],
                "generated_sources": [
                    "/home/username/code/scipy/build/scipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so.p/_decomp_update.c"
                ],
                "unity_sources": []
            },
            {
                "linker": [
                    "/home/username/mambaforge/envs/scipy-dev/bin/x86_64-conda-linux-gnu-cc"
                ],
                "parameters": [
                    "-L/home/username/mambaforge/envs/scipy-dev/lib",
                    "-Wl,--as-needed",
                    "-Wl,--allow-shlib-undefined",
                    "-shared",
                    "-fPIC",
                    "-lm",
                    "-Wl,-O2",
                    "-Wl,--sort-common",
                    "-Wl,--as-needed",
                    "-Wl,-z,relro",
                    "-Wl,-z,now",
                    "-Wl,--allow-shlib-undefined",
                    "-Wl,-rpath,/home/username/mambaforge/envs/scipy-dev/lib",
                    "-Wl,-rpath-link,/home/username/mambaforge/envs/scipy-dev/lib",
                    "-Wl,--version-script=/home/username/code/scipy/tools/building/link-version-pyinit.map"
                ]
            }
        ],
        "extra_files": [],
        "subproject": null,
        "dependencies": [
            "dep17362272261413320709411132636391318505",
            "numpy",
            "python-3.14"
        ],
        "depends": [],
        "installed": true,
        "install_filename": [
            "/home/username/code/scipy/build-install/lib/python3.14/site-packages/scipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so"
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

      "/home/username/code/scipy/build/scipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so":{
         "destination":"{py_platlib}/scipy/linalg/_decomp_update.cpython-314-x86_64-linux-gnu.so",
         "tag":"runtime",
         "subproject":null,
         "install_rpath":null,
         "build_rpaths":[]
      },

We may also be interested in knowing what dependencies were detected by the
configure stage of the build. So we look in ``intro-dependencies.json``
(abbreviated here - the full file also lists ``threads``, ``cblas`` and
``lapack``, and has a few more keys per entry):

.. code-block:: json

    [
       {
          "name":"python-3.14",
          "type":"pkgconfig",
          "version":"3.14",
          "compile_args":[
             "-I/home/username/mambaforge/envs/scipy-dev/include/python3.14"
          ],
          "link_args":[

          ],
          "meson_variables":[]
       },
       {
          "name":"numpy",
          "type":"config-tool",
          "version":"2.4.4",
          "compile_args":[
             "-I/home/username/mambaforge/envs/scipy-dev/lib/python3.14/site-packages/numpy/_core/include"
          ],
          "link_args":[

          ],
          "meson_variables":[
             "_numpy_dep"
          ]
       },
       {
          "name":"blas",
          "type":"pkgconfig",
          "version":"0.3.33",
          "compile_args":[
             "-I/home/username/mambaforge/envs/scipy-dev/include"
          ],
          "link_args":[
             "/home/username/mambaforge/envs/scipy-dev/lib/libopenblas.so"
          ],
          "meson_variables":[
             "blas"
          ]
       }
    ]

This tells us which dependencies were found, and how. Note the ``type`` field:
``numpy`` is found through the ``numpy-config`` executable ("config-tool"),
while ``blas`` and ``python`` are found through pkg-config. The
``meson_variables`` field is handy for mapping an entry back to the
``meson.build`` code that looked it up.
