Building SciPy with CLion Nova
==============================

CLion Nova now supports meson based projects. This guide shows how to setup CLion Nova to work with SciPy. After
following this guide you will have access to CLion Nova's syntax highlighting, code analysis and refactoring tools.

.. warning::
   Jetbrains offers two versions of CLion: CLion and CLion Nova. At the moment only CLion Nova offers the required
   features that allow us to setup a mixed Python project with meson.

.. warning::
   This guide is still a work in progress. It has only been tested on Linux and may not work in other environments.


Installing CLion Nova
---------------------

The recommended way to manage Jetbrains (the creator of CLion Nova) products is to install `Jetbrains Toolbox`_. Once
the toolbox is installed select "CLion Nova" and click the "install" button.


Assumed build setup
-------------------

This guide makes assumptions about how you are building SciPy. In particular, it assumes that you have followed the
Basic Workflow in the "dev_quickstart.rst" section. As a summary this is::

    git clone git@github.com:YOURUSERNAME/scipy.git scipy
    cd scipy
    git remote add upstream https://github.com/scipy/scipy.git
    git submodule update --init

Finding the meson build configuration settings
----------------------------------------------

This guide will require changing several meson variables within CLion Nova. These variables can be configured at::

    File > Settings > Build, Execution, Deployment > Meson

Pointing CLion Nova at meson
----------------------------

The conda environment we have installed above installs it's own version of meson. We need to point CLion Nova at this
version. To find where meson is located execute the command::

    (scipy-dev) kai@kai:~/Projects/scipy$ which meson
    /home/kai/miniconda3/envs/scipy-dev/bin/meson

Now that we know where meson is we need to change the "Meson executable" variable to be the meson path specified above.

Pointing meson at the venv's compilers
--------------------------------------

CLion Nova executes the meson command in the global environment. This means meson cannot find the correct executables
from our virtual environment. Instead we have to manually point to them. This can be done in two ways:

    1. By setting an environment variable for each missing executable
    2. By setting a `native-file` for meson

In my case I chose to go with the native file. This file requires::

    (scipy-dev) kai@kai:~/Projects/scipy$ cat local-config.ini
    [binaries]
    c = '/home' / 'kai' / 'miniconda3' / 'envs' / 'scipy-dev' / 'bin' / 'gcc'
    cpp = '/home' / 'kai' / 'miniconda3' / 'envs' / 'scipy-dev' / 'bin' / 'g++'
    cython = '/home' / 'kai' / 'miniconda3' / 'envs' / 'scipy-dev' / 'bin' / 'cython'
    fortran = '/home' / 'kai' / 'miniconda3' / 'envs' / 'scipy-dev' / 'bin' / 'gfortran'
    pythran = '/home' / 'kai' / 'miniconda3' / 'envs' / 'scipy-dev' / 'bin' / 'pythran'
    pkg-config = '/home' / 'kai' / 'miniconda3' / 'envs' / 'scipy-dev' / 'bin' / 'pkg-config'

Now we need to add the file to our meson settings. To the `setup options` box add::

    --native-file local-config.ini

Testing that meson is configured correctly
------------------------------------------

meson should now be configured correctly. Open the top level `meson.build` file in SciPy. CLion Nova will now load the
meson project and show the results in the build window.

References
----------

.. _Jetbrains Toolbox: https://www.jetbrains.com/toolbox-app/
