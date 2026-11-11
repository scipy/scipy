# HiGHS - Linear optimization software

<!-- ![Build Status](https://github.com/ERGO-Code/HiGHS/actions/workflows/build.yml/badge.svg) -->

[![Build Status][fast_build_svg]][fast_build_link]
[![Build Status][linux_build_svg]][linux_build_link]
[![Build Status][macos_build_svg]][macos_build_link]
[![Build Status][windows_build_svg]][windows_build_link]
\
[![Conan Center](https://img.shields.io/conan/v/highs)](https://conan.io/center/recipes/highs)
\
[![PyPi](https://img.shields.io/pypi/v/highspy.svg)](https://pypi.python.org/pypi/highspy)
[![PyPi](https://img.shields.io/pypi/dm/highspy.svg)](https://pypi.python.org/pypi/highspy)
\
[![NuGet version](https://img.shields.io/nuget/v/Highs.Native.svg)](https://www.nuget.org/packages/Highs.Native)
[![NuGet download](https://img.shields.io/nuget/dt/Highs.Native.svg)](https://www.nuget.org/packages/Highs.Native)

[fast_build_svg]: https://github.com/ERGO-Code/HiGHS/actions/workflows/build-fast.yml/badge.svg
[fast_build_link]: https://github.com/ERGO-Code/HiGHS/actions/workflows/build-fast.yml
[linux_build_svg]: https://github.com/ERGO-Code/HiGHS/actions/workflows/build-linux.yml/badge.svg
[linux_build_link]: https://github.com/ERGO-Code/HiGHS/actions/workflows/build-linux.yml
[macos_build_svg]: https://github.com/ERGO-Code/HiGHS/actions/workflows/build-macos.yml/badge.svg
[macos_build_link]: https://github.com/ERGO-Code/HiGHS/actions/workflows/build-macos.yml
[windows_build_svg]: https://github.com/ERGO-Code/HiGHS/actions/workflows/build-windows.yml/badge.svg
[windows_build_link]: https://github.com/ERGO-Code/HiGHS/actions/workflows/build-windows.yml

- [About HiGHS](#about-highs)
- [Documentation](#documentation)
- [Installation](#installation)
  - [Build from source using CMake](#build-from-source-using-cmake)
  - [Build with Meson*](#build-with-meson)
  - [Build with Nix*](#build-with-nix)
  - [Precompiled binaries](#precompiled-binaries)
- [Running HiGHS](#running-highs)
- [Interfaces](#interfaces)
  - [Python](#python)
  - [C](#c)
  - [CSharp](#csharp)
  - [Fortran](#fortran)
- [Reference](#reference)

<a id="about-highs"></a>

## About HiGHS

HiGHS is a high performance serial and parallel solver for large scale sparse
linear optimization problems of the form

$$ \min \quad \dfrac{1}{2}x^TQx + c^Tx \qquad \textrm{s.t.}~ \quad L \leq Ax \leq U; \quad l \leq x \leq u $$

where $Q$ must be positive semi-definite and, if $Q$ is zero, there may be a requirement that some of the variables take integer values. Thus HiGHS can solve linear programming (LP) problems, convex quadratic programming (QP) problems, and mixed integer programming (MIP) problems. It is mainly written in C++, but also has some C. It has been developed and tested on various Linux, MacOS and Windows installations. No third-party dependencies are required.

HiGHS has primal and dual revised simplex solvers, originally written by Qi Huangfu and further developed by Julian Hall. It also has an interior point solver for LP written by Lukas Schork, an active set solver for QP written by Michael Feldmeier, and a MIP solver written by Leona Gottwald. Other features have been added by Julian Hall and Ivet Galabova, who manages the software engineering of HiGHS and interfaces to C, C#, FORTRAN, Julia and Python.

Find out more about HiGHS at https://www.highs.dev.

Although HiGHS is freely available under the MIT license, we would be pleased to learn about users' experience and give advice via email sent to highsopt@gmail.com.

<a id="documentation"></a>

## Documentation

Documentation is available at https://ergo-code.github.io/HiGHS/.

<a id="installation"></a>

## Installation

<a id="build-from-source-using-cmake"></a>

### Build from source using CMake

HiGHS uses CMake as build system, and requires at least version 3.15. To generate build files in a new subdirectory called 'build', run:

```shell
    cmake -S . -B build
    cmake --build build
```
This installs the executable `bin/highs` and the library `lib/highs`.

To test whether the compilation was successful, change into the build directory and run

```shell
    ctest
```
More details on building with CMake can be found in `HiGHS/cmake/README.md`.

<a id="build-with-meson"></a>

#### Build with Meson

As an alternative, HiGHS can be installed using the `meson` build interface:
``` sh
meson setup bbdir -Dwith_tests=True
meson test -C bbdir
```
_The meson build files are provided by the community and are not officially supported by the HiGHS development team._ **If you use this method and encounter issues, please consider contributing fixes or updates by checking the [HiGHS Contribution Guide](https://github.com/ERGO-Code/HiGHS/blob/master/CONTRIBUTING.md).**

<a id="build-with-nix"></a>

#### Build with Nix

There is a nix flake that provides the `highs` binary:

```shell
nix run .
```

You can even run [without installing
anything](https://determinate.systems/posts/nix-run/), supposing you have
installed [nix](https://nixos.org/download.html):

```shell
nix run github:ERGO-Code/HiGHS
```

The nix flake also provides the python package:

```shell
nix build .#highspy
tree result/
```

And a devShell for testing it:

```shell
nix develop .#highspy
python
>>> import highspy
>>> highspy.Highs()
```

_The nix build files are provided by the community and are not officially supported by the HiGHS development team._

<a id="precompiled-binaries"></a>

### Precompiled binaries

Precompiled static binaries are available at https://github.com/ERGO-Code/HiGHS/releases.

Additionally, there is one package containing shared libraries for Windows x64.

The `*-mit` binary packages contain HiGHS and are MIT-licenced.
The `*-apache` binary packages contain HiGHS with HiPO and are Apache-licenced, due to the licensing of the dependencies of HiPO. For more information, see [THIRD_PARTY_NOTICES.md](https://github.com/ERGO-Code/HiGHS/blob/master/THIRD_PARTY_NOTICES.md).

If you have any questions or requests for more platforms and binaries, please get in touch with us at hello@highs.dev.

<a id="running-highs"></a>

## Running HiGHS

HiGHS can read MPS files and (CPLEX) LP files, and the following command
solves the model in `ml.mps`

```shell
    highs ml.mps
```
#### Command line options

When HiGHS is run from the command line, some fundamental option values may be
specified directly. Many more may be specified via a file. Formally, the usage
is:

```shell
$ bin/highs --help
usage:
      ./bin/highs [options] [file]

options:
      --model_file file          File of model to solve.
      --options_file file        File containing HiGHS options.
      --read_solution_file file  File of solution to read.
      --read_basis_file text     File of initial basis to read.
      --write_model_file text    File for writing out model.
      --solution_file text       File for writing out solution.
      --write_basis_file text    File for writing out final basis.
      --presolve text            Set presolve option to:
                                   "choose" * default
                                   "on"
                                   "off"
      --solver text              Set solver option to:
                                   "choose" * default
                                   "simplex"
                                   "ipm"
      --parallel text            Set parallel option to:
                                   "choose" * default
                                   "on"
                                   "off"
      --threads int              Set maximum number of threads to use:
                                 0: automatic * default
      --run_crossover text       Set run_crossover option to:
                                   "choose"
                                   "on" * default
                                   "off"
      --time_limit float         Run time limit (seconds - double).
      --random_seed int          Seed to initialize random number
                                 generation.
      --ranging text             Compute cost, bound, RHS and basic
                                 solution ranging:
                                   "on"
                                   "off" * default
  -v, --version                  Print version.
  -h, --help                     Print help.

```
For a full list of options, see the [options page](https://ergo-code.github.io/HiGHS/stable/options/definitions/) of the documentation website.

<a id="interfaces"></a>

## Interfaces

There are HiGHS interfaces for C, C#, FORTRAN, and Python in `HiGHS/highs/interfaces`, with example driver files in `HiGHS/examples/`. More on language and modelling interfaces can be found at https://ergo-code.github.io/HiGHS/stable/interfaces/other/.

We are happy to give a reasonable level of support via email sent to highsopt@gmail.com.

<a id="python"></a>

### Python

The python package `highspy` is a thin wrapper around HiGHS and is available on [PyPi](https://pypi.org/project/highspy/). It can be easily installed via `pip` by running

```shell
$ pip install highspy
```

Alternatively, `highspy` can be built from source.  Download the HiGHS source code and run

```shell
pip install .
```
from the root directory.

The python package `highspy` depends on the `numpy` package and `numpy` will be installed as well, if it is not already present.

The installation can be tested using the small example `HiGHS/examples/call_highs_from_python_highspy.py`.

The [Google Colab Example Notebook](https://colab.research.google.com/drive/1JmHF53OYfU-0Sp9bzLw-D2TQyRABSjHb?usp=sharing) also demonstrates how to call `highspy`.

To use HiPO in Python, the additional `highspy-extras` is required. It contains the HiPO dependencies and is Apache 2.0 licensed. It can be installed with

```shell
$ pip install highspy[extras]
```

For more details, see https://ergo-code.github.io/HiGHS/ at Interfaces -> Python -> HiPO in Python.

<a id="c"></a>

### C
The C API is in `HiGHS/highs/interfaces/highs_c_api.h`. It is included in the default build. For more details, check out the documentation website https://ergo-code.github.io/HiGHS/.

<a id="csharp"></a>

### CSharp

The nuget package Highs.Native is on https://www.nuget.org, at https://www.nuget.org/packages/Highs.Native/.

It can be added to your C# project with `dotnet`

```shell
dotnet add package Highs.Native --version 1.16.0
```

The nuget package contains runtime libraries for

* `win-x64`
* `win-x86`
* `linux-x64`
* `linux-arm64`
* `macos-x64`
* `macos-arm64`

Details for building locally can be found in `nuget/README.md`.

<a id="fortran"></a>

### Fortran

The Fortran API is in `HiGHS/highs/interfaces/highs_fortran_api.f90`. It is *not* included in the default build. For more details, check out the documentation website https://ergo-code.github.io/HiGHS/.


<a id="reference"></a>

## Reference

If you use HiGHS in an academic context, please acknowledge this and cite the following article.

Parallelizing the dual revised simplex method
Q. Huangfu and J. A. J. Hall
Mathematical Programming Computation, 10 (1), 119-142, 2018.
DOI: [10.1007/s12532-017-0130-5](https://link.springer.com/article/10.1007/s12532-017-0130-5)
