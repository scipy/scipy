# Third Party Licenses and Acknowledgements

The majority of the HiGHS source code is available under the [MIT license](https://opensource.org/license/MIT).

Code in the `/extern` directory was originally developed by third-parties and is
licensed under additional licenses.

## amd

The source code in `/extern/amd` is distributed under the [BSD-3 license](https://opensource.org/license/bsd-3-clause)
at `/extern/amd/License.txt`.

It was originally developed by Timothy Davis.

The upstream source code is available at:

 * https://github.com/DrTimothyAldenDavis/SuiteSparse

To avoid compiling this code into HiGHS, use `-DHIPO=OFF` (the default value).

## cli

The source code in `/extern/CLI11.hpp` is distributed under a license stated in the source code. Although the license is not named, the text appears to correspond to the [BSD-3 license](https://opensource.org/license/bsd-3-clause).

The upstream source code is available at:

* https://github.com/CLIUtils/CLI11

CLI11 is only used to parse command line input for the HiGHS executable, so does not affect the license status of the HiGHS library or language interfaces to it.

## metis

The source code in `/extern/metis` is distributed under the [Apache 2.0 license](https://opensource.org/license/apache-2-0)
at `/extern/metis/LICENSE.txt`.

It was originally developed by George Karypis.

The upstream source code is available at:

 * https://github.com/KarypisLab/METIS
 * https://github.com/KarypisLab/GKlib

To avoid compiling this code into HiGHS, use `-DHIPO=OFF` (the default value).

## pdqsort

The source code in `/extern/pdqsort` is distributed under the [zlib license](https://opensource.org/license/zlib)
at `/extern/pdqsort/license.txt`.

It was originally developed by Orson Peters.

The upstream source code is available at:

 * https://github.com/orlp/pdqsort

## rcm

The source code in `/extern/rcm` is distributed under the [MIT license](https://opensource.org/license/MIT)
at `/extern/rcm/LICENSE`.

It was originally developed by Alan George, Joseph Liu, and John Burkardt.

The upstream source code is available at:

 * https://people.sc.fsu.edu/~jburkardt/cpp_src/rcm/rcm.html

To avoid compiling this code into HiGHS, use `-DHIPO=OFF` (the default value).

## zstr

The source code in `/extern/zstr` is distributed under the [MIT license](https://opensource.org/license/MIT)
at `/extern/zstr/LICENSE`.

It was originally developed by Matei David.

The upstream source code is available at:

 * https://github.com/mateidavid/zstr
