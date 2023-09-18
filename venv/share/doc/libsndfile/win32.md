---
layout: page
---

# Building libsndfile on Win32

**Note : For pre-compiled binaries for windows, both for win32 and win64, see
the main web page.**

There are currently two build systems; the official GNU autotool based one and
a more limited and experimental CMake based build system.

libsndfile is written to be compiled by a compiler which supports large chunks
of the 1999 ISO C Standard (tested with GCC, Clang and Visual Studio 2015).

It is recommended to use CMake and Visual Studio to build libsndfile on Windows
but you can try the [MinGW](http://www.mingw.org/) compiler suite with Autotools
or CMake buildsystem.
