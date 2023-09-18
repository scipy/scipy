---
layout: home
title: The libsndfile Home Page
---

Libsndfile is a C library for reading and writing files containing sampled sound
(such as MS Windows WAV and the Apple/SGI AIFF format) through one standard
library interface. It is released in source code format under the
[Gnu Lesser General Public License](http://www.gnu.org/copyleft/lesser.html).

The library was written to compile and run on a Linux system but should compile
and run on just about any Unix (including MacOS X).
There are also pre-compiled binaries available for 32 and 64 bit windows.

It was designed to handle both little-endian (such as WAV) and big-endian (such
as AIFF) data, and to compile and run correctly on little-endian (such as Intel
and DEC/Compaq Alpha) processor systems as well as big-endian processor systems
such as Motorola 68k, Power PC, MIPS and Sparc. Hopefully the design of the
library will also make it easy to extend for reading and writing new sound file
formats.

It has been compiled and tested (at one time or another) on the following
systems:

* Every platform supported by Debian GNU/Linux including x86_64-linux-gnu,
  i486-linux-gnu, powerpc-linux-gnu, sparc-linux-gnu, alpha-linux-gnu,
  mips-linux-gnu and armel-linux-gnu.
* powerpc-apple-darwin7.0 (Mac OS X 10.3)
* sparc-sun-solaris2.8 (using gcc)
* mips-sgi-irix5.3 (using gcc)
* QNX 6.0
* i386-unknown-openbsd2.9
* Microsoft Windows

At the moment, each new release is being tested on i386 Linux, x86_64 Linux,
PowerPC Linux, Win32 and Win64.

## Features

libsndfile has the following main features :

* Ability to read and write a large number of [file formats](formats.html).
* A simple, elegant and easy to use Applications Programming Interface.
* Usable on Unix, Win32, MacOS and others.
* On the fly format conversion, including endian-ness swapping, type conversion
  and bitwidth scaling.
* Optional normalisation when reading floating point data from files containing
  integer data.
* Ability to open files in read/write mode.
* The ability to write the file header without closing the file (only on files
  open for write or read/write).
* Ability to query the library about all supported formats and retrieve text
  strings describing each format.

libsndfile has a comprehensive test suite so that each release is as bug free
as possible.
When new bugs are found, new tests are added to the test suite to ensure that
these bugs don't creep back into the code.
When new features are added, tests are added to the test suite to make sure that
these features continue to work correctly even when they are old features.

## History

My first attempt at reading and writing WAV files was in 1990 or so under
Windows 3.1. I started using Linux in early 1995 and contributed some code to
the [wavplay](http://www.vaxxine.com/ve3wwg/gnuwave.html) program. That
contributed code would eventually mutate into this library. As one of my
interests is Digital Signal Processing (DSP) I decided that as well as reading
data from an audio file in the native format (typically 16 bit short integers)
it would also be useful to be able to have the library do the conversion to
floating point numbers for DSP applications. It then dawned on me that whatever
file format (anything from 8 bit unsigned chars, to 32 bit floating point
numbers) the library should be able to convert the data to whatever format the
library user wishes to use it in. For example, in a sound playback program, the
library caller typically wants the sound data in 16 bit short integers to dump
into a sound card even though the data in the file may be 32 bit floating point
numbers (ie Microsoft's WAVE_FORMAT_IEEE_FLOAT format). Another example would be
someone doing speech recognition research who has recorded some speech as a 16
bit WAV file but wants to process it as double precision floating point numbers.

Here is the release history for libsndfile:

* Version 0.0.8 (Feb 15 1999) First official release.
* Version 0.0.28 (Apr 26 2002) Final release of version 0 of libsndfile.
* Version 1.0.0rc1 (Jun 24 2002) Release candidate 1 of version 1 of libsndfile.
* Version 1.0.0rc6 (Aug 14 2002) MacOS 9 fixes.
* Version 1.0.0 (Aug 16 2002) First 1.0.X release.
* Version 1.0.1 (Sep 14 2002) Added MAT4 and MAT5 file formats.
* Version 1.0.2 (Nov 24 2002) Added VOX ADPCM format.
* Version 1.0.3 (Dec 09 2002) Fixes for Linux on ia64 CPUs.
* Version 1.0.4 (Feb 02 2003) New file formats and functionality.
* Version 1.0.5 (May 03 2003) One new file format and new functionality.
* Version 1.0.6 (Feb 08 2004) Large file fix for Linux/Solaris, new
  functionality   and Win32 improvements.
* Version 1.0.7 (Feb 24 2004) Fix build problems on MacOS X and fix ia64/MIPS
  etc clip mode detection.
* Version 1.0.8 (Mar 14 2004) Minor bug fixes.
* Version 1.0.9 (Mar 30 2004) Add AVR format. Improve handling of some WAV
  files.
* Version 1.0.10 (Jun 15 2004) Minor bug fixes. Fix support for Win32 MinGW
  compiler.
* Version 1.0.11 (Nov 15 2004) Add SD2 file support, reading of loop data in WAV
  and AIFF. Minor bug fixes.
* Version 1.0.12 (Sep 30 2005) Add FLAC and CAF file support, virtual I/O
  interface. Minor bug fixes and cleanups.
* Version 1.0.13 (Jan 21 2006) Add read/write of instrument chunks. Minor bug
  fixes.
* Version 1.0.14 (Feb 19 2006) Minor bug fixes. Start shipping windows
  binary/source ZIP.
* Version 1.0.15 (Mar 16 2006) Minor bug fixes.
* Version 1.0.16 (Apr 30 2006) Add support for RIFX. Other minor feature
  enhancements and bug fixes.
* Version 1.0.17 (Aug 31 2006) Add C++ wrapper sndfile.hh. Minor bug fixes and
  cleanups.
* Version 1.0.18 (Feb 07 2009) Add Ogg/Vorbis support, remove captive
  libraries, many new features and bug fixes. Generate Win32 and Win64
  pre-compiled binaries.
* Version 1.0.19 (Mar 02 2009) Fix for CVE-2009-0186. Huge number of minor fixes
  as a result of static analysis.
* Version 1.0.20 (May 14 2009) Fix for potential heap overflow.
* Version 1.0.21 (December 13 2009) Bunch of minor bug fixes.
* Version 1.0.22 (October 04 2010) Bunch of minor bug fixes.
* Version 1.0.23 (October 10 2010) Minor bug fixes.
* Version 1.0.24 (March 23 2011) Minor bug fixes.
* Version 1.0.25 (July 13 2011) Fix for Secunia Advisory SA45125. Minor bug
  fixes and improvements.
* Version 1.0.26 (November 22 2015) Fix for CVE-2014-9496, CVE-2014-9756 and
  CVE-2015-7805. Add ALAC/CAF support. Minor bug fixes and improvements.
* Version 1.0.27 (June 19 2016) Fix a seek regression in 1.0.26. Add metadata
  read/write for CAF and RF64. FIx PAF endian-ness issue.
* Version 1.0.28 (April 2 2017) Fix buffer overruns in FLAC and ID3 handling
  code. Reduce default header memory requirements. Fix detection of Large File
  Support for 32 bit systems.
* Version 1.0.29 (August 15 2020) Opus support, build system improvements and
  bug fixes.
* Version 1.0.30 (September 19 2020) Bugfix release. Fix file descriptor leaks
  in sf_open_fd () function. Fix critical CMake bug leading to broken ABI on
  Linux platforms. Other numerous fixes to CMake build system, consider it
  stable now. Fix some memory leaks. Fix handling of some SD2 files. Update
  documentation. Integrate GitHub Actions for faster test builds and Oss-Fuzz
  for fuzzing tests. Move sndfile.h.in from src/ to include/ directory. To avoid
  problems, delete old generated sndfile.h from $(top_builddir)/src.
* Version 1.0.31 (January 24 2021) Bugfix release. Fix multiple memory leaks
  reported by OSS-Fuzz. More SSE2-optimized functions for x86 and amd64.
* Version 1.1.0 (March 27 2022) Minor release, backward compatible with previous
  releases. Added long-awaited MP3 support. Numerous improvements and bugfixes.
* Version 1.2.0 (December 25 2022) Various bugfixes,
  removed artificial samplerate limit
* Version 1.2.1 (August 12 2023) Patch release, various bugfixes.
* Version 1.2.2 (August 13 2023) Patch release, various bugfixes.

## Similar or Related Projects

* [SoX](http://sox.sourceforge.net/) is a program for converting between sound
  file formats.
* [Wavplay](http://www.hitsquad.com/smm/programs/WavPlay/) started out as a
  minimal WAV file player under Linux and has mutated into Gnuwave, a
  client/server application for more general multimedia and games sound
  playback.
* [Audiofile](http://www.68k.org/~michael/audiofile/) (libaudiofile) is a
  library similar to libsndfile but with a different programming interface. The
  author Michael Pruett has set out to clone (and fix some bugs in) the
  libaudiofile library which ships with SGI's IRIX OS.
* [sndlib.tar.gz](ftp://ccrma-ftp.stanford.edu/pub/Lisp/sndlib.tar.gz) is
  another library written by Bill Schottstaedt of CCRMA.

## Licensing

libsndfile is released under the terms of the GNU Lesser General Public License,
of which there are two versions;
[version 2.1](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
and
[version 3](http://www.gnu.org/copyleft/lesser.html).
To maximise the compatibility of libsndfile, the user may choose to use
libsndfile under either of the above two licenses.
You can also read a simple explanation of the ideas behind the GPL and the LGPL
[here](http://www.gnu.org/copyleft/lesser.html).

You can use libsndfile with
[Free Software](http://www.gnu.org/),
[Open Source](http://www.opensource.org/),
proprietary, shareware or other closed source applications as long as libsndfile
is used as a dynamically loaded library and you abide by a small number of other
conditions (read the LGPL for more info).
With applications released under the GNU GPL you can also use libsndfile
statically linked to your application.

I would like to see libsndfile used as widely as possible but I would prefer it
if you released software that uses libsndfile as
[Free Software](http://www.gnu.org/)
or
[Open Source](http://www.opensource.org/).
However, if you put in a great deal of effort building a significant application
which simply uses libsndfile for file I/O, then I have no problem with you
releasing that as closed source and charging as much money as you want for it as
long as you abide by [the license](http://www.gnu.org/copyleft/lesser.html).

## Download

Check latest version on
[GitHub Releases page](https://github.com/libsndfile/libsndfile/releases/).

Binaries and source packages are signed by current releaser David Seifert aka
@SoapGentoo. You can verify signatures with his public GPG key:

```
-----BEGIN PGP PUBLIC KEY BLOCK-----
Version: GnuPG v2

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
=K3x+
-----END PGP PUBLIC KEY BLOCK-----
```

## See Also

* [sndfile-tools](https://github.com/libsndfile/sndfile-tools): a small
collection of programs which use libsndfile.
