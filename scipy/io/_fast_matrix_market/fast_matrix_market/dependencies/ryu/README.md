By Ulf Adams

See https://github.com/ulfjack/ryu

# Ryu

This project contains routines to convert IEEE-754 floating-point numbers to decimal strings using shortest, fixed %f, and scientific %e formatting. The primary implementation is in C, and there is a port of the shortest conversion to Java. All algorithms have been published in peer-reviewed publications. At the time of this writing, these are the fastest known float-to-string conversion algorithms. The fixed, and scientific conversion routines are several times faster than the usual implementations of sprintf (we compared against glibc, Apple's libc, MSVC, and others).

