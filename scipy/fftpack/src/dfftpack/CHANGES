This documents the changes done by HCP to make fftpack into dfftpack

(1) Renamed all files corresponding to subroutines in the API
    i.e. ones documented as callable by the luser. Names chosen to match
    the ones in libsunperf.

(2) Inserted  IMPLICIT DOUBLE PRECISION (A-H,O-Z) after every
    subroutine statement. This makes everything that used to be a real 
    into a double.

(3) Replaced floating constants with Double Prec. constants. All
    0. become 0.D0 etc and PI, SQRT(2) etc. expanded to dble prec.
 
(4) Replaced DIMENSION FOO(1) with DIMENSION FOO(*) where foo
    is an array argument of a subroutine. I only did this in the places 
    where g77 notices it, so the compile looks cleaner.

(5) Replaced COMPLEX with DOUBLE COMPLEX. Now, this is not standard 
    fortran 77, so the whole thing may fall apart if you have a VERY
    vanilla Fortran 77 compiler. On the other hand, the only place a 
    complex is _declared_ as such is in the test program. If you don't have 
    DOUBLE COMPLEX my guess is that the library will work, except for 
    the routines  ZFFTI, ZFFTB and ZFFTF.

(6) Updated the file doc 