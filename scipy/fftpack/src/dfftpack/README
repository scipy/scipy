                      DFFTPACK V1.0
*****************************************************************

        A Double precision clone by Hugh C. Pumphrey  of:

                      FFTPACK
               version 4  april 1985

     a package of fortran subprograms for the fast fourier
      transform of periodic and other symmetric sequences

                         by

                  paul n swarztrauber

  national center for atmospheric research  boulder,colorado 80307

   which is sponsored by the national science foundation

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


this package consists of programs which perform fast fourier
transforms for both double complex and (double precision) real
periodic sequences and certain other symmetric sequences that are
listed below.

1.   dffti     initialize  dfftf and dfftb
2.   dfftf     forward transform of a real periodic sequence
3.   dfftb     backward transform of a real coefficient array

4.   dzffti    initialize dzfftf and dzfftb
5.   dzfftf    a simplified real periodic forward transform
6.   dzfftb    a simplified real periodic backward transform

7.   dsinti     initialize dsint
8.   dsint      sine transform of a real odd sequence

9.   dcosti     initialize dcost
10.  dcost      cosine transform of a real even sequence

11.  dsinqi     initialize dsinqf and dsinqb
12.  dsinqf     forward sine transform with odd wave numbers
13.  dsinqb     unnormalized inverse of dsinqf

14.  dcosqi     initialize dcosqf and dcosqb
15.  dcosqf     forward cosine transform with odd wave numbers
16.  dcosqb     unnormalized inverse of dcosqf

17.  zffti     initialize zfftf and zfftb
18.  zfftf     forward transform of a double complex periodic sequence
19.  zfftb     unnormalized inverse of zfftf


******************************************************************

subroutine dffti(n,wsave)

  ****************************************************************

subroutine dffti initializes the array wsave which is used in
both dfftf and dfftb. the prime factorization of n together with
a tabulation of the trigonometric functions are computed and
stored in wsave.

input parameter

n       the length of the sequence to be transformed.

output parameter

wsave   a work array which must be dimensioned at least 2*n+15.
        the same work array can be used for both dfftf and dfftb
        as long as n remains unchanged. different wsave arrays
        are required for different values of n. the contents of
        wsave must not be changed between calls of dfftf or dfftb.

******************************************************************

subroutine dfftf(n,r,wsave)

******************************************************************

subroutine dfftf computes the fourier coefficients of a real
perodic sequence (fourier analysis). the transform is defined
below at output parameter r.

input parameters

n       the length of the array r to be transformed.  the method
        is most efficient when n is a product of small primes.
        n may change so long as different work arrays are provided

r       a real array of length n which contains the sequence
        to be transformed

wsave   a work array which must be dimensioned at least 2*n+15.
        in the program that calls dfftf. the wsave array must be
        initialized by calling subroutine dffti(n,wsave) and a
        different wsave array must be used for each different
        value of n. this initialization does not have to be
        repeated so long as n remains unchanged thus subsequent
        transforms can be obtained faster than the first.
        the same wsave array can be used by dfftf and dfftb.


output parameters

r       r(1) = the sum from i=1 to i=n of r(i)

        if n is even set l =n/2   , if n is odd set l = (n+1)/2

          then for k = 2,...,l

             r(2*k-2) = the sum from i = 1 to i = n of

                  r(i)*cos((k-1)*(i-1)*2*pi/n)

             r(2*k-1) = the sum from i = 1 to i = n of

                 -r(i)*sin((k-1)*(i-1)*2*pi/n)

        if n is even

             r(n) = the sum from i = 1 to i = n of

                  (-1)**(i-1)*r(i)

 *****  note
             this transform is unnormalized since a call of dfftf
             followed by a call of dfftb will multiply the input
             sequence by n.

wsave   contains results which must not be destroyed between
        calls of dfftf or dfftb.


******************************************************************

subroutine dfftb(n,r,wsave)

******************************************************************

subroutine dfftb computes the real perodic sequence from its
fourier coefficients (fourier synthesis). the transform is defined
below at output parameter r.

input parameters

n       the length of the array r to be transformed.  the method
        is most efficient when n is a product of small primes.
        n may change so long as different work arrays are provided

r       a real array of length n which contains the sequence
        to be transformed

wsave   a work array which must be dimensioned at least 2*n+15.
        in the program that calls dfftb. the wsave array must be
        initialized by calling subroutine dffti(n,wsave) and a
        different wsave array must be used for each different
        value of n. this initialization does not have to be
        repeated so long as n remains unchanged thus subsequent
        transforms can be obtained faster than the first.
        the same wsave array can be used by dfftf and dfftb.


output parameters

r       for n even and for i = 1,...,n

             r(i) = r(1)+(-1)**(i-1)*r(n)

                  plus the sum from k=2 to k=n/2 of

                   2.*r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)

                  -2.*r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)

        for n odd and for i = 1,...,n

             r(i) = r(1) plus the sum from k=2 to k=(n+1)/2 of

                  2.*r(2*k-2)*cos((k-1)*(i-1)*2*pi/n)

                 -2.*r(2*k-1)*sin((k-1)*(i-1)*2*pi/n)

 *****  note
             this transform is unnormalized since a call of dfftf
             followed by a call of dfftb will multiply the input
             sequence by n.

wsave   contains results which must not be destroyed between
        calls of dfftb or dfftf.


******************************************************************

subroutine dzffti(n,wsave)

******************************************************************

subroutine dzffti initializes the array wsave which is used in
both dzfftf and dzfftb. the prime factorization of n together with
a tabulation of the trigonometric functions are computed and
stored in wsave.

input parameter

n       the length of the sequence to be transformed.

output parameter

wsave   a work array which must be dimensioned at least 3*n+15.
        the same work array can be used for both dzfftf and dzfftb
        as long as n remains unchanged. different wsave arrays
        are required for different values of n.


******************************************************************

subroutine dzfftf(n,r,azero,a,b,wsave)

******************************************************************

subroutine dzfftf computes the fourier coefficients of a real
perodic sequence (fourier analysis). the transform is defined
below at output parameters azero,a and b. dzfftf is a simplified
but slower version of dfftf.

input parameters

n       the length of the array r to be transformed.  the method
        is must efficient when n is the product of small primes.

r       a real array of length n which contains the sequence
        to be transformed. r is not destroyed.


wsave   a work array which must be dimensioned at least 3*n+15.
        in the program that calls dzfftf. the wsave array must be
        initialized by calling subroutine dzffti(n,wsave) and a
        different wsave array must be used for each different
        value of n. this initialization does not have to be
        repeated so long as n remains unchanged thus subsequent
        transforms can be obtained faster than the first.
        the same wsave array can be used by dzfftf and dzfftb.

output parameters

azero   the sum from i=1 to i=n of r(i)/n

a,b     for n even b(n/2)=0. and a(n/2) is the sum from i=1 to
        i=n of (-1)**(i-1)*r(i)/n

        for n even define kmax=n/2-1
        for n odd  define kmax=(n-1)/2

        then for  k=1,...,kmax

             a(k) equals the sum from i=1 to i=n of

                  2./n*r(i)*cos(k*(i-1)*2*pi/n)

             b(k) equals the sum from i=1 to i=n of

                  2./n*r(i)*sin(k*(i-1)*2*pi/n)


******************************************************************

subroutine dzfftb(n,r,azero,a,b,wsave)

******************************************************************

subroutine dzfftb computes a real perodic sequence from its
fourier coefficients (fourier synthesis). the transform is
defined below at output parameter r. dzfftb is a simplified
but slower version of dfftb.

input parameters

n       the length of the output array r.  the method is most
        efficient when n is the product of small primes.

azero   the constant fourier coefficient

a,b     arrays which contain the remaining fourier coefficients
        these arrays are not destroyed.

        the length of these arrays depends on whether n is even or
        odd.

        if n is even n/2    locations are required
        if n is odd (n-1)/2 locations are required

wsave   a work array which must be dimensioned at least 3*n+15.
        in the program that calls dzfftb. the wsave array must be
        initialized by calling subroutine dzffti(n,wsave) and a
        different wsave array must be used for each different
        value of n. this initialization does not have to be
        repeated so long as n remains unchanged thus subsequent
        transforms can be obtained faster than the first.
        the same wsave array can be used by dzfftf and dzfftb.


output parameters

r       if n is even define kmax=n/2
        if n is odd  define kmax=(n-1)/2

        then for i=1,...,n

             r(i)=azero plus the sum from k=1 to k=kmax of

             a(k)*cos(k*(i-1)*2*pi/n)+b(k)*sin(k*(i-1)*2*pi/n)

********************* complex notation **************************

        for j=1,...,n

        r(j) equals the sum from k=-kmax to k=kmax of

             c(k)*exp(i*k*(j-1)*2*pi/n)

        where

             c(k) = .5*cmplx(a(k),-b(k))   for k=1,...,kmax

             c(-k) = conjg(c(k))

             c(0) = azero

                  and i=sqrt(-1)

*************** amplitude - phase notation ***********************

        for i=1,...,n

        r(i) equals azero plus the sum from k=1 to k=kmax of

             alpha(k)*cos(k*(i-1)*2*pi/n+beta(k))

        where

             alpha(k) = sqrt(a(k)*a(k)+b(k)*b(k))

             cos(beta(k))=a(k)/alpha(k)

             sin(beta(k))=-b(k)/alpha(k)

******************************************************************

subroutine dsinti(n,wsave)

******************************************************************

subroutine dsinti initializes the array wsave which is used in
subroutine dsint. the prime factorization of n together with
a tabulation of the trigonometric functions are computed and
stored in wsave.

input parameter

n       the length of the sequence to be transformed.  the method
        is most efficient when n+1 is a product of small primes.

output parameter

wsave   a work array with at least int(2.5*n+15) locations.
        different wsave arrays are required for different values
        of n. the contents of wsave must not be changed between
        calls of dsint.

******************************************************************

subroutine dsint(n,x,wsave)

******************************************************************

subroutine dsint computes the discrete fourier sine transform
of an odd sequence x(i). the transform is defined below at
output parameter x.

dsint is the unnormalized inverse of itself since a call of dsint
followed by another call of dsint will multiply the input sequence
x by 2*(n+1).

the array wsave which is used by subroutine dsint must be
initialized by calling subroutine dsinti(n,wsave).

input parameters

n       the length of the sequence to be transformed.  the method
        is most efficient when n+1 is the product of small primes.

x       an array which contains the sequence to be transformed


wsave   a work array with dimension at least int(2.5*n+15)
        in the program that calls dsint. the wsave array must be
        initialized by calling subroutine dsinti(n,wsave) and a
        different wsave array must be used for each different
        value of n. this initialization does not have to be
        repeated so long as n remains unchanged thus subsequent
        transforms can be obtained faster than the first.

output parameters

x       for i=1,...,n

             x(i)= the sum from k=1 to k=n

                  2*x(k)*sin(k*i*pi/(n+1))

             a call of dsint followed by another call of
             dsint will multiply the sequence x by 2*(n+1).
             hence dsint is the unnormalized inverse
             of itself.

wsave   contains initialization calculations which must not be
        destroyed between calls of dsint.

******************************************************************

subroutine dcosti(n,wsave)

******************************************************************

subroutine dcosti initializes the array wsave which is used in
subroutine dcost. the prime factorization of n together with
a tabulation of the trigonometric functions are computed and
stored in wsave.

input parameter

n       the length of the sequence to be transformed.  the method
        is most efficient when n-1 is a product of small primes.

output parameter

wsave   a work array which must be dimensioned at least 3*n+15.
        different wsave arrays are required for different values
        of n. the contents of wsave must not be changed between
        calls of dcost.

******************************************************************

subroutine dcost(n,x,wsave)

******************************************************************

subroutine dcost computes the discrete fourier cosine transform
of an even sequence x(i). the transform is defined below at output
parameter x.

dcost is the unnormalized inverse of itself since a call of dcost
followed by another call of dcost will multiply the input sequence
x by 2*(n-1). the transform is defined below at output parameter x

the array wsave which is used by subroutine dcost must be
initialized by calling subroutine dcosti(n,wsave).

input parameters

n       the length of the sequence x. n must be greater than 1.
        the method is most efficient when n-1 is a product of
        small primes.

x       an array which contains the sequence to be transformed

wsave   a work array which must be dimensioned at least 3*n+15
        in the program that calls dcost. the wsave array must be
        initialized by calling subroutine dcosti(n,wsave) and a
        different wsave array must be used for each different
        value of n. this initialization does not have to be
        repeated so long as n remains unchanged thus subsequent
        transforms can be obtained faster than the first.

output parameters

x       for i=1,...,n

            x(i) = x(1)+(-1)**(i-1)*x(n)

             + the sum from k=2 to k=n-1

                 2*x(k)*cos((k-1)*(i-1)*pi/(n-1))

             a call of dcost followed by another call of
             dcost will multiply the sequence x by 2*(n-1)
             hence dcost is the unnormalized inverse
             of itself.

wsave   contains initialization calculations which must not be
        destroyed between calls of dcost.

******************************************************************

subroutine dsinqi(n,wsave)

******************************************************************

subroutine dsinqi initializes the array wsave which is used in
both dsinqf and dsinqb. the prime factorization of n together with
a tabulation of the trigonometric functions are computed and
stored in wsave.

input parameter

n       the length of the sequence to be transformed. the method
        is most efficient when n is a product of small primes.

output parameter

wsave   a work array which must be dimensioned at least 3*n+15.
        the same work array can be used for both dsinqf and dsinqb
        as long as n remains unchanged. different wsave arrays
        are required for different values of n. the contents of
        wsave must not be changed between calls of dsinqf or dsinqb.

******************************************************************

subroutine dsinqf(n,x,wsave)

******************************************************************

subroutine dsinqf computes the fast fourier transform of quarter
wave data. that is , dsinqf computes the coefficients in a sine
series representation with only odd wave numbers. the transform
is defined below at output parameter x.

dsinqb is the unnormalized inverse of dsinqf since a call of dsinqf
followed by a call of dsinqb will multiply the input sequence x
by 4*n.

the array wsave which is used by subroutine dsinqf must be
initialized by calling subroutine dsinqi(n,wsave).


input parameters

n       the length of the array x to be transformed.  the method
        is most efficient when n is a product of small primes.

x       an array which contains the sequence to be transformed

wsave   a work array which must be dimensioned at least 3*n+15.
        in the program that calls dsinqf. the wsave array must be
        initialized by calling subroutine dsinqi(n,wsave) and a
        different wsave array must be used for each different
        value of n. this initialization does not have to be
        repeated so long as n remains unchanged thus subsequent
        transforms can be obtained faster than the first.

output parameters

x       for i=1,...,n

             x(i) = (-1)**(i-1)*x(n)

                + the sum from k=1 to k=n-1 of

                2*x(k)*sin((2*i-1)*k*pi/(2*n))

             a call of dsinqf followed by a call of
             dsinqb will multiply the sequence x by 4*n.
             therefore dsinqb is the unnormalized inverse
             of dsinqf.

wsave   contains initialization calculations which must not
        be destroyed between calls of dsinqf or dsinqb.

******************************************************************

subroutine dsinqb(n,x,wsave)

******************************************************************

subroutine dsinqb computes the fast fourier transform of quarter
wave data. that is , dsinqb computes a sequence from its
representation in terms of a sine series with odd wave numbers.
the transform is defined below at output parameter x.

dsinqf is the unnormalized inverse of dsinqb since a call of dsinqb
followed by a call of dsinqf will multiply the input sequence x
by 4*n.

the array wsave which is used by subroutine dsinqb must be
initialized by calling subroutine dsinqi(n,wsave).


input parameters

n       the length of the array x to be transformed.  the method
        is most efficient when n is a product of small primes.

x       an array which contains the sequence to be transformed

wsave   a work array which must be dimensioned at least 3*n+15.
        in the program that calls dsinqb. the wsave array must be
        initialized by calling subroutine dsinqi(n,wsave) and a
        different wsave array must be used for each different
        value of n. this initialization does not have to be
        repeated so long as n remains unchanged thus subsequent
        transforms can be obtained faster than the first.

output parameters

x       for i=1,...,n

             x(i)= the sum from k=1 to k=n of

               4*x(k)*sin((2k-1)*i*pi/(2*n))

             a call of dsinqb followed by a call of
             dsinqf will multiply the sequence x by 4*n.
             therefore dsinqf is the unnormalized inverse
             of dsinqb.

wsave   contains initialization calculations which must not
        be destroyed between calls of dsinqb or dsinqf.

******************************************************************

subroutine dcosqi(n,wsave)

******************************************************************

subroutine dcosqi initializes the array wsave which is used in
both dcosqf and dcosqb. the prime factorization of n together with
a tabulation of the trigonometric functions are computed and
stored in wsave.

input parameter

n       the length of the array to be transformed.  the method
        is most efficient when n is a product of small primes.

output parameter

wsave   a work array which must be dimensioned at least 3*n+15.
        the same work array can be used for both dcosqf and dcosqb
        as long as n remains unchanged. different wsave arrays
        are required for different values of n. the contents of
        wsave must not be changed between calls of dcosqf or dcosqb.

******************************************************************

subroutine dcosqf(n,x,wsave)

******************************************************************

subroutine dcosqf computes the fast fourier transform of quarter
wave data. that is , dcosqf computes the coefficients in a cosine
series representation with only odd wave numbers. the transform
is defined below at output parameter x

dcosqf is the unnormalized inverse of dcosqb since a call of dcosqf
followed by a call of dcosqb will multiply the input sequence x
by 4*n.

the array wsave which is used by subroutine dcosqf must be
initialized by calling subroutine dcosqi(n,wsave).


input parameters

n       the length of the array x to be transformed.  the method
        is most efficient when n is a product of small primes.

x       an array which contains the sequence to be transformed

wsave   a work array which must be dimensioned at least 3*n+15
        in the program that calls dcosqf. the wsave array must be
        initialized by calling subroutine dcosqi(n,wsave) and a
        different wsave array must be used for each different
        value of n. this initialization does not have to be
        repeated so long as n remains unchanged thus subsequent
        transforms can be obtained faster than the first.

output parameters

x       for i=1,...,n

             x(i) = x(1) plus the sum from k=2 to k=n of

                2*x(k)*cos((2*i-1)*(k-1)*pi/(2*n))

             a call of dcosqf followed by a call of
             cosqb will multiply the sequence x by 4*n.
             therefore dcosqb is the unnormalized inverse
             of dcosqf.

wsave   contains initialization calculations which must not
        be destroyed between calls of dcosqf or dcosqb.

******************************************************************

subroutine dcosqb(n,x,wsave)

******************************************************************

subroutine dcosqb computes the fast fourier transform of quarter
wave data. that is , dcosqb computes a sequence from its
representation in terms of a cosine series with odd wave numbers.
the transform is defined below at output parameter x.

dcosqb is the unnormalized inverse of dcosqf since a call of dcosqb
followed by a call of dcosqf will multiply the input sequence x
by 4*n.

the array wsave which is used by subroutine dcosqb must be
initialized by calling subroutine dcosqi(n,wsave).


input parameters

n       the length of the array x to be transformed.  the method
        is most efficient when n is a product of small primes.

x       an array which contains the sequence to be transformed

wsave   a work array that must be dimensioned at least 3*n+15
        in the program that calls dcosqb. the wsave array must be
        initialized by calling subroutine dcosqi(n,wsave) and a
        different wsave array must be used for each different
        value of n. this initialization does not have to be
        repeated so long as n remains unchanged thus subsequent
        transforms can be obtained faster than the first.

output parameters

x       for i=1,...,n

             x(i)= the sum from k=1 to k=n of

               4*x(k)*cos((2*k-1)*(i-1)*pi/(2*n))

             a call of dcosqb followed by a call of
             dcosqf will multiply the sequence x by 4*n.
             therefore dcosqf is the unnormalized inverse
             of dcosqb.

wsave   contains initialization calculations which must not
        be destroyed between calls of dcosqb or dcosqf.

******************************************************************

subroutine zffti(n,wsave)

******************************************************************

subroutine zffti initializes the array wsave which is used in
both zfftf and zfftb. the prime factorization of n together with
a tabulation of the trigonometric functions are computed and
stored in wsave.

input parameter

n       the length of the sequence to be transformed

output parameter

wsave   a work array which must be dimensioned at least 4*n+15
        the same work array can be used for both zfftf and zfftb
        as long as n remains unchanged. different wsave arrays
        are required for different values of n. the contents of
        wsave must not be changed between calls of zfftf or zfftb.

******************************************************************

subroutine zfftf(n,c,wsave)

******************************************************************

subroutine zfftf computes the forward complex discrete fourier
transform (the fourier analysis). equivalently , zfftf computes
the fourier coefficients of a complex periodic sequence.
the transform is defined below at output parameter c.

the transform is not normalized. to obtain a normalized transform
the output must be divided by n. otherwise a call of zfftf
followed by a call of zfftb will multiply the sequence by n.

the array wsave which is used by subroutine zfftf must be
initialized by calling subroutine zffti(n,wsave).

input parameters


n      the length of the complex sequence c. the method is
       more efficient when n is the product of small primes. n

c      a complex array of length n which contains the sequence

wsave   a real work array which must be dimensioned at least 4n+15
        in the program that calls zfftf. the wsave array must be
        initialized by calling subroutine zffti(n,wsave) and a
        different wsave array must be used for each different
        value of n. this initialization does not have to be
        repeated so long as n remains unchanged thus subsequent
        transforms can be obtained faster than the first.
        the same wsave array can be used by zfftf and zfftb.

output parameters

c      for j=1,...,n

           c(j)=the sum from k=1,...,n of

                 c(k)*exp(-i*(j-1)*(k-1)*2*pi/n)

                       where i=sqrt(-1)

wsave   contains initialization calculations which must not be
        destroyed between calls of subroutine zfftf or zfftb

******************************************************************

subroutine zfftb(n,c,wsave)

******************************************************************

subroutine zfftb computes the backward complex discrete fourier
transform (the fourier synthesis). equivalently , zfftb computes
a complex periodic sequence from its fourier coefficients.
the transform is defined below at output parameter c.

a call of zfftf followed by a call of zfftb will multiply the
sequence by n.

the array wsave which is used by subroutine zfftb must be
initialized by calling subroutine zffti(n,wsave).

input parameters


n      the length of the complex sequence c. the method is
       more efficient when n is the product of small primes.

c      a complex array of length n which contains the sequence

wsave   a real work array which must be dimensioned at least 4n+15
        in the program that calls zfftb. the wsave array must be
        initialized by calling subroutine zffti(n,wsave) and a
        different wsave array must be used for each different
        value of n. this initialization does not have to be
        repeated so long as n remains unchanged thus subsequent
        transforms can be obtained faster than the first.
        the same wsave array can be used by zfftf and zfftb.

output parameters

c      for j=1,...,n

           c(j)=the sum from k=1,...,n of

                 c(k)*exp(i*(j-1)*(k-1)*2*pi/n)

                       where i=sqrt(-1)

wsave   contains initialization calculations which must not be
        destroyed between calls of subroutine zfftf or zfftb



["send index for vfftpk" describes a vectorized version of fftpack]

