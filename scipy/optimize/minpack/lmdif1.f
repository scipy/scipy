      recursive
     *subroutine lmdif1(fcn,m,n,x,fvec,tol,info,iwa,wa,lwa)
      integer m,n,info,lwa
      integer iwa(n)
      double precision tol
      double precision x(n),fvec(m),wa(lwa)
      external fcn
c     **********
c
c     subroutine lmdif1
c
c     the purpose of lmdif1 is to minimize the sum of the squares of
c     m nonlinear functions in n variables by a modification of the
c     levenberg-marquardt algorithm. this is done by using the more
c     general least-squares solver lmdif. the user must provide a
c     subroutine which calculates the functions. the jacobian is
c     then calculated by a forward-difference approximation.
c
c     the subroutine statement is
c
c       subroutine lmdif1(fcn,m,n,x,fvec,tol,info,iwa,wa,lwa)
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions. fcn must be declared
c         in an external statement in the user calling
c         program, and should be written as follows.
c
c         subroutine fcn(m,n,x,fvec,iflag)
c         integer m,n,iflag
c         double precision x(n),fvec(m)
c         ----------
c         calculate the functions at x and
c         return this vector in fvec.
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of lmdif1.
c         in this case set iflag to a negative integer.
c
c       m is a positive integer input variable set to the number
c         of functions.
c
c       n is a positive integer input variable set to the number
c         of variables. n must not exceed m.
c
c       x is an array of length n. on input x must contain
c         an initial estimate of the solution vector. on output x
c         contains the final estimate of the solution vector.
c
c       fvec is an output array of length m which contains
c         the functions evaluated at the output x.
c
c       tol is a nonnegative input variable. termination occurs
c         when the algorithm estimates either that the relative
c         error in the sum of squares is at most tol or that
c         the relative error between x and the solution is at
c         most tol.
c
c       info is an integer output variable. if the user has
c         terminated execution, info is set to the (negative)
c         value of iflag. see description of fcn. otherwise,
c         info is set as follows.
c
c         info = 0  improper input parameters.
c
c         info = 1  algorithm estimates that the relative error
c                   in the sum of squares is at most tol.
c
c         info = 2  algorithm estimates that the relative error
c                   between x and the solution is at most tol.
c
c         info = 3  conditions for info = 1 and info = 2 both hold.
c
c         info = 4  fvec is orthogonal to the columns of the
c                   jacobian to machine precision.
c
c         info = 5  number of calls to fcn has reached or
c                   exceeded 200*(n+1).
c
c         info = 6  tol is too small. no further reduction in
c                   the sum of squares is possible.
c
c         info = 7  tol is too small. no further improvement in
c                   the approximate solution x is possible.
c
c       iwa is an integer work array of length n.
c
c       wa is a work array of length lwa.
c
c       lwa is a positive integer input variable not less than
c         m*n+5*n+m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... lmdif
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer maxfev,mode,mp5n,nfev,nprint
      double precision epsfcn,factor,ftol,gtol,xtol,zero
      data factor,zero /1.0d2,0.0d0/
      info = 0
c
c     check the input parameters for errors.
c
      if (n .le. 0 .or. m .lt. n .or. tol .lt. zero
     *    .or. lwa .lt. m*n + 5*n + m) go to 10
c
c     call lmdif.
c
      maxfev = 200*(n + 1)
      ftol = tol
      xtol = tol
      gtol = zero
      epsfcn = zero
      mode = 1
      nprint = 0
      mp5n = m + 5*n
      call lmdif(fcn,m,n,x,fvec,ftol,xtol,gtol,maxfev,epsfcn,wa(1),
     *           mode,factor,nprint,info,nfev,wa(mp5n+1),m,iwa,
     *           wa(n+1),wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
      if (info .eq. 8) info = 4
   10 continue
      return
c
c     last card of subroutine lmdif1.
c
      end
