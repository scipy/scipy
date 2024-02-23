      recursive
     *subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
      integer n,info,lwa
      double precision tol
      double precision x(n),fvec(n),wa(lwa)
      external fcn
c     **********
c
c     subroutine hybrd1
c
c     the purpose of hybrd1 is to find a zero of a system of
c     n nonlinear functions in n variables by a modification
c     of the powell hybrid method. this is done by using the
c     more general nonlinear equation solver hybrd. the user
c     must provide a subroutine which calculates the functions.
c     the jacobian is then calculated by a forward-difference
c     approximation.
c
c     the subroutine statement is
c
c       subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions. fcn must be declared
c         in an external statement in the user calling
c         program, and should be written as follows.
c
c         subroutine fcn(n,x,fvec,iflag)
c         integer n,iflag
c         double precision x(n),fvec(n)
c         ----------
c         calculate the functions at x and
c         return this vector in fvec.
c         ---------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of hybrd1.
c         in this case set iflag to a negative integer.
c
c       n is a positive integer input variable set to the number
c         of functions and variables.
c
c       x is an array of length n. on input x must contain
c         an initial estimate of the solution vector. on output x
c         contains the final estimate of the solution vector.
c
c       fvec is an output array of length n which contains
c         the functions evaluated at the output x.
c
c       tol is a nonnegative input variable. termination occurs
c         when the algorithm estimates that the relative error
c         between x and the solution is at most tol.
c
c       info is an integer output variable. if the user has
c         terminated execution, info is set to the (negative)
c         value of iflag. see description of fcn. otherwise,
c         info is set as follows.
c
c         info = 0   improper input parameters.
c
c         info = 1   algorithm estimates that the relative error
c                    between x and the solution is at most tol.
c
c         info = 2   number of calls to fcn has reached or exceeded
c                    200*(n+1).
c
c         info = 3   tol is too small. no further improvement in
c                    the approximate solution x is possible.
c
c         info = 4   iteration is not making good progress.
c
c       wa is a work array of length lwa.
c
c       lwa is a positive integer input variable not less than
c         (n*(3*n+13))/2.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... hybrd
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer index,j,lr,maxfev,ml,mode,mu,nfev,nprint
      double precision epsfcn,factor,one,xtol,zero
      data factor,one,zero /1.0d2,1.0d0,0.0d0/
      info = 0
c
c     check the input parameters for errors.
c
      if (n .le. 0 .or. tol .lt. zero .or. lwa .lt. (n*(3*n + 13))/2)
     *   go to 20
c
c     call hybrd.
c
      maxfev = 200*(n + 1)
      xtol = tol
      ml = n - 1
      mu = n - 1
      epsfcn = zero
      mode = 2
      do 10 j = 1, n
         wa(j) = one
   10    continue
      nprint = 0
      lr = (n*(n + 1))/2
      index = 6*n + lr
      call hybrd(fcn,n,x,fvec,xtol,maxfev,ml,mu,epsfcn,wa(1),mode,
     *           factor,nprint,info,nfev,wa(index+1),n,wa(6*n+1),lr,
     *           wa(n+1),wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
      if (info .eq. 5) info = 4
   20 continue
      return
c
c     last card of subroutine hybrd1.
c
      end
