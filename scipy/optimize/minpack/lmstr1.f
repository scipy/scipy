      recursive
     *subroutine lmstr1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,ipvt,wa,
     *                  lwa)
      integer m,n,ldfjac,info,lwa
      integer ipvt(n)
      double precision tol
      double precision x(n),fvec(m),fjac(ldfjac,n),wa(lwa)
      external fcn
c     **********
c
c     subroutine lmstr1
c
c     the purpose of lmstr1 is to minimize the sum of the squares of
c     m nonlinear functions in n variables by a modification of
c     the levenberg-marquardt algorithm which uses minimal storage.
c     this is done by using the more general least-squares solver
c     lmstr. the user must provide a subroutine which calculates
c     the functions and the rows of the jacobian.
c
c     the subroutine statement is
c
c       subroutine lmstr1(fcn,m,n,x,fvec,fjac,ldfjac,tol,info,
c                         ipvt,wa,lwa)
c
c     where
c
c       fcn is the name of the user-supplied subroutine which
c         calculates the functions and the rows of the jacobian.
c         fcn must be declared in an external statement in the
c         user calling program, and should be written as follows.
c
c         subroutine fcn(m,n,x,fvec,fjrow,iflag)
c         integer m,n,iflag
c         double precision x(n),fvec(m),fjrow(n)
c         ----------
c         if iflag = 1 calculate the functions at x and
c         return this vector in fvec.
c         if iflag = i calculate the (i-1)-st row of the
c         jacobian at x and return this vector in fjrow.
c         ----------
c         return
c         end
c
c         the value of iflag should not be changed by fcn unless
c         the user wants to terminate execution of lmstr1.
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
c       fjac is an output n by n array. the upper triangle of fjac
c         contains an upper triangular matrix r such that
c
c                t     t           t
c               p *(jac *jac)*p = r *r,
c
c         where p is a permutation matrix and jac is the final
c         calculated jacobian. column j of p is column ipvt(j)
c         (see below) of the identity matrix. the lower triangular
c         part of fjac contains information generated during
c         the computation of r.
c
c       ldfjac is a positive integer input variable not less than n
c         which specifies the leading dimension of the array fjac.
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
c         info = 5  number of calls to fcn with iflag = 1 has
c                   reached 100*(n+1).
c
c         info = 6  tol is too small. no further reduction in
c                   the sum of squares is possible.
c
c         info = 7  tol is too small. no further improvement in
c                   the approximate solution x is possible.
c
c       ipvt is an integer output array of length n. ipvt
c         defines a permutation matrix p such that jac*p = q*r,
c         where jac is the final calculated jacobian, q is
c         orthogonal (not stored), and r is upper triangular.
c         column j of p is column ipvt(j) of the identity matrix.
c
c       wa is a work array of length lwa.
c
c       lwa is a positive integer input variable not less than 5*n+m.
c
c     subprograms called
c
c       user-supplied ...... fcn
c
c       minpack-supplied ... lmstr
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, dudley v. goetschel, kenneth e. hillstrom,
c     jorge j. more
c
c     **********
      integer maxfev,mode,nfev,njev,nprint
      double precision factor,ftol,gtol,xtol,zero
      data factor,zero /1.0d2,0.0d0/
      info = 0
c
c     check the input parameters for errors.
c
      if (n .le. 0 .or. m .lt. n .or. ldfjac .lt. n .or. tol .lt. zero
     *    .or. lwa .lt. 5*n + m) go to 10
c
c     call lmstr.
c
      maxfev = 100*(n + 1)
      ftol = tol
      xtol = tol
      gtol = zero
      mode = 1
      nprint = 0
      call lmstr(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,maxfev,
     *           wa(1),mode,factor,nprint,info,nfev,njev,ipvt,wa(n+1),
     *           wa(2*n+1),wa(3*n+1),wa(4*n+1),wa(5*n+1))
      if (info .eq. 8) info = 4
   10 continue
      return
c
c     last card of subroutine lmstr1.
c
      end
