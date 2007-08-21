* wsc@research.bell-labs.com Mon Dec 30 16:55 EST 1985
* W. S. Cleveland
* Bell Laboratories
* Murray Hill NJ 07974
* 
* outline of this file:
*    lines 1-72   introduction
*        73-177   documentation for lowess
*       178-238   ratfor version of lowess
*       239-301   documentation for lowest
*       302-350   ratfor version of lowest
*       351-end   test driver and fortran version of lowess and lowest
* 
*   a multivariate version is available by "send dloess from a"
* 
*              COMPUTER PROGRAMS FOR LOCALLY WEIGHTED REGRESSION
* 
*             This package consists  of  two  FORTRAN  programs  for
*        smoothing    scatterplots   by   robust   locally   weighted
*        regression, or lowess.   The  principal  routine  is  LOWESS
*        which   computes   the  smoothed  values  using  the  method
*        described in The Elements of Graphing Data, by William S.
*        Cleveland    (Wadsworth,    555 Morego   Street,   Monterey,
*        California 93940).
* 
*             LOWESS calls a support routine, LOWEST, the code for
*        which is included. LOWESS also calls a routine  SORT,  which
*        the user must provide.
* 
*             To reduce the computations, LOWESS  requires  that  the
*        arrays  X  and  Y,  which  are  the  horizontal and vertical
*        coordinates, respectively, of the scatterplot, be such  that
*        X  is  sorted  from  smallest  to  largest.   The  user must
*        therefore use another sort routine which will sort X  and  Y
*        according  to X.
*             To summarize the scatterplot, YS,  the  fitted  values,
*        should  be  plotted  against X.   No  graphics  routines are
*        available in the package and must be supplied by the user.
* 
*             The FORTRAN code for the routines LOWESS and LOWEST has
*        been   generated   from   higher   level   RATFOR   programs
*        (B. W. Kernighan, ``RATFOR:  A Preprocessor for  a  Rational
*        Fortran,''  Software Practice and Experience, Vol. 5 (1975),
*        which are also included.
* 
*             The following are data and output from LOWESS that  can
*        be  used  to check your implementation of the routines.  The
*        notation (10)v means 10 values of v.
* 
* 
* 
* 
*        X values:
*          1  2  3  4  5  (10)6  8  10  12  14  50
* 
*        Y values:
*           18  2  15  6  10  4  16  11  7  3  14  17  20  12  9  13  1  8  5  19
* 
* 
*        YS values with F = .25, NSTEPS = 0, DELTA = 0.0
*         13.659  11.145  8.701  9.722  10.000  (10)11.300  13.000  6.440  5.596
*           5.456  18.998
* 
*        YS values with F = .25, NSTEPS = 0 ,  DELTA = 3.0
*          13.659  12.347  11.034  9.722  10.511  (10)11.300  13.000  6.440  5.596
*            5.456  18.998
* 
*        YS values with F = .25, NSTEPS = 2, DELTA = 0.0
*          14.811  12.115  8.984  9.676  10.000  (10)11.346  13.000  6.734  5.744
*            5.415  18.998
* 
* 
* 
* 
*                                   LOWESS
* 
* 
* 
*        Calling sequence
* 
*        CALL LOWESS(X,Y,N,F,NSTEPS,DELTA,YS,RW,RES)
* 
*        Purpose
* 
*        LOWESS computes the smooth of a scatterplot of Y  against  X
*        using  robust  locally  weighted regression.  Fitted values,
*        YS, are computed at each of the  values  of  the  horizontal
*        axis in X.
* 
*        Argument description
* 
*              X = Input; abscissas of the points on the
*                  scatterplot; the values in X must be ordered
*                  from smallest to largest.
*              Y = Input; ordinates of the points on the
*                  scatterplot.
*              N = Input; dimension of X,Y,YS,RW, and RES.
*              F = Input; specifies the amount of smoothing; F is
*                  the fraction of points used to compute each
*                  fitted value; as F increases the smoothed values
*                  become smoother; choosing F in the range .2 to
*                  .8 usually results in a good fit; if you have no
*                  idea which value to use, try F = .5.
*         NSTEPS = Input; the number of iterations in the robust
*                  fit; if NSTEPS = 0, the nonrobust fit is
*                  returned; setting NSTEPS equal to 2 should serve
*                  most purposes.
*          DELTA = input; nonnegative parameter which may be used
*                  to save computations; if N is less than 100, set
*                  DELTA equal to 0.0; if N is greater than 100 you
*                  should find out how DELTA works by reading the
*                  additional instructions section.
*             YS = Output; fitted values; YS(I) is the fitted value
*                  at X(I); to summarize the scatterplot, YS(I)
*                  should be plotted against X(I).
*             RW = Output; robustness weights; RW(I) is the weight
*                  given to the point (X(I),Y(I)); if NSTEPS = 0,
*                  RW is not used.
*            RES = Output; residuals; RES(I) = Y(I)-YS(I).
* 
* 
*        Other programs called
* 
*               LOWEST
*               SSORT
* 
*        Additional instructions
* 
*        DELTA can be used to save computations.   Very  roughly  the
*        algorithm  is  this:   on the initial fit and on each of the
*        NSTEPS iterations locally weighted regression fitted  values
*        are computed at points in X which are spaced, roughly, DELTA
*        apart; then the fitted values at the  remaining  points  are
*        computed  using  linear  interpolation.   The  first locally
*        weighted regression (l.w.r.) computation is carried  out  at
*        X(1)  and  the  last  is  carried  out at X(N).  Suppose the
*        l.w.r. computation is carried out at  X(I).   If  X(I+1)  is
*        greater  than  or  equal  to  X(I)+DELTA,  the  next  l.w.r.
*        computation is carried out at X(I+1).   If  X(I+1)  is  less
*        than X(I)+DELTA, the next l.w.r.  computation is carried out
*        at the largest X(J) which is greater than or equal  to  X(I)
*        but  is not greater than X(I)+DELTA.  Then the fitted values
*        for X(K) between X(I)  and  X(J),  if  there  are  any,  are
*        computed  by  linear  interpolation  of the fitted values at
*        X(I) and X(J).  If N is less than 100 then DELTA can be  set
*        to  0.0  since  the  computation time will not be too great.
*        For larger N it is typically not necessary to carry out  the
*        l.w.r.  computation for all points, so that much computation
*        time can be saved by taking DELTA to be  greater  than  0.0.
*        If  DELTA =  Range  (X)/k  then,  if  the  values  in X were
*        uniformly  scattered  over  the  range,  the   full   l.w.r.
*        computation  would be carried out at approximately k points.
*        Taking k to be 50 often works well.
* 
*        Method
* 
*        The fitted values are computed by using the nearest neighbor
*        routine  and  robust locally weighted regression of degree 1
*        with the tricube weight function.  A few additional features
*        have  been  added.  Suppose r is FN truncated to an integer.
*        Let  h  be  the  distance  to  the  r-th  nearest   neighbor
*        from X(I).   All  points within h of X(I) are used.  Thus if
*        the r-th nearest neighbor is exactly the  same  distance  as
*        other  points,  more  than r points can possibly be used for
*        the smooth at  X(I).   There  are  two  cases  where  robust
*        locally  weighted regression of degree 0 is actually used at
*        X(I).  One case occurs when  h  is  0.0.   The  second  case
*        occurs  when  the  weighted  standard error of the X(I) with
*        respect to the weights w(j) is  less  than  .001  times  the
*        range  of the X(I), where w(j) is the weight assigned to the
*        j-th point of X (the tricube  weight  times  the  robustness
*        weight)  divided by the sum of all of the weights.  Finally,
*        if the w(j) are all zero for the smooth at X(I), the  fitted
*        value is taken to be Y(I).
* 
* 
* 
* 
*  subroutine lowess(x,y,n,f,nsteps,delta,ys,rw,res)
*  real x(n),y(n),ys(n),rw(n),res(n)
*  logical ok
*  if (n<2){ ys(1) = y(1); return }
*  ns = max0(min0(ifix(f*float(n)),n),2)  # at least two, at most n points
*  for(iter=1; iter<=nsteps+1; iter=iter+1){      # robustness iterations
*         nleft = 1; nright = ns
*         last = 0        # index of prev estimated point
*         i = 1   # index of current point
*         repeat{
*                 while(nright<n){
*  # move nleft, nright to right if radius decreases
*                         d1 = x(i)-x(nleft)
*                         d2 = x(nright+1)-x(i)
*  # if d1<=d2 with x(nright+1)==x(nright), lowest fixes
*                         if (d1<=d2) break
*  # radius will not decrease by move right
*                         nleft = nleft+1
*                         nright = nright+1
*                         }
*                 call lowest(x,y,n,x(i),ys(i),nleft,nright,res,iter>1,rw,ok)
*  # fitted value at x(i)
*                 if (!ok) ys(i) = y(i)
*  # all weights zero - copy over value (all rw==0)
*                 if (last<i-1) { # skipped points -- interpolate
*                         denom = x(i)-x(last)    # non-zero - proof?
*                         for(j=last+1; j<i; j=j+1){
*                                 alpha = (x(j)-x(last))/denom
*                                 ys(j) = alpha*ys(i)+(1.0-alpha)*ys(last)
*                                 }
*                         }
*                 last = i        # last point actually estimated
*                 cut = x(last)+delta     # x coord of close points
*                 for(i=last+1; i<=n; i=i+1){     # find close points
*                         if (x(i)>cut) break     # i one beyond last pt within cut
*                         if(x(i)==x(last)){      # exact match in x
*                                 ys(i) = ys(last)
*                                 last = i
*                                 }
*                         }
*                 i=max0(last+1,i-1)
*  # back 1 point so interpolation within delta, but always go forward
*                 } until(last>=n)
*         do i = 1,n      # residuals
*                 res(i) = y(i)-ys(i)
*         if (iter>nsteps) break  # compute robustness weights except last time
*         do i = 1,n
*                 rw(i) = abs(res(i))
*         call sort(rw,n)
*         m1 = 1+n/2; m2 = n-m1+1
*         cmad = 3.0*(rw(m1)+rw(m2))      # 6 median abs resid
*         c9 = .999*cmad; c1 = .001*cmad
*         do i = 1,n {
*                 r = abs(res(i))
*                 if(r<=c1) rw(i)=1.      # near 0, avoid underflow
*                 else if(r>c9) rw(i)=0.  # near 1, avoid underflow
*                 else rw(i) = (1.0-(r/cmad)**2)**2
*                 }
*         }
*  return
*  end
* 
* 
* 
* 
*                                   LOWEST
* 
* 
* 
*        Calling sequence
* 
*        CALL LOWEST(X,Y,N,XS,YS,NLEFT,NRIGHT,W,USERW,RW,OK)
* 
*        Purpose
* 
*        LOWEST is a support routine for LOWESS and  ordinarily  will
*        not  be  called  by  the  user.   The  fitted  value, YS, is
*        computed  at  the  value,  XS,  of  the   horizontal   axis.
*        Robustness  weights,  RW,  can  be employed in computing the
*        fit.
* 
*        Argument description
* 
* 
*              X = Input; abscissas of the points on the
*                  scatterplot; the values in X must be ordered
*                  from smallest to largest.
*              Y = Input; ordinates of the points on the
*                  scatterplot.
*              N = Input; dimension of X,Y,W, and RW.
*             XS = Input; value of the horizontal axis at which the
*                  smooth is computed.
*             YS = Output; fitted value at XS.
*          NLEFT = Input; index of the first point which should be
*                  considered in computing the fitted value.
*         NRIGHT = Input; index of the last point which should be
*                  considered in computing the fitted value.
*              W = Output; W(I) is the weight for Y(I) used in the
*                  expression for YS, which is the sum from
*                  I = NLEFT to NRIGHT of W(I)*Y(I); W(I) is
*                  defined only at locations NLEFT to NRIGHT.
*          USERW = Input; logical variable; if USERW is .TRUE., a
*                  robust fit is carried out using the weights in
*                  RW; if USERW is .FALSE., the values in RW are
*                  not used.
*             RW = Input; robustness weights.
*             OK = Output; logical variable; if the weights for the
*                  smooth are all 0.0, the fitted value, YS, is not
*                  computed and OK is set equal to .FALSE.; if the
*                  fitted value is computed OK is set equal to
* 
* 
*        Method
* 
*        The smooth at XS is computed using (robust) locally weighted
*        regression of degree 1.  The tricube weight function is used
*        with h equal to the maximum of XS-X(NLEFT) and X(NRIGHT)-XS.
*        Two  cases  where  the  program  reverts to locally weighted
*        regression of degree 0 are described  in  the  documentation
*        for LOWESS.
* 
* 
* 
* 
*  subroutine lowest(x,y,n,xs,ys,nleft,nright,w,userw,rw,ok)
*  real x(n),y(n),w(n),rw(n)
*  logical userw,ok
*  range = x(n)-x(1)
*  h = amax1(xs-x(nleft),x(nright)-xs)
*  h9 = .999*h
*  h1 = .001*h
*  a = 0.0        # sum of weights
*  for(j=nleft; j<=n; j=j+1){     # compute weights (pick up all ties on right)
*         w(j)=0.
*         r = abs(x(j)-xs)
*         if (r<=h9) {    # small enough for non-zero weight
*                 if (r>h1) w(j) = (1.0-(r/h)**3)**3
*                 else      w(j) = 1.
*                 if (userw) w(j) = rw(j)*w(j)
*                 a = a+w(j)
*                 }
*         else if(x(j)>xs)break   # get out at first zero wt on right
*         }
*  nrt=j-1        # rightmost pt (may be greater than nright because of ties)
*  if (a<=0.0) ok = FALSE
*  else { # weighted least squares
*         ok = TRUE
*         do j = nleft,nrt
*                 w(j) = w(j)/a   # make sum of w(j) == 1
*         if (h>0.) {     # use linear fit
*                 a = 0.0
*                 do j = nleft,nrt
*                         a = a+w(j)*x(j) # weighted center of x values
*                 b = xs-a
*                 c = 0.0
*                 do j = nleft,nrt
*                         c = c+w(j)*(x(j)-a)**2
*                 if(sqrt(c)>.001*range) {
*  # points are spread out enough to compute slope
*                         b = b/c
*                         do j = nleft,nrt
*                                 w(j) = w(j)*(1.0+b*(x(j)-a))
*                         }
*                 }
*         ys = 0.0
*         do j = nleft,nrt
*                 ys = ys+w(j)*y(j)
*         }
*  return
*  end
* 
* 
* 
c  test driver for lowess
c  for expected output, see introduction
      double precision x(20), y(20), ys(20), rw(20), res(20)
      data x /1,2,3,4,5,10*6,8,10,12,14,50/
      data y /18,2,15,6,10,4,16,11,7,3,14,17,20,12,9,13,1,8,5,19/
      call lowess(x,y,20,.25,0,0.,ys,rw,res)
      write(6,*) ys
      call lowess(x,y,20,.25,0,3.,ys,rw,res)
      write(6,*) ys
      call lowess(x,y,20,.25,2,0.,ys,rw,res)
      write(6,*) ys
      end
c**************************************************************
c  Fortran output from ratfor
c
      subroutine lowess(x, y, n, f, nsteps, delta, ys, rw, res)
      integer n, nsteps
      double precision x(n), y(n), f, delta, ys(n), rw(n), res(n)
      integer nright, i, j, iter, last, mid(2), ns, nleft
      double precision cut, cmad, r, d1, d2
      double precision c1, c9, alpha, denom, dabs
      logical ok
      if (n .ge. 2) goto 1
         ys(1) = y(1)
         return
c at least two, at most n points
   1  ns = max(min(int(f*dble(n)), n), 2)
      iter = 1
         goto  3
   2     iter = iter+1
   3     if (iter .gt. nsteps+1) goto  22
c robustness iterations
         nleft = 1
         nright = ns
c index of prev estimated point
         last = 0
c index of current point
         i = 1
   4        if (nright .ge. n) goto  5
c move nleft, nright to right if radius decreases
               d1 = x(i)-x(nleft)
c if d1<=d2 with x(nright+1)==x(nright), lowest fixes
               d2 = x(nright+1)-x(i)
               if (d1 .le. d2) goto  5
c radius will not decrease by move right
               nleft = nleft+1
               nright = nright+1
               goto  4
c fitted value at x(i)
   5        call lowest(x, y, n, x(i), ys(i), nleft, nright, res, iter
     +     .gt. 1, rw, ok)
            if (.not. ok) ys(i) = y(i)
c all weights zero - copy over value (all rw==0)
            if (last .ge. i-1) goto 9
               denom = x(i)-x(last)
c skipped points -- interpolate
c non-zero - proof?
               j = last+1
                  goto  7
   6              j = j+1
   7              if (j .ge. i) goto  8
                  alpha = (x(j)-x(last))/denom
                  ys(j) = alpha*ys(i)+(1.D0-alpha)*ys(last)
                  goto  6
   8           continue
c last point actually estimated
   9        last = i
c x coord of close points
            cut = x(last)+delta
            i = last+1
               goto  11
  10           i = i+1
  11           if (i .gt. n) goto  13
c find close points
               if (x(i) .gt. cut) goto  13
c i one beyond last pt within cut
               if (x(i) .ne. x(last)) goto 12
                  ys(i) = ys(last)
c exact match in x
                  last = i
  12           continue
               goto  10
c back 1 point so interpolation within delta, but always go forward
  13        i = max(last+1, i-1)
  14        if (last .lt. n) goto  4
c residuals
         do  15 i = 1, n
            res(i) = y(i)-ys(i)
  15        continue
         if (iter .gt. nsteps) goto  22
c compute robustness weights except last time
         do  16 i = 1, n
            rw(i) = dabs(res(i))
  16        continue
         call ssort(rw,n)
         mid(1) = n/2+1
         mid(2) = n-mid(1)+1
c 6 median abs resid
         cmad = 3.D0*(rw(mid(1))+rw(mid(2)))
         c9 = .999999D0*cmad
         c1 = .000001D0*cmad
         do  21 i = 1, n
            r = dabs(res(i))
            if (r .gt. c1) goto 17
               rw(i) = 1.D0
c near 0, avoid underflow
               goto  20
  17           if (r .le. c9) goto 18
                  rw(i) = 0.D0
c near 1, avoid underflow
                  goto  19
  18              rw(i) = (1.D0-(r/cmad)**2.D0)**2.D0
  19        continue
  20        continue
  21        continue
         goto  2
  22  return
      end
      
      
      subroutine lowest(x, y, n, xs, ys, nleft, nright, w, userw
     +, rw, ok)
      integer n
      integer nleft, nright
      double precision x(n), y(n), xs, ys, w(n), rw(n)
      logical userw, ok
      integer nrt, j
      double precision dabs, a, b, c, h, r
      double precision h1, dsqrt, h9, max, range
      range = x(n)-x(1)
      h = max(xs-x(nleft), x(nright)-xs)
      h9 = .999999D0*h
      h1 = .000001D0*h
c sum of weights
      a = 0.D0
      j = nleft
         goto  2
   1     j = j+1
   2     if (j .gt. n) goto  7
c compute weights (pick up all ties on right)
         w(j) = 0.D0
         r = dabs(x(j)-xs)
         if (r .gt. h9) goto 5
            if (r .le. h1) goto 3
               w(j) = (1.D0-(r/h)**3.D0)**3.D0
c small enough for non-zero weight
               goto  4
   3           w(j) = 1.D0
   4        if (userw) w(j) = rw(j)*w(j)
            a = a+w(j)
            goto  6
   5        if (x(j) .gt. xs) goto  7
c get out at first zero wt on right
   6     continue
         goto  1
c rightmost pt (may be greater than nright because of ties)
   7  nrt = j-1
      if (a .gt. 0.D0) goto 8
         ok = .false.
         goto  16
   8     ok = .true.
c weighted least squares
         do  9 j = nleft, nrt
c make sum of w(j) == 1
            w(j) = w(j)/a
   9        continue
         if (h .le. 0.D0) goto 14
            a = 0.D0
c use linear fit
            do  10 j = nleft, nrt
c weighted center of x values
               a = a+w(j)*x(j)
  10           continue
            b = xs-a
            c = 0.D0
            do  11 j = nleft, nrt
               c = c+w(j)*(x(j)-a)**2
  11           continue
            if (dsqrt(c) .le. .0000001D0*range) goto 13
               b = b/c
c points are spread out enough to compute slope
               do  12 j = nleft, nrt
                  w(j) = w(j)*(b*(x(j)-a)+1.D0)
  12              continue
  13        continue
  14     ys = 0.D0
         do  15 j = nleft, nrt
            ys = ys+w(j)*y(j)
  15        continue
  16  return
      end

      
      subroutine ssort(a,n)

C Sorting by Hoare method, C.A.C.M. (1961) 321, modified by Singleton
C C.A.C.M. (1969) 185.
	  double precision a(n)
	  integer iu(16), il(16)
      integer p

      i =1
      j = n
      m = 1
  5   if (i.ge.j) goto 70
c first order a(i),a(j),a((i+j)/2), and use median to split the data
 10   k=i
      ij=(i+j)/2
      t=a(ij)
      if(a(i) .le. t) goto 20
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
 20   l=j
      if(a(j).ge.t) goto 40
      a(ij)=a(j)
      a(j)=t
      t=a(ij)
      if(a(i).le.t) goto 40
      a(ij)=a(i)
      a(i)=t
      t=a(ij)
      goto 40
 30   a(l)=a(k)
      a(k)=tt
 40   l=l-1
      if(a(l) .gt. t) goto 40
      tt=a(l)
c split the data into a(i to l) .lt. t, a(k to j) .gt. t
 50   k=k+1
      if(a(k) .lt. t) goto 50
      if(k .le. l) goto 30
      p=m
      m=m+1
c split the larger of the segments
      if (l-i .le. j-k) goto 60
      il(p)=i
      iu(p)=l
      i=k
      goto 80
 60   il(p)=k
      iu(p)=j
      j=l
      goto 80
 70   m=m-1
      if(m .eq. 0) return
      i =il(m)
      j=iu(m)
c short sections are sorted by bubble sort
 80   if (j-i .gt. 10) goto 10
      if (i .eq. 1) goto 5
      i=i-1
 90   i=i+1
      if(i .eq. j) goto 70
      t=a(i+1)
      if(a(i) .le. t) goto 90
      k=i
 100  a(k+1)=a(k)
      k=k-1
      if(t .lt. a(k)) goto 100
      a(k+1)=t
      goto 90

      end      