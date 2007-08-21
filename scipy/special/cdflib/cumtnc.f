      SUBROUTINE cumtnc(t,df,pnonc,cum,ccum)
C**********************************************************************
C
C     SUBROUTINE CUMTNC(T,DF,PNONC,CUM,CCUM)
C
C                 CUMulative Non-Central T-distribution
C
C
C                              Function
C
C
C     Computes the integral from -infinity to T of the non-central
C     t-density.
C
C
C                              Arguments
C
C
C     T --> Upper limit of integration of the non-central t-density.
C                                                  T is DOUBLE PRECISION
C
C     DF --> Degrees of freedom of the non-central t-distribution.
C                                                  DF is DOUBLE PRECISIO
C
C     PNONC --> Non-centrality parameter of the non-central t distibutio
C                                                  PNONC is DOUBLE PRECI
C
C     CUM <-- Cumulative t-distribution.
C                                                  CCUM is DOUBLE PRECIS
C
C     CCUM <-- Compliment of Cumulative t-distribution.
C                                                  CCUM is DOUBLE PRECIS
C
C
C                              Method
C
C     Upper tail    of  the  cumulative  noncentral t   using
C     formulae from page 532  of Johnson, Kotz,  Balakrishnan, Coninuous
C     Univariate Distributions, Vol 2, 2nd Edition.  Wiley (1995)
C
C     This implementation starts the calculation at i = lambda,
C     which is near the largest Di.  It then sums forward and backward.
C***********************************************************************
C     .. Parameters ..

      DOUBLE PRECISION one,zero,half,two,onep5
      PARAMETER (one=1.0d0,zero=0.0d0,half=0.5d0,two=2.0d0,onep5=1.5d0)
      DOUBLE PRECISION conv
      PARAMETER (conv=1.0d-7)
      DOUBLE PRECISION tiny
      PARAMETER (tiny=1.0d-10)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION ccum,cum,df,pnonc,t
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION alghdf,b,bb,bbcent,bcent,cent,d,dcent,dpnonc,
     +                 dum1,dum2,e,ecent,halfdf,lambda,lnomx,lnx,omx,
     +                 pnonc2,s,scent,ss,sscent,t2,term,tt,twoi,x,
     +                 xi,xlnd,xlne
      INTEGER ierr
      LOGICAL qrevs
C     ..
C     .. External Functions ..
      DOUBLE PRECISION gamln
      EXTERNAL gamln
C     ..
C     .. External Subroutines ..
      EXTERNAL bratio,cumnor,cumt
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC abs,exp,int,log,max,min
C     ..

C     Case pnonc essentially zero

      IF (abs(pnonc).LE.tiny) THEN
          CALL cumt(t,df,cum,ccum)
          RETURN

      END IF

      qrevs = t .LT. zero
      IF (qrevs) THEN
          tt = -t
          dpnonc = -pnonc

      ELSE
          tt = t
          dpnonc = pnonc
      END IF

      pnonc2 = dpnonc*dpnonc
      t2 = tt*tt

      IF (abs(tt).LE.tiny) THEN
          CALL cumnor(-pnonc,cum,ccum)
          RETURN

      END IF

      lambda = half*pnonc2
      x = df/ (df+t2)
      omx = one - x

      lnx = log(x)
      lnomx = log(omx)

      halfdf = half*df
      alghdf = gamln(halfdf)

C     ******************** Case i = lambda

      cent = int(lambda)

      IF (cent.LT.one) cent = one

C     Compute d=T(2i) in log space and offset by exp(-lambda)

      xlnd = cent*log(lambda) - gamln(cent+one) - lambda

      dcent = exp(xlnd)

C     Compute e=t(2i+1) in log space offset by exp(-lambda)

      xlne = (cent+half)*log(lambda) - gamln(cent+onep5) - lambda
      ecent = exp(xlne)

      IF (dpnonc.LT.zero) ecent = -ecent

C     Compute bcent=B(2*cent)

      CALL bratio(halfdf,cent+half,x,omx,bcent,dum1,ierr)

C     compute bbcent=B(2*cent+1)

      CALL bratio(halfdf,cent+one,x,omx,bbcent,dum2,ierr)

C     Case bcent and bbcent are essentially zero
C     Thus t is effectively infinite

      IF ((bcent+bbcent).LT.tiny) THEN
          IF (qrevs) THEN
              cum = zero
              ccum = one

          ELSE
              cum = one
              ccum = zero
          END IF

          RETURN

      END IF

C     Case bcent and bbcent are essentially one
C     Thus t is effectively zero

      IF ((dum1+dum2).LT.tiny) THEN
          CALL cumnor(-pnonc,cum,ccum)
          RETURN

      END IF

C     First term in ccum is D*B + E*BB

      ccum = dcent*bcent + ecent*bbcent

C     compute s(cent) = B(2*(cent+1)) - B(2*cent))

      scent = gamln(halfdf+cent+half) - gamln(cent+onep5) - alghdf +
     +        halfdf*lnx + (cent+half)*lnomx
      scent = exp(scent)

C     compute ss(cent) = B(2*cent+3) - B(2*cent+1)

      sscent = gamln(halfdf+cent+one) - gamln(cent+two) - alghdf +
     +         halfdf*lnx + (cent+one)*lnomx
      sscent = exp(sscent)

C     ******************** Sum Forward

      xi = cent + one
      twoi = two*xi

      d = dcent

      e = ecent

      b = bcent

      bb = bbcent

      s = scent

      ss = sscent

   10 b = b + s
      bb = bb + ss

      d = (lambda/xi)*d
      e = (lambda/ (xi+half))*e

      term = d*b + e*bb

      ccum = ccum + term

      s = s*omx* (df+twoi-one)/ (twoi+one)

      ss = ss*omx* (df+twoi)/ (twoi+two)

      xi = xi + one
      twoi = two*xi

      IF (abs(term).GT.conv*ccum) GO TO 10

C     ******************** Sum Backward

      xi = cent
      twoi = two*xi

      d = dcent

      e = ecent

      b = bcent

      bb = bbcent

      s = scent* (one+twoi)/ ((df+twoi-one)*omx)

      ss = sscent* (two+twoi)/ ((df+twoi)*omx)

   20 b = b - s
      bb = bb - ss

      d = d* (xi/lambda)

      e = e* ((xi+half)/lambda)

      term = d*b + e*bb

      ccum = ccum + term

      xi = xi - one

      IF (xi.LT.half) GO TO 30

      twoi = two*xi

      s = s* (one+twoi)/ ((df+twoi-one)*omx)

      ss = ss* (two+twoi)/ ((df+twoi)*omx)

      IF (abs(term).GT.conv*ccum) GO TO 20

   30 CONTINUE

      IF (qrevs) THEN
          cum = half*ccum
          ccum = one - cum

      ELSE
          ccum = half*ccum
          cum = one - ccum
      END IF

C     Due to roundoff error the answer may not lie between zero and one
C     Force it to do so

      cum = max(min(cum,one),zero)
      ccum = max(min(ccum,one),zero)

      RETURN

      END
