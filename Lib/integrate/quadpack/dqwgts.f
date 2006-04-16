      double precision function dqwgts(x,a,b,alfa,beta,integr)
c***begin prologue  dqwgts
c***refer to dqk15w
c***routines called  (none)
c***revision date  810101   (yymmdd)
c***keywords  weight function, algebraico-logarithmic
c             end-point singularities
c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
c***purpose  this function subprogram is used together with the
c            routine dqaws and defines the weight function.
c***end prologue  dqwgts
c
      double precision a,alfa,b,beta,bmx,dlog,x,xma
      integer integr
c***first executable statement  dqwgts
      xma = x-a
      bmx = b-x
      dqwgts = xma**alfa*bmx**beta
      go to (40,10,20,30),integr
   10 dqwgts = dqwgts*dlog(xma)
      go to 40
   20 dqwgts = dqwgts*dlog(bmx)
      go to 40
   30 dqwgts = dqwgts*dlog(xma)*dlog(bmx)
   40 return
      end
