      double precision function dqwgtf(x,omega,p2,p3,p4,integr)
c***begin prologue  dqwgtf
c***refer to   dqk15w
c***routines called  (none)
c***revision date 810101   (yymmdd)
c***keywords  cos or sin in weight function
c***author  piessens,robert, appl. math. & progr. div. - k.u.leuven
c           de doncker,elise,appl. math. * progr. div. - k.u.leuven
c***end prologue  dqwgtf
c
      double precision dcos,dsin,omega,omx,p2,p3,p4,x
      integer integr
c***first executable statement  dqwgtf
      omx = omega*x
      go to(10,20),integr
   10 dqwgtf = dcos(omx)
      go to 30
   20 dqwgtf = dsin(omx)
   30 return
      end
