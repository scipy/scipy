      subroutine timer(ttime)
      double precision ttime
c
      real temp
c
c     This routine computes cpu time in double precision; it makes use of 
c     the intrinsic f90 cpu_time therefore a conversion type is
c     needed.
c
c           J.L Morales  Departamento de Matematicas, 
c                        Instituto Tecnologico Autonomo de Mexico
c                        Mexico D.F.
c
c           J.L Nocedal  Department of Electrical Engineering and
c                        Computer Science.
c                        Northwestern University. Evanston, IL. USA
c                         
c                        January 21, 2011
c

c     Don't actually return the time --- in case of old Fortran compiler

      ttime = 0

      return

      end
      
