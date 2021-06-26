C+-----------------------------------------------------------------------+
C| Program       : Direct.f (subfile DIRserial.f)                        |
C| Last modified : 04-12-2001                                            |
C| Written by    : Joerg Gablonsky                                       |
C| SUBROUTINEs, which differ depENDing on the serial or parallel version.|
C+-----------------------------------------------------------------------+

C+-----------------------------------------------------------------------+
C| SUBROUTINE for sampling.                                              |
C+-----------------------------------------------------------------------+
      SUBROUTINE DIRSamplef(c,ArrayI,delta,sample,new,length,
     +           dwrit,logfile,f,free,maxI,point,fcn,x,l,fmin,
     +           minpos,u,n,maxfunc,maxdeep,oops,fmax,
     +            IFeasiblef,IInfesiblef,
     +           iidata, iisize, ddata, idsize, cdata, icsize) 
      IMPLICIT NONE

C+-----------------------------------------------------------------------+
C| JG 07/16/01 fcn must be declared external.                            |
C+-----------------------------------------------------------------------+
      EXTERNAL fcn

      INTEGER n,maxfunc,maxdeep,oops
      INTEGER maxI,ArrayI(n),sample,new
      INTEGER length(maxfunc,n),free,point(maxfunc),i
C+-----------------------------------------------------------------------+
C| JG 07/16/01 Removed fcn.                                              |
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION c(maxfunc,n),delta,x(n)
      DOUBLE PRECISION l(n),u(n),f(maxfunc,2)
      DOUBLE PRECISION fmin
      INTEGER pos,j,dwrit,logfile,minpos
      INTEGER helppoint,kret
C+-----------------------------------------------------------------------+
C| JG 01/22/01 Added variable to keep track of the maximum value found.  |
C|             Added variable to keep track IF feasible point was found. |
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION fmax
      INTEGER  IFeasiblef,IInfesiblef
C+-----------------------------------------------------------------------+
C| Variables to pass user defined data to the function to be optimized.  |
C+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      DOUBLE PRECISION ddata(idsize)
      CHARACTER*40 cdata(icsize)
      
C+-----------------------------------------------------------------------+
C| Set the pointer to the first function to be evaluated,                |
C| store this position also in helppoint.                                |
C+-----------------------------------------------------------------------+
      pos = new
      helppoint = pos
C+-----------------------------------------------------------------------+
C| Iterate over all points, where the function should be                 |
C| evaluated.                                                            |
C+-----------------------------------------------------------------------+
      DO 40,j=1,maxI + maxI
C+-----------------------------------------------------------------------+
C| Copy the position into the helparrayy x.                              |
C+-----------------------------------------------------------------------+
         DO 60,i=1,n
           x(i) = c(pos,i)
60       CONTINUE
C+-----------------------------------------------------------------------+
C| Call the function.                                                    |
C+-----------------------------------------------------------------------+
         CALL DIRinfcn(fcn,x,l,u,n,f(pos,1),kret,
     +                 iidata, iisize, ddata, idsize, cdata, icsize)
C+-----------------------------------------------------------------------+
C| Remember IF an infeasible point has been found.                       |
C+-----------------------------------------------------------------------+
         IInfesiblef = max(IInfesiblef,kret)
         IF (kret .eq. 0) then
C+-----------------------------------------------------------------------+
C| IF the function evaluation was O.K., set the flag in                  |
C| f(pos,2). Also mark that a feasible point has been found.             |
C+-----------------------------------------------------------------------+
           f(pos,2) = 0.D0
            IFeasiblef = 0
C+-----------------------------------------------------------------------+
C| JG 01/22/01 Added variable to keep track of the maximum value found.  |
C+-----------------------------------------------------------------------+
           fmax = max(f(pos,1),fmax)
         END if
         IF (kret .ge. 1) then
C+-----------------------------------------------------------------------+
C|  IF the function could not be evaluated at the given point,            |
C| set flag to mark this (f(pos,2) and store the maximum                 |
C| box-sidelength in f(pos,1).                                           |
C+-----------------------------------------------------------------------+
           f(pos,2) = 2.D0
           f(pos,1) = fmax
         END if
C+-----------------------------------------------------------------------+
C|  IF the function could not be evaluated due to a failure in            |
C| the setup, mark this.                                                 |
C+-----------------------------------------------------------------------+
         IF (kret .eq. -1) then
           f(pos,2) = -1.D0
         END if
C+-----------------------------------------------------------------------+
C| Set the position to the next point, at which the function             |
C| should e evaluated.                                                   |
C+-----------------------------------------------------------------------+
         pos = point(pos)
40    CONTINUE
      pos = helppoint
C+-----------------------------------------------------------------------+
C| Iterate over all evaluated points and see, IF the minimal             |
C| value of the function has changed.  IF this has happEND,               |
C| store the minimal value and its position in the array.                |
C| Attention: Only valied values are checked!!                           |
C+-----------------------------------------------------------------------+
      DO 50,j=1,maxI + maxI
        IF ((f(pos,1) .LT. fmin) .and. (f(pos,2) .eq. 0)) THEN
          fmin = f(pos,1) 
          minpos = pos
        END IF
        pos = point(pos)
50    CONTINUE
      END

C+-----------------------------------------------------------------------+
C| Problem-specific Initialisation                                       |
C+-----------------------------------------------------------------------+
      SUBROUTINE DIRInitSpecific(x,n,
     +   iidata, iisize, ddata, idsize, cdata, icsize)
      IMPLICIT NONE
      INTEGER n
      DOUBLE PRECISION x(n)
C+-----------------------------------------------------------------------+
C| Variables to pass user defined data to the function to be optimized.  |
C+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      DOUBLE PRECISION ddata(idsize)
      CHARACTER*40 cdata(icsize)
C+-----------------------------------------------------------------------+
C| Problem - specific variables !                                        |
C+-----------------------------------------------------------------------+

C+-----------------------------------------------------------------------+
C| END of problem - specific variables !                                 |
C+-----------------------------------------------------------------------+

C+-----------------------------------------------------------------------+
C| Start of problem-specific initialisation                              |
C+-----------------------------------------------------------------------+

C+-----------------------------------------------------------------------+
C| END of problem-specific initialisation                                |
C+-----------------------------------------------------------------------+
      END
