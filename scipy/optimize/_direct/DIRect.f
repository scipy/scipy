C+-----------------------------------------------------------------------+
C| Program       : Direct.f                                              |
C| Last modified : 07-16-2001                                            |
C| Written by    : Joerg Gablonsky (jmgablon@unity.ncsu.edu)             |
C|                 North Carolina State University                       |
C|                 Dept. of Mathematics                                  |
C| DIRECT is a method to solve problems of the form:                     |
C|              min f: Q --> R,                                          |
C| where f is the function to be minimized and Q is an n-dimensional     |
C| hyperrectangle given by the the following equation:                   |
C|                                                                       |
C|       Q={ x : l(i) <= x(i) <= u(i), i = 1,...,n }.                    |
C| Note: This version of DIRECT can also handle hidden constraints. By   |
C|       this we mean that the function may not be defined over the whole|
C|       hyperrectangle Q, but only over a subset, which is not given    |
C|       analytically.                                                   |
C|                                                                       |
C| We now give a brief outline of the algorithm:                         |
C|                                                                       |
C|   The algorithm starts with mapping the hyperrectangle Q to the       |
C|   n-dimensional unit hypercube. DIRECT then samples the function at   |
C|   the center of this hypercube and at 2n more points, 2 in each       |
C|   coordinate direction. Uisng these function values, DIRECT then      |
C|   divides the domain into hyperrectangles, each having exactly one of |
C|   the sampling points as its center. In each iteration, DIRECT chooses|
C|   some of the existing hyperrectangles to be further divided.         |
C|   We provide two different strategies of how to decide which          |
C|   hyperrectangles DIRECT divides and several different convergence    |
C|   criteria.                                                           |
C|                                                                       |
C|   DIRECT was designed to solve problems where the function f is       |
C|   Lipschitz continues. However, DIRECT has proven to be effective on  |
C|   more complex problems than these.                                   |
C+-----------------------------------------------------------------------+

      SUBROUTINE Direct(fcn, x, n, eps, maxf, maxT, fmin, l, u,
     +                  algmethod, Ierror, logfilename,  
     +                  fglobal, fglper, volper, sigmaper,
     +                  iidata, iisize, ddata, idsize, cdata, icsize,
     +                  disp)

C+-----------------------------------------------------------------------+
C|    SUBROUTINE Direct                                                  |
C| On entry                                                              |
C|     fcn -- The argument containing the name of the user-supplied      |
C|            SUBROUTINE that returns values for the function to be      |
C|            minimized.                                                 |
C|       n -- The dimension of the problem.                              |
C|     eps -- Exceeding value. If eps > 0, we use the same epsilon for   |
C|            all iterations. If eps < 0, we use the update formula from |
C|            Jones:                                                     |
C|               eps = max(1.D-4*abs(fmin),epsfix),                      |
C|            where epsfix = abs(eps), the absolute value of eps which is|
C|            passed to the function.                                    |
C|    maxf -- The maximum number of function evaluations.                |
C|    maxT -- The maximum number of iterations.                          |
C|            Direct stops when either the maximum number of iterations  |
C|            is reached or more than maxf function-evalutions were made.|
C|       l -- The lower bounds of the hyperbox.                          |
C|       u -- The upper bounds of the hyperbox.                          |
C|algmethod-- Choose the method, that is either use the original method  |
C|            as described by Jones et.al. (0) or use our modification(1)|
C| logfile -- File-Handle for the logfile. DIRECT expects this file to be|
C|            opened and closed by the user outside of DIRECT. We moved  |
C|            this to the outside so the user can add extra informations |
C|            to this file before and after the call to DIRECT.          |
C| fglobal -- Function value of the global optimum. If this value is not |
C|            known (that is, we solve a real problem, not a testproblem)|
C|            set this value to -1.D100 and fglper (see below) to 0.D0.  |
C|  fglper -- Terminate the optimization when the percent error          |
C|                100(f_min - fglobal)/max(1,abs(fglobal)) < fglper.     |
C|  volper -- Terminate the optimization when the volume of the          |
C|            hyperrectangle S with f(c(S)) = fmin is less then volper   |
C|            percent of the volume of the original hyperrectangle.      |
C|sigmaper -- Terminate the optimization when the measure of the         |
C|            hyperrectangle S with f(c(S)) = fmin is less then sigmaper.|
C|                                                                       |
C| User data that is passed through without being changed:               |
C|  iidata -- Integer Array of user data. This array of size iisize is   |
C|            passed to the function to be optimized and can be used to  |
C|            transfer data to this function. The contents are not       |
C|            changed by DIRECT.                                         |
C|  iisize -- Size of array iidata.                                      |
C|   ddata -- DOUBLE PRECISION array of user data. See iidata.           |
C|  idsize -- Size of array ddata.                                       |
C|   cdata -- Character array. See iidata.                               |
C|  icsize -- Size of array ddata.                                       |
C|                                                                       |
C|  disp -- display progress or not                                      |
C|                                                                       |
C| On return                                                             |
C|                                                                       |
C|       x -- The final point obtained in the optimization process.      |
C|            X should be a good approximation to the global minimum     |
C|            for the function within the hyper-box.                     |
C|                                                                       |
C|    fmin -- The value of the function at x.                            |
C|  Ierror -- Error flag. If Ierror is lower 0, an error has occured. The|
C|            values of Ierror mean                                      |
C|            Fatal errors :                                             |
C|             -1   u(i) <= l(i) for some i.                             |
C|             -2   maxf is too large.                                   |
C|             -3   Initialization in DIRpreprc failed.                  |
C|             -4   Error in DIRSamplepoints, that is there was an error |
C|                  in the creation of the sample points.                |
C|             -5   Error in DIRSamplef, that is an error occured while  |
C|                  the function was sampled.                            |
C|             -6   Error in DIRDoubleInsert, that is an error occured   |
C|                  DIRECT tried to add all hyperrectangles with the same|
C|                  size and function value at the center. Either        |
C|                  increase maxdiv or use our modification (Jones = 1). |
C|            Termination values :                                       |
C|              1   Number of function evaluations done is larger then   |
C|                  maxf.                                                |
C|              2   Number of iterations is equal to maxT.               |
C|              3   The best function value found is within fglper of    |
C|                  the (known) global optimum, that is                  |
C|                   100(fmin - fglobal/max(1,|fglobal|))  < fglper.     |
C|                  Note that this termination signal only occurs when   |
C|                  the global optimal value is known, that is, a test   |
C|                  function is optimized.                               |
C|              4   The volume of the hyperrectangle with fmin at its    |
C|                  center is less than volper percent of the volume of  |
C|                  the original hyperrectangle.                         |
C|              5   The measure of the hyperrectangle with fmin at its   |
C|                  center is less than sigmaper.                        |
C|                                                                       |
C| SUBROUTINEs used :                                                    |
C|                                                                       |
C| DIRheader, DIRInitSpecific, DIRInitList, DIRpreprc, DIRInit, DIRChoose|
C| DIRDoubleInsert, DIRGet_I, DIRSamplepoints, DIRSamplef, DIRDivide     |
C| DIRInsertList, DIRreplaceInf, DIRWritehistbox, DIRsummary, Findareas  |
C|                                                                       |
C| Functions used :                                                      |
C|                                                                       |
C| DIRgetMaxdeep, DIRgetlevel                                            |
C+-----------------------------------------------------------------------+

      IMPLICIT NONE
C+-----------------------------------------------------------------------+
C| Parameters                                                            |
C+-----------------------------------------------------------------------+

C+-----------------------------------------------------------------------+
C| The maximum of function evaluations allowed.                          |
C| The maximum dept of the algorithm.                                    |
C| The maximum number of divisions allowed.                              |
C| The maximal dimension of the problem.                                 |
C+-----------------------------------------------------------------------+
      INTEGER maxfunc, maxdeep, maxdiv, MaxDim, mdeep
      PARAMETER (Maxfunc = 90000)
      PARAMETER (maxdeep = 600)
      PARAMETER (maxdiv = 3000)
      PARAMETER (MaxDim = 64)

C+-----------------------------------------------------------------------+
C| Global Variables.                                                     |
C+-----------------------------------------------------------------------+
      INTEGER JONES
      COMMON /directcontrol/ JONES

C+-----------------------------------------------------------------------+
C| EXTERNAL Variables.                                                   |
C+-----------------------------------------------------------------------+
      EXTERNAL fcn
      INTEGER n, maxf, maxT, algmethod, Ierror, logfile, dwrit, disp
      CHARACTER*(*) logfilename
Cf2py intent(in) logfilename
      DOUBLE PRECISION  x(n),fmin,eps,l(n),u(n)
      DOUBLE PRECISION fglobal, fglper, volper, sigmaper

C+-----------------------------------------------------------------------+
C| User Variables.                                                       |
C| These can be used to pass user defined data to the function to be     |
C| optimized.                                                            |
C+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      DOUBLE PRECISION ddata(idsize)
      Character*40 cdata(icsize)

C+-----------------------------------------------------------------------+
C| Place to define, if needed, some application-specific variables.      |
C| Note: You should try to use the arrays defined above for this.        |
C+-----------------------------------------------------------------------+

C+-----------------------------------------------------------------------+
C| End of application - specific variables !                             |
C+-----------------------------------------------------------------------+


C+-----------------------------------------------------------------------+
C| Internal variables :                                                  |
C|       f -- values of functions.                                       |
C|divfactor-- Factor used for termination with known global minimum.     |
C|  anchor -- anchors of lists with deepness i, -1 is anchor for list of |
C|            NaN - values.                                              |
C|       S -- List of potentially optimal points.                        |
C|   point -- lists                                                      |
C|    free -- first free position                                        |
C|       c -- midpoints of arrays                                        |
C|  thirds -- Precalculated values of 1/3^i.                             |
C|  levels -- Length of intervals.                                       |
C|  length -- Length of intervall (index)                                |
C|       t -- actual iteration                                           |
C|       j -- loop-variable                                              |
C| actdeep -- the actual minimal interval-length index                   |
C|  Minpos -- position of the actual minimum                             |
C|    file -- The filehandle for a datafile.                             |
C|  maxpos -- The number of intervalls, which are truncated.             |
C|    help -- A help variable.                                           |
C| numfunc -- The actual number of function evaluations.                 |
C|   file2 -- The filehandle for an other datafile.                      |
C|  ArrayI -- Array with the indexes of the sides with maximum length.   |
C|    maxi -- Number of directions with maximal side length.             |
C|    oops -- Flag which shows if anything went wrong in the             |
C|            initialisation.                                            |
C|   cheat -- Obsolete. If equal 1, we don't allow Ktilde > kmax.        |
C|  writed -- If writed=1, store final division to plot with Matlab.     |
C|   List2 -- List of indicies of intervalls, which are to be truncated. |
C|       i -- Another loop-variable.                                     |
C|actmaxdeep-- The actual maximum (minimum) of possible Interval length. |
C|  oldpos -- The old index of the minimum. Used to print only, if there |
C|            is a new minimum found.                                    |
C|  tstart -- The start of the outer loop.                               |
C|   start -- The postion of the starting point in the inner loop.       |
C| Newtosample -- The total number of points to sample in the inner loop.|
C|       w -- Array used to divide the intervalls                        |
C|    kmax -- Obsolete. If cheat = 1, Ktilde was not allowed to be larger|
C|            than kmax. If Ktilde > kmax, we set ktilde = kmax.         |
C|   delta -- The distance to new points from center of old hyperrec.    |
C|    pos1 -- Help variable used as an index.                            |
C| version -- Store the version number of DIRECT.                        |
C| oldmaxf -- Store the original function budget.                        |
C|increase -- Flag used to keep track if function budget was increased   |
C|            because no feasible point was found.                       |
C| freeold -- Keep track which index was free before. Used with          |
C|            SUBROUTINE DIRReplaceInf.                                  |
C|actdeep_div-- Keep track of the current depths for divisions.          |
C|    oldl -- Array used to store the original bounds of the domain.     |
C|    oldu -- Array used to store the original bounds of the domain.     |
C|  epsfix -- If eps < 0, we use Jones update formula. epsfix stores the |
C|            absolute value of epsilon.                                 |
C|iepschange-- flag iepschange to store if epsilon stays fixed or is     |
C|             changed.                                                  |
C|DIRgetMaxdeep-- Function to calculate the level of a hyperrectangle.   |
C|DIRgetlevel-- Function to calculate the level and stage of a hyperrec. |
C|    fmax -- Keep track of the maximum value of the function found.     |
C|Ifeasiblef-- Keep track if a feasible point has  been found so far.    |
C|             Ifeasiblef = 0 means a feasible point has been found,     |
C|             Ifeasiblef = 1 no feasible point has been found.          |
C|   dwrit -- Controls the level of output. So far not used a lot, set to|
C|            0. If its value is set to 2, more output, exspecially from |
C|            Subroutine DIRChoose.                                      |
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION  f(maxfunc,2), divfactor
      INTEGER anchor(-1:maxdeep),S(maxdiv,2)
      INTEGER point(maxfunc), free
      DOUBLE PRECISION  c(maxfunc,MaxDim)
      DOUBLE PRECISION  thirds(0:maxdeep),levels(0:maxdeep)
      INTEGER length(maxfunc,MaxDim),t,j,actdeep
      INTEGER Minpos,file,maxpos,help,numfunc,file2
      INTEGER ArrayI(MaxDim),maxi,oops,cheat,writed
      INTEGER List2(MaxDim,2),i,actmaxdeep,oldpos
      INTEGER tstart,start,Newtosample
      DOUBLE PRECISION  w(MaxDim),kmax, delta
      INTEGER pos1
C+-----------------------------------------------------------------------+
C| JG 09/25/00 Version counter.                                          |
C+-----------------------------------------------------------------------+
      INTEGER version
      INTEGER oldmaxf,increase, freeold
C+-----------------------------------------------------------------------+
C| JG 09/24/00 Add another actdeep to keep track of the current depths   |
C|             for divisions.                                            |
C+-----------------------------------------------------------------------+
      INTEGER actdeep_div
      DOUBLE PRECISION  oldl(MaxDim), oldu(MaxDim)
C+-----------------------------------------------------------------------+
C|JG 01/13/01 Added epsfix for epsilon update. If eps < 0, we use Jones  |
C|            update formula. epsfix stores the absolute value of epsilon|
C|            then. Also added flag iepschange to store if epsilon stays |
C|            fixed or is changed.                                       |
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION epsfix
      INTEGER iepschange

      INTEGER DIRGetMaxdeep, DIRgetlevel
C+-----------------------------------------------------------------------+
C| JG 01/22/01 fmax is used to keep track of the maximum value found.    |
C+-----------------------------------------------------------------------+
      DOUBLE PRECISION fmax
C+-----------------------------------------------------------------------+
C| JG 01/22/01 Ifeasiblef is used to keep track if a feasible point has  |
C|             been found so far. Ifeasiblef = 0 means a feasible point  |
C|             has been found, Ifeasiblef = 1 if not.                    |
C| JG 03/09/01 IInfeasible is used to keep track if an infeasible point  |
C|             has been found. IInfeasible > 0 means a infeasible point  |
C|             has been found, IInfeasible = 0 if not.                   |
C+-----------------------------------------------------------------------+
      INTEGER Ifeasiblef, IInfesiblef

C+-----------------------------------------------------------------------+
C+-----------------------------------------------------------------------+
C|                            Start of code.                             |
C+-----------------------------------------------------------------------+
C+-----------------------------------------------------------------------+
C|  Define and open the logfile                                          |
C+-----------------------------------------------------------------------+
      logfile    = 2
CC      open(logfile, file=logfilename)

C+-----------------------------------------------------------------------+
      writed = 0
      dwrit  = disp
      JONES  = algmethod
C+-----------------------------------------------------------------------+
C| Save the upper and lower bounds.                                      |
C+-----------------------------------------------------------------------+
      DO 150,i=1,n
        oldu(i) = u(i)
        oldl(i) = l(i)
150   CONTINUE

C+-----------------------------------------------------------------------+
C| Set version.                                                          |
C+-----------------------------------------------------------------------+
      version = 204
C+-----------------------------------------------------------------------+
C| Set parameters.                                                       |
C|    If cheat > 0, we do not allow \tilde{K} to be larger than kmax, and|
C|    set \tilde{K} to set value if necessary. Not used anymore.         |
C+-----------------------------------------------------------------------+
      cheat = 0
      kmax  = 1.D10
      mdeep = maxdeep
C+-----------------------------------------------------------------------+
C| Write the header of the logfile.                                      |
C+-----------------------------------------------------------------------+
CC      CALL DIRheader(logfile, version, x, n, eps, maxf, maxT, l, u,
CC     +               algmethod, maxfunc, maxdeep, fglobal, fglper,
CC     +               Ierror, epsfix, iepschange, volper, sigmaper,
CC     +               iidata, iisize, ddata, idsize, cdata,
CC     +               icsize)
C+-----------------------------------------------------------------------+
C| If an error has occured while writing the header (we do some checking |
C| of variables there), return to the main program.                      |
C+-----------------------------------------------------------------------+
CC      IF (Ierror .lt. 0) then
CC         RETURN
CC      END IF
C+-----------------------------------------------------------------------+
C| If the known global minimum is equal 0, we cannot divide by it.       |
C| Therefore we set it to 1. If not, we set the divisionfactor to the    |
C| absolute value of the global minimum.                                 |
C+-----------------------------------------------------------------------+
      IF (fglobal .EQ. 0.D0) then
         divfactor = 1.D0
      ELSE
         divfactor = abs(fglobal)
      END IF

C+-----------------------------------------------------------------------+
C| Start of application-specific initialisation.                         |
C+-----------------------------------------------------------------------+
      CALL DIRInitSpecific(x,n,
     +   iidata, iisize, ddata, idsize, cdata, icsize)
C+-----------------------------------------------------------------------+
C| End of application-specific initialisation.                           |
C+-----------------------------------------------------------------------+

C+-----------------------------------------------------------------------+
C| Save the budget given by the user. The variable maxf will be changed  |
C| if in the beginning no feasible points are found.                     |
C+-----------------------------------------------------------------------+
      oldmaxf = maxf
      increase = 0
C+-----------------------------------------------------------------------+
C| Initialiase the lists.                                                |
C+-----------------------------------------------------------------------+
      CALL DIRInitList(anchor,free,point,f,maxfunc,maxdeep)
C+-----------------------------------------------------------------------+
C| Call the routine to initialise the mapping of x from the n-dimensional|
C| unit cube to the hypercube given by u and l. If an error occured,     |
C| give out a error message and return to the main program with the error|
C| flag set.                                                             |
C| JG 07/16/01 Changed call to remove unused data.                       |
C+-----------------------------------------------------------------------+
      CALL DIRpreprc(u,l,n,l,u,oops)
      IF (oops .GT. 0) THEN
        Write(*,10005)
CC        Write(logfile,10005)
        IError = -3
        Return
      END IF
      tstart = 2
C+-----------------------------------------------------------------------+
C| Initialise the algorithm DIRECT.                                      |
C+-----------------------------------------------------------------------+
C+-----------------------------------------------------------------------+
C| Added variable to keep track of the maximum value found.              |
C+-----------------------------------------------------------------------+
      CALL DIRInit(f,fcn,c,length,actdeep,point,anchor,free,
     +   dwrit,logfile,ArrayI,maxI,List2,w,x,l,u,fmin,minpos,
     +   thirds,levels,maxfunc,maxdeep,n,MaxDim,fmax,Ifeasiblef,
     +   IInfesiblef, Ierror,
     +   iidata, iisize, ddata, idsize, cdata, icsize)
C+-----------------------------------------------------------------------+
C| Added error checking.                                                 |
C+-----------------------------------------------------------------------+
      IF (Ierror .lt. 0) then
         IF (Ierror .eq. -4) THEN
             IF (dwrit .gt. 0) THEN
                Write(*,10006)
             ENDIF
CC            Write(logfile,10006)
            return
         END IF
         IF (Ierror .eq. -5) THEN
            Write(*,10007)
CC            Write(logfile,10007)
            return
         END IF
      END IF
      
      numfunc = 1 + maxI + maxI
      actmaxdeep = 1
      oldpos = 0
      tstart = 2
C+-----------------------------------------------------------------------+
C| If no feasible point has been found, give out the iteration, the      |
C| number of function evaluations and a warning. Otherwise, give out     |
C| the iteration, the number of function evaluations done and fmin.      |
C+-----------------------------------------------------------------------+
      IF (Ifeasiblef .gt. 0) then
        IF (dwrit .gt. 0) THEN
            write(*,10012) tstart-1,numfunc
        ENDIF
CC        write(logfile,10012) t,numfunc
      ELSE
        IF (dwrit .gt. 0) THEN
            Write(*,10002) numfunc, fmin, fmax
        ENDIF
CC        write(logfile,10012) t,numfunc
CC        Write(logfile,10003) tstart-1,numfunc,fmin
      END IF
C+-----------------------------------------------------------------------+
C+-----------------------------------------------------------------------+
C| Main loop!                                                            |
C+-----------------------------------------------------------------------+
C+-----------------------------------------------------------------------+
      DO 10,t=tstart,MaxT
C+-----------------------------------------------------------------------+
C| Choose the sample points. The indices of the sample points are stored |
C| in the list S.                                                        |
C+-----------------------------------------------------------------------+
        actdeep = actmaxdeep
        CALL DIRChoose(anchor,S,maxdeep,f,fmin,eps,levels,maxpos,
     +     length,maxfunc,maxdeep,maxdiv,n,logfile,dwrit,cheat,kmax,
     +     Ifeasiblef)
C+-----------------------------------------------------------------------+
C| Add other hyperrectangles to S, which have the same level and the same|
C| function value at the center as the ones found above (that are stored |
C| in S). This is only done if we use the original DIRECT algorithm.     |
C| JG 07/16/01 Added Errorflag.                                          |
C+-----------------------------------------------------------------------+
        IF (algmethod .EQ. 0) THEN
          CALL DIRDoubleInsert(anchor, S, maxpos, point, f,
     +           maxdeep, maxfunc, maxdiv, Ierror)
          IF (Ierror .eq. -6) THEN
              Write(*,10020)
              Write(*,10021)
              Write(*,10022)
              Write(*,10023)
CC              Write(logfile,10020)
CC              Write(logfile,10021)
CC              Write(logfile,10022)
CC              Write(logfile,10023)
              return
          END IF
        ENDIF
        oldpos = minpos
C+-----------------------------------------------------------------------+
C| Initialise the number of sample points in this outer loop.            |
C+-----------------------------------------------------------------------+
        Newtosample = 0
        DO 20,j=1,maxpos
           actdeep = S(j,2)
C+-----------------------------------------------------------------------+
C| If the actual index is a point to sample, do it.                      |
C+-----------------------------------------------------------------------+
           IF (S(j,1) .GT. 0) THEN
C+-----------------------------------------------------------------------+
C| JG 09/24/00 Calculate the value delta used for sampling points.       |
C+-----------------------------------------------------------------------+
              actdeep_div = DIRGetmaxdeep(S(j,1),length,maxfunc,n)
              delta = thirds(actdeep_div+1)
              actdeep = S(j,2)
C+-----------------------------------------------------------------------+
C| If the current dept of division is only one under the maximal allowed |
C| dept, stop the computation.                                           |
C+-----------------------------------------------------------------------+
              IF (actdeep+1 .GE. mdeep) THEN
                 Write(*,10004)
CC                 write(logfile,10004)
                 Ierror = -6
                 GOTO 100
              END IF
              actmaxdeep = max(actdeep,actmaxdeep)
              help = S(j,1)
              IF (.NOT. (anchor(actdeep) .EQ. help)) THEN
                pos1 = anchor(actdeep)
                DO WHILE (.NOT. (point(pos1) .EQ. help))
                  pos1 = point(pos1)
                END DO
                point(pos1) = point(help)
              ELSE
                anchor(actdeep) = point(help)
              END IF
              IF (actdeep .lt. 0) then
                actdeep = f(help,1)
              END IF
C+-----------------------------------------------------------------------+
C| Get the Directions in which to decrease the intervall-length.         |
C+-----------------------------------------------------------------------+
              CALL DIRGet_I(length,help,ArrayI,maxI,n,maxfunc)
C+-----------------------------------------------------------------------+
C| Sample the function. To do this, we first calculate the points where  |
C| we need to sample the function. After checking for errors, we then do |
C| the actual evaluation of the function, again followed by checking for |
C| errors.                                                               |
C+-----------------------------------------------------------------------+
              CALL DIRSamplepoints(c,ArrayI,delta,help,
     +             start,length,dwrit,logfile,f,free,maxI,point,
     +             fcn,x,l,fmin,minpos,u,n,maxfunc,maxdeep,oops)
              IF (oops .GT. 0) THEN
                Write(*,10006)
CC                Write(logfile,10006)
                IError = -4
                return
              END IF
              Newtosample = newtosample + maxI
C+-----------------------------------------------------------------------+
C| JG 01/22/01 Added variable to keep track of the maximum value found.  |
C+-----------------------------------------------------------------------+
              CALL DIRSamplef(c,ArrayI,delta,help,start,
     +            length,dwrit,logfile,f,free,maxI,point,fcn,x,l,
     +            fmin,minpos,u,n,maxfunc,maxdeep,oops,fmax,
     +            Ifeasiblef,IInfesiblef,
     +            iidata, iisize, ddata, idsize, cdata, icsize)
              IF (oops .GT. 0) THEN
                Write(*,10007)
CC                Write(logfile,10007)
                IError = -5
                return
              END IF
C+-----------------------------------------------------------------------+
C| Divide the intervalls.                                                |
C+-----------------------------------------------------------------------+
              CALL DIRDivide(start,actdeep_div,length,point,
     +             ArrayI,help,List2,w,maxI,f,maxfunc,maxdeep,n)
C+-----------------------------------------------------------------------+
C| Insert the new intervalls into the list (sorted).                     |
C+-----------------------------------------------------------------------+
              CALL DIRInsertList(start,anchor,point,f,maxI,length,
     +                    maxfunc,maxdeep,n,help)
C+-----------------------------------------------------------------------+
C| Increase the number of function evaluations.                          |
C+-----------------------------------------------------------------------+
              numfunc = numfunc + maxI + maxI
           END IF
C+-----------------------------------------------------------------------+
C| End of main loop.                                                     |
C+-----------------------------------------------------------------------+
20      CONTINUE
C+-----------------------------------------------------------------------+
C| If there is a new minimum, show the actual iteration, the number of   |
C| function evaluations, the minimum value of f (so far) and the position|
C| in the array.                                                         |
C+-----------------------------------------------------------------------+
        IF (oldpos .LT. minpos) THEN
            IF (dwrit .gt. 0) THEN
                Write(*,10002) numfunc,fmin, fmax
            END IF
CC          Write(logfile,10003) t,numfunc,fmin
        END IF
C+-----------------------------------------------------------------------+
C| If no feasible point has been found, give out the iteration, the      |
C| number of function evaluations and a warning.                         |
C+-----------------------------------------------------------------------+
        IF (Ifeasiblef .gt. 0) then
          write(*,10012) t,numfunc
CC          write(logfile,10012) t,numfunc
        END IF
C+-----------------------------------------------------------------------+
C+-----------------------------------------------------------------------+
C|                       Termination Checks                              |
C+-----------------------------------------------------------------------+

C+-----------------------------------------------------------------------+
C| JG 01/22/01 Calculate the index for the hyperrectangle at which       |
C|             fmin is assumed. We then calculate the volume of this     |
C|             hyperrectangle and store it in delta. This delta can be   |
C|             used to stop DIRECT once the volume is below a certain    |
C|             percentage of the original volume. Since the original     |
C|             is 1 (scaled), we can stop once delta is below a certain  |
C|             percentage, given by volper.                              |
C+-----------------------------------------------------------------------+
        Ierror = Jones
        Jones = 0
        actdeep_div = DIRGetlevel(minpos,length,maxfunc,n)
        Jones = Ierror
C+-----------------------------------------------------------------------+
C| JG 07/16/01 Use precalculated values to calculate volume.             |
C+-----------------------------------------------------------------------+
        delta = thirds(actdeep_div)*100
        IF (delta .LE. volper) THEN
           Ierror = 4
           Write(*,10011) delta, volper
CC           Write(logfile,10011) delta, volper
           GOTO 100
        END IF
C+-----------------------------------------------------------------------+
C| JG 01/23/01 Calculate the measure for the hyperrectangle at which     |
C|             fmin is assumed. If this measure is smaller then sigmaper,|
C|             we stop DIRECT.                                           |
C+-----------------------------------------------------------------------+
        actdeep_div = DIRGetlevel(minpos,length,maxfunc,n)
        delta = levels(actdeep_div)
        IF (delta .LE. sigmaper) THEN
           Ierror = 5
           Write(*,10013) delta, sigmaper
CC           Write(logfile,10013) delta, sigmaper
           GOTO 100
        END IF
C+-----------------------------------------------------------------------+
C| If the best found function value is within fglper of the (known)      |
C| global minimum value, terminate. This only makes sense if this optimal|
C| value is known, that is, in test problems.                            |
C+-----------------------------------------------------------------------+
        IF ((100*(fmin - fglobal)/divfactor) .LE. fglper)
     +   THEN
           Ierror = 3
           Write(*,10010)
C           Write(logfile,10010)
           GOTO 100
        END IF
C+-----------------------------------------------------------------------+
C| Find out if there are infeasible points which are near feasible ones. |
C| If this is the case, replace the function value at the center of the  |
C| hyper rectangle by the lowest function value of a nearby function.    |
C| If no infeasible points exist (IInfesiblef = 0), skip this.           |
C+-----------------------------------------------------------------------+
        IF (IInfesiblef .gt. 0) THEN
           CALL DIRreplaceInf(free,freeold,f,c,thirds,length,anchor,
     +       point,u,l,maxfunc,maxdeep,maxdim,n,logfile, fmax)
        ENDIF
        freeold = free
C+-----------------------------------------------------------------------+
C| If iepschange = 1, we use the epsilon change formula from Jones.      |
C+-----------------------------------------------------------------------+
        IF (iepschange .eq. 1) then
           eps = max(1.D-4*abs(fmin),epsfix)
        END IF
C+-----------------------------------------------------------------------+
C| If no feasible point has been found yet, set the maximum number of    |
C| function evaluations to the number of evaluations already done plus   |
C| the budget given by the user.                                         |
C| If the budget has already be increased, increase it again. If a       |
C| feasible point has been found, remark that and reset flag. No further |
C| increase is needed.                                                   |
C+-----------------------------------------------------------------------+
        IF (increase .eq. 1) then        
           maxf = numfunc + oldmaxf
           IF (Ifeasiblef .eq. 0) then
CC             write(logfile,10031) maxf
             increase = 0
           END IF
        END IF
C+-----------------------------------------------------------------------+
C| Check if the number of function evaluations done is larger than the   |
C| allocated budget. If this is the case, check if a feasible point was  |
C| found. If this is a case, terminate. If no feasible point was found,  |
C| increase the budget and set flag increase.                            |
C+-----------------------------------------------------------------------+
        IF (numfunc .GT. maxf) THEN
           IF (Ifeasiblef .eq. 0) then
              Ierror = 1
              IF (dwrit .gt. 0) THEN
                  Write(*,10008)
              END IF              
CC              Write(logfile,10008)
              GOTO 100
           ELSE
              increase = 1
CC              write(logfile,10030) numfunc
              maxf = numfunc+ oldmaxf
           END IF
        END IF
10    CONTINUE
C+-----------------------------------------------------------------------+
C+-----------------------------------------------------------------------+
C| End of main loop.                                                     |
C+-----------------------------------------------------------------------+
C+-----------------------------------------------------------------------+

C+-----------------------------------------------------------------------+
C| The algorithm stopped after maxT iterations.                          |
C+-----------------------------------------------------------------------+
      Ierror = 2
      IF (dwrit .gt. 0) THEN
          Write(*,10009)
      END IF
CC      Write(logfile,10009)

100   CONTINUE
C+-----------------------------------------------------------------------+
C| Store the position of the minimum in x.                               |
C+-----------------------------------------------------------------------+
      DO 50,i=1,n
         x(i) = c(Minpos,i)*l(i)+l(i)*u(i)
         u(i) = oldu(i)
         l(i) = oldl(i)
50    CONTINUE
C+-----------------------------------------------------------------------+
C| Store the number of function evaluations in maxf.                     |
C+-----------------------------------------------------------------------+
      maxf = numfunc

C+-----------------------------------------------------------------------+
C| If needed, save the final division in a file for use with Matlab.     |
C+-----------------------------------------------------------------------+
      writed = 0
      IF (writed .EQ. 1) THEN
        file  = 12
        file2 =  0
        CALL DIRWritehistbox(point,f,thirds,c,anchor,actdeep,file,
     +               l,u,file2,maxfunc,maxdeep,n,MaxDim,length)
      END IF

C+-----------------------------------------------------------------------+
C| Give out a summary of the run.                                        |
C+-----------------------------------------------------------------------+
      CALL DIRsummary(logfile,x,l,u,n,fmin,fglobal, numfunc, Ierror)

C+-----------------------------------------------------------------------+
C| Close the logfile.                                                    |
C+-----------------------------------------------------------------------+
CC      close(logfile)
      
C+-----------------------------------------------------------------------+
C| Format statements.                                                    |
C+-----------------------------------------------------------------------+
10002  FORMAT(i5," & ",f18.10," & ",f18.10," \\\\ ")
10003  FORMAT(i5,'       ', i5,'     ',f18.10)
10004  FORMAT('WARNING : Maximum number of levels reached. Increase
     +        maxdeep.')
10005  FORMAT('WARNING : Initialisation in DIRpreprc failed.')
10006  FORMAT('WARNING : Error occured in routine DIRsamplepoints.')
10007  FORMAT('WARNING : Error occured in routine DIRsamplef.')
10008  FORMAT('DIRECT stopped: numfunc >= maxf.')
10009  FORMAT('DIRECT stopped: maxT iterations.')
10010  FORMAT('DIRECT stopped: fmin within fglper of global minimum.')
10011  FORMAT('DIRECT stopped: Volume of S_min is ',d8.2,
     +       '% < ',d8.2,'% of the original volume.')
10012  FORMAT('No feasible point found in ',I4,' iterations ',
     +       'and ',I5,' function evaluations.')
10013  FORMAT('DIRECT stopped: Measure of S_min = ',d8.2,' < '
     +       ,d8.2,'.')
10020  FORMAT('WARNING : Capacity of array S in DIRDoubleInsert'
     +        ' reached. Increase maxdiv.')
10021  FORMAT('This means that there are a lot of hyperrectangles')
10022  FORMAT('with the same function value at the center. We')
10023  FORMAT('suggest to use our modification instead (Jones = 1)')

10030 FORMAT('DIRECT could not find a feasible ',
     +        'point after ', I5, ' function evaluations. ',
     +        'DIRECT continues until a feasible point is found.')
10031 FORMAT('DIRECT found a feasible point. ',
     +         'The adjusted budget is now set to ',I5,'.')
      END
