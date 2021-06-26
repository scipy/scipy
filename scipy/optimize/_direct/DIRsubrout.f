C+-----------------------------------------------------------------------+
C| Program       : Direct.f (subfile DIRsubrout.f)                       |
C| Last modified : 07-16-2001                                            |
C| Written by    : Joerg Gablonsky                                       |
C| Subroutines used by the algorithm DIRECT.                             |
C+-----------------------------------------------------------------------+

C+-----------------------------------------------------------------------+
C|    SUBROUTINE DIRChoose                                               |
C|    Decide, which is the next sampling point.                          |
C|    Changed 09/25/00 JG                                                |
C|         Added maxdiv to call and changed S to size maxdiv.            |
C|    Changed 01/22/01 JG                                                |
C|         Added Ifeasiblef to call to keep track if a feasible point has|
C|         been found.                                                   |
C|    Changed 07/16/01 JG                                                |
C|         Changed if statement to prevent run-time errors.              |                                                   |
C+-----------------------------------------------------------------------+
      SUBROUTINE DIRChoose(anchor,S,actdeep,f,fmin,eps,thirds,maxpos,
     +      length,maxfunc,maxdeep,maxdiv,n,logfile,dwrit,cheat,kmax,
     +      Ifeasiblef)

      IMPLICIT None
      Double Precision MaxLower
      Parameter (MaxLower = 1.D20)
      Integer maxfunc,maxdeep,n, maxdiv
      Integer Anchor(-1:maxdeep),S(maxdiv,2),length(maxfunc,n)
      Integer actdeep,dwrit
      Integer DIRGetlevel
      Double Precision f(maxfunc,2),fmin,eps,thirds(0:maxfunc)
      Integer maxpos,i,j,k,i_,j_,logfile,cheat
      Double Precision help2,helplower,helpgreater,kmax
      Integer novalue,novaluedeep
      Integer Ifeasiblef

        helplower = MaxLower
        helpgreater = 0.D0
        k = 1
        if (Ifeasiblef .ge. 1) then
           DO 1001,j=0,actdeep
              IF (anchor(j) .GT. 0) THEN
                 S(k,1) = anchor(j)
                 S(k,2) = DIRGetlevel(S(k,1),length,maxfunc,n)
                 goto 12
              END IF
1001       CONTINUE
12         k = k + 1
           maxpos = 1
           return
        else
           DO 10,j=0,actdeep
              IF (anchor(j) .GT. 0) THEN
                 S(k,1) = anchor(j)
                 S(k,2) = DIRGetlevel(S(k,1),length,maxfunc,n)
                 k = k + 1
              END IF
10         CONTINUE
        END IF
        novalue = 0
        if (anchor(-1) .gt. 0) then
          novalue = anchor(-1)
          novaluedeep = DIRGetlevel(novalue,length,maxfunc,n)
        end if
        maxpos = k - 1
       DO 11,j=k-1,maxdeep
         S(k,1) = 0
11     CONTINUE
       DO 40,j=maxpos,1,-1
         helplower = Maxlower
         helpgreater = 0.D0
         j_ = S(j,1)
         DO 30,i=1,j-1
           i_ = S(i,1)
C+-----------------------------------------------------------------------+
C| JG 07/16/01 Changed IF statement into two to prevent run-time errors  |
C|             which could occur if the compiler checks the second       |
C|             expression in an .AND. statement although the first       |
C|             statement is already not true.                            |
C+-----------------------------------------------------------------------+
           IF ((i_ .GT. 0) .AND. .NOT. (i .EQ. j)) THEN
             IF (f(i_,2) .le. 1.D0) THEN
               help2 = thirds(S(i,2)) - thirds(S(j,2))
               help2 = (f(i_,1) - f(j_,1))/help2
               IF (help2 .LE. 0.D0) THEN
CC                 IF (dwrit .EQ. 2) THEN
CC                   Write(logfile,*) "thirds > 0,help2 <= 0"
CC                 END IF
                 GOTO 60
               END IF
               IF (help2 .LT. helplower) THEN
CC                 IF (dwrit .EQ. 2) THEN
CC                   Write(logfile,*) "helplower = ",help2
CC                 END IF
                 helplower = help2
               END IF
               END IF
           END IF
30       CONTINUE
         DO 31,i=j+1,maxpos
           i_ = S(i,1)
C+-----------------------------------------------------------------------+
C| JG 07/16/01 Changed IF statement into two to prevent run-time errors  |
C|             which could occur if the compiler checks the second       |
C|             expression in an .AND. statement although the first       |
C|             statement is already not true.                            |
C+-----------------------------------------------------------------------+
           IF ((i_ .GT. 0) .AND. .NOT. (i .EQ. j)) THEN
             IF (f(i_,2) .le. 1.D0) THEN 
               help2 = thirds(S(i,2)) - thirds(S(j,2))
               help2 = (f(i_,1) - f(j_,1))/help2
               IF (help2 .LE. 0.D0) THEN
CC                 IF (dwrit .EQ. 2) THEN
CC                   Write(logfile,*) "thirds < 0,help2 <= 0"
CC                 END IF
                 GOTO 60
               END IF
               IF (help2 .GT. helpgreater) THEN
CC                 IF (dwrit .EQ. 2) THEN
CC                   Write(logfile,*) "helpgreater = ",help2
CC                 END IF
                 helpgreater = help2
               END IF
             END IF
           END IF
31       CONTINUE
         IF ((helplower .GT. Maxlower) .AND. 
     +       (helpgreater .gt. 0)) THEN
           helplower = helpgreater
           helpgreater = helpgreater - 1.D0
         END IF
         IF (helpgreater .LE. helplower) THEN
           IF ((cheat .EQ. 1) .AND. (helplower .GT. Kmax)) THEN
             helplower = Kmax
           END IF
           IF ((f(j_,1) - helplower * thirds(S(j,2))) .GT. 
     +        (fmin - eps*abs(fmin))) THEN
CC              IF (dwrit .EQ. 2) THEN
CC                Write(logfile,*) "> fmin - eps|fmin|"
CC              END IF
              GOTO 60
            END IF
         ELSE
CC           IF (dwrit .EQ. 2) THEN
CC           Write(logfile,*) "helpgreater > helplower",helpgreater,
CC     +            helplower,helpgreater - helplower
CC           END IF
           GOTO 60
         END IF
         GOTO 40
60       S(j,1) = 0
40     CONTINUE
       if (novalue .gt. 0) then
          maxpos = maxpos + 1
          S(maxpos,1) = novalue
          S(maxpos,2) = novaluedeep
       end if
      END

C+-----------------------------------------------------------------------+
C| INTEGER Function GetmaxDeep                                           |
C| function to get the maximal length (1/length) of the n-dimensional    |
C| rectangle with midpoint pos.                                          |
C|                                                                       |
C| On Return :                                                           |
C|    the maximal length                                                 |
C|                                                                       |
C| pos     -- the position of the midpoint in the array length           |
C| length  -- the array with the dimensions                              |
C| maxfunc -- the leading dimension of length                            |
C| n       -- the dimension of the problem                               |
C|                                                                       |
C+-----------------------------------------------------------------------+

      INTEGER Function DIRGetmaxDeep(pos,length,maxfunc,n)
      IMPLICIT None
      INTEGER pos,maxfunc,n,length(maxfunc,n),help,i

      help = length(pos,1)
      DO 10,i = 2,n
        help = min(help, length(pos,i))
10    CONTINUE            
      DIRGetMaxDeep = help
      END

C+-----------------------------------------------------------------------+
C| INTEGER Function DIRGetlevel                                          |
C| Returns the level of the hyperrectangle. Depending on the value of the|
C| global variable JONES. IF JONES equals 0, the level is given by       |
C|               kN + p, where the rectangle has p sides with a length of|
C|             1/3^(k+1), and N-p sides with a length of 1/3^k.          |
C| If JONES equals 1, the level is the power of 1/3 of the length of the |
C| longest side hyperrectangle.                                          |
C|                                                                       |
C| On Return :                                                           |
C|    the maximal length                                                 |
C|                                                                       |
C| pos     -- the position of the midpoint in the array length           |
C| length  -- the array with the dimensions                              |
C| maxfunc -- the leading dimension of length                            |
C| n       -- the dimension of the problem                               |
C|                                                                       |
C+-----------------------------------------------------------------------+

      INTEGER Function DIRGetlevel(pos,length,maxfunc,n)
      IMPLICIT None
      INTEGER pos,maxfunc,n,length(maxfunc,n),help,i,p,k
C JG 09/15/00 Added variable JONES (see above)
      Integer JONES
      COMMON /directcontrol/ JONES

      IF (JONES .eq. 0) THEN
        help = length(pos,1)
        k = help
        p = 1
        DO 100,i = 2,n
          IF (length(pos,i) .LT. k) THEN
            k = length(pos,i)
          END IF
          IF (length(pos,i) .EQ. help) THEN
             p = p + 1
          END IF
100     CONTINUE
        IF (k .EQ. help) THEN            
           DIRGetLevel = k*n + n-p
        ELSE
           DIRGetLevel = k*n + p
        END IF
      ELSE
        help = length(pos,1)
        DO 10,i = 2,n
          IF (length(pos,i) .LT. help) THEN
            help = length(pos,i)
          END IF
10      CONTINUE            
        DIRGetLevel = help
      END IF
      END

C+-----------------------------------------------------------------------+
C|    SUBROUTINE DIRDoubleInsert                                         |
C|      Routine to make sure that if there are several potential optimal |
C|      hyperrectangles of the same level (i.e. hyperrectangles that have|
C|      the same level and the same function value at the center), all of|
C|      them are divided. This is the way as originally described in     |
C|      Jones et.al.                                                     |
C| JG 07/16/01 Added errorflag to calling sequence. We check if more     |
C|             we reach the capacity of the array S. If this happens, we |
C|             return to the main program with an error.                 |
C+-----------------------------------------------------------------------+
      SUBROUTINE DIRDoubleInsert(anchor, S, maxpos, point, f, 
     +           maxdeep, maxfunc, maxdiv, IError)
      IMPLICIT None
      Integer maxpos, maxdeep, maxfunc, maxdiv
      Integer anchor(-1:maxdeep),S(maxdiv,2)
      Integer point(maxfunc)
      Double Precision  f(maxfunc,2)
      
C+-----------------------------------------------------------------------+
C| JG 07/16/01 Added flag to prevent run time-errors on some systems.    |
C+-----------------------------------------------------------------------+
      Integer iflag, IError
      Integer oldmaxpos, i, pos, help, actdeep
      
      oldmaxpos = maxpos
      DO 10, i = 1,oldmaxpos
        IF (S(i,1) .GT. 0) THEN
          actdeep = S(i,2)
          help = anchor(actdeep)
          pos = point(help)
          iflag = 0
C+-----------------------------------------------------------------------+
C| JG 07/16/01 Added flag to prevent run time-errors on some systems. On |
C|             some systems the second conditions in an AND statement is |
C|             evaluated even if the first one is already not true.      |
C+-----------------------------------------------------------------------+
          DO WHILE ((pos .GT. 0) .AND. (iflag .eq. 0))
             if (f(pos,1) - f(help,1) .LE. 1.D-13) then
               if (maxpos .lt. maxdiv) then
                 maxpos = maxpos + 1
                 S(maxpos,1) = pos
                 S(maxpos,2) = actdeep
                 pos = point(pos)
               else
C+-----------------------------------------------------------------------+
C| JG 07/16/01 Maximum number of elements possible in S has been reached!|
C+-----------------------------------------------------------------------+
                 IError = -6
                 return
               end if
             else
               iflag = 1
             end if
          END DO
        END IF 
10    CONTINUE
      END

C+-----------------------------------------------------------------------+
C| JG Added 09/25/00                                                     |
C|                       SUBROUTINE DIRreplaceInf                        |
C|                                                                       |
C| Find out if there are infeasible points which are near feasible ones. |
C| If this is the case, replace the function value at the center of the  |
C| hyper rectangle by the lowest function value of a nearby function.    |
C+-----------------------------------------------------------------------+
      SUBROUTINE DIRreplaceInf(free,freeold,f,c,thirds,length,anchor,
     +       point,c1,c2,maxfunc,maxdeep,maxdim,n,logfile, fmax)
      Implicit None
      Integer maxfunc, maxdeep, maxdim, n, free, freeold, logfile
      Double Precision  f(maxfunc,2)
      Integer anchor(-1:maxdeep)
      Integer point(maxfunc)
      Double Precision  c1(n),c2(n)
      Double Precision  c(maxfunc,MaxDim)
      Double Precision  thirds(0:maxdeep)
      Integer length(maxfunc,MaxDim)
      
      Integer LMaxDim
      PARAMETER (LMaxDim = 32)
      
      Double Precision sidelength
      Double Precision a(LmaxDim),b(LmaxDim),x(LmaxDim)
      Integer i,j,k,l, help, Isinbox
      Integer DIRgetmaxdeep
C+-----------------------------------------------------------------------+
C| JG 01/22/01 Added variable to keep track of the maximum value found.  |
C+-----------------------------------------------------------------------+
      Double Precision fmax
      
      DO 10, i = 1,free-1
         if (f(i,2) .gt. 0) then
C+-----------------------------------------------------------------------+
C| Get the maximum side length of the hyper rectangle and then set the   |
C| new side length to this lengths times the growth factor.              |
C+-----------------------------------------------------------------------+
            help = DIRgetmaxdeep(i,length,maxfunc,n)
            sidelength = thirds(help)*2.D0
C+-----------------------------------------------------------------------+
C| Set the Center and the upper and lower bounds of the rectangles.      |  
C+-----------------------------------------------------------------------+
            do 20, j = 1,n
               sidelength = thirds(length(i,j))               
               a(j) = c(i,j) - sidelength
               b(j) = c(i,j) + sidelength
20          continue
C+-----------------------------------------------------------------------+
C| The function value is reset to 'Inf', since it may have been changed  |
C| in an earlier iteration and now the feasible point which was close    |
C| is not close anymore (since the hyper rectangle surrounding the       |
C| current point may have shrunk).                                       |
C+-----------------------------------------------------------------------+            
            f(i,1) = 1.0E+6
            f(i,2) = 2.D0
C+-----------------------------------------------------------------------+
C| Check if any feasible point is near this infeasible point.            |
C+-----------------------------------------------------------------------+
            DO 30, k = 1,free-1
C+-----------------------------------------------------------------------+
C| If the point k is feasible, check if it is near.                      |
C+-----------------------------------------------------------------------+
               if (f(k,2) .eq. 0) then
C+-----------------------------------------------------------------------+
C| Copy the coordinates of the point k into x.                           |
C+-----------------------------------------------------------------------+
                  DO 40, l = 1,n
                     x(l) = c(k,l)
40                continue
C+-----------------------------------------------------------------------+
C| Check if the point k is near the infeasible point, if so, replace the |
C| value at 
C+-----------------------------------------------------------------------+
                  if (Isinbox(x,a,b,n,Lmaxdim) .eq. 1) then
                     f(i,1) = min(f(i,1), f(k,1))
                     f(i,2) = 1.D0
                  end if
               end if               
30          continue
            if (f(i,2) .eq. 1.0D0) then
               f(i,1) = f(i,1) + 1.0E-6*abs(f(i,1))
               do 200,l=1,n
                  x(l)=c(i,l)*c1(l)+c(i,l)*c2(l)
200            continue
               CALL DIRResortlist(i,anchor,f,point,length,n,maxfunc,
     +              maxdim,maxdeep,logfile)
            else
C+-----------------------------------------------------------------------+
C| JG 01/22/01                                                           |
C| Replaced fixed value for infeasible points with maximum value found,  |
C| increased by 1.                                                       |
C+-----------------------------------------------------------------------+
              if (.NOT. (fmax .eq. f(i,1))) then
                f(i,1) = max(fmax + 1.0D0,f(i,1))
              end if
            end if
         end if
10    continue
1000  format(20f18.8) 
      END 
      
C+-----------------------------------------------------------------------+
C| JG Added 09/25/00                                                     |
C|                                                                       |
C|                       SUBROUTINE DIRResortlist                        |
C|                                                                       |
C| Resort the list so that the infeasible point is in the list with the  |
C| replaced value.                                                       |
C+-----------------------------------------------------------------------+
      SUBROUTINE DIRResortlist(replace,anchor,f,point,length,n,maxfunc,
     +           maxdim,maxdeep,logfile)
      Implicit None
      Integer maxfunc, maxdim, maxdeep, n, logfile
      Integer replace
      Double Precision  f(maxfunc,2)
      Integer anchor(-1:maxdeep)
      Integer point(maxfunc)
      Integer length(maxfunc,MaxDim)
      Integer start, l, i, pos
      Integer DIRgetlevel
      
C+-----------------------------------------------------------------------+
C| Get the length of the hyper rectangle with infeasible mid point and   |
C| Index of the corresponding list.                                      |
C+-----------------------------------------------------------------------+
C JG 09/25/00 Replaced with DIRgetlevel
C      l = DIRgetmaxDeep(replace,length,maxfunc,n)
      l = DIRgetlevel(replace,length,maxfunc,n)
      start = anchor(l)
C+-----------------------------------------------------------------------+
C| If the hyper rectangle with infeasibel midpoint is already the start  |
C| of the list, give out message, nothing to do.                         |
C+-----------------------------------------------------------------------+
      if (replace .eq. start) then
CC         write(logfile,*) 'No resorting of list necessarry, since new ',
CC     + 'point is already anchor of list .',l
      else
C+-----------------------------------------------------------------------+
C| Take the hyper rectangle with infeasible midpoint out of the list.    |
C+-----------------------------------------------------------------------+
         pos = start
         do 10, i = 1,maxfunc
            if (point(pos) .eq. replace) then
               point(pos) = point(replace)
               goto 20
            else
               pos = point(pos)
            end if
            if (pos .eq. 0) then
               write(logfile,*) 'Error in DIRREsortlist: We went ',
     + 'through the whole list and could not find the point to 
     +  replace!!'
               goto 20
            end if
10       continue
C+-----------------------------------------------------------------------+
C| If the anchor of the list has a higher value than the value of a      |
C| nearby point, put the infeasible point at the beginning of the list.  |
C+-----------------------------------------------------------------------+
20       if (f(start,1) .gt. f(replace,1)) then
            anchor(l) = replace
            point(replace) = start
C            write(logfile,*) 'Point is replacing current anchor for '
C     +             , 'this list ',l,replace,start
         else
C+-----------------------------------------------------------------------+
C| Insert the point into the list according to its (replaced) function   |
C| value.                                                                |
C+-----------------------------------------------------------------------+
            pos = start
            do 30, i = 1,maxfunc
C+-----------------------------------------------------------------------+
C| The point has to be added at the end of the list.                     |
C+-----------------------------------------------------------------------+
               if (point(pos) .eq. 0) then
                  point(replace) = point(pos)
                  point(pos) = replace
C                  write(logfile,*) 'Point is added at the end of the '
C     +             , 'list ',l, replace
                  goto 40
               else
                  if (f(point(pos),1) .gt. f(replace,1)) then
                     point(replace) = point(pos)
                     point(pos) = replace
C                     write(logfile,*) 'There are points with a higher '
C     +               ,'f-value in the list ',l,replace, pos
                     goto 40
                  end if
                  pos = point(pos)
               end if
30          continue
40          pos = pos
         end if
      end if   
      end
C+-----------------------------------------------------------------------+
C|    SUBROUTINE DIRInsertList                                           |
C|    Changed 02-24-2000                                                 |
C|      Got rid of the distinction between feasible and infeasible points|
C|      I could do this since infeasible points get set to a high        |
C|      function value, which may be replaced by a function value of a   |
C|      nearby function at the end of the main loop.                     |
C+-----------------------------------------------------------------------+
      SUBROUTINE DIRInsertList(new,anchor,point,f,maxI,
     +                      length,maxfunc,maxdeep,n,samp)
      IMPLICIT None
      INTEGER maxfunc,maxdeep,j,maxI,n,samp
      INTEGER pos1,pos2,pos,new,deep,anchor(-1:maxdeep)
      INTEGER point(maxfunc),length(maxfunc,n)
C JG 09/24/00 Changed this to Getlevel
      INTEGER DIRGetlevel
      Double Precision  f(maxfunc,2)

      DO 10,j = 1,maxI
        pos1 = new
        pos2 = point(pos1)
        new = point(pos2)
C JG 09/24/00 Changed this to Getlevel
C        deep = DIRGetMaxdeep(pos1,length,maxfunc,n)
        deep = DIRGetlevel(pos1,length,maxfunc,n)
        IF (anchor(deep) .EQ. 0) THEN
           IF (f(pos2,1) .LT. f(pos1,1)) THEN
              anchor(deep) = pos2
              point(pos2) = pos1
              point(pos1) = 0
           ELSE
              anchor(deep) = pos1
              point(pos2) = 0
           END IF
        ELSE
           pos = anchor(deep)
           IF (f(pos2,1) .LT. f(pos1,1)) THEN
              IF (f(pos2,1) .LT. f(pos,1)) THEN
                 anchor(deep) = pos2
C JG 08/30/00 Fixed bug. Sorting was not correct when 
C      f(pos2,1) < f(pos1,1) < f(pos,1)
                 IF (f(pos1,1) .LT. f(pos,1)) THEN
                    point(pos2) = pos1
                    point(pos1) = pos
                 ELSE
                    point(pos2) = pos
                    CALL DIRInsert(pos,pos1,point,f,maxfunc)
                 END IF
              ELSE
                 CALL DIRInsert(pos,pos2,point,f,maxfunc)
                 CALL DIRInsert(pos,pos1,point,f,maxfunc)
              END IF
           ELSE
              IF (f(pos1,1) .LT. f(pos,1)) THEN
C JG 08/30/00 Fixed bug. Sorting was not correct when 
C      f(pos1,1) < f(pos2,1) < f(pos,1)
                 anchor(deep) = pos1
                 IF (f(pos,1) .LT. f(pos2,1)) THEN
                    point(pos1) = pos
                    CALL DIRInsert(pos,pos2,point,f,maxfunc)
                 ELSE
                    point(pos1) = pos2
                    point(pos2) = pos
                 END IF
              ELSE
                 CALL DIRInsert(pos,pos1,point,f,maxfunc)
                 CALL DIRInsert(pos,pos2,point,f,maxfunc)
              END IF              
           END IF
        END IF
10    CONTINUE
C JG 09/24/00 Changed this to Getlevel
C      deep = DIRGetMaxdeep(samp,length,maxfunc,n)
      deep = DIRGetlevel(samp,length,maxfunc,n)
      pos = anchor(deep)
      IF (f(samp,1) .LT. f(pos,1)) THEN
         anchor(deep) = samp
         point(samp) = pos
      ELSE
         CALL DIRInsert(pos,samp,point,f,maxfunc)
      END IF
      END

C+-----------------------------------------------------------------------+
C|    SUBROUTINE DIRInsertList2  (Old way to do it.)                     |
C+-----------------------------------------------------------------------+
      SUBROUTINE DIRInsertList_2(start,j,k,List2,w,maxI,n)
      IMPLICIT None
      INTEGER start,n,j,k
      INTEGER List2(n,2)
      Double Precision w(n)
      INTEGER pos,i,maxI

        pos = start
        IF (start .EQ. 0) THEN
          List2(j,1) = 0
          start = j
          GOTO 50
        END IF
        IF (w(start) .GT. w(j)) THEN
          List2(j,1) = start
          start = j
        ELSE
          DO 10,i=1,maxI
            IF (List2(pos,1) .EQ. 0) THEN
              List2(j,1) = 0
              List2(pos,1) =  j
              GOTO 50
            ELSE
              IF (w(j) .LT. w(List2(pos,1))) THEN
                List2(j,1) = List2(pos,1)
                List2(pos,1) = j
                GOTO 50
              END IF
            END IF
            pos = List2(pos,1)
10        CONTINUE
        END IF
50     List2(j,2) = k 

      END

C+-----------------------------------------------------------------------+
C|    SUBROUTINE DIRSearchmin                                            |
C|    Search for the minimum in the list.                                !
C+-----------------------------------------------------------------------+
      SUBROUTINE DIRSearchmin(start,List2,pos,k,n)
      IMPLICIT None
      Integer start,pos,k,n
      INTEGER List2(n,2)

        k = start
        pos = List2(start,2)
        start = List2(start,1)
      END

C+-----------------------------------------------------------------------+
C|    SUBROUTINE DIRInit                                                 |
C|    Initialise all needed variables and do the first run of the        |
C|    algorithm.                                                         |
C|    Changed 02/24/2000                                                 |
C|       Changed fcn Double precision to fcn external!                   |
C|    Changed 09/15/2000                                                 |
C|       Added distinction between Jones way to characterize rectangles  |
C|       and our way. Common variable JONES controls which way we use.   |
C|          JONES = 0    Jones way (Distance from midpoint to corner)    |
C|          JONES = 1    Our way (Length of longest side)                |
C|    Changed 09/24/00                                                   |
C|       Added array levels. Levels contain the values to characterize   |
C|       the hyperrectangles.                                            |
C|    Changed 01/22/01                                                   |
C|       Added variable fmax to keep track of maximum value found.       |
C|       Added variable Ifeasiblef to keep track if feasibel point has   |
C|       been found.                                                     |
C|    Changed 01/23/01                                                   |
C|       Added variable Ierror to keep track of errors.                  |
C+-----------------------------------------------------------------------+
      SUBROUTINE DIRInit(f,fcn,c,length,actdeep,point,anchor,free,
     + dwrit,logfile,ArrayI,maxI,List2,w,x,l,u,fmin,minpos,thirds,
     + levels,maxfunc,maxdeep,n,maxor,fmax,Ifeasiblef,IInfeasible, 
     + Ierror,
     + iidata, iisize, ddata, idsize, cdata, icsize)
      IMPLICIT None
      Integer maxfunc,maxdeep,n,maxor
      Double Precision  f(maxfunc,2),c(maxfunc,maxor),fmin
      Double Precision  x(n),delta, thirds(0:maxdeep)
      Double Precision  levels(0:maxdeep)
      Integer length(maxfunc,maxor),actdeep,minpos,i,oops,j
      Integer point(maxfunc),anchor(-1:maxdeep),free
      Integer ArrayI(maxor),maxI,new,dwrit,logfile,List2(maxor,2)
      Double Precision  w(maxor)
      External fcn
      Double Precision help2,l(n),u(n)
      Double Precision costmin
C+-----------------------------------------------------------------------+
C| JG 01/22/01 Added variable to keep track of the maximum value found.  |
C+-----------------------------------------------------------------------+
      Double Precision fmax
C+-----------------------------------------------------------------------+
C| JG 01/22/01 Added variable Ifeasiblef to keep track if feasibel point |
C|             has been found.                                           |
C| JG 01/23/01 Added variable Ierror to keep track of errors.            |
C| JG 03/09/01 Added IInfeasible to keep track if an infeasible point has|
C|             been found.                                               |
C+-----------------------------------------------------------------------+
      Integer Ifeasiblef, Ierror, IInfeasible
C JG 09/15/00 Added variable JONES (see above)
      Integer JONES, help
      COMMON /directcontrol/ JONES
C+-----------------------------------------------------------------------+
C| Variables to pass user defined data to the function to be optimized.  |
C+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      Double Precision ddata(idsize)
      Character*40 cdata(icsize)
        
      fmin = 1.D20
      costmin = fmin
C JG 09/15/00 If Jones way of characterising rectangles is used, 
C             initialise thirds to reflect this.
      IF (JONES .eq. 0) THEN
        DO 5,j = 0,n-1
          w(j+1) = 0.5D0 * dsqrt(n - j + j/9.D0)
5       CONTINUE
        help2 = 1.D0
        DO 10,i = 1,maxdeep/n
          DO 8, j = 0, n-1
            levels((i-1)*n+j) = w(j+1) / help2
8         CONTINUE
          help2 = help2 * 3.D0
10      CONTINUE
      ELSE
C JG 09/15/00 Initialiase levels to contain 1/j
        help2 = 3.D0
        DO 11,i = 1,maxdeep
          levels(i) = 1.D0 / help2
          help2 = help2 * 3.D0
11      CONTINUE
        levels(0) = 1.D0
      ENDIF
      help2 = 3.D0
      DO 21,i = 1,maxdeep
        thirds(i) = 1.D0 / help2
        help2 = help2 * 3.D0
21    CONTINUE
      thirds(0) = 1.D0        
      DO 20,i=1,n
        c(1,i) = 0.5D0
        x(i) = 0.5D0
        length(1,i) = 0
20    CONTINUE
      CALL DIRinfcn(fcn,x,l,u,n,f(1,1),help,
     +                iidata, iisize, ddata, idsize, cdata, icsize)
      f(1,2) = help
      IInfeasible = help
      fmax = f(1,1)
C 09/25/00 Added this
C      if (f(1,1) .ge. 1.E+6) then
      if (f(1,2) .gt. 0.D0) then
        f(1,1) = 1.D6
        fmax = f(1,1)
        Ifeasiblef = 1
      else
        Ifeasiblef = 0
      end if

C JG 09/25/00 Remove IF
      fmin = f(1,1)
      costmin = f(1,1)
      minpos = 1
      actdeep = 2
      point(1) = 0
      free = 2
      delta = thirds(1)
      CALL DIRGet_I(length,1,ArrayI,maxI,n,maxfunc)
      new = free
      CALL DIRSamplepoints(c,ArrayI,delta,1,new,length,
     +           dwrit,logfile,f,free,maxI,point,fcn,x,l,
     +           fmin,minpos,u,n,
     +           maxfunc,maxdeep,oops)
C+-----------------------------------------------------------------------+
C| JG 01/23/01 Added error checking.                                     |
C+-----------------------------------------------------------------------+
      IF (oops .GT. 0) THEN
         IError = -4
         return
      END IF
C+-----------------------------------------------------------------------+
C| JG 01/22/01 Added variable to keep track of the maximum value found.  |
C|             Added variable to keep track if feasible point was found. |
C+-----------------------------------------------------------------------+
      CALL DIRSamplef(c,ArrayI,delta,1,new,length,
     +         dwrit,logfile,f,free,maxI,point,fcn,x,l,
     +         fmin,minpos,u,n,maxfunc,maxdeep,oops,fmax,Ifeasiblef,
     +         IInfeasible,
     +         iidata, iisize, ddata, idsize, cdata, icsize)
C+-----------------------------------------------------------------------+
C| JG 01/23/01 Added error checking.                                     |
C+-----------------------------------------------------------------------+
      IF (oops .GT. 0) THEN
         IError = -5
         return
      END IF
      CALL DIRDivide(new,0,length,point,
     +  ArrayI,1,List2,w,maxI,f,maxfunc,maxdeep,n)
      CALL DIRInsertList(new,anchor,point,f,maxI,length,
     +                    maxfunc,maxdeep,n,1)
      END

C+-----------------------------------------------------------------------+
C|    SUBROUTINE DIRDivide                                               |
C|    Subroutine to divide the hyper rectangles according to the rules.  |
C|    Changed 02-24-2000                                                 |
C|      Replaced if statement by min (line 367)                          |
C+-----------------------------------------------------------------------+
      SUBROUTINE DIRDivide(new,currentlength,length,point,
     +	   ArrayI,sample,List2,w,maxI,f,maxfunc,maxdeep,n)
      IMPLICIT None
      INTEGER start,new,maxfunc,maxdeep,n,sample
      INTEGER currentlength,length(maxfunc,n)
      INTEGER point(maxfunc)
      INTEGER List2(n,2),maxI,ArrayI(n)
      Double Precision f(maxfunc,2),w(n)
      INTEGER pos,i,j,k,pos2

        start = 0
        pos = new
        DO 10,i=1,maxI
          j = ArrayI(i)
          w(j) = f(pos,1)
          k = pos
          pos = point(pos)
          w(j) = min(f(pos,1),w(j))
          pos = point(pos)
          CALL DIRInsertList_2(start,j,k,list2,w,maxI,n)
10     CONTINUE
       IF (pos .GT. 0) THEN
           Write(*,*) "Error Divide"
           STOP
       END IF
       DO 20,j=1,maxI
         CALL DIRSearchmin(start,List2,pos,k,n)
         pos2 = start
         length(sample,k) = currentlength + 1
         DO 30,i=1,maxI-j+1
           length(pos,k) = currentlength + 1
           pos = point(pos)
           length(pos,k) = currentlength + 1
C JG 07/10/01 pos2 = 0 at the end of the 30-loop. Since we end
C             the loop now, we do not need to reassign pos and pos2.
           if (pos2 .gt. 0) then
             pos = List2(pos2,2)
             pos2 = List2(pos2,1)
           end if
30       CONTINUE
20     CONTINUE
      END

C+-----------------------------------------------------------------------+
C|    SUBROUTINE DIRSamplepoints                                         |
C|    Subroutine to sample the new points.                               |
C+-----------------------------------------------------------------------+
      SUBROUTINE DIRSamplepoints(c,ArrayI,delta,sample,start,length,
     +           dwrit,logfile,f,free,maxI,point,fcn,x,l,fmin,minpos,
     +           u,n,
     +           maxfunc,maxdeep,oops) 
      IMPLICIT None
      INTEGER n,maxfunc,maxdeep,oops
      INTEGER maxI,ArrayI(n),sample
      INTEGER length(maxfunc,n),free,point(maxfunc)
      Double Precision c(maxfunc,n),delta,x(n),l(n)
      Double Precision u(n),f(maxfunc,2)
      Double Precision fmin
      INTEGER pos,j,k,dwrit,logfile,minpos
      Integer start
      External fcn

      oops = 0
      pos = free
      start = free
      DO 10,k=1,maxI+maxI
        DO 20,j=1,n
          length(free,j) = length(sample,j)
          c(free,j) = c(sample,j)
20      CONTINUE
        pos = free
        free = point(free)
        IF (free .EQ. 0) THEN
           Write(*,1000) 
           Write(*,1001) 
CC           IF (dwrit .EQ. 2) THEN
CC             Write(logfile,1000)
CC             Write(logfile,1001)
CC           END IF
           oops = 1
           RETURN
1000  FORMAT("Error, no more free positions !")
1001  FORMAT("Increase maxfunc !")
        END IF
10    CONTINUE
      point(pos) = 0
      pos = start
      DO 30,j=1,maxI
         c(pos,ArrayI(j)) = c(sample,ArrayI(j)) + delta
         pos = point(pos)
         c(pos,ArrayI(j)) = c(sample,ArrayI(j)) - delta
         pos = point(pos)
30    CONTINUE
      IF (pos .GT. 0) THEN
          Write(*,2000)
CC          IF (dwrit .EQ. 2) THEN
CC             Write(logfile,2000)
CC           END IF
          STOP
2000      FORMAT("Error ! ") 
      END IF
      END


C+-----------------------------------------------------------------------+
C|    SUBROUTINE DIRGet_I                                                |
C+-----------------------------------------------------------------------+
      SUBROUTINE DIRGet_I(length,pos,ArrayI,maxi,n,maxfunc)
      IMPLICIT None
      Integer maxfunc,n,maxi,pos
      Integer length(maxfunc,n),ArrayI(n),i,help,j

      j = 1
      help = length(pos,1)
      DO 10,i = 2,n
        IF (length(pos,i) .LT. help) THEN
           help = length(pos,i)
        END IF
10    CONTINUE
      DO 20,i = 1,n
        IF (length(pos,i) .EQ. help) THEN
           ArrayI(j) = i
           j = j + 1
        END IF
20    CONTINUE
      maxi = j - 1
      END

C+-----------------------------------------------------------------------+
C|    SUBROUTINE DIRInitList                                             |
C|    Initialise the list.                                               |
C+-----------------------------------------------------------------------+
      SUBROUTINE DIRInitList(anchor,free,point,f,maxfunc,maxdeep)
      IMPLICIT None
      Integer maxdeep,maxfunc
      Double Precision f(maxfunc,2)
C   f -- values of functions.
      Integer anchor(-1:maxdeep)
C   anchor -- anchors of lists with deep i
      Integer point(maxfunc), free
C   point -- lists
C   free  -- first free position
      Integer i

      DO 10,i = -1,maxdeep
        anchor(i) = 0
10    CONTINUE
      DO 20,i = 1,maxfunc
        f(i,1) = 0.D0
        f(i,2) = 0
        point(i) = i + 1
C       point(i) = 0
20    CONTINUE
      point(maxfunc) = 0
      free = 1
      END


C+-----------------------------------------------------------------------+
C|    SUBROUTINE DIRInsert3                                              |
C+-----------------------------------------------------------------------+
         SUBROUTINE DIRInsert3(pos1,pos2,pos3,deep,anchor,point,free,
     +                   f,fmin,minpos,maxfunc,maxdeep)
      IMPLICIT None
      INTEGER maxfunc,maxdeep
      INTEGER deep,free,pos1,pos2,pos3
      INTEGER anchor(-1:maxdeep),point(maxfunc)
      Double Precision f(maxfunc,2),fmin
      INTEGER pos,minpos

      CALL DIRSort3(pos1,pos2,pos3,f,maxfunc)
      IF (anchor(deep) .EQ. 0) THEN
        anchor(deep) = pos1
        point(pos1) = pos2
        point(pos2) = pos3
        point(pos3) = 0
      ELSE
        pos = anchor(deep)
        IF (f(pos1,1) .LT. f(pos,1)) THEN
          anchor(deep) = pos1
          point(pos1) = pos
        ELSE
          CALL DIRInsert(pos,pos1,point,f,maxfunc)
        END IF
        CALL DIRInsert(pos,pos2,point,f,maxfunc)
        CALL DIRInsert(pos,pos3,point,f,maxfunc)	  	
      END IF
        IF ((f(pos1,1) .LT. fmin) .and. (f(pos1,2) .eq. 0)) THEN
           fmin = f(pos1,1)
           minpos = pos1
        END IF
      END

C+-----------------------------------------------------------------------+
C|    SUBROUTINE DIRInsert                                               |
C+-----------------------------------------------------------------------+
      SUBROUTINE DIRInsert(start,ins,point,f,maxfunc)
      IMPLICIT None
      INTEGER maxfunc,start,ins
      INTEGER point(maxfunc)
      Double Precision f(maxfunc,2)
      INTEGER i,help

C JG 09/17/00 Rewrote this routine.      
C      DO 10,i = 1,maxfunc
C        IF (f(ins,1) .LT. f(point(start),1)) THEN
C          help = point(start)
C          point(start) = ins
C          point(ins) = help
C          GOTO 20
C        END IF
C        IF (point(start) .EQ. 0) THEN
C           point(start) = ins
C           point(ins) = 0
C           GOTO 20
C        END IF
C        start = point(start)   
C10    CONTINUE
C20    END
      DO 10,i = 1,maxfunc
        IF (point(start) .EQ. 0) THEN
           point(start) = ins
           point(ins) = 0
           GOTO 20
        ELSE
          IF (f(ins,1) .LT. f(point(start),1)) THEN
            help = point(start)
            point(start) = ins
            point(ins) = help
            GOTO 20
          END IF
        END IF
        start = point(start)   
10    CONTINUE
20    END

C+-----------------------------------------------------------------------+
C|    SUBROUTINE DIRSort3                                                |
C+-----------------------------------------------------------------------+
      SUBROUTINE DIRSort3(pos1,pos2,pos3,f,maxfunc)
      IMPLICIT None
      INTEGER maxfunc
      Double Precision f(maxfunc,2)
      INTEGER pos1,pos2,pos3,help

      IF (f(pos1,1) .LT. f(pos2,1)) THEN
         IF (f(pos1,1) .LT. f(pos3,1)) THEN
           IF (f(pos3,1) .LT. f(pos2,1)) THEN
             help = pos2
             pos2 = pos3
             pos3 = help
           END IF
         ELSE
           help = pos1
           pos1 = pos3
           pos3 = pos2
           pos2 = help
         END IF 
      ELSE 
         IF (f(pos2,1) .LT. f(pos3,1)) THEN
           IF (f(pos3,1) .LT. f(pos1,1)) THEN
             help = pos1
             pos1 = pos2
             pos2 = pos3
             pos3 = help
           ELSE
             help = pos1
             pos1 = pos2
             pos2 = help
           END IF
         ELSE
           help = pos1
           pos1 = pos3
           pos3 = help
         END IF
      END IF
      END

C+-----------------------------------------------------------------------+
C|                                                                       |
C|                       SUBROUTINE DIRPREPRC                            |
C|                                                                       |
C| Subroutine DIRpreprc uses an afine mapping to map the hyper-box given |
C| by the constraints on the variable x onto the n-dimensional unit cube.|
C| This mapping is done using the following equation:                    |
C|                                                                       |
C|               x(i)=x(i)/(u(i)-l(i))-l(i)/(u(i)-l(i)).                 |
C|                                                                       |
C| DIRpreprc checks if the bounds l and u are well-defined. That is, if  |
C|                                                                       |
C|               l(i) < u(i) forevery i.                                 |
C|                                                                       |
C| On entry                                                              |
C|                                                                       |
C|          u -- A double-precision vector of length n. The vector       |
C|               containing the upper bounds for the n independent       |
C|               variables.                                              |
C|                                                                       |
C|          l -- A double-precision vector of length n. The vector       |
C|               containing the lower bounds for the n independent       |
C|               variables.                                              |
C|                                                                       |
C|          n -- An integer. The dimension of the problem.               |
C|                                                                       |
C| On return                                                             |
C|                                                                       |
C|        xs1 -- A double-precision vector of length n, used for scaling |
C|               and unscaling the vector x.                             |
C|                                                                       |
C|        xs2 -- A double-precision vector of length n, used for scaling |
C|               and unscaling the vector x.                             |
C|                                                                       |
C|                                                                       |
C|       oops -- An integer. If an upper bound is less than a lower      |
C|               bound or if the initial point is not in the             |
C|               hyper-box oops is set to 1 and iffco terminates.        |
C|                                                                       |
C+-----------------------------------------------------------------------+
      subroutine DIRpreprc(u,l,n,xs1,xs2,oops)

        integer n,i,oops

        Double Precision u(n),l(n),xs1(n),xs2(n)
        Double Precision help


        oops=0

        do 20 i=1,n

C+-----------------------------------------------------------------------+
C| Check if the hyper-box is well-defined.                               |
C+-----------------------------------------------------------------------+
            if(u(i).le.l(i))then
               oops=1
               return
            end if

20    continue

C+-----------------------------------------------------------------------+
C| Scale the initial iterate so that it is in the unit cube.             |
C+-----------------------------------------------------------------------+
      do 50 i=1,n
            help=(u(i)-l(i))
            xs2(i)=l(i)/help
            xs1(i)=help
 50     continue
      
        return
        end

C+-----------------------------------------------------------------------+
C|                                                                       |
C|                       SUBROUTINE DIRINFCN                             |
C|                                                                       |
C| Subroutine DIRinfcn unscales the variable x for use in the            |
C| user-supplied function evaluation subroutine fcn. After fcn returns   |
C| to DIRinfcn, DIRinfcn then rescales x for use by DIRECT.              |
C|                                                                       |
C| On entry                                                              |
C|                                                                       |
C|        fcn -- The argument containing the name of the user-supplied   |
C|               subroutine that returns values for the function to be   |
C|               minimized.                                              |
C|                                                                       |
C|          x -- A double-precision vector of length n. The point at     |
C|               which the derivative is to be evaluated.                |
C|                                                                       |
C|        xs1 -- A double-precision vector of length n. Used for         |
C|               scaling and unscaling the vector x by DIRinfcn.         |
C|                                                                       |
C|        xs2 -- A double-precision vector of length n. Used for         |
C|               scaling and unscaling the vector x by DIRinfcn.         |
C|                                                                       |
C|          n -- An integer. The dimension of the problem.               |
C|       kret -- An Integer. If kret =  1, the point is infeasible,      |
C|                              kret = -1, bad problem set up,           |
C|                              kret =  0, feasible.                     |
C|                                                                       |
C| On return                                                             |
C|                                                                       |
C|          f -- A double-precision scalar.                              |
C|                                                                       |
C| Subroutines and Functions                                             |
C|                                                                       |
C| The subroutine whose name is passed through the argument fcn.         |
C|                                                                       |
C+-----------------------------------------------------------------------+
      subroutine DIRinfcn(fcn,x,c1,c2,n,f,flag,
     +                  iidata, iisize, ddata, idsize, cdata, icsize)

      implicit none
      integer n,i, flag

      Double Precision x(n),c1(n),c2(n),f
      EXTERNAL fcn
C+-----------------------------------------------------------------------+
C| Variables to pass user defined data to the function to be optimized.  |
C+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      Double Precision ddata(idsize)
      Character*40 cdata(icsize)


C+-----------------------------------------------------------------------+
C| Unscale the variable x.                                               |
C+-----------------------------------------------------------------------+
      do 20 i=1,n
        x(i)=(x(i)+c2(i))*c1(i)
20    continue
C+-----------------------------------------------------------------------+
C| Call the function-evaluation subroutine fcn.                          |
C+-----------------------------------------------------------------------+
      f = 0.D0
      CALL fcn(n,x,f,flag,iidata, iisize, ddata, idsize, cdata, icsize)

C+-----------------------------------------------------------------------+
C| Rescale the variable x.                                               |
C+-----------------------------------------------------------------------+
      do 30 i=1,n
        x(i)=x(i)/c1(i)-c2(i)
 30   continue
      return

      end

C+-----------------------------------------------------------------------+
C| Subroutines to output the iteration history (with boxes) to a file.   |
C+-----------------------------------------------------------------------+
      SUBROUTINE DIRWriteHistBox(point,f,thirds,c,anchor,actdeep,file,
     +               l,u,file2,maxfunc,maxdeep,n,maxor,length)
      IMPLICIT None
      Integer maxfunc,maxdeep,actdeep,file,file2,n,maxor
      Integer length(maxfunc,maxor)
      Double Precision f(maxfunc,2),thirds(maxfunc)
      Double Precision c(maxfunc,n),l(n),u(n)
      Integer anchor(-1:maxdeep),i,point(maxfunc)
      Double Precision ufact(128),help
      Integer j

      do 40,i=1,n
        ufact(i) = (u(i) - l(i))
40    continue
      OPEN(12,FILE="matlab/DIRECT_histbox.dat",STATUS='UNKNOWN')
      DO 10,i = 1,maxfunc
         help = 0.D0
         Do 50, j=1,n
            help = max(help,c(i,j))
50       CONTINUE
         IF (help .GT. 0) THEN
          write(12,1000) f(i,1), f(i,2), (c(i,j)*ufact(j)+l(j),j=1,n), 
     +         (thirds(length(i,j)+1)*ufact(j),j=1,n)
        END IF
10    CONTINUE
      CLOSE(12)
1000  FORMAT(40E18.10)
      END



      SUBROUTINE DIRheader(logfile,version,x,n,eps,maxf,maxT,l,u,
     +               algmethod, maxfunc, maxdeep,    
     +               fglobal, fglper, Ierror,epsfix, iepschange,
     +               volper, sigmaper, 
     +               iidata, iisize, ddata, idsize, cdata, icsize)
      IMPLICIT None
      Integer logfile, version,n, maxf, maxT
      Integer algmethod, Ierror, i, maxfunc, maxdeep
      Double Precision x(n), l(n), u(n), eps    
      Double Precision fglobal, fglper, volper, sigmaper
      Integer iepschange
      Double Precision epsfix
C+-----------------------------------------------------------------------+
C| Variables to pass user defined data to the function to be optimized.  |
C+-----------------------------------------------------------------------+
      INTEGER iisize, idsize, icsize
      INTEGER iidata(iisize)
      Double Precision ddata(idsize)
      Character*40 cdata(icsize)
      
      Integer Imainver, Isubver, Isubsubver, Ihelp, numerrors

CC      Write(logfile,900)
      numerrors = 0
      IError = 0
      Imainver = INT(version/100)
      Ihelp = version - Imainver*100
      Isubver = INT(Ihelp /10)
      Ihelp = Ihelp - Isubver*10
      Isubsubver = Ihelp
C+-----------------------------------------------------------------------+
C| JG 01/13/01 Added check for epsilon. If epsilon is smaller 0, we use  |
C|             the update formula from Jones. We then set the flag       |
C|             iepschange to 1, and store the absolute value of eps in   |
C|             epsfix. epsilon is then changed after each iteration.     |
C+-----------------------------------------------------------------------+
      If (eps .lt. 0.D0) then
        iepschange = 1
        epsfix =  -eps
        eps = -eps
      else
        iepschange = 0
        epsfix = 1.D100
      endif
      
CC      write(logfile,100) Imainver, Isubver, Isubsubver
      write(*,100) Imainver, Isubver, Isubsubver
C+-----------------------------------------------------------------------+
C| JG 07/16/01 Removed printout of contents in cdata(1).                 |
C+-----------------------------------------------------------------------+
C      write(*,*) cdata(1)
      write(*,200) n
      write(*,201) eps
      if (iepschange .eq. 1) then
         write(*,206)
      else
         write(*,207)
      end if
      write(*,202) maxf
      write(*,203) maxT
      write(*,204) fglobal
      write(*,205) fglper
      write(*,208) volper
      write(*,209) sigmaper
      
C+-----------------------------------------------------------------------+
C| JG 07/16/01 Removed printout of contents in cdata(1).                 |
C+-----------------------------------------------------------------------+
C      write(logfile,*) cdata(1)
CC      write(logfile,200) n
CC      write(logfile,201) eps
CC      write(logfile,202) maxf
CC      write(logfile,203) maxT
CC      write(logfile,204) fglobal
CC      write(logfile,205) fglper
CC      write(logfile,208) volper
CC      write(logfile,209) sigmaper
      if (iepschange .eq. 1) then
CC         write(logfile,206)
      else
CC         write(logfile,207)
      end if
      if (algmethod .eq. 0) then
         write(*,*) 'Jones original DIRECT algorithm is used.'
CC         write(logfile,*) 'Jones original DIRECT algorithm is used.'
      else
         write(*,*) 'Our modification of the DIRECT algorithm is used.'
CC         write(logfile,*) 'Our modification of the DIRECT algorithm',
CC     +                    ' is used.'
      end if
      do 1010, i = 1,n
         IF (u(i) .le. l(i)) then
            Ierror = -1
            write(*,153) i,l(i), u(i)
CC            write(logfile,153) i,l(i), u(i)
            numerrors = numerrors + 1
         else
            write(*,152) i,l(i), u(i)
CC            write(logfile,152) i,l(i), u(i)
         end if
1010  continue
C+-----------------------------------------------------------------------+
C| If there are to many function evaluations or to many iteration, note  |
C| this and set the error flag accordingly. Note: If more than one error |
C| occurred, we give out an extra message.                               |
C+-----------------------------------------------------------------------+
      IF ((maxf+20) .GT. maxfunc) THEN
         Write(*,10001) maxf, maxfunc
CC         Write(logfile,10001) maxf, maxfunc
         numerrors = numerrors + 1
         IError = -2
      END IF
      if (IError .lt. 0) then
CC         write(logfile,120)
         write(*,120)
         if (numerrors .eq. 1) then
            write(*,105)
CC            write(logfile,105)
         else
            write(*,110) numerrors
CC            write(logfile,110) numerrors
         end if
      end if
CC      write(logfile,120)
      write(*,120)
      if (IError .ge. 0) then
CC         write(logfile,*) 'Iteration   # of f-eval.   fmin'
      end if

10001 FORMAT("WARNING : The maximum number of function evaluations (",
     +       I6,") is higher then the constant maxfunc (",I6,
     +       "). Increase maxfunc in subroutine DIRECT or ",
     + "decrease the maximum number of function evaluations.")
10005 FORMAT("WARNING : The maximum number of iterations (",I5,
     +        ") is higher then the constant maxdeep (",I5,
     +  "). Increase maxdeep or decrease the number of iterations. ")
      
100   FORMAT('DIRECT Version ',I1,'.',I1,'.',I1)
105   FORMAT('WARNING : There was one error in the input!')
110   FORMAT('WARNING : There were ',I2,' errors in the input!')
900   FORMAT('--------------------------------- Log file -------------'
     +       ,'-------------------')
120   FORMAT('--------------------------------------------------------'
     +       ,'-------------------')
152   format('Bounds on variable x',i2,'    : ', F12.5,
     + ' <= xi <= ', F12.5)
153   format('WARNING : Bounds on variable x',i2,'    : ', F12.5,
     + ' <= xi <= ', F12.5)
200   format(' Problem Dimension n                    : ',I6) 
201   format(' Eps value                              : ',E12.4)
202   format(' Maximum number of f-evaluations (maxf) : ',I6)
203   format(' Maximum number of iterations (MaxT)    : ',I6)
204   format(' Value of f_global                      : ',E12.4)
205   format(' Global percentage wanted               : ',E12.4)
206   format(' Epsilon is changed using the Jones formula.')
207   format(' Epsilon is constant.')
208   format(' Volume percentage wanted               : ',E12.4)
209   format(' Measure percentage wanted              : ',E12.4)

      END

      SUBROUTINE DIRsummary(logfile,x,l,u,n,fmin,
     +               fglobal, numfunc, Ierror)
      IMPLICIT None
      Integer logfile, n
      Integer Ierror, numfunc, i
      Double Precision x(n), l(n), u(n)
      Double Precision fglobal , fmin

CC      Write(logfile,900)
CC      Write(logfile,1000) fmin
CC      Write(logfile,1010) numfunc
CC      if (fglobal .gt. -1.D99) then
CC         write(logfile,1001) 100*(fmin-fglobal)/max(1.D0,abs(fglobal))
CC      end if
CC      Write(logfile,1002) 
CC      do 100, i = 1,n
CC         write(logfile,1003) i, x(i), x(i)-l(i), u(i) - x(i)
CC100   continue
CC      write(logfile,1200)
      
900   FORMAT('--------------------------------- Summary -------------'
     +       ,'-------------------')
1000  FORMAT('Final function value           : ',F12.7) 
1010  FORMAT('Number of function evaluations : ',I12) 
1001  FORMAT('Final function value is within ',F10.5, 
     +        ' percent of global optimum.')
1002  FORMAT('Index  Final solution   x(i) - l(i)   ',
     +       'u(i) - x(i) ') 
1003  FORMAT(i5,' ',F12.7,'    ',F12.7,'   ',F12.7)
1200  FORMAT('--------------------------------------------------------'
     +       ,'-------------------')

      END
      
      SUBROUTINE DIRMaxf_to_high1(maxf,maxfunc,dwrit,logfile)
      IMPLICIT None
      INTEGER maxf,maxfunc,dwrit,logfile

      IF (dwrit .gt. 0) THEN
          Write(*,10001) maxf
          Write(*,10002) maxfunc
          Write(*,10003)
          Write(*,10004)
      END IF
      IF (dwrit .EQ. 2) THEN
         Write(logfile,10001) maxf
         Write(logfile,10002) maxfunc
         Write(logfile,10003)
         Write(logfile,10004)
      END IF

10001 FORMAT("The maximum number of function evaluations (",
     +       I6,") is ")
10002 FORMAT("higher then the constant maxfunc (",I6,
     +       "). Increase ")
10003 FORMAT("maxfunc in the SUBROUTINE DIRECT or decrease ")
10004 FORMAT("the maximum number of function evaluations.")
      END

      SUBROUTINE DIRMaxT_to_high1(maxT,maxdeep,dwrit,logfile)
      IMPLICIT None
      INTEGER maxT,maxdeep,dwrit,logfile

      Write(*,10001) maxT
      Write(*,10002) maxdeep
      Write(*,10003)
      IF (dwrit .EQ. 2) THEN
         Write(logfile,10001) maxT
         Write(logfile,10002) maxdeep
         Write(logfile,10003)
      END IF

10001 FORMAT("The maximum number of iterations (",I5,
     +        ") is higher ")
10002 FORMAT("then the constant maxdeep (",I5,
     +       "). Increase maxdeep ")
10003 FORMAT("or decrease the number of iterations. ")
      END


      INTEGER Function Isinbox(x,a,b,n,Lmaxdim)
      IMPLICIT None

      Integer n, Lmaxdim
      Double Precision a(Lmaxdim),b(Lmaxdim),x(Lmaxdim)
      
      Integer i,outofbox
      
      outofbox = 1
      DO 1000, i = 1,n
        IF ((a(i) .gt. x(i)) .or. (b(i) .lt. x(i))) then
          outofbox = 0
          goto 1010
        end if
1000  continue
1010  Isinbox = outofbox
      end
