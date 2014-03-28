      subroutine gehrd(min_lwork,max_lwork,prefix,n,lo,hi)
      integer min_lwork,max_lwork,n,lo,hi
      character prefix
c
c     Returned maxwrk is acctually optimal lwork.
c
cf2py intent(out,out=minwrk) :: min_lwork
cf2py intent(out,out=maxwrk) :: max_lwork
cf2py intent(in) :: prefix
cf2py intent(in) :: n,lo,hi

      INTEGER NB
      EXTERNAL ILAENV
      INTRINSIC MIN

      NB = MIN( 64, ILAENV( 1, prefix // 'GEHRD', ' ', n, lo, hi, -1 ) )
      max_lwork = n * NB
      min_lwork = MIN(max_lwork,MAX(1,n))

      end

      subroutine orghr(min_lwork,max_lwork,prefix,n,lo,hi)
      integer min_lwork,max_lwork,n,lo,hi
      character prefix
c
c     Returned maxwrk is acctually optimal lwork.
c
cf2py intent(out,out=minwrk) :: min_lwork
cf2py intent(out,out=maxwrk) :: max_lwork
cf2py intent(in) :: prefix
cf2py intent(in) :: n,lo,hi

      INTEGER NB
      EXTERNAL ILAENV
      INTRINSIC MIN, MAX

      NB = MIN( 64, ILAENV( 1, prefix // 'ORGHR', ' ', n, lo, hi, -1 ) )
      max_lwork = MAX(1, (hi - lo) * NB)
      min_lwork = MAX(1, hi - lo)

      end

      subroutine gesdd(min_lwork,max_lwork,prefix,m,n,compute_uv)
      integer min_lwork,max_lwork,m,n,compute_uv
      character prefix

cf2py callstatement (*f2py_func)(&min_lwork,&max_lwork,prefix,&m,&n,&compute_uv)
cf2py callprotoargument int*,int*,char*,int*,int*,int*
cf2py intent(out,out=minwrk) :: min_lwork
cf2py intent(out,out=maxwrk) :: max_lwork
cf2py intent(in) :: prefix
cf2py intent(in) :: m,n,compute_uv

      INTEGER MINMN, MNTHR, MINWRK, MAXWRK, SMLSIZ, BDSPAC, BDSPAN
      INTEGER            ILAENV, WRKBL
      EXTERNAL           ILAENV
      INTRINSIC          INT, MAX, MIN

      MINMN = MIN( M, N )
      MNTHR = INT( MINMN*11.0D0 / 6.0D0 )
      MINWRK = 1
      MAXWRK = 1
      SMLSIZ = ILAENV( 9, prefix // 'GESDD', ' ', 0, 0, 0, 0 )
      IF( M.GE.N ) THEN
*
*           Compute space needed for DBDSDC
*
         BDSPAC = 3*N*N + 7*N
         BDSPAN = MAX( 12*N+4, 8*N+2+SMLSIZ*( SMLSIZ+8 ) )
         IF( M.GE.MNTHR ) THEN
            IF (compute_uv.eq.0) THEN
*
*     Path 1 (M much larger than N, JOBZ='N')
*
               MAXWRK = N + N*ILAENV( 1, prefix // 'GEQRF', ' ', 
     $              M, N, -1,
     $              -1 )
               MAXWRK = MAX( MAXWRK, 3*N+2*N*
     $              ILAENV( 1, prefix // 'GEBRD', ' ',
     $              N, N, -1, -1 ) )
               MAXWRK = MAX( MAXWRK, BDSPAC )
               MINWRK = BDSPAC
            ELSE
*     
*     Path 4 (M much larger than N, JOBZ='A')
*
               WRKBL = N + N*ILAENV( 1, prefix // 'GEQRF', ' ',
     $              M, N, -1, -1 )
               WRKBL = MAX( WRKBL, N+M*ILAENV( 1, prefix // 'ORGQR',
     $              ' ', M,
     $              M, N, -1 ) )
               WRKBL = MAX( WRKBL, 3*N+2*N*
     $              ILAENV( 1, prefix // 'GEBRD', ' ',
     $              N, N, -1, -1 ) )
               WRKBL = MAX( WRKBL, 3*N+N*
     $              ILAENV( 1, prefix // 'ORMBR', 'QLN', 
     $              N, N, N, -1 ) )
               WRKBL = MAX( WRKBL, 3*N+N*
     $              ILAENV( 1, prefix // 'ORMBR', 'PRT', 
     $              N, N, N, -1 ) )
               WRKBL = MAX( WRKBL, BDSPAC+2*N )
               MAXWRK = N*N + WRKBL
               MINWRK = BDSPAC + N*N + M + N
            ENDIF
         ELSE
*
*     Path 5 (M at least N, but not much larger)
*     
            WRKBL = 3*N + ( M+N )*ILAENV( 1, prefix // 'GEBRD', ' ',
     $           M, N, -1, -1)
            IF (compute_uv.eq.0) THEN
               MAXWRK = MAX(WRKBL,BDSPAC + 3*N)
               MINWRK = 3*N + MAX(M,BDSPAC)
            ELSE
               MAXWRK = MAX( MAXWRK, 3*N+M*
     $              ILAENV( 1, prefix // 'ORMBR', 'QLN', 
     $              M, M, N, -1 ) )
               MAXWRK = MAX( MAXWRK, 3*N+N*
     $              ILAENV( 1, prefix // 'ORMBR', 'PRT', 
     $              N, N, N, -1 ) )
               MAXWRK = MAX( MAXWRK, BDSPAC+2*N+M )
               MINWRK = BDSPAC + 2*N + M
            ENDIF
         ENDIF
      ELSE
*
*     Compute space needed for DBDSDC
*     
         BDSPAC = 3*M*M + 7*M
         BDSPAN = MAX( 12*M+4, 8*M+2+SMLSIZ*( SMLSIZ+8 ) )
         IF( N.GE.MNTHR ) THEN
            IF( compute_uv.eq.0 ) THEN
*     
*     Path 1t (N much larger than M, JOBZ='N')
*     
               MAXWRK = M + M*ILAENV( 1, prefix // 'GELQF', ' ',
     $              M, N, -1,
     $              -1 )
               MAXWRK = MAX( MAXWRK, 3*M+2*M*
     $              ILAENV( 1, prefix // 'GEBRD', ' ',
     $              M, M, -1, -1 ) )
               MAXWRK = MAX( MAXWRK, BDSPAC )
               MINWRK = BDSPAC
            ELSE
*
*     Path 4t (N much larger than M, JOBZ='A')
*     
               WRKBL = M + M*ILAENV( 1, prefix // 'GELQF', ' ',
     $              M, N, -1, -1 )
               WRKBL = MAX( WRKBL, M+N*ILAENV( 1, prefix // 'ORGLQ',
     $              ' ', N,
     $              N, M, -1 ) )
               WRKBL = MAX( WRKBL, 3*M+2*M*
     $              ILAENV( 1, prefix // 'GEBRD', ' ', 
     $              M, M, -1, -1 ) )
               WRKBL = MAX( WRKBL, 3*M+M*
     $              ILAENV( 1, prefix // 'ORMBR', 'QLN', 
     $              M, M, M, -1 ) )
               WRKBL = MAX( WRKBL, 3*M+M*
     $              ILAENV( 1, prefix // 'ORMBR', 'PRT', 
     $              M, M, M, -1 ) )
               WRKBL = MAX( WRKBL, BDSPAC+2*M )
               MAXWRK = WRKBL + M*M
               MINWRK = BDSPAC + M*M + M + N
            ENDIF
         ELSE
            WRKBL = 3*M + ( M+N )*ILAENV( 1, prefix // 'GEBRD', ' ',
     $           M, N, -1,
     $           -1 )
            IF (compute_uv.eq.0) THEN
               MAXWRK = MAX(WRKBL,BDSPAC + 3*M)
               MINWRK = 3*M + MAX(N,BDSPAC)               
            ELSE
               MAXWRK = MAX( MAXWRK, 3*M+M*
     $              ILAENV( 1, prefix // 'ORMBR', 'QLN', 
     $              M, M, N, -1 ) )
               MAXWRK = MAX( MAXWRK, 3*M+N*
     $              ILAENV( 1, prefix // 'ORMBR', 'PRT', 
     $              N, N, M, -1 ) )
               MAXWRK = MAX( MAXWRK, BDSPAC+2*M )
               MINWRK = BDSPAC + 2*M + N
            ENDIF
         ENDIF
      ENDIF
      min_lwork = MINWRK
      max_lwork = MAX(MINWRK,MAXWRK)
      end

      subroutine gelss(min_lwork,max_lwork,prefix,m,n,nrhs)

      integer min_lwork,max_lwork,m,n,nrhs
      character prefix      

cf2py callstatement (*f2py_func)(&min_lwork,&max_lwork,prefix,&m,&n,&nrhs)
cf2py callprotoargument int*,int*,char*,int*,int*,int*
cf2py intent(out,out=minwrk) :: min_lwork
cf2py intent(out,out=maxwrk) :: max_lwork
cf2py intent(in) :: prefix
cf2py intent(in) :: m,n,nrhs

      INTEGER MAXWRK, MINMN, MINWRK, MM, MNTHR
      INTEGER ILAENV, BDSPAC, MAXMN
      EXTERNAL ILAENV
      INTRINSIC          MAX, MIN

      MINMN = MIN( M, N )
      MAXMN = MAX( M, N )
      MNTHR = ILAENV( 6, prefix // 'GELSS', ' ', M, N, NRHS, -1 )
      MINWRK = 1
      MAXWRK = 0
      MM = M
      IF( M.GE.N .AND. M.GE.MNTHR ) THEN
*     
*     Path 1a - overdetermined, with many more rows than columns
*     
         MM = N
         MAXWRK = MAX( MAXWRK, N+N*ILAENV( 1, prefix //  'GEQRF', ' ',
     $        M, N, -1, -1 ) )
         MAXWRK = MAX( MAXWRK, N+NRHS*
     $        ILAENV( 1, prefix //  'ORMQR', 'LT', M, NRHS, N, -1 ) )
      END IF
      IF( M.GE.N ) THEN
*     
*     Path 1 - overdetermined or exactly determined
*     
*     Compute workspace neede for BDSQR
*     
         BDSPAC = MAX( 1, 5*N )
         MAXWRK = MAX( MAXWRK, 3*N+( MM+N )*
     $        ILAENV( 1, prefix // 'GEBRD', ' ', MM, N, -1, -1 ) )
         MAXWRK = MAX( MAXWRK, 3*N+NRHS*
     $        ILAENV( 1, prefix // 'ORMBR', 'QLT', MM, NRHS, N, -1 ) )
         MAXWRK = MAX( MAXWRK, 3*N+( N-1 )*
     $        ILAENV( 1, prefix // 'ORGBR', 'P', N, N, N, -1 ) )
         MAXWRK = MAX( MAXWRK, BDSPAC )
         MAXWRK = MAX( MAXWRK, N*NRHS )
         MINWRK = MAX( 3*N+MM, 3*N+NRHS, BDSPAC )
         MAXWRK = MAX( MINWRK, MAXWRK )
      END IF
         
      IF( N.GT.M ) THEN
*     
*     Compute workspace neede for DBDSQR
*
         BDSPAC = MAX( 1, 5*M )
         MINWRK = MAX( 3*M+NRHS, 3*M+N, BDSPAC )
         IF( N.GE.MNTHR ) THEN
*     
*     Path 2a - underdetermined, with many more columns
*     than rows
*     
            MAXWRK = M + M*ILAENV( 1, prefix // 'GELQF', ' ', 
     $           M, N, -1, -1 )
            MAXWRK = MAX( MAXWRK, M*M+4*M+2*M*
     $           ILAENV( 1, prefix //  'GEBRD', ' ', M, M, -1, -1 ) )
            MAXWRK = MAX( MAXWRK, M*M+4*M+NRHS*
     $           ILAENV( 1, prefix // 'ORMBR', 'QLT', M, NRHS, M, -1 ))
            MAXWRK = MAX( MAXWRK, M*M+4*M+( M-1 )*
     $           ILAENV( 1, prefix // 'ORGBR', 'P', M, M, M, -1 ) )
            MAXWRK = MAX( MAXWRK, M*M+M+BDSPAC )
            IF( NRHS.GT.1 ) THEN
               MAXWRK = MAX( MAXWRK, M*M+M+M*NRHS )
            ELSE
               MAXWRK = MAX( MAXWRK, M*M+2*M )
            END IF
            MAXWRK = MAX( MAXWRK, M+NRHS*
     $           ILAENV( 1, prefix // 'ORMLQ', 'LT', N, NRHS, M, -1 ) )
            
         ELSE
*     
*     Path 2 - underdetermined
*     
            MAXWRK = 3*M + ( N+M )*ILAENV( 1, prefix // 'GEBRD', ' ',
     $           M, N, -1, -1 )
            MAXWRK = MAX( MAXWRK, 3*M+NRHS*
     $           ILAENV( 1, prefix // 'ORMBR', 'QLT', M, NRHS, M, -1 ) )
            MAXWRK = MAX( MAXWRK, 3*M+M*
     $           ILAENV( 1, prefix // 'ORGBR', 'P', M, N, M, -1 ) )
            MAXWRK = MAX( MAXWRK, BDSPAC )
            MAXWRK = MAX( MAXWRK, N*NRHS )
         END IF
      END IF
      MAXWRK = MAX( MINWRK, MAXWRK )
      MINWRK = MAX( MINWRK, 1 )

      min_lwork = MINWRK
      max_lwork = MAXWRK
      end

      subroutine getri(min_lwork,max_lwork,prefix,n)
      integer min_lwork,max_lwork,n
      character prefix
cf2py callstatement (*f2py_func)(&min_lwork,&max_lwork,prefix,&n)
cf2py callprotoargument int*,int*,char*,int*
cf2py intent(out,out=minwrk) :: min_lwork
cf2py intent(out,out=maxwrk) :: max_lwork
cf2py intent(in) :: prefix
cf2py intent(in) :: n
      INTEGER ILAENV, NB
      EXTERNAL ILAENV
      NB = ILAENV( 1, prefix // 'GETRI', ' ', N, -1, -1, -1 )
      min_lwork = N
      max_lwork = N*NB
      end

      subroutine geev(min_lwork,max_lwork,prefix,n,
     $     compute_vl,compute_vr)

      integer min_lwork,max_lwork,n,compute_vl,compute_vr
      character prefix
cf2py callstatement (*f2py_func)(&min_lwork,&max_lwork,prefix,&n,&compute_vl,&compute_vr)
cf2py callprotoargument int*,int*,char*,int*,int*,int*
cf2py intent(out,out=minwrk) :: min_lwork
cf2py intent(out,out=maxwrk) :: max_lwork
cf2py integer optional,intent(in) :: compute_vl = 1,compute_vr = 1
cf2py intent(in) :: prefix
cf2py intent(in) :: n

      LOGICAL WANTVL, WANTVR
      INTEGER ILAENV, MINWRK, MAXWRK, MAXB, HSWORK, K
      EXTERNAL ILAENV
      INTRINSIC          MAX, MIN

      WANTVL = compute_vl.eq.1
      WANTVR = compute_vr.eq.1

      MINWRK = 1
      MAXWRK = 2*N + N*ILAENV( 1, prefix // 'GEHRD', ' ', N, 1, N, 0 )
      IF( ( .NOT.WANTVL ) .AND. ( .NOT.WANTVR ) ) THEN
         MINWRK = MAX( 1, 3*N )
         MAXB = MAX( ILAENV( 8, prefix // 'HSEQR', 'EN', N, 1, N, -1 )
     $        , 2 )
         K = MIN( MAXB, N, MAX( 2, ILAENV( 4, prefix // 'HSEQR', 'EN', N
     $        , 1,  N, -1 ) ) )
         HSWORK = MAX( K*( K+2 ), 2*N )
         MAXWRK = MAX( MAXWRK, N+1, N+HSWORK )
      ELSE
         MINWRK = MAX( 1, 4*N )
         MAXWRK = MAX( MAXWRK, 2*N+( N-1 )*
     $        ILAENV( 1, prefix // 'ORGHR', ' ', N, 1, N, -1 ) )
         MAXB = MAX( ILAENV( 8, prefix // 'HSEQR', 'SV', N, 1, N, -1 ), 
     $        2  )
         K = MIN( MAXB, N, MAX( 2, ILAENV( 4, prefix // 'HSEQR', 'SV', N
     $        , 1, N, -1 ) ) )
         HSWORK = MAX( K*( K+2 ), 2*N )
         MAXWRK = MAX( MAXWRK, N+1, N+HSWORK )
         MAXWRK = MAX( MAXWRK, 4*N )
      END IF
      min_lwork = MINWRK
      max_lwork = MAXWRK
      end

      subroutine heev(min_lwork,max_lwork,prefix,n,lower)

      integer min_lwork,max_lwork,n,lower
      character prefix
cf2py callstatement (*f2py_func)(&min_lwork,&max_lwork,prefix,&n,&lower)
cf2py callprotoargument int*,int*,char*,int*,int*
cf2py intent(out,out=minwrk) :: min_lwork
cf2py intent(out,out=maxwrk) :: max_lwork
cf2py integer optional,intent(in) :: lower = 0
cf2py intent(in) :: prefix
cf2py intent(in) :: n

      CHARACTER UPLO
      INTEGER ILAENV, NB
      EXTERNAL ILAENV
      INTRINSIC MAX

      UPLO = 'L'
      if (lower.eq.0) then
        UPLO = 'U'
      endif

      NB = ILAENV( 1, prefix // 'HETRD', UPLO, N, -1, -1, -1 )

      min_lwork = MAX(1,2*N-1)
      max_lwork = MAX( 1, ( NB+1 )*N )

      end

      subroutine syev(min_lwork,max_lwork,prefix,n,lower)

      integer min_lwork,max_lwork,n,lower
      character prefix
cf2py callstatement (*f2py_func)(&min_lwork,&max_lwork,prefix,&n,&lower)
cf2py callprotoargument int*,int*,char*,int*,int*
cf2py intent(out,out=minwrk) :: min_lwork
cf2py intent(out,out=maxwrk) :: max_lwork
cf2py integer optional,intent(in) :: lower = 0
cf2py intent(in) :: prefix
cf2py intent(in) :: n

      CHARACTER UPLO
      INTEGER ILAENV, NB
      EXTERNAL ILAENV
      INTRINSIC MAX

      UPLO = 'L'
      if (lower.eq.0) then
        UPLO = 'U'
      end if

      NB = ILAENV( 1, prefix // 'SYTRD', UPLO, N, -1, -1, -1 )

      min_lwork = MAX(1,3*N-1)
      max_lwork = MAX( 1, ( NB+2 )*N )

      end

      subroutine gees(min_lwork,max_lwork,prefix,n,compute_v)

      integer min_lwork,max_lwork,n,compute_v
      character prefix

cf2py callstatement (*f2py_func)(&min_lwork,&max_lwork,prefix,&n,&compute_v)
cf2py callprotoargument int*,int*,char*,int*,int*
cf2py intent(out,out=minwrk) :: min_lwork
cf2py intent(out,out=maxwrk) :: max_lwork
cf2py integer optional,intent(in) :: compute_v = 1
cf2py intent(in) :: prefix
cf2py intent(in) :: n

      INTEGER            HSWORK, MAXWRK, MINWRK, MAXB, K
      INTEGER            ILAENV
      EXTERNAL           ILAENV
      INTRINSIC          MAX, MIN

      MAXWRK = N + N*ILAENV( 1, prefix // 'GEHRD', ' ', N, 1, N, 0 )
      MINWRK = MAX( 1, 2*N )
      IF( compute_v.eq.0 ) THEN
         MAXB = MAX( ILAENV( 8, prefix // 'HSEQR',
     $        'SN', N, 1, N, -1 ), 2 )
         K = MIN( MAXB, N, MAX( 2, ILAENV( 4, prefix // 'HSEQR', 
     $    'SN', N, 1, N, -1 ) ) )
         HSWORK = MAX( K*( K+2 ), 2*N )
         MAXWRK = MAX( MAXWRK, HSWORK, 1 )
      ELSE
         MAXWRK = MAX( MAXWRK, N+( N-1 )*
     $        ILAENV( 1, prefix // 'UNGHR', ' ', N, 1, N, -1 ) )
         MAXB = MAX( ILAENV( 8, prefix // 'HSEQR',
     $        'EN', N, 1, N, -1 ), 2 )
         K = MIN( MAXB, N, MAX( 2, ILAENV( 4, prefix // 'HSEQR',
     $        'EN', N, 1, N, -1 ) ) )
         HSWORK = MAX( K*( K+2 ), 2*N )
         MAXWRK = MAX( MAXWRK, HSWORK, 1 )
      END IF

      min_lwork = MINWRK
      max_lwork = MAXWRK

      end

      subroutine geqrf(min_lwork,max_lwork,prefix,m,n)

      integer min_lwork,max_lwork,m,n
      character prefix      

cf2py callstatement (*f2py_func)(&min_lwork,&max_lwork,prefix,&m,&n)
cf2py callprotoargument int*,int*,char*,int*,int*
cf2py intent(out,out=minwrk) :: min_lwork
cf2py intent(out,out=maxwrk) :: max_lwork
cf2py intent(in) :: prefix
cf2py intent(in) :: m,n

      INTEGER NB
      INTEGER ILAENV
      EXTERNAL ILAENV
      INTRINSIC MAX

      NB = ILAENV( 1, prefix // 'GEQRF', ' ', M, N, -1, -1 )

      min_lwork = MAX(1,N)
      max_lwork = MAX(1,N*NB)
      end

      subroutine gqr(min_lwork,max_lwork,prefix,m,n)

      integer min_lwork,max_lwork,m,n
      character prefix

cf2py callstatement (*f2py_func)(&min_lwork,&max_lwork,prefix,&m,&n)
cf2py callprotoargument int*,int*,char*,int*,int*
cf2py intent(out,out=minwrk) :: min_lwork
cf2py intent(out,out=maxwrk) :: max_lwork
cf2py intent(in) :: prefix
cf2py intent(in) :: m,n

      INTEGER NB
      INTEGER ILAENV
      EXTERNAL ILAENV
      INTRINSIC MAX

      if ((prefix.eq.'d').or.(prefix.eq.'s')
     $     .or.(prefix.eq.'D').or.(prefix.eq.'S')) then
         NB = ILAENV( 1, prefix // 'ORGQR', ' ', M, N, -1, -1 )
      else
         NB = ILAENV( 1, prefix // 'UNGQR', ' ', M, N, -1, -1 )
      endif
      min_lwork = MAX(1,N)
      max_lwork = MAX(1,N*NB)
      end
