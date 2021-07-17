ccccccccccccccccccccc Read ncol, nnz from HB file cccccccccccccccccccccc

      subroutine readhb_hdr(filename,nc,nz)
      implicit none

      character*256 filename
      integer nc, nz
      integer nrow, ncol, nnzero

      integer lunit
      parameter(lunit=10)
      CHARACTER      TITLE*72 , KEY*8    , MXTYPE*3 , RHSTYP*3,
     1               PTRFMT*16, INDFMT*16, VALFMT*20, RHSFMT*20

      INTEGER        TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD,
     1               NELTVL, I,  NRHS  , NRHSIX

      open(lunit,file=filename,status='old')

C     ------------------------
C     ... READ IN HEADER BLOCK
C     ------------------------
      read ( lunit, 1000 ) title , key   ,
     1                     totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     2                     mxtype, nrow  , ncol  , nnzero, neltvl,
     3     ptrfmt, indfmt, valfmt, rhsfmt
 1000 format ( a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20 )
      close(lunit)
      nc = ncol
      nz = nnzero
      end


ccccccccccccccccccccc Harwell-Boeing format cccccccccccccccccccccccccccc
c Note: this is also known as compressed column storage (CCS) format.

      subroutine readhb(filename,nc,nz,values,colptr,rowind)
      implicit none
c
c     Read matrix in Harwell-Boeing format
c

      integer mmax,nmax,kmax,nnzmax
      integer densemmax, densenmax
      parameter(densemmax=1000,densenmax=1000)
      parameter(mmax=5000,nmax=5000,kmax=1000)
      parameter(nnzmax=densemmax*densenmax)
      integer nrow, ncol, nnzero

      character*256 filename
      integer nz, nc
      integer rowind(nz), colptr(nc+1)
      complex values(nz)

      integer lunit
      parameter(lunit=10)
      CHARACTER      TITLE*72 , KEY*8    , MXTYPE*3 , RHSTYP*3,
     1               PTRFMT*16, INDFMT*16, VALFMT*20, RHSFMT*20

      INTEGER        TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD,
     1               NELTVL, I,  NRHS  , NRHSIX
      logical lsame
      external lsame


      open(lunit,file=filename,status='old')

C    ------------------------
C     ... READ IN HEADER BLOCK
C     ------------------------
      read ( lunit, 1000 ) title , key   ,
     1                     totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     2                     mxtype, nrow  , ncol  , nnzero, neltvl,
     3                     ptrfmt, indfmt, valfmt, rhsfmt
      if  ( rhscrd .gt. 0 )
     1    read ( lunit, 1001 ) rhstyp, nrhs, nrhsix
 1000 format ( a72, a8 / 5i14 / a3, 11x, 4i14 / 2a16, 2a20 )
 1001 format ( a3, 11x, 2i14 )


c      write (*,1000) title , key   ,
c     1                     totcrd, ptrcrd, indcrd, valcrd, rhscrd,
c     2                     mxtype, nrow  , ncol  , nnzero, neltvl,
c     3                     ptrfmt, indfmt, valfmt, rhsfmt
c     -------------------------
c     ... READ MATRIX STRUCTURE
c     -------------------------

      if (ncol.gt.nmax) stop 'ERROR in readHB: ncol > nmax'
      if (nrow.gt.mmax) stop 'ERROR in readHB: nrow > mmax'
      if (nnzero.gt.nnzmax) stop 'ERROR in readHB: nnzero > nnzmax'

      read ( lunit, ptrfmt ) ( colptr (i), i = 1, ncol+1 )
      read ( lunit, indfmt ) ( rowind (i), i = 1, nnzero )

      if  ( valcrd .gt. 0 )  then

c         ----------------------
c         ... read matrix values
c         ----------------------

          read ( lunit, valfmt ) ( values (i), i = 1, nnzero )

      endif

      close(lunit)

      end
