*--------------------------------------------------------------------
*\Documentation
*
*\Name: ICOPY
*
*\Description:
*     ICOPY copies an integer vector lx to an integer vector ly.
*
*\Usage:
*     call icopy ( n, lx, inc, ly, incy )
*
*\Arguments:
*    n        integer (input)
*             On entry, n is the number of elements of lx to be
c             copied to ly.
*
*    lx       integer array (input)
*             On entry, lx is the integer vector to be copied.
*
*    incx     integer (input)
*             On entry, incx is the increment between elements of lx.
*
*    ly       integer array (input)
*             On exit, ly is the integer vector that contains the
*             copy of lx.
*
*    incy     integer (input)
*             On entry, incy is the increment between elements of ly.
*
*\Enddoc
*
*--------------------------------------------------------------------
*
      subroutine icopy( n, lx, incx, ly, incy )
*
*     ----------------------------
*     Specifications for arguments
*     ----------------------------
      integer    incx, incy, n
      integer    lx( 1 ), ly( 1 )
*
*     ----------------------------------
*     Specifications for local variables
*     ----------------------------------
      integer           i, ix, iy
*
*     --------------------------
*     First executable statement
*     --------------------------
      if( n.le.0 )
     $   return
      if( incx.eq.1 .and. incy.eq.1 )
     $   go to 20
c
c.....code for unequal increments or equal increments
c     not equal to 1
      ix = 1
      iy = 1
      if( incx.lt.0 )
     $   ix = ( -n+1 )*incx + 1
      if( incy.lt.0 )
     $   iy = ( -n+1 )*incy + 1
      do 10 i = 1, n
         ly( iy ) = lx( ix )
         ix = ix + incx
         iy = iy + incy
   10 continue
      return
c
c.....code for both increments equal to 1
c
   20 continue
      do 30 i = 1, n
         ly( i ) = lx( i )
   30 continue
      return
      end
