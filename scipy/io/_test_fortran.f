      subroutine read_unformatted_double(m, n, k, a, filename)
      implicit none
      integer :: m, n, k
      double precision :: a(m,n,k)
      character*4096 :: filename
      open(10, file=trim(filename), form='unformatted')
      read(10) a
      close(10)
      end subroutine

      subroutine read_unformatted_int(m, n, k, a, filename)
      implicit none
      integer :: m, n, k
      integer :: a(m,n,k)
      character*4096 :: filename
      open(10, file=trim(filename), form='unformatted')
      read(10) a
      close(10)
      end subroutine

      subroutine read_unformatted_mixed(m, n, k, a, b, filename)
      implicit none
      integer :: m, n, k
      double precision :: a(m,n)
      integer :: b(k)
      character*4096 :: filename
      open(10, file=trim(filename), form='unformatted')
      read(10) a, b
      close(10)
      end subroutine
