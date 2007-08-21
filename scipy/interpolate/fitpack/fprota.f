      subroutine fprota(cos,sin,a,b)
c  subroutine fprota applies a givens rotation to a and b.
c  ..
c  ..scalar arguments..
      real*8 cos,sin,a,b
c ..local scalars..
      real*8 stor1,stor2
c  ..
      stor1 = a
      stor2 = b
      b = cos*stor2+sin*stor1
      a = cos*stor1-sin*stor2
      return
      end
