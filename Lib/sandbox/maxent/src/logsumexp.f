c File: logsumexp.f
c -----------------
c Functions for computing logs of sums of exponentials without numerical
c overflow.
c 
c Copyright: Ed Schofield, 2003-2005
c License: BSD-style (see LICENSE.txt in root SciPy directory)
c -----------------

        double precision function logsumexp(a,n)
c       Compute log(e^a_1 + e^a_2 + ... + e^a_n)
c
        double precision :: a(n)
        integer :: n
cf2py   intent(hide),depend(a) :: n = len(a)
cf2py   double precision intent(in),dimension(n) :: a    

        double precision :: a_i, b_i
        double precision neginfinity
        integer :: i
    
        neginfinity = log(0.0)

c       Loop until we have a value > -inf
        j = 1
        do while (a(j).eq.neginfinity)
           if (j.lt.n) then
              j = j + 1
           else
c       If all inputs are -inf, return -inf as the result
              logsumexp = neginfinity
              return
           end if
        end do

        b_i = a(j)
        do i = j+1,n
            a_i = a(i)
            if (a_i.eq.neginfinity) then
               continue
            end if
            b_i = max(b_i, a_i) + log(1+exp(-abs(b_i-a_i)))
        enddo
    
        logsumexp = b_i
        end function logsumexp
        

        double complex function logsumexpcomplex(a,n)
c       Compute log(e^a_1 + e^a_2 + ... + e^a_n) when the a_i are
c       possibly complex, and return the complex logarithm.

c       A version of logsumexp for complex inputs, such as the output of
c       robustarraylog().  We expect
c  
c       cmath.exp(logsumexpcomplex(robustarraylog(values))) ~= sum(values)
c
c       except for a small rounding error in both real and imag components.
c       To recover just the real component, use A.real, where A is the
c       complex return value.
c
        double complex :: a(n)
        integer :: n
cf2py   intent(hide),depend(a) :: n = len(a)
cf2py   double complex intent(in),dimension(n) :: a    
          
        double complex :: a_i, b_i
        double complex neginfinity
        integer :: i, j
        
        neginfinity = log(0.0)

c       Loop until we have a value > -inf
        j = 1
        do while (a(j).eq.neginfinity)
           if (j.lt.n) then
              j = j + 1
           else
c       If all inputs are -inf, return -inf as the result
              logsumexpcomplex = neginfinity
              return
           end if
        end do
        
        b_i = a(j)
        do i = j+1,n
           a_i = a(i)
           if (a_i.eq.neginfinity) then
              continue
           end if
           if (DBLE(b_i) > DBLE(a_i)) then
               b_i = b_i + log(1.+exp(a_i-b_i))
           else
               b_i = a_i + log(1.+exp(b_i-a_i))
           end if
        enddo
        logsumexpcomplex = b_i
        end function logsumexpcomplex
        
          

        subroutine robustarraylog(input, m, output)
c       An array version of a 'robust logarithm'.  Operates on a real
c       array {x}, for each element returning:
c           log(x) if x > 0,
c           the complex log of x if x < 0, or
c           float('-inf') if x == 0.
c
cf2py   intent(hide),depend(input)  :: m = len(input)
cf2py   double precision intent(in),dimension(m) :: input
cf2py   double complex intent(out),dimension(m) :: output    
        double precision :: input(m)
        double complex :: output(m)
        integer m, i
        do i=1,m
           if (input(i) .gt. 0.0) then
              output(i) = log(input(i))
           else if (input(i) .lt. 0.0) then
              output(i) = log(DCMPLX(input(i),0.0))
           else if (input(i) .ne. input(i)) then
c           This is only true if input(i) is nan.  In this case raise an
c           exception
cf2py   callstatement {f2py_success = 0;}
              output(i) = input(i)
           else
              output(i) = log(0.0)
           end if
        enddo
        end subroutine
