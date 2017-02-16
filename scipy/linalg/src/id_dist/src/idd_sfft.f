c       this file contains the following user-callable routines:
c
c
c       routine idd_sffti initializes routine idd_sfft.
c
c       routine idd_sfft rapidly computes a subset of the entries
c       of the DFT of a vector, composed with permutation matrices
c       both on input and on output.
c
c       routine idd_ldiv finds the greatest integer less than or equal
c       to a specified integer, that is divisible by another (larger)
c       specified integer.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine idd_ldiv(l,n,m)
c
c       finds the greatest integer less than or equal to l
c       that divides n.
c
c       input:
c       l -- integer at least as great as m
c       n -- integer divisible by m
c
c       output:
c       m -- greatest integer less than or equal to l that divides n
c
        implicit none
        integer n,l,m
c
c
        m = l
c
 1000   continue
        if(m*(n/m) .eq. n) goto 2000
c
          m = m-1
          goto 1000
c
 2000   continue
c
c
        return
        end
c
c
c
c
        subroutine idd_sffti(l,ind,n,wsave)
c
c       initializes wsave for using routine idd_sfft.
c
c       input:
c       l -- number of pairs of entries in the output of idd_sfft
c            to compute
c       ind -- indices of the pairs of entries in the output
c              of idd_sfft to compute; the indices must be chosen
c              in the range from 1 to n/2
c       n -- length of the vector to be transformed
c
c       output:
c       wsave -- array needed by routine idd_sfft for processing
c                (the present routine does not use the last n elements
c                 of wsave, but routine idd_sfft does)
c
        implicit none
        integer l,ind(l),n
        complex*16 wsave(2*l+15+4*n)
c
c
        if(l .eq. 1) call idd_sffti1(ind,n,wsave)
        if(l .gt. 1) call idd_sffti2(l,ind,n,wsave)
c
c
        return
        end
c
c
c
c
        subroutine idd_sffti1(ind,n,wsave)
c
c       routine idd_sffti serves as a wrapper around
c       the present routine; please see routine idd_sffti
c       for documentation.
c
        implicit none
        integer ind,n,k
        real*8 r1,twopi,wsave(2*(2+15+4*n)),fact
c
        r1 = 1
        twopi = 2*4*atan(r1)
c
c
        fact = 1/sqrt(r1*n)
c
c
        do k = 1,n
          wsave(k) = cos(twopi*(k-1)*ind/(r1*n))*fact
        enddo ! k
c
        do k = 1,n
          wsave(n+k) = -sin(twopi*(k-1)*ind/(r1*n))*fact
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine idd_sffti2(l,ind,n,wsave)
c
c       routine idd_sffti serves as a wrapper around
c       the present routine; please see routine idd_sffti
c       for documentation.
c
        implicit none
        integer l,ind(l),n,nblock,ii,m,idivm,imodm,i,j,k
        real*8 r1,twopi,fact
        complex*16 wsave(2*l+15+4*n),ci,twopii
c
        ci = (0,1)
        r1 = 1
        twopi = 2*4*atan(r1)
        twopii = twopi*ci
c
c
c       Determine the block lengths for the FFTs.
c
        call idd_ldiv(l,n,nblock)
        m = n/nblock
c
c
c       Initialize wsave for using routine dfftf.
c
        call dffti(nblock,wsave)
c
c
c       Calculate the coefficients in the linear combinations
c       needed for the direct portion of the calculation.
c
        fact = 1/sqrt(r1*n)
c
        ii = 2*l+15
c
        do j = 1,l
c
c
          i = ind(j)
c
c
          if(i .le. n/2-m/2) then
c
            idivm = (i-1)/m
            imodm = (i-1)-m*idivm
c
            do k = 1,m
              wsave(ii+m*(j-1)+k) = exp(-twopii*(k-1)*imodm/(r1*m))
     1         * exp(-twopii*(k-1)*(idivm+1)/(r1*n)) * fact
            enddo ! k
c
          endif ! i .le. n/2-m/2
c
c
          if(i .gt. n/2-m/2) then
c
            idivm = i/(m/2)
            imodm = i-(m/2)*idivm
c
            do k = 1,m
              wsave(ii+m*(j-1)+k) = exp(-twopii*(k-1)*imodm/(r1*m))
     1                            * fact
            enddo ! k
c
          endif ! i .gt. n/2-m/2
c
c
        enddo ! j
c
c
        return
        end
c
c
c
c
        subroutine idd_sfft(l,ind,n,wsave,v)
c
c       computes a subset of the entries of the DFT of v,
c       composed with permutation matrices both on input and on output,
c       via a two-stage procedure (debugging code routine dfftf2 above
c       is supposed to calculate the full vector from which idd_sfft
c       returns a subset of the entries, when dfftf2 has
c       the same parameter nblock as in the present routine).
c
c       input:
c       l -- number of pairs of entries in the output to compute
c       ind -- indices of the pairs of entries in the output
c              to compute; the indices must be chosen
c              in the range from 1 to n/2
c       n -- length of v; n must be a positive integer power of 2
c       v -- vector to be transformed
c       wsave -- processing array initialized by routine idd_sffti
c
c       output:
c       v -- pairs of entries indexed by ind are given
c            their appropriately transformed values
c
c       _N.B._: n must be a positive integer power of 2.
c
c       references:
c       Sorensen and Burrus, "Efficient computation of the DFT with
c            only a subset of input or output points,"
c            IEEE Transactions on Signal Processing, 41 (3): 1184-1200,
c            1993.
c       Woolfe, Liberty, Rokhlin, Tygert, "A fast randomized algorithm
c            for the approximation of matrices," Applied and
c            Computational Harmonic Analysis, 25 (3): 335-366, 2008;
c            Section 3.3.
c
        implicit none
        integer l,ind(l),n
        real*8 v(n)
        complex*16 wsave(2*l+15+4*n)
c
c
        if(l .eq. 1) call idd_sfft1(ind,n,v,wsave)
        if(l .gt. 1) call idd_sfft2(l,ind,n,v,wsave)
c
c
        return
        end
c
c
c
c
        subroutine idd_sfft1(ind,n,v,wsave)
c
c       routine idd_sfft serves as a wrapper around
c       the present routine; please see routine idd_sfft
c       for documentation.
c
        implicit none
        integer ind,n,k
        real*8 v(n),r1,twopi,sumr,sumi,fact,wsave(2*(2+15+4*n))
c
        r1 = 1
        twopi = 2*4*atan(r1)
c
c
        if(ind .lt. n/2) then
c
c
          sumr = 0
c
          do k = 1,n
            sumr = sumr+wsave(k)*v(k)
          enddo ! k
c
c
          sumi = 0
c
          do k = 1,n
            sumi = sumi+wsave(n+k)*v(k)
          enddo ! k
c
c
        endif ! ind .lt. n/2
c
c
        if(ind .eq. n/2) then
c
c
          fact = 1/sqrt(r1*n)
c
c
          sumr = 0
c
          do k = 1,n
            sumr = sumr+v(k)
          enddo ! k
c
          sumr = sumr*fact
c
c
          sumi = 0
c
          do k = 1,n/2
            sumi = sumi+v(2*k-1)
            sumi = sumi-v(2*k)
          enddo ! k
c
          sumi = sumi*fact
c
c
        endif ! ind .eq. n/2
c
c
        v(2*ind-1) = sumr
        v(2*ind) = sumi
c
c
        return
        end
c
c
c
c
        subroutine idd_sfft2(l,ind,n,v,wsave)
c
c       routine idd_sfft serves as a wrapper around
c       the present routine; please see routine idd_sfft
c       for documentation.
c
        implicit none
        integer n,m,l,k,j,ind(l),i,idivm,nblock,ii,iii,imodm
        real*8 r1,twopi,v(n),rsum,fact
        complex*16 wsave(2*l+15+4*n),ci,sum
c
        ci = (0,1)
        r1 = 1
        twopi = 2*4*atan(r1)
c
c
c       Determine the block lengths for the FFTs.
c
        call idd_ldiv(l,n,nblock)
c
c
        m = n/nblock
c
c
c       FFT each block of length nblock of v.
c
        do k = 1,m
          call dfftf(nblock,v(nblock*(k-1)+1),wsave)
        enddo ! k
c
c
c       Transpose v to obtain wsave(2*l+15+2*n+1 : 2*l+15+3*n).
c
        iii = 2*l+15+2*n
c
        do k = 1,m
          do j = 1,nblock/2-1
            wsave(iii+m*(j-1)+k) = v(nblock*(k-1)+2*j)
     1                           + ci*v(nblock*(k-1)+2*j+1)
          enddo ! j
        enddo ! k
c
c       Handle the purely real frequency components separately.
c
        do k = 1,m
          wsave(iii+m*(nblock/2-1)+k) = v(nblock*(k-1)+nblock)
          wsave(iii+m*(nblock/2)+k) = v(nblock*(k-1)+1)
        enddo ! k
c
c
c       Directly calculate the desired entries of v.
c
        ii = 2*l+15
c
        do j = 1,l
c
c
          i = ind(j)
c
c
          if(i .le. n/2-m/2) then
c
            idivm = (i-1)/m
            imodm = (i-1)-m*idivm
c
            sum = 0
c
            do k = 1,m
              sum = sum + wsave(iii+m*idivm+k) * wsave(ii+m*(j-1)+k)
            enddo ! k
c
            v(2*i-1) = sum
            v(2*i) = -ci*sum
c
          endif ! i .le. n/2-m/2
c
c
          if(i .gt. n/2-m/2) then
c
            if(i .lt. n/2) then
c
              idivm = i/(m/2)
              imodm = i-(m/2)*idivm
c
              sum = 0
c
              do k = 1,m
                sum = sum + wsave(iii+m*(nblock/2)+k)
     1              * wsave(ii+m*(j-1)+k)
              enddo ! k
c
              v(2*i-1) = sum
              v(2*i) = -ci*sum
c
            endif
c
            if(i .eq. n/2) then
c
              fact = 1/sqrt(r1*n)
c
c
              rsum = 0
c
              do k = 1,m
                rsum = rsum + wsave(iii+m*(nblock/2)+k)
              enddo ! k
c
              v(n-1) = rsum*fact
c
c
              rsum = 0
c
              do k = 1,m/2
                rsum = rsum + wsave(iii+m*(nblock/2)+2*k-1)
                rsum = rsum - wsave(iii+m*(nblock/2)+2*k)
              enddo ! k
c
              v(n) = rsum*fact
c
            endif
c
          endif ! i .gt. n/2-m/2
c
c
        enddo ! j
c
c
        return
        end
