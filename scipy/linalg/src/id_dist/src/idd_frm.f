c       this file contains the following user-callable routines:
c
c
c       routine idd_frm transforms a vector via a composition
c       of Rokhlin's random transform, random subselection, and an FFT.
c
c       routine idd_sfrm transforms a vector into a vector
c       of specified length via a composition
c       of Rokhlin's random transform, random subselection, and an FFT.
c
c       routine idd_frmi initializes routine idd_frm.
c
c       routine idd_sfrmi initializes routine idd_sfrm.
c
c       routine idd_pairsamps calculates the indices of the pairs
c       of integers to which the individual integers
c       in a specified set belong.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine idd_frm(m,n,w,x,y)
c
c       transforms x into y via a composition
c       of Rokhlin's random transform, random subselection, and an FFT.
c       In contrast to routine idd_sfrm, the present routine works best
c       when the length of the transformed vector is the integer n
c       output by routine idd_frmi, or when the length
c       is not specified, but instead determined a posteriori
c       using the output of the present routine. The transformed vector
c       output by the present routine is randomly permuted.
c
c       input:
c       m -- length of x
c       n -- greatest integer expressible as a positive integer power
c            of 2 that is less than or equal to m, as obtained
c            from the routine idd_frmi; n is the length of y
c       w -- initialization array constructed by routine idd_frmi
c       x -- vector to be transformed
c
c       output:
c       y -- transform of x
c
c       reference:
c       Halko, Martinsson, Tropp, "Finding structure with randomness:
c            probabilistic algorithms for constructing approximate
c            matrix decompositions," SIAM Review, 53 (2): 217-288,
c            2011.
c
        implicit none
        integer m,iw,n,k
        real*8 w(17*m+70),x(m),y(n)
c
c
c       Apply Rokhlin's random transformation to x, obtaining
c       w(16*m+71 : 17*m+70).
c
        iw = w(3+m+n)
        call idd_random_transf(x,w(16*m+70+1),w(iw))
c
c
c       Subselect from  w(16*m+71 : 17*m+70)  to obtain y.
c
        call idd_subselect(n,w(3),m,w(16*m+70+1),y)
c
c
c       Copy y into  w(16*m+71 : 16*m+n+70).
c
        do k = 1,n
          w(16*m+70+k) = y(k)
        enddo ! k
c
c
c       Fourier transform  w(16*m+71 : 16*m+n+70).
c
        call dfftf(n,w(16*m+70+1),w(4+m+n))
c
c
c       Permute  w(16*m+71 : 16*m+n+70)  to obtain y.
c
        call idd_permute(n,w(3+m),w(16*m+70+1),y)
c
c
        return
        end
c
c
c
c
        subroutine idd_sfrm(l,m,n,w,x,y)
c
c       transforms x into y via a composition
c       of Rokhlin's random transform, random subselection, and an FFT.
c       In contrast to routine idd_frm, the present routine works best
c       when the length l of the transformed vector is known a priori.
c
c       input:
c       l -- length of y; l must be less than or equal to n
c       m -- length of x
c       n -- greatest integer expressible as a positive integer power
c            of 2 that is less than or equal to m, as obtained
c            from the routine idd_sfrmi
c       w -- initialization array constructed by routine idd_sfrmi
c       x -- vector to be transformed
c
c       output:
c       y -- transform of x
c
c       _N.B._: l must be less than or equal to n.
c
c       reference:
c       Halko, Martinsson, Tropp, "Finding structure with randomness:
c            probabilistic algorithms for constructing approximate
c            matrix decompositions," SIAM Review, 53 (2): 217-288,
c            2011.
c
        implicit none
        integer m,iw,n,l,l2
        real*8 w(27*m+90),x(m),y(l)
c
c
c       Retrieve the number of pairs of outputs to be calculated
c       via sfft.
c
        l2 = w(3)
c
c
c       Apply Rokhlin's random transformation to x, obtaining
c       w(25*m+91 : 26*m+90).
c
        iw = w(4+m+l+l2)
        call idd_random_transf(x,w(25*m+90+1),w(iw))
c
c
c       Subselect from  w(25*m+91 : 26*m+90)  to obtain
c       w(26*m+91 : 26*m+n+90).
c
        call idd_subselect(n,w(4),m,w(25*m+90+1),w(26*m+90+1))
c
c
c       Fourier transform  w(26*m+91 : 26*m+n+90).
c
        call idd_sfft(l2,w(4+m+l),n,w(5+m+l+l2),w(26*m+90+1))
c
c
c       Copy the desired entries from  w(26*m+91 : 26*m+n+90)
c       to y.
c
        call idd_subselect(l,w(4+m),n,w(26*m+90+1),y)
c
c
        return
        end
c
c
c
c
        subroutine idd_pairsamps(n,l,ind,l2,ind2,marker)
c
c       calculates the indices of the l2 pairs of integers
c       to which the l individual integers from ind belong.
c       The integers in ind may range from 1 to n.
c
c       input:
c       n -- upper bound on the integers in ind
c            (the number 1 must be a lower bound);
c            n must be even
c       l -- length of ind
c       ind -- integers selected from 1 to n
c
c       output:
c       l2 -- length of ind2
c       ind2 -- indices in the range from 1 to n/2 of the pairs
c               of integers to which the entries of ind belong
c
c       work:
c       marker -- must be at least n/2 integer elements long
c
c       _N.B._: n must be even.
c
        implicit none
        integer l,n,ind(l),ind2(l),marker(n/2),l2,k
c
c
c       Unmark all pairs.
c
        do k = 1,n/2
          marker(k) = 0
        enddo ! k
c
c
c       Mark the required pairs.
c
        do k = 1,l
          marker((ind(k)+1)/2) = marker((ind(k)+1)/2)+1
        enddo ! k
c
c
c       Record the required pairs in indpair.
c
        l2 = 0
c
        do k = 1,n/2
c
          if(marker(k) .ne. 0) then
            l2 = l2+1
            ind2(l2) = k
          endif
c
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine idd_permute(n,ind,x,y)
c
c       copy the entries of x into y, rearranged according
c       to the permutation specified by ind.
c
c       input:
c       n -- length of ind, x, and y
c       ind -- permutation of n objects
c       x -- vector to be permuted
c
c       output:
c       y -- permutation of x
c
        implicit none
        integer n,ind(n),k
        real*8 x(n),y(n)
c
c
        do k = 1,n
          y(k) = x(ind(k))
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine idd_subselect(n,ind,m,x,y)
c
c       copies into y the entries of x indicated by ind.
c
c       input:
c       n -- number of entries of x to copy into y
c       ind -- indices of the entries in x to copy into y
c       m -- length of x
c       x -- vector whose entries are to be copied
c
c       output:
c       y -- collection of entries of x specified by ind
c
        implicit none
        integer n,ind(n),m,k
        real*8 x(m),y(n)
c
c
        do k = 1,n
          y(k) = x(ind(k))
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine idd_frmi(m,n,w)
c
c       initializes data for the routine idd_frm.
c
c       input:
c       m -- length of the vector to be transformed
c
c       output:
c       n -- greatest integer expressible as a positive integer power
c            of 2 that is less than or equal to m
c       w -- initialization array to be used by routine idd_frm
c
c
c       glossary for the fully initialized w:
c
c       w(1) = m
c       w(2) = n
c       w(3:2+m) stores a permutation of m objects
c       w(3+m:2+m+n) stores a permutation of n objects
c       w(3+m+n) = address in w of the initialization array
c                  for idd_random_transf
c       w(4+m+n:int(w(3+m+n))-1) stores the initialization array
c                                for dfft
c       w(int(w(3+m+n)):16*m+70) stores the initialization array
c                                for idd_random_transf
c
c
c       _N.B._: n is an output of the present routine;
c               this routine changes n.
c
c
        implicit none
        integer m,n,l,nsteps,keep,lw,ia
        real*8 w(17*m+70)
c
c
c       Find the greatest integer less than or equal to m
c       which is a power of two.
c
        call idd_poweroftwo(m,l,n)
c
c
c       Store m and n in w.
c
        w(1) = m
        w(2) = n
c
c
c       Store random permutations of m and n objects in w.
c
        call id_randperm(m,w(3))
        call id_randperm(n,w(3+m))
c
c
c       Store the address within w of the idd_random_transf_init
c       initialization data.
c
        ia = 4+m+n+2*n+15
        w(3+m+n) = ia
c
c
c       Store the initialization data for dfft in w.
c
        call dffti(n,w(4+m+n))
c
c
c       Store the initialization data for idd_random_transf_init in w.
c
        nsteps = 3
        call idd_random_transf_init(nsteps,m,w(ia),keep)
c
c
c       Calculate the total number of elements used in w.
c
        lw = 3+m+n+2*n+15 + 3*nsteps*m+2*m+m/4+50
c
        if(16*m+70 .lt. lw) then
          call prinf('lw = *',lw,1)
          call prinf('16m+70 = *',16*m+70,1)
          stop
        endif
c
c
        return
        end
c
c
c
c
        subroutine idd_sfrmi(l,m,n,w)
c
c       initializes data for the routine idd_sfrm.
c
c       input:
c       l -- length of the transformed (output) vector
c       m -- length of the vector to be transformed
c
c       output:
c       n -- greatest integer expressible as a positive integer power
c            of 2 that is less than or equal to m
c       w -- initialization array to be used by routine idd_sfrm
c
c
c       glossary for the fully initialized w:
c
c       w(1) = m
c       w(2) = n
c       w(3) = l2
c       w(4:3+m) stores a permutation of m objects
c       w(4+m:3+m+l) stores the indices of the l outputs which idd_sfft
c                    calculates
c       w(4+m+l:3+m+l+l2) stores the indices of the l2 pairs of outputs
c                         which idd_sfft calculates
c       w(4+m+l+l2) = address in w of the initialization array
c                     for idd_random_transf
c       w(5+m+l+l2:int(w(4+m+l+l2))-1) stores the initialization array
c                                      for idd_sfft
c       w(int(w(4+m+l+l2)):25*m+90) stores the initialization array
c                                   for idd_random_transf
c
c
c       _N.B._: n is an output of the present routine;
c               this routine changes n.
c
c
        implicit none
        integer l,m,n,idummy,nsteps,keep,lw,l2,ia
        real*8 w(27*m+90)
c
c
c       Find the greatest integer less than or equal to m
c       which is a power of two.
c
        call idd_poweroftwo(m,idummy,n)
c
c
c       Store m and n in w.
c
        w(1) = m
        w(2) = n
c
c
c       Store random permutations of m and n objects in w.
c
        call id_randperm(m,w(4))
        call id_randperm(n,w(4+m))
c
c
c       Find the pairs of integers covering the integers in
c       w(4+m : 3+m+(l+1)/2).
c
        call idd_pairsamps(n,l,w(4+m),l2,w(4+m+2*l),w(4+m+3*l))
        w(3) = l2
        call idd_copyints(l2,w(4+m+2*l),w(4+m+l))
c
c
c       Store the address within w of the idd_random_transf_init
c       initialization data.
c
        ia = 5+m+l+l2+4*l2+30+8*n
        w(4+m+l+l2) = ia
c
c
c       Store the initialization data for idd_sfft in w.
c
        call idd_sffti(l2,w(4+m+l),n,w(5+m+l+l2))
c
c
c       Store the initialization data for idd_random_transf_init in w.
c
        nsteps = 3
        call idd_random_transf_init(nsteps,m,w(ia),keep)
c
c
c       Calculate the total number of elements used in w.
c
        lw = 4+m+l+l2+4*l2+30+8*n + 3*nsteps*m+2*m+m/4+50
c
        if(25*m+90 .lt. lw) then
          call prinf('lw = *',lw,1)
          call prinf('25m+90 = *',25*m+90,1)
          stop
        endif
c
c
        return
        end
c
c
c
c
        subroutine idd_copyints(n,ia,ib)
c
c       copies ia into ib.
c
c       input:
c       n -- length of ia and ib
c       ia -- array to be copied
c
c       output:
c       ib -- copy of ia
c
        implicit none
        integer n,ia(n),ib(n),k
c
c
        do k = 1,n
          ib(k) = ia(k)
        enddo ! k
c
c
        return
        end
c
c
c
c
        subroutine idd_poweroftwo(m,l,n)
c
c       computes l = floor(log_2(m)) and n = 2**l.
c
c       input:
c       m -- integer whose log_2 is to be taken
c
c       output:
c       l -- floor(log_2(m))
c       n -- 2**l
c
        implicit none
        integer l,m,n
c
c
        l = 0
        n = 1
c
 1000   continue
          l = l+1
          n = n*2
        if(n .le. m) goto 1000
c
        l = l-1
        n = n/2
c
c
        return
        end
