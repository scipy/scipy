c       this file contains the following user-callable routines:
c
c 
c       routine idd_random_transf applies rapidly
c       a random orthogonal matrix to a user-supplied vector.
c   
c       routine idd_random_transf_inverse applies rapidly
c       the inverse of the operator applied
c       by routine idd_random_transf.
c
c       routine idz_random_transf applies rapidly
c       a random unitary matrix to a user-supplied vector.
c
c       routine idz_random_transf_inverse applies rapidly
c       the inverse of the operator applied
c       by routine idz_random_transf.
c
c       routine idd_random_transf_init initializes data
c       for routines idd_random_transf and idd_random_transf_inverse.
c
c       routine idz_random_transf_init initializes data
c       for routines idz_random_transf and idz_random_transf_inverse.
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine idd_random_transf_init(nsteps,n,w,keep)
        implicit real *8 (a-h,o-z)
        save
        dimension w(*)
c
c       prepares and stores in array w the data used
c       by the routines idd_random_transf and idd_random_transf_inverse
c       to apply rapidly a random orthogonal matrix
c       to an arbitrary user-specified vector.
c
c       input:
c       nsteps -- the degree of randomness of the operator
c                 to be applied
c       n -- the size of the matrix to be applied     
c
c       output:
c       w -- the first keep elements of w contain all the data
c            to be used by routines idd_random_tranf
c            and idd_random_transf_inverse. Please note that
c            the number of elements used by the present routine
c            is also equal to keep. This array should be at least
c            3*nsteps*n + 2*n + n/4 + 50 real*8 elements long.
c       keep - the number of elements in w actually used 
c              by the present routine; keep is also the number
c              of elements that must not be changed between the call
c              to this routine and subsequent calls to routines
c              idd_random_transf and idd_random_transf_inverse.
c
c
c        . . . allocate memory 
c
        ninire=2
c
        ialbetas=10
        lalbetas=2*n*nsteps+10 
c
        iixs=ialbetas+lalbetas
        lixs=n*nsteps/ninire+10
c
        iww=iixs+lixs
        lww=2*n+n/4+20
c
        keep=iww+lww
c
        w(1)=ialbetas+0.1
        w(2)=iixs+0.1
        w(3)=nsteps+0.1
        w(4)=iww+0.1        
        w(5)=n+0.1
c
        call idd_random_transf_init0(nsteps,n,w(ialbetas),w(iixs))
c
        return
        end
c
c 
c 
c
c 
        subroutine idz_random_transf_init(nsteps,n,w,keep)
        implicit real *8 (a-h,o-z)
        save
        dimension w(*)
c
c       prepares and stores in array w the data used
c       by routines idz_random_transf and idz_random_transf_inverse
c       to apply rapidly a random unitary matrix
c       to an arbitrary user-specified vector.
c
c       input:
c       nsteps -- the degree of randomness of the operator
c                 to be applied
c       n -- the size of the matrix to be applied     
c
c       output:
c       w -- the first keep elements of w contain all the data
c            to be used by routines idz_random_transf
c            and idz_random_transf_inverse. Please note that
c            the number of elements used by the present routine
c            is also equal to keep. This array should be at least
c            5*nsteps*n + 2*n + n/4 + 60 real*8 elements long.
c       keep - the number of elements in w actually used
c              by the present routine; keep is also the number
c              of elements that must not be changed between the call
c              to this routine and subsequent calls to routines
c              idz_random_transf and idz_random_transf_inverse.
c
c
c        . . . allocate memory 
c
        ninire=2
c
        ialbetas=10
        lalbetas=2*n*nsteps+10 
c
        igammas=ialbetas+lalbetas
        lgammas=2*n*nsteps+10
c
        iixs=igammas+lgammas
        lixs=n*nsteps/ninire+10
c
        iww=iixs+lixs
        lww=2*n+n/4+20
c
        keep=iww+lww
c
        w(1)=ialbetas+0.1
        w(2)=iixs+0.1
        w(3)=nsteps+0.1
        w(4)=iww+0.1        
        w(5)=n+0.1
        w(6)=igammas+0.1
c
        call idz_random_transf_init0(nsteps,n,w(ialbetas),
     1      w(igammas),w(iixs))
c
        return
        end
c
c 
c 
c 
c 
        subroutine idd_random_transf(x,y,w)
        implicit real *8 (a-h,o-z)
        save
        dimension x(*),y(*),w(*)
c
c       applies rapidly a random orthogonal matrix
c       to the user-specified real vector x, 
c       using the data in array w stored there by a preceding 
c       call to routine idd_random_transf_init.
c
c       input:
c       x -- the vector of length n to which the random matrix is
c            to be applied
c       w -- array containing all initialization data
c
c       output:
c       y -- the result of applying the random matrix to x
c
c
c        . . . allocate memory
c
        ialbetas=w(1)
        iixs=w(2)
        nsteps=w(3)
        iww=w(4)
        n=w(5)
c
        call idd_random_transf0(nsteps,x,y,n,w(iww),
     1      w(ialbetas),w(iixs))
c
        return
        end
c
c 
c 
c 
c 
        subroutine idd_random_transf_inverse(x,y,w)
        implicit real *8 (a-h,o-z)
        save
        dimension x(*),y(*),w(*)
c
c       applies rapidly a random orthogonal matrix
c       to the user-specified real vector x, 
c       using the data in array w stored there by a preceding 
c       call to routine idd_random_transf_init.
c       The transformation applied by the present routine is
c       the inverse of the transformation applied
c       by routine idd_random_transf.
c
c       input:
c       x -- the vector of length n to which the random matrix is
c            to be applied
c       w -- array containing all initialization data
c
c       output:
c       y -- the result of applying the random matrix to x
c
c
c        . . . allocate memory
c
        ialbetas=w(1)
        iixs=w(2)
        nsteps=w(3)
        iww=w(4)
        n=w(5)
c
        call idd_random_transf0_inv(nsteps,x,y,n,w(iww),
     1      w(ialbetas),w(iixs))
c
        return
        end
c
c 
c 
c 
c 
        subroutine idz_random_transf(x,y,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(*),y(*)
        dimension w(*)
c
c       applies rapidly a random unitary matrix
c       to the user-specified vector x, 
c       using the data in array w stored there by a preceding 
c       call to routine idz_random_transf_init.
c
c       input:
c       x -- the vector of length n to which the random matrix is
c            to be applied
c       w -- array containing all initialization data
c
c       output:
c       y -- the result of applying the random matrix to x
c
c
c        . . . allocate memory
c
        ialbetas=w(1)
        iixs=w(2)
        nsteps=w(3)
        iww=w(4)
        n=w(5)
        igammas=w(6)
c
        call idz_random_transf0(nsteps,x,y,n,w(iww),w(ialbetas),
     1      w(igammas),w(iixs))
c
        return
        end
c
c 
c 
c 
c 
        subroutine idz_random_transf_inverse(x,y,w)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(*),y(*)
        dimension w(*)
c
c       applies rapidly a random unitary matrix
c       to the user-specified vector x,
c       using the data in array w stored there by a preceding 
c       call to routine idz_random_transf_init.
c       The transformation applied by the present routine is
c       the inverse of the transformation applied
c       by routine idz_random_transf.
c
c       input:
c       x -- the vector of length n to which the random matrix is
c            to be applied
c       w -- array containing all initialization data
c
c       output:
c       y -- the result of applying the random matrix to x
c
c
c        . . . allocate memory
c
        ialbetas=w(1)
        iixs=w(2)
        nsteps=w(3)
        iww=w(4)
        n=w(5)
        igammas=w(6)
c
        call idz_random_transf0_inv(nsteps,x,y,n,w(iww),
     1      w(ialbetas),w(igammas),w(iixs))
c
        return
        end
c
c 
c 
c 
c 
        subroutine idd_random_transf0_inv(nsteps,x,y,n,w2,albetas,iixs)
        implicit real *8 (a-h,o-z)
        save
        dimension x(*),y(*),w2(*),albetas(2,n,*),iixs(n,*)
c
c       routine idd_random_transf_inverse serves as a memory wrapper
c       for the present routine; see routine idd_random_transf_inverse
c       for documentation.
c
        do 1200 i=1,n
c
        w2(i)=x(i)
 1200 continue
c
        do 2000 ijk=nsteps,1,-1
c
        call idd_random_transf00_inv(w2,y,n,albetas(1,1,ijk),
     1      iixs(1,ijk) )
c
        do 1400 j=1,n
c
        w2(j)=y(j)
 1400 continue
 2000 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine idd_random_transf00_inv(x,y,n,albetas,ixs)
        implicit real *8 (a-h,o-z)
        save
        dimension x(*),y(*),albetas(2,*),ixs(*)
c
c       implements one step of the random transform required
c       by routine idd_random_transf0_inv (please see the latter).
c
c
c        implement 2 \times 2 matrices
c
        do 1600 i=1,n
        y(i)=x(i)
 1600 continue
c
        do 1800 i=n-1,1,-1
c
        alpha=albetas(1,i)
        beta=albetas(2,i)
c
        a=y(i)
        b=y(i+1)
c
        y(i)=alpha*a-beta*b
        y(i+1)=beta*a+alpha*b
 1800 continue
c
c        implement the permutation
c
        do 2600 i=1,n
c
        j=ixs(i)
        x(j)=y(i)
 2600 continue
c
        do 2800 i=1,n
c
        y(i)=x(i)
 2800 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine idz_random_transf0_inv(nsteps,x,y,n,w2,albetas,
     1      gammas,iixs)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(*),y(*),w2(*),gammas(n,*)
        dimension albetas(2,n,*),iixs(n,*)
c
c       routine idz_random_transf_inverse serves as a memory wrapper
c       for the present routine; please see routine
c       idz_random_transf_inverse for documentation.
c
        do 1200 i=1,n
c
        w2(i)=x(i)
 1200 continue
c
        do 2000 ijk=nsteps,1,-1
c
        call idz_random_transf00_inv(w2,y,n,albetas(1,1,ijk),
     1      gammas(1,ijk),iixs(1,ijk) )
c
        do 1400 j=1,n
c
        w2(j)=y(j)
 1400 continue
 2000 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine idz_random_transf00_inv(x,y,n,albetas,gammas,ixs)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(*),y(*),gammas(*),a,b
        dimension albetas(2,*),ixs(*)
c
c       implements one step of the random transform
c       required by routine idz_random_transf0_inv
c       (please see the latter).
c
c        implement 2 \times 2 matrices
c
        do 1600 i=n-1,1,-1
c
        alpha=albetas(1,i)
        beta=albetas(2,i)
c
        a=x(i)
        b=x(i+1)
c
        x(i)=alpha*a-beta*b
        x(i+1)=beta*a+alpha*b
 1600 continue
c
c        implement the permutation
c        and divide by the random numbers on the unit circle
c        (or, equivalently, multiply by their conjugates)
c
        do 1800 i=1,n
c
        j=ixs(i)
        y(j)=x(i)*conjg(gammas(i))
 1800 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine idd_random_transf0(nsteps,x,y,n,w2,albetas,iixs)
        implicit real *8 (a-h,o-z)
        save
        dimension x(*),y(*),w2(*),albetas(2,n,*),iixs(n,*)
c
c       routine idd_random_transf serves as a memory wrapper
c       for the present routine; please see routine idd_random_transf
c       for documentation.
c
        do 1200 i=1,n
c
        w2(i)=x(i)
 1200 continue
c
        do 2000 ijk=1,nsteps
c
        call idd_random_transf00(w2,y,n,albetas(1,1,ijk),iixs(1,ijk) )
c
        do 1400 j=1,n
c
        w2(j)=y(j)
 1400 continue
 2000 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine idd_random_transf00(x,y,n,albetas,ixs)
        implicit real *8 (a-h,o-z)
        save
        dimension x(*),y(*),albetas(2,*),ixs(*)
c
c       implements one step of the random transform
c       required by routine idd_random_transf0 (please see the latter).
c
c        implement the permutation
c
        do 1600 i=1,n
c
        j=ixs(i)
        y(i)=x(j)
 1600 continue
c
c        implement 2 \times 2 matrices
c
        do 1800 i=1,n-1
c
        alpha=albetas(1,i)
        beta=albetas(2,i)
c
        a=y(i)
        b=y(i+1)
c
        y(i)=alpha*a+beta*b
        y(i+1)=-beta*a+alpha*b
 1800 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine idz_random_transf_init0(nsteps,n,albetas,gammas,ixs)
        implicit real *8 (a-h,o-z)
        save
        dimension albetas(2,n,*),ixs(n,*)
        complex *16 gammas(n,*)
c
c       routine idz_random_transf_init serves as a memory wrapper
c       for the present routine; please see routine
c       idz_random_transf_init for documentation.
c
        do 2000 ijk=1,nsteps
c
        call idz_random_transf_init00(n,albetas(1,1,ijk),
     1      gammas(1,ijk),ixs(1,ijk) )
 2000 continue
        return
        end
c
c 
c 
c 
c 
        subroutine idz_random_transf_init00(n,albetas,gammas,ixs)
        implicit real *8 (a-h,o-z)
        save
        dimension albetas(2,*),gammas(*),ixs(*)
c
c       constructs one stage of the random transform
c       initialized by routine idz_random_transf_init0
c       (please see the latter).
c
        done=1
        twopi=2*4*atan(done)
c
c        construct the random permutation
c
        ifrepeat=0
        call id_randperm(n,ixs)
c
c        construct the random variables
c
        call id_srand(2*n,albetas)
        call id_srand(2*n,gammas)
c
        do 1300 i=1,n
c
        albetas(1,i)=2*albetas(1,i)-1
        albetas(2,i)=2*albetas(2,i)-1
        gammas(2*i-1)=2*gammas(2*i-1)-1
        gammas(2*i)=2*gammas(2*i)-1
 1300 continue
c
c        construct the random 2 \times 2 transformations
c
        do 1400 i=1,n
c
        d=albetas(1,i)**2+albetas(2,i)**2
        d=1/sqrt(d)
        albetas(1,i)=albetas(1,i)*d
        albetas(2,i)=albetas(2,i)*d
 1400 continue
c
c        construct the random multipliers on the unit circle
c
        do 1500 i=1,n
c
        d=gammas(2*i-1)**2+gammas(2*i)**2
        d=1/sqrt(d)
c
c        fill the real part
c
        gammas(2*i-1)=gammas(2*i-1)*d
c
c        fill the imaginary part
c
        gammas(2*i)=gammas(2*i)*d
 1500 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine idz_random_transf0(nsteps,x,y,n,w2,albetas,
     1      gammas,iixs)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(*),y(*),w2(*),gammas(n,*)
        dimension albetas(2,n,*),iixs(n,*)
c
c       routine idz_random_transf serves as a memory wrapper
c       for the present routine; please see routine idz_random_transf
c       for documentation.
c
        do 1200 i=1,n
c
        w2(i)=x(i)
 1200 continue
c
        do 2000 ijk=1,nsteps
c
        call idz_random_transf00(w2,y,n,albetas(1,1,ijk),
     1      gammas(1,ijk),iixs(1,ijk) )
        do 1400 j=1,n
c
        w2(j)=y(j)
 1400 continue
 2000 continue
c
        return
        end
c
c 
c 
c 
c 
        subroutine idz_random_transf00(x,y,n,albetas,gammas,ixs)
        implicit real *8 (a-h,o-z)
        save
        complex *16 x(*),y(*),gammas(*),a,b
        dimension albetas(2,*),ixs(*)
c
c       implements one step of the random transform
c       required by routine idz_random_transf0 (please see the latter).
c
c        implement the permutation
c        and multiply by the random numbers
c        on the unit circle
c
        do 1600 i=1,n
c
        j=ixs(i)
        y(i)=x(j)*gammas(i)
 1600 continue
c
c        implement 2 \times 2 matrices
c
        do 2600 i=1,n-1
c
        alpha=albetas(1,i)
        beta=albetas(2,i)
c
        a=y(i)
        b=y(i+1)
c
        y(i)=alpha*a+beta*b
        y(i+1)=-beta*a+alpha*b
 2600 continue
c
        return
        end
c
c
c
c
c
        subroutine idd_random_transf_init0(nsteps,n,albetas,ixs)
        implicit real *8 (a-h,o-z)
        save
        dimension albetas(2,n,*),ixs(n,*)
c
c       routine idd_random_transf_init serves as a memory wrapper
c       for the present routine; please see routine
c       idd_random_transf_init for documentation.
c
        do 2000 ijk=1,nsteps
c
        call idd_random_transf_init00(n,albetas(1,1,ijk),ixs(1,ijk) )
 2000 continue
        return
        end
c
c 
c 
c 
c 
        subroutine idd_random_transf_init00(n,albetas,ixs)
        implicit real *8 (a-h,o-z)
        save
        dimension albetas(2,*),ixs(*)
c
c       constructs one stage of the random transform
c       initialized by routine idd_random_transf_init0
c       (please see the latter).
c
c        construct the random permutation
c
        ifrepeat=0
        call id_randperm(n,ixs)
c
c        construct the random variables
c
        call id_srand(2*n,albetas)
c
        do 1300 i=1,n
c
        albetas(1,i)=2*albetas(1,i)-1
        albetas(2,i)=2*albetas(2,i)-1
 1300 continue
c
c        construct the random 2 \times 2 transformations
c
        do 1400 i=1,n
c
        d=albetas(1,i)**2+albetas(2,i)**2
        d=1/sqrt(d)
        albetas(1,i)=albetas(1,i)*d
        albetas(2,i)=albetas(2,i)*d
 1400 continue
        return
        end
