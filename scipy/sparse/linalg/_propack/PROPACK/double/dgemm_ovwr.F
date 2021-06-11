c
c     (C) Rasmus Munk Larsen, Stanford University, 2004
c

      subroutine dgemm_ovwr(transa,m,n,k,alpha,A,lda,beta,B,ldb,
     c     dwork,ldwork)
c     
c     compute B <- alpha*op(A)*B + beta*B
c
      implicit none
      character*1 transa
      integer m,n,k,lda,ldb,ldwork
      double precision alpha,beta,A(lda,*),B(ldb,*),dwork(ldwork)
      integer i,j,l,blocksize

      if((m.le.0).or.(n.le.0).or.(k.le.0)) return
      if (ldwork.lt.m) stop 'Too little workspace in DGEMM_OVWR'
      if (m.gt.ldb) stop 'm>ldb in DGEMM_OVWR'
      blocksize = int(ldwork/m)
      do i=1,n-blocksize+1,blocksize
         call dgemm(transa,'N',m,blocksize,k,alpha,A,lda,
     c              B(1,i),ldb,0D0,dwork,m)
         if (beta.eq.0D0) then
            do j=0,blocksize-1
               do l=1,m
                  B(l,i+j)  = dwork(j*m+l)
               enddo
            enddo
         else
            do j=0,blocksize-1
               do l=1,m
                  B(l,i+j)  = dwork(j*m+l) + beta*B(l,i+j)
               enddo
            enddo
         endif
      enddo
      call dgemm(transa,'N',m,n-i+1,k,alpha,A,lda,
     c           B(1,i),ldb,0D0,dwork,m)
      if (beta.eq.0D0) then
         do j=0,n-i
            do l=1,m
               B(l,i+j)  = dwork(j*m+l)
            enddo
         enddo
      else
         do j=0,n-i
            do l=1,m
               B(l,i+j)  = dwork(j*m+l) + beta*B(l,i+j)
            enddo            
         enddo
      endif      
      return
      end


      subroutine dgemm_ovwr_left(transb,m,n,k,alpha,A,lda,beta,B,ldb,
     c     dwork,ldwork)
c     
c     compute  A <- alpha*A*op(B)
c
      implicit none
      character*1 transb
      integer m,n,k,lda,ldb,ldwork
      double precision alpha,beta,A(lda,*),B(ldb,*),dwork(ldwork)
      integer i,j,l,blocksize

      if((m.le.0).or.(n.le.0).or.(k.le.0)) return
      if (ldwork.lt.n) stop 'Too little workspace in DGEMM_OVWR_LEFT'
      blocksize = int(ldwork/n)
      do i=1,m-blocksize+1,blocksize
         call dgemm('n',transb,blocksize,n,k,alpha,A(i,1),lda,
     c              B,ldb,0d0,dwork,blocksize)
         do j=0,n-1
            do l=0,blocksize-1
               A(i+l,j+1) = dwork(j*blocksize+1+l)
            enddo
         enddo
      enddo
      call dgemm('n',transb,m-i+1,n,k,alpha,A(i,1),lda,
     c           B,ldb,0d0,dwork,m-i+1)
      do j=0,n-1
         do l=0,m-i
            A(i+l,j+1) = dwork(j*(m-i+1)+1+l)
         enddo
      enddo
      return
      end


