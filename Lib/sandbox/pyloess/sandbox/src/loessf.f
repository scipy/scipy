      subroutine ehg126(d,n,vc,x,v,nvmax)
      integer d,execnt,i,j,k,n,nv,nvmax,vc
      DOUBLE PRECISION machin,alpha,beta,mu,t
      DOUBLE PRECISION v(nvmax,d),x(n,d)
      DOUBLE PRECISION D1MACH
      external D1MACH
      save machin,execnt
      data execnt /0/
c     MachInf -> machin
      execnt=execnt+1
      if(execnt.eq.1)then
         machin=D1MACH(2)
      end if
c     fill in vertices for bounding box of $x$
c     lower left, upper right
      do 3 k=1,d
         alpha=machin
         beta=-machin
         do 4 i=1,n
            t=x(i,k)
            alpha=min(alpha,t)
            beta=max(beta,t)
    4    continue
c        expand the box a little
         mu=0.005D0*max(beta-alpha,1.d-10*max(DABS(alpha),DABS(beta))+1.
     +d-30)
         alpha=alpha-mu
         beta=beta+mu
         v(1,k)=alpha
         v(vc,k)=beta
    3 continue
c     remaining vertices
      do 5 i=2,vc-1
         j=i-1
         do 6 k=1,d
            v(i,k)=v(1+mod(j,2)*(vc-1),k)
            j=DFLOAT(j)/2.D0
    6    continue
    5 continue
      return
      end
C----------------------------------------------------------------------C
C     cpvert
      subroutine ehg125(p,nv,v,vhit,nvmax,d,k,t,r,s,f,l,u)
      logical i1,i2,match
      integer d,execnt,h,i,i3,j,k,m,mm,nv,nvmax,p,r,s
      integer f(r,0:1,s),l(r,0:1,s),u(r,0:1,s),vhit(nvmax)
      DOUBLE PRECISION t
      DOUBLE PRECISION v(nvmax,d)
      external ehg182
      save execnt
      data execnt /0/
      execnt=execnt+1
      h=nv
      do 3 i=1,r
         do 4 j=1,s
            h=h+1
            do 5 i3=1,d
               v(h,i3)=v(f(i,0,j),i3)
    5       continue
            v(h,k)=t
c           check for redundant vertex
            match=.false.
            m=1
c           top of while loop
    6       if(.not.match)then
               i1=(m.le.nv)
            else
               i1=.false.
            end if
            if(.not.(i1))goto 7
               match=(v(m,1).eq.v(h,1))
               mm=2
c              top of while loop
    8          if(match)then
                  i2=(mm.le.d)
               else
                  i2=.false.
               end if
               if(.not.(i2))goto 9
                  match=(v(m,mm).eq.v(h,mm))
                  mm=mm+1
                  goto 8
c              bottom of while loop
    9          m=m+1
               goto 6
c           bottom of while loop
    7       m=m-1
            if(match)then
               h=h-1
            else
               m=h
               if(vhit(1).ge.0)then
                  vhit(m)=p
               end if
            end if
            l(i,0,j)=f(i,0,j)
            l(i,1,j)=m
            u(i,0,j)=m
            u(i,1,j)=f(i,1,j)
    4    continue
    3 continue
      nv=h
      if(.not.(nv.le.nvmax))then
         call ehg182(180)
      end if
      return
      end
C-----------------------------------------------------------------------
C     descend
      integer function ehg138(i,z,a,xi,lo,hi,ncmax)
      logical i1
      integer d,execnt,i,j,nc,ncmax
      integer a(ncmax),hi(ncmax),lo(ncmax)
      DOUBLE PRECISION xi(ncmax),z(8)
      save execnt
      data execnt /0/
      execnt=execnt+1
c     descend tree until leaf or ambiguous
      j=i
c     top of while loop
    3 if(a(j).ne.0)then
         i1=(z(a(j)).ne.xi(j))
      else
         i1=.false.
      end if
      if(.not.(i1))goto 4
         if(z(a(j)).le.xi(j))then
            j=lo(j)
         else
            j=hi(j)
         end if
         goto 3
c     bottom of while loop
    4 ehg138=j
      return
      end

C----------------------------------------------------------------------C
C     select q-th smallest by partial sorting
      subroutine ehg106(il,ir,k,nk,p,pi,n)
      integer execnt,i,ii,il,ir,j,k,l,n,nk,r
      integer pi(n)
      DOUBLE PRECISION t
      DOUBLE PRECISION p(nk,n)
      save execnt
      data execnt /0/
      execnt=execnt+1
c     find the $k$-th smallest of $n$ elements
c     Floyd+Rivest, CACM Mar '75, Algorithm 489
      l=il
      r=ir
c     top of while loop
    3 if(.not.(l.lt.r))goto 4
c        to avoid recursion, sophisticated partition deleted
c        partition $x sub {l..r}$ about $t$
         t=p(1,pi(k))
         i=l
         j=r
         ii=pi(l)
         pi(l)=pi(k)
         pi(k)=ii
         if(t.lt.p(1,pi(r)))then
            ii=pi(l)
            pi(l)=pi(r)
            pi(r)=ii
         end if
c        top of while loop
    5    if(.not.(i.lt.j))goto 6
            ii=pi(i)
            pi(i)=pi(j)
            pi(j)=ii
            i=i+1
            j=j-1
c           top of while loop
    7       if(.not.(p(1,pi(i)).lt.t))goto 8
               i=i+1
               goto 7
c           bottom of while loop
    8       continue
c           top of while loop
    9       if(.not.(t.lt.p(1,pi(j))))goto 10
               j=j-1
               goto 9
c           bottom of while loop
   10       goto 5
c        bottom of while loop
    6    if(p(1,pi(l)).eq.t)then
            ii=pi(l)
            pi(l)=pi(j)
            pi(j)=ii
         else
            j=j+1
            ii=pi(r)
            pi(r)=pi(j)
            pi(j)=ii
         end if
         if(j.le.k)then
            l=j+1
         end if
         if(k.le.j)then
            r=j-1
         end if
         goto 3
c     bottom of while loop
    4 return
      end

C----------------------------------------------------------------------C
C     l2fit,l2tr computational kernel
      subroutine ehg127(q,n,d,nf,f,x,psi,y,rw,kernel,k,dist,eta,b,od,w,r
     +cond,sing,sigma,u,e,dgamma,qraux,work,tol,dd,tdeg,cdeg,s)
      integer column,d,dd,execnt,i,i3,i9,info,inorm2,j,jj,jpvt,k,kernel,
     +n,nf,od,sing,tdeg
      integer cdeg(8),psi(n)
      double precision machep,f,i1,i10,i2,i4,i5,i6,i7,i8,rcond,rho,scal,
     +tol
      double precision g(15),sigma(15),u(15,15),e(15,15),b(nf,k),colnor(
     +15),dist(n),eta(nf),dgamma(15),q(d),qraux(15),rw(n),s(0:od),w(nf),
     +work(15),x(n,d),y(n)
      external ehg106,ehg182,ehg184,dqrdc,dqrsl,dsvdc
      integer idamax
      external idamax
      double precision d1mach
      external d1mach
      double precision ddot
      external ddot
      save machep,execnt
      data execnt /0/
c     colnorm -> colnor
c     E -> g
c     MachEps -> machep
c     V -> e
c     X -> b
      execnt=execnt+1
      if(execnt.eq.1)then
         machep=d1mach(4)
      end if
c     sort by distance
      do 3 i3=1,n
         dist(i3)=0
    3 continue
      do 4 j=1,dd
         i4=q(j)
         do 5 i3=1,n
            dist(i3)=dist(i3)+(x(i3,j)-i4)**2
    5    continue
    4 continue
      call ehg106(1,n,nf,1,dist,psi,n)
      rho=dist(psi(nf))*max(1.d0,f)
      if(.not.(0.lt.rho))then
         call ehg182(120)
      end if
c     compute neighborhood weights
      if(kernel.eq.2)then
         do 6 i=1,nf
            if(dist(psi(i)).lt.rho)then
               i1=dsqrt(rw(psi(i)))
            else
               i1=0
            end if
            w(i)=i1
    6    continue
      else
         do 7 i3=1,nf
            w(i3)=dsqrt(dist(psi(i3))/rho)
    7    continue
         do 8 i3=1,nf
            w(i3)=dsqrt(rw(psi(i3))*(1-w(i3)**3)**3)
    8    continue
      end if
      if(dabs(w(idamax(nf,w,1))).eq.0)then
         call ehg184('at ',q,dd,1)
         call ehg184('radius ',rho,1,1)
         if(.not..false.)then
            call ehg182(121)
         end if
      end if
c     fill design matrix
      column=1
      do 9 i3=1,nf
         b(i3,column)=w(i3)
    9 continue
      if(tdeg.ge.1)then
         do 10 j=1,d
            if(cdeg(j).ge.1)then
               column=column+1
               i5=q(j)
               do 11 i3=1,nf
                  b(i3,column)=w(i3)*(x(psi(i3),j)-i5)
   11          continue
            end if
   10    continue
      end if
      if(tdeg.ge.2)then
         do 12 j=1,d
            if(cdeg(j).ge.1)then
               if(cdeg(j).ge.2)then
                  column=column+1
                  i6=q(j)
                  do 13 i3=1,nf
                     b(i3,column)=w(i3)*(x(psi(i3),j)-i6)**2
   13             continue
               end if
               do 14 jj=j+1,d
                  if(cdeg(jj).ge.1)then
                     column=column+1
                     i7=q(j)
                     i8=q(jj)
                     do 15 i3=1,nf
                        b(i3,column)=w(i3)*(x(psi(i3),j)-i7)*(x(psi(i3),
     +jj)-i8)
   15                continue
                  end if
   14          continue
            end if
   12    continue
         k=column
      end if
      do 16 i3=1,nf
         eta(i3)=w(i3)*y(psi(i3))
   16 continue
c     equilibrate columns
      do 17 j=1,k
         scal=0
         do 18 inorm2=1,nf
            scal=scal+b(inorm2,j)**2
   18    continue
         scal=dsqrt(scal)
         if(0.lt.scal)then
            do 19 i3=1,nf
               b(i3,j)=b(i3,j)/scal
   19       continue
            colnor(j)=scal
         else
            colnor(j)=1
         end if
   17 continue
c     singular value decomposition
      call dqrdc(b,nf,nf,k,qraux,jpvt,work,0)
      call dqrsl(b,nf,nf,k,qraux,eta,work,eta,eta,work,work,1000,info)
      do 20 i9=1,k
         do 21 i3=1,k
            u(i3,i9)=0
   21    continue
   20 continue
      do 22 i=1,k
         do 23 j=i,k
            u(i,j)=b(i,j)
   23    continue
   22 continue
      call dsvdc(u,15,k,k,sigma,g,u,15,e,15,work,21,info)
      if(.not.(info.eq.0))then
         call ehg182(182)
      end if
      tol=sigma(1)*(100*machep)
      rcond=min(rcond,sigma(k)/sigma(1))
      if(sigma(k).le.tol)then
         sing=sing+1
         if(sing.eq.1)then
            call ehg184('Warning. pseudoinverse used at',q,d,1)
            call ehg184('neighborhood radius',dsqrt(rho),1,1)
            call ehg184('reciprocal condition number ',rcond,1,1)
         else
            if(sing.eq.2)then
               call ehg184('There are other near singularities as well.'
     +,rho,1,1)
            end if
         end if
      end if
c     compensate for equilibration
      do 24 j=1,k
         i10=colnor(j)
         do 25 i3=1,k
            e(j,i3)=e(j,i3)/i10
   25    continue
   24 continue
c     solve least squares problem
      do 26 j=1,k
         if(tol.lt.sigma(j))then
            i2=ddot(k,u(1,j),1,eta,1)/sigma(j)
         else
            i2=0.d0
         end if
         dgamma(j)=i2
   26 continue
      do 27 j=0,od
c        bug fix 2006-07-04 for k=1, od>1.   (thanks btyner@gmail.com)
         if(j.lt.k)then
            s(j)=ddot(k,e(j+1,1),15,dgamma,1)
         else
            s(j)=0.d0
         end if
   27 continue
      return
      end

C----------------------------------------------------------------------C
C     lowesb after workspace expansion
      subroutine ehg131(x,y,rw,trl,diagl,kernel,k,n,d,nc,ncmax,vc,nv,nvm
     +ax,nf,f,a,c,hi,lo,pi,psi,v,vhit,vval,xi,dist,eta,b,ntol,fd,w,vval2
     +,rcond,sing,dd,tdeg,cdeg,lq,lf,setlf)
      logical setlf
      integer identi,d,dd,execnt,i1,i2,j,k,kernel,n,nc,ncmax,nf,ntol,nv,
     +nvmax,sing,tdeg,vc
      integer lq(nvmax,nf),a(ncmax),c(vc,ncmax),cdeg(8),hi(ncmax),lo(ncm
     +ax),pi(n),psi(n),vhit(nvmax)
      double precision f,fd,rcond,trl
      double precision lf(0:d,nvmax,nf),b(*),delta(8),diagl(n),dist(n),e
     +ta(nf),rw(n),v(nvmax,d),vval(0:d,nvmax),vval2(0:d,nvmax),w(nf),x(n
     +,d),xi(ncmax),y(n)
      external ehg126,ehg182,ehg139,ehg124
      double precision dnrm2
      external dnrm2
      save execnt
      data execnt /0/
c     Identity -> identi
c     X -> b
      execnt=execnt+1
      if(.not.(d.le.8))then
         call ehg182(101)
      end if
c     build $k$-d tree
      call ehg126(d,n,vc,x,v,nvmax)
      nv=vc
      nc=1
      do 3 j=1,vc
         c(j,nc)=j
         vhit(j)=0
    3 continue
      do 4 i1=1,d
         delta(i1)=v(vc,i1)-v(1,i1)
    4 continue
      fd=fd*dnrm2(d,delta,1)
      do 5 identi=1,n
         pi(identi)=identi
    5 continue
      call ehg124(1,n,d,n,nv,nc,ncmax,vc,x,pi,a,xi,lo,hi,c,v,vhit,nvmax,
     +ntol,fd,dd)
c     smooth
      if(trl.ne.0)then
         do 6 i2=1,nv
            do 7 i1=0,d
               vval2(i1,i2)=0
    7       continue
    6    continue
      end if
      call ehg139(v,nvmax,nv,n,d,nf,f,x,pi,psi,y,rw,trl,kernel,k,dist,di
     +st,eta,b,d,w,diagl,vval2,nc,vc,a,xi,lo,hi,c,vhit,rcond,sing,dd,tde
     +g,cdeg,lq,lf,setlf,vval)
      return
      end

C----------------------------------------------------------------------C
C     lowese after workspace expansion
      subroutine ehg133(n,d,vc,nvmax,nc,ncmax,a,c,hi,lo,v,vval,xi,m,z,s)
      integer d,execnt,i,i1,m,nc,ncmax,nv,nvmax,vc
      integer a(ncmax),c(vc,ncmax),hi(ncmax),lo(ncmax)
      double precision delta(8),s(m),v(nvmax,d),vval(0:d,nvmax),xi(ncmax
     +),z(m,d)
      double precision ehg128
      external ehg128
      save execnt
      data execnt /0/
      execnt=execnt+1
      do 3 i=1,m
         do 4 i1=1,d
            delta(i1)=z(i,i1)
    4    continue
         s(i)=ehg128(delta,d,ncmax,vc,a,xi,lo,hi,c,v,nvmax,vval)
    3 continue
      return
      end
      subroutine ehg140(iw,i,j)
      integer execnt,i,j
      integer iw(i)
      save execnt
      data execnt /0/
      execnt=execnt+1
      iw(i)=j
      return
      end

C----------------------------------------------------------------------C
C     delta1,2 from trL
      subroutine ehg141(trl,n,deg,k,d,nsing,dk,delta1,delta2)
      integer d,deg,dk,k,n,nsing
      external ehg176
      double precision ehg176
      double precision corx,delta1,delta2,trl,z
      double precision c(48), c1, c2, c3, c4
c     coef, d, deg, del
      data c / .2971620d0,.3802660d0,.5886043d0,.4263766d0,.3346498d0,.6
     +271053d0,.5241198d0,.3484836d0,.6687687d0,.6338795d0,.4076457d0,.7
     +207693d0,.1611761d0,.3091323d0,.4401023d0,.2939609d0,.3580278d0,.5
     +555741d0,.3972390d0,.4171278d0,.6293196d0,.4675173d0,.4699070d0,.6
     +674802d0,.2848308d0,.2254512d0,.2914126d0,.5393624d0,.2517230d0,.3
     +898970d0,.7603231d0,.2969113d0,.4740130d0,.9664956d0,.3629838d0,.5
     +348889d0,.2075670d0,.2822574d0,.2369957d0,.3911566d0,.2981154d0,.3
     +623232d0,.5508869d0,.3501989d0,.4371032d0,.7002667d0,.4291632d0,.4
     +930370d0 /
      if(deg.eq.0) dk=1
      if(deg.eq.1) dk=d+1
      if(deg.eq.2) dk=dfloat((d+2)*(d+1))/2.d0
      corx=dsqrt(k/dfloat(n))
      z=(dsqrt(k/trl)-corx)/(1-corx)
      if(nsing .eq. 0 .and. 1 .lt. z)   call ehg184('Chernobyl! trL<k',t
     +rl,1,1)
      if(z .lt. 0) call ehg184('Chernobyl! trL>n',trl,1,1)
      z=min(1.0d0,max(0.0d0,z))
      c4=dexp(ehg176(z))
      i=1+3*(min(d,4)-1+4*(deg-1))
      if(d.le.4)then
         c1=c(i)
         c2=c(i+1)
         c3=c(i+2)
      else
         c1=c(i)+(d-4)*(c(i)-c(i-3))
         c2=c(i+1)+(d-4)*(c(i+1)-c(i-2))
         c3=c(i+2)+(d-4)*(c(i+2)-c(i-1))
      endif
      delta1=n-trl*dexp(c1*z**c2*(1-z)**c3*c4)
      i=i+24
      if(d.le.4)then
         c1=c(i)
         c2=c(i+1)
         c3=c(i+2)
      else
         c1=c(i)+(d-4)*(c(i)-c(i-3))
         c2=c(i+1)+(d-4)*(c(i+1)-c(i-2))
         c3=c(i+2)+(d-4)*(c(i+2)-c(i-1))
      endif
      delta2=n-trl*dexp(c1*z**c2*(1-z)**c3*c4)
      return
      end

C----------------------------------------------------------------------C
C     exact delta
      subroutine lowesc(n,l,ll,trl,delta1,delta2)
      integer execnt,i,j,n
      double precision delta1,delta2,trl
      double precision l(n,n),ll(n,n)
      double precision ddot
      external ddot
      save execnt
      data execnt /0/
      execnt=execnt+1
c     compute $LL~=~(I-L)(I-L)'$
      do 3 i=1,n
         l(i,i)=l(i,i)-1
    3 continue
      do 4 i=1,n
         do 5 j=1,i
            ll(i,j)=ddot(n,l(i,1),n,l(j,1),n)
    5    continue
    4 continue
      do 6 i=1,n
         do 7 j=i+1,n
            ll(i,j)=ll(j,i)
    7    continue
    6 continue
      do 8 i=1,n
         l(i,i)=l(i,i)+1
    8 continue
c     accumulate first two traces
      trl=0
      delta1=0
      do 9 i=1,n
         trl=trl+l(i,i)
         delta1=delta1+ll(i,i)
    9 continue
c     $delta sub 2 = "tr" LL sup 2$
      delta2=0
      do 10 i=1,n
         delta2=delta2+ddot(n,ll(i,1),n,ll(1,i),1)
   10 continue
      return
      end

C----------------------------------------------------------------------C
C     compute derived k-d tree information
      subroutine ehg169(d,vc,nc,ncmax,nv,nvmax,v,a,xi,c,hi,lo)
      integer d,execnt,i,j,k,mc,mv,nc,ncmax,nv,nvmax,p,vc
      integer a(ncmax),c(vc,ncmax),hi(ncmax),lo(ncmax),novhit(1)
      DOUBLE PRECISION v(nvmax,d),xi(ncmax)
      external ehg125,ehg182
      integer ifloor
      external ifloor
      save execnt
      data execnt /0/
      execnt=execnt+1
c     as in bbox
c     remaining vertices
      do 3 i=2,vc-1
         j=i-1
         do 4 k=1,d
            v(i,k)=v(1+mod(j,2)*(vc-1),k)
            j=ifloor(DFLOAT(j)/2.D0)
    4    continue
    3 continue
c     as in ehg131
      mc=1
      mv=vc
      novhit(1)=-1
      do 5 j=1,vc
         c(j,mc)=j
    5 continue
c     as in rbuild
      p=1
c     top of while loop
    6 if(.not.(p.le.nc))goto 7
         if(a(p).ne.0)then
            k=a(p)
c           left son
            mc=mc+1
            lo(p)=mc
c           right son
            mc=mc+1
            hi(p)=mc
            call ehg125(p,mv,v,novhit,nvmax,d,k,xi(p),2**(k-1),2**(d-k),
     +c(1,p),c(1,lo(p)),c(1,hi(p)))
         end if
         p=p+1
         goto 6
c     bottom of while loop
    7 if(.not.(mc.eq.nc))then
         call ehg182(193)
      end if
      if(.not.(mv.eq.nv))then
         call ehg182(193)
      end if
      return
      end

C----------------------------------------------------------------------C
C     loeval for delta
       DOUBLE PRECISION function ehg176(z)
       DOUBLE PRECISION z(*)
       integer d,vc,nv,nc
       integer a(17), c(2,17)
       integer hi(17), lo(17)
       DOUBLE PRECISION v(10,1)
       DOUBLE PRECISION vval(0:1,10)
       DOUBLE PRECISION xi(17)
       DOUBLE PRECISION ehg128
       data d,vc,nv,nc /1,2,10,17/
       data a(1) /1/
       data hi(1),lo(1),xi(1) /3,2,0.3705D0/
       data c(1,1) /1/
       data c(2,1) /2/
       data a(2) /1/
       data hi(2),lo(2),xi(2) /5,4,0.2017D0/
       data c(1,2) /1/
       data c(2,2) /3/
       data a(3) /1/
       data hi(3),lo(3),xi(3) /7,6,0.5591D0/
       data c(1,3) /3/
       data c(2,3) /2/
       data a(4) /1/
       data hi(4),lo(4),xi(4) /9,8,0.1204D0/
       data c(1,4) /1/
       data c(2,4) /4/
       data a(5) /1/
       data hi(5),lo(5),xi(5) /11,10,0.2815D0/
       data c(1,5) /4/
       data c(2,5) /3/
       data a(6) /1/
       data hi(6),lo(6),xi(6) /13,12,0.4536D0/
       data c(1,6) /3/
       data c(2,6) /5/
       data a(7) /1/
       data hi(7),lo(7),xi(7) /15,14,0.7132D0/
       data c(1,7) /5/
       data c(2,7) /2/
       data a(8) /0/
       data c(1,8) /1/
       data c(2,8) /6/
       data a(9) /0/
       data c(1,9) /6/
       data c(2,9) /4/
       data a(10) /0/
       data c(1,10) /4/
       data c(2,10) /7/
       data a(11) /0/
       data c(1,11) /7/
       data c(2,11) /3/
       data a(12) /0/
       data c(1,12) /3/
       data c(2,12) /8/
       data a(13) /0/
       data c(1,13) /8/
       data c(2,13) /5/
       data a(14) /0/
       data c(1,14) /5/
       data c(2,14) /9/
       data a(15) /1/
       data hi(15),lo(15),xi(15) /17,16,0.8751D0/
       data c(1,15) /9/
       data c(2,15) /2/
       data a(16) /0/
       data c(1,16) /9/
       data c(2,16) /10/
       data a(17) /0/
       data c(1,17) /10/
       data c(2,17) /2/
       data vval(0,1) /-9.0572D-2/
       data v(1,1) /-5.D-3/
       data vval(1,1) /4.4844D0/
       data vval(0,2) /-1.0856D-2/
       data v(2,1) /1.005D0/
       data vval(1,2) /-0.7736D0/
       data vval(0,3) /-5.3718D-2/
       data v(3,1) /0.3705D0/
       data vval(1,3) /-0.3495D0/
       data vval(0,4) /2.6152D-2/
       data v(4,1) /0.2017D0/
       data vval(1,4) /-0.7286D0/
       data vval(0,5) /-5.8387D-2/
       data v(5,1) /0.5591D0/
       data vval(1,5) /0.1611D0/
       data vval(0,6) /9.5807D-2/
       data v(6,1) /0.1204D0/
       data vval(1,6) /-0.7978D0/
       data vval(0,7) /-3.1926D-2/
       data v(7,1) /0.2815D0/
       data vval(1,7) /-0.4457D0/
       data vval(0,8) /-6.4170D-2/
       data v(8,1) /0.4536D0/
       data vval(1,8) /3.2813D-2/
       data vval(0,9) /-2.0636D-2/
       data v(9,1) /0.7132D0/
       data vval(1,9) /0.3350D0/
       data vval(0,10) /4.0172D-2/
       data v(10,1) /0.8751D0/
       data vval(1,10) /-4.1032D-2/
       ehg176=ehg128(z,d,nc,vc,a,xi,lo,hi,c,v,nv,vval)
       end

C----------------------------------------------------------------------C
C     approximate delta
      subroutine lowesa(trl,n,d,tau,nsing,delta1,delta2)
      integer d,dka,dkb,execnt,n,nsing,tau
      double precision alpha,d1a,d1b,d2a,d2b,delta1,delta2,trl
      external ehg141
      save execnt
      data execnt /0/
      execnt=execnt+1
      call ehg141(trl,n,1,tau,d,nsing,dka,d1a,d2a)
      call ehg141(trl,n,2,tau,d,nsing,dkb,d1b,d2b)
      alpha=dfloat(tau-dka)/dfloat(dkb-dka)
      delta1=(1-alpha)*d1a+alpha*d1b
      delta2=(1-alpha)*d2a+alpha*d2b
      return
      end

C----------------------------------------------------------------------C
C     lowesl after workspace expansion
      subroutine ehg191(m,z,l,d,n,nf,nv,ncmax,vc,a,xi,lo,hi,c,v,nvmax,vv
     +al2,lf,lq)
      integer lq1,d,execnt,i,i1,i2,j,m,n,nc,ncmax,nf,nv,nvmax,p,vc
      integer lq(nvmax,nf),a(ncmax),c(vc,ncmax),hi(ncmax),lo(ncmax)
      double precision l(m,n),lf(0:d,nvmax,nf),v(nvmax,d),vval2(0:d,nvma
     +x),xi(ncmax),z(m,d),zi(8)
      double precision ehg128
      external ehg128
      save execnt
      data execnt /0/
      execnt=execnt+1
      do 3 j=1,n
         do 4 i2=1,nv
            do 5 i1=0,d
               vval2(i1,i2)=0
    5       continue
    4    continue
         do 6 i=1,nv
c           linear search for i in Lq
            lq1=lq(i,1)
            lq(i,1)=j
            p=nf
c           top of while loop
    7       if(.not.(lq(i,p).ne.j))goto 8
               p=p-1
               goto 7
c           bottom of while loop
    8       lq(i,1)=lq1
            if(lq(i,p).eq.j)then
               do 9 i1=0,d
                  vval2(i1,i)=lf(i1,i,p)
    9          continue
            end if
    6    continue
         do 10 i=1,m
            do 11 i1=1,d
               zi(i1)=z(i,i1)
   11       continue
            l(i,j)=ehg128(zi,d,ncmax,vc,a,xi,lo,hi,c,v,nvmax,vval2)
   10    continue
    3 continue
      return
      end

C----------------------------------------------------------------------C
C     trL approximation
      subroutine ehg196(tau,d,f,trl)
      integer d,dka,dkb,execnt,tau
      double precision alpha,f,trl,trla,trlb
      external ehg197
      save execnt
      data execnt /0/
      execnt=execnt+1
      call ehg197(1,tau,d,f,dka,trla)
      call ehg197(2,tau,d,f,dkb,trlb)
      alpha=dfloat(tau-dka)/dfloat(dkb-dka)
      trl=(1-alpha)*trla+alpha*trlb
      return
      end

C----------------------------------------------------------------------C
C     for deg 1,2
      subroutine ehg197(deg,tau,d,f,dk,trl)
      integer d,deg,dk,tau
      double precision trl, f
      dk = 0
      if(deg.eq.1) dk=d+1
      if(deg.eq.2) dk=dfloat((d+2)*(d+1))/2.d0
      g1 = (-0.08125d0*d+0.13d0)*d+1.05d0
      trl = dk*(1+max(0.d0,(g1-f)/f))
      return
      end

C----------------------------------------------------------------------C
C     lowesr after workspace expansion
      subroutine ehg192(y,d,n,nf,nv,nvmax,vval,lf,lq)
      integer d,execnt,i,i1,i2,j,n,nf,nv,nvmax
      integer lq(nvmax,nf)
      DOUBLE PRECISION i3
      DOUBLE PRECISION lf(0:d,nvmax,nf),vval(0:d,nvmax),y(n)
      save execnt
      data execnt /0/
      execnt=execnt+1
      do 3 i2=1,nv
         do 4 i1=0,d
            vval(i1,i2)=0
    4    continue
    3 continue
      do 5 i=1,nv
         do 6 j=1,nf
            i3=y(lq(i,j))
            do 7 i1=0,d
               vval(i1,i)=vval(i1,i)+i3*lf(i1,i,j)
    7       continue
    6    continue
    5 continue
      return
      end


C----------------------------------------------------------------------C
C     eval
      DOUBLE PRECISION function ehg128(z,d,ncmax,vc,a,xi,lo,hi,c,v,nvmax
     +,vval)
      logical i10,i2,i3,i4,i5,i6,i7,i8,i9
      integer d,execnt,i,i1,i11,i12,ig,ii,j,lg,ll,m,nc,ncmax,nt,nv,nvmax
     +,ur,vc
      integer a(ncmax),c(vc,ncmax),hi(ncmax),lo(ncmax),t(20)
      DOUBLE PRECISION ge,gn,gs,gw,gpe,gpn,gps,gpw,h,phi0,phi1,psi0,psi1
     +,s,sew,sns,v0,v1,xibar
      DOUBLE PRECISION g(0:8,256),g0(0:8),g1(0:8),v(nvmax,d),vval(0:d,nv
     +max),xi(ncmax),z(d)
      external ehg182,ehg184
      save execnt
      data execnt /0/
      execnt=execnt+1
c     locate enclosing cell
      nt=1
      t(nt)=1
      j=1
c     top of while loop
    3 if(.not.(a(j).ne.0))goto 4
         nt=nt+1
c     bug fix 2006-07-18 (thanks, btyner@gmail.com)
         if(z(a(j)).le.xi(j))then
            i1=lo(j)
         else
            i1=hi(j)
         end if
         t(nt)=i1
         if(.not.(nt.lt.20))then
            call ehg182(181)
         end if
         j=t(nt)
         goto 3
c     bottom of while loop
    4 continue
c     tensor
      do 5 i12=1,vc
         do 6 i11=0,d
            g(i11,i12)=vval(i11,c(i12,j))
    6    continue
    5 continue
      lg=vc
      ll=c(1,j)
      ur=c(vc,j)
      do 7 i=d,1,-1
         h=(z(i)-v(ll,i))/(v(ur,i)-v(ll,i))
         if(h.lt.-.001D0)then
            call ehg184('eval ',z,d,1)
            call ehg184('lowerlimit ',v(ll,1),d,nvmax)
         else
            if(1.001D0.lt.h)then
               call ehg184('eval ',z,d,1)
               call ehg184('upperlimit ',v(ur,1),d,nvmax)
            end if
         end if
         if(-.001D0.le.h)then
            i2=(h.le.1.001D0)
         else
            i2=.false.
         end if
         if(.not.i2)then
            call ehg182(122)
         end if
         lg=DFLOAT(lg)/2.D0
         do 8 ig=1,lg
c           Hermite basis
            phi0=(1-h)**2*(1+2*h)
            phi1=h**2*(3-2*h)
            psi0=h*(1-h)**2
            psi1=h**2*(h-1)
            g(0,ig)=phi0*g(0,ig)+phi1*g(0,ig+lg)+(psi0*g(i,ig)+psi1*g(i,
     +ig+lg))*(v(ur,i)-v(ll,i))
            do 9 ii=1,i-1
               g(ii,ig)=phi0*g(ii,ig)+phi1*g(ii,ig+lg)
    9       continue
    8    continue
    7 continue
      s=g(0,1)
c     blending
      if(d.eq.2)then
c        ----- North -----
         v0=v(ll,1)
         v1=v(ur,1)
         do 10 i11=0,d
            g0(i11)=vval(i11,c(3,j))
   10    continue
         do 11 i11=0,d
            g1(i11)=vval(i11,c(4,j))
   11    continue
         xibar=v(ur,2)
         m=nt-1
c        top of while loop
   12    if(m.eq.0)then
            i4=.true.
         else
            if(a(t(m)).eq.2)then
               i3=(xi(t(m)).eq.xibar)
            else
               i3=.false.
            end if
            i4=i3
         end if
         if(.not.(.not.i4))goto 13
            m=m-1
c           voidp junk
            goto 12
c        bottom of while loop
   13    if(m.ge.1)then
            m=hi(t(m))
c           top of while loop
   14       if(.not.(a(m).ne.0))goto 15
               if(z(a(m)).le.xi(m))then
                  m=lo(m)
               else
                  m=hi(m)
               end if
               goto 14
c           bottom of while loop
   15       if(v0.lt.v(c(1,m),1))then
               v0=v(c(1,m),1)
               do 16 i11=0,d
                  g0(i11)=vval(i11,c(1,m))
   16          continue
            end if
            if(v(c(2,m),1).lt.v1)then
               v1=v(c(2,m),1)
               do 17 i11=0,d
                  g1(i11)=vval(i11,c(2,m))
   17          continue
            end if
         end if
         h=(z(1)-v0)/(v1-v0)
c        Hermite basis
         phi0=(1-h)**2*(1+2*h)
         phi1=h**2*(3-2*h)
         psi0=h*(1-h)**2
         psi1=h**2*(h-1)
         gn=phi0*g0(0)+phi1*g1(0)+(psi0*g0(1)+psi1*g1(1))*(v1-v0)
         gpn=phi0*g0(2)+phi1*g1(2)
c        ----- South -----
         v0=v(ll,1)
         v1=v(ur,1)
         do 18 i11=0,d
            g0(i11)=vval(i11,c(1,j))
   18    continue
         do 19 i11=0,d
            g1(i11)=vval(i11,c(2,j))
   19    continue
         xibar=v(ll,2)
         m=nt-1
c        top of while loop
   20    if(m.eq.0)then
            i6=.true.
         else
            if(a(t(m)).eq.2)then
               i5=(xi(t(m)).eq.xibar)
            else
               i5=.false.
            end if
            i6=i5
         end if
         if(.not.(.not.i6))goto 21
            m=m-1
c           voidp junk
            goto 20
c        bottom of while loop
   21    if(m.ge.1)then
            m=lo(t(m))
c           top of while loop
   22       if(.not.(a(m).ne.0))goto 23
               if(z(a(m)).le.xi(m))then
                  m=lo(m)
               else
                  m=hi(m)
               end if
               goto 22
c           bottom of while loop
   23       if(v0.lt.v(c(3,m),1))then
               v0=v(c(3,m),1)
               do 24 i11=0,d
                  g0(i11)=vval(i11,c(3,m))
   24          continue
            end if
            if(v(c(4,m),1).lt.v1)then
               v1=v(c(4,m),1)
               do 25 i11=0,d
                  g1(i11)=vval(i11,c(4,m))
   25          continue
            end if
         end if
         h=(z(1)-v0)/(v1-v0)
c        Hermite basis
         phi0=(1-h)**2*(1+2*h)
         phi1=h**2*(3-2*h)
         psi0=h*(1-h)**2
         psi1=h**2*(h-1)
         gs=phi0*g0(0)+phi1*g1(0)+(psi0*g0(1)+psi1*g1(1))*(v1-v0)
         gps=phi0*g0(2)+phi1*g1(2)
c        ----- East -----
         v0=v(ll,2)
         v1=v(ur,2)
         do 26 i11=0,d
            g0(i11)=vval(i11,c(2,j))
   26    continue
         do 27 i11=0,d
            g1(i11)=vval(i11,c(4,j))
   27    continue
         xibar=v(ur,1)
         m=nt-1
c        top of while loop
   28    if(m.eq.0)then
            i8=.true.
         else
            if(a(t(m)).eq.1)then
               i7=(xi(t(m)).eq.xibar)
            else
               i7=.false.
            end if
            i8=i7
         end if
         if(.not.(.not.i8))goto 29
            m=m-1
c           voidp junk
            goto 28
c        bottom of while loop
   29    if(m.ge.1)then
            m=hi(t(m))
c           top of while loop
   30       if(.not.(a(m).ne.0))goto 31
               if(z(a(m)).le.xi(m))then
                  m=lo(m)
               else
                  m=hi(m)
               end if
               goto 30
c           bottom of while loop
   31       if(v0.lt.v(c(1,m),2))then
               v0=v(c(1,m),2)
               do 32 i11=0,d
                  g0(i11)=vval(i11,c(1,m))
   32          continue
            end if
            if(v(c(3,m),2).lt.v1)then
               v1=v(c(3,m),2)
               do 33 i11=0,d
                  g1(i11)=vval(i11,c(3,m))
   33          continue
            end if
         end if
         h=(z(2)-v0)/(v1-v0)
c        Hermite basis
         phi0=(1-h)**2*(1+2*h)
         phi1=h**2*(3-2*h)
         psi0=h*(1-h)**2
         psi1=h**2*(h-1)
         ge=phi0*g0(0)+phi1*g1(0)+(psi0*g0(2)+psi1*g1(2))*(v1-v0)
         gpe=phi0*g0(1)+phi1*g1(1)
c        ----- West -----
         v0=v(ll,2)
         v1=v(ur,2)
         do 34 i11=0,d
            g0(i11)=vval(i11,c(1,j))
   34    continue
         do 35 i11=0,d
            g1(i11)=vval(i11,c(3,j))
   35    continue
         xibar=v(ll,1)
         m=nt-1
c        top of while loop
   36    if(m.eq.0)then
            i10=.true.
         else
            if(a(t(m)).eq.1)then
               i9=(xi(t(m)).eq.xibar)
            else
               i9=.false.
            end if
            i10=i9
         end if
         if(.not.(.not.i10))goto 37
            m=m-1
c           voidp junk
            goto 36
c        bottom of while loop
   37    if(m.ge.1)then
            m=lo(t(m))
c           top of while loop
   38       if(.not.(a(m).ne.0))goto 39
               if(z(a(m)).le.xi(m))then
                  m=lo(m)
               else
                  m=hi(m)
               end if
               goto 38
c           bottom of while loop
   39       if(v0.lt.v(c(2,m),2))then
               v0=v(c(2,m),2)
               do 40 i11=0,d
                  g0(i11)=vval(i11,c(2,m))
   40          continue
            end if
            if(v(c(4,m),2).lt.v1)then
               v1=v(c(4,m),2)
               do 41 i11=0,d
                  g1(i11)=vval(i11,c(4,m))
   41          continue
            end if
         end if
         h=(z(2)-v0)/(v1-v0)
c        Hermite basis
         phi0=(1-h)**2*(1+2*h)
         phi1=h**2*(3-2*h)
         psi0=h*(1-h)**2
         psi1=h**2*(h-1)
         gw=phi0*g0(0)+phi1*g1(0)+(psi0*g0(2)+psi1*g1(2))*(v1-v0)
         gpw=phi0*g0(1)+phi1*g1(1)
c        NS
         h=(z(2)-v(ll,2))/(v(ur,2)-v(ll,2))
c        Hermite basis
         phi0=(1-h)**2*(1+2*h)
         phi1=h**2*(3-2*h)
         psi0=h*(1-h)**2
         psi1=h**2*(h-1)
         sns=phi0*gs+phi1*gn+(psi0*gps+psi1*gpn)*(v(ur,2)-v(ll,2))
c        EW
         h=(z(1)-v(ll,1))/(v(ur,1)-v(ll,1))
c        Hermite basis
         phi0=(1-h)**2*(1+2*h)
         phi1=h**2*(3-2*h)
         psi0=h*(1-h)**2
         psi1=h**2*(h-1)
         sew=phi0*gw+phi1*ge+(psi0*gpw+psi1*gpe)*(v(ur,1)-v(ll,1))
         s=(sns+sew)-s
      end if
      ehg128=s
      return
      end


C----------------------------------------------------------------------C
C
      integer function ifloor(x)
      DOUBLE PRECISION x
      ifloor=x
      if(ifloor.gt.x) ifloor=ifloor-1
      end
      DOUBLE PRECISION functionDSIGN(a1,a2)
      DOUBLE PRECISION a1, a2
      DSIGN=DABS(a1)
      if(a2.ge.0)DSIGN=-DSIGN
      end
      subroutine ehg136(u,lm,m,n,d,nf,f,x,psi,y,rw,kernel,k,dist,eta,b,o
     +d,o,ihat,w,rcond,sing,dd,tdeg,cdeg,s)
      integer identi,d,dd,execnt,i,i1,ihat,info,j,k,kernel,l,lm,m,n,nf,o
     +d,sing,tdeg
      integer cdeg(8),psi(n)
      double precision f,i2,rcond,scale,tol
      double precision o(m,n),sigma(15),e(15,15),g(15,15),b(nf,k),dist(n
     +),eta(nf),dgamma(15),q(8),qraux(15),rw(n),s(0:od,m),u(lm,d),w(nf),
     +work(15),x(n,d),y(n)
      external ehg127,ehg182,dqrsl
      double precision ddot
      external ddot
      save execnt
      data execnt /0/
c     V -> g
c     U -> e
c     Identity -> identi
c     L -> o
c     X -> b
      execnt=execnt+1
      if(.not.(k.le.nf-1))then
         call ehg182(104)
      end if
      if(.not.(k.le.15))then
         call ehg182(105)
      end if
      do 3 identi=1,n
         psi(identi)=identi
    3 continue
      do 4 l=1,m
         do 5 i1=1,d
            q(i1)=u(l,i1)
    5    continue
         call ehg127(q,n,d,nf,f,x,psi,y,rw,kernel,k,dist,eta,b,od,w,rcon
     +d,sing,sigma,e,g,dgamma,qraux,work,tol,dd,tdeg,cdeg,s(0,l))
         if(ihat.eq.1)then
c           $L sub {l,l} =
c           V sub {1,:} SIGMA sup {+} U sup T
c           (Q sup T W e sub i )$
            if(.not.(m.eq.n))then
               call ehg182(123)
            end if
c           find $i$ such that $l = psi sub i$
            i=1
c           top of while loop
    6       if(.not.(l.ne.psi(i)))goto 7
               i=i+1
               if(.not.(i.lt.nf))then
                  call ehg182(123)
               end if
               goto 6
c           bottom of while loop
    7       do 8 i1=1,nf
               eta(i1)=0
    8       continue
            eta(i)=w(i)
c           $eta = Q sup T W e sub i$
            call dqrsl(b,nf,nf,k,qraux,eta,eta,eta,eta,eta,eta,1000,info
     +)
c           $gamma = U sup T eta sub {1:k}$
            do 9 i1=1,k
               dgamma(i1)=0
    9       continue
            do 10 j=1,k
               i2=eta(j)
               do 11 i1=1,k
                  dgamma(i1)=dgamma(i1)+i2*e(j,i1)
   11          continue
   10       continue
c           $gamma = SIGMA sup {+} gamma$
            do 12 j=1,k
               if(tol.lt.sigma(j))then
                  dgamma(j)=dgamma(j)/sigma(j)
               else
                  dgamma(j)=0.d0
               end if
   12       continue
c           voidp junk
c           voidp junk
            o(l,1)=ddot(k,g(1,1),15,dgamma,1)
         else
            if(ihat.eq.2)then
c              $L sub {l,:} =
c              V sub {1,:} SIGMA sup {+}
c              ( U sup T Q sup T ) W $
               do 13 i1=1,n
                  o(l,i1)=0
   13          continue
               do 14 j=1,k
                  do 15 i1=1,nf
                     eta(i1)=0
   15             continue
                  do 16 i1=1,k
                     eta(i1)=e(i1,j)
   16             continue
                  call dqrsl(b,nf,nf,k,qraux,eta,eta,work,work,work,work
     +,10000,info)
                  if(tol.lt.sigma(j))then
                     scale=1.d0/sigma(j)
                  else
                     scale=0.d0
                  end if
                  do 17 i1=1,nf
                     eta(i1)=eta(i1)*(scale*w(i1))
   17             continue
                  do 18 i=1,nf
                     o(l,psi(i))=o(l,psi(i))+g(1,j)*eta(i)
   18             continue
   14          continue
            end if
         end if
    4 continue
      return
      end

C----------------------------------------------------------------------C
C     l2tr
      subroutine ehg139(v,nvmax,nv,n,d,nf,f,x,pi,psi,y,rw,trl,kernel,k,d
     +ist,phi,eta,b,od,w,diagl,vval2,ncmax,vc,a,xi,lo,hi,c,vhit,rcond,si
     +ng,dd,tdeg,cdeg,lq,lf,setlf,s)
      logical setlf
      integer identi,d,dd,execnt,i,i2,i3,i5,i6,ii,ileaf,info,j,k,kernel,
     +l,n,nc,ncmax,nf,nleaf,nv,nvmax,od,sing,tdeg,vc
      integer lq(nvmax,nf),a(ncmax),c(vc,ncmax),cdeg(8),hi(ncmax),leaf(2
     +56),lo(ncmax),phi(n),pi(n),psi(n),vhit(nvmax)
      DOUBLE PRECISION f,i1,i4,i7,rcond,scale,term,tol,trl
      DOUBLE PRECISION lf(0:d,nvmax,nf),sigma(15),u(15,15),e(15,15),b(nf
     +,k),diagl(n),dist(n),eta(nf),DGAMMA(15),q(8),qraux(15),rw(n),s(0:o
     +d,nv),v(nvmax,d),vval2(0:d,nv),w(nf),work(15),x(n,d),xi(ncmax),y(n
     +),z(8)
      external ehg127,ehg182,DQRSL,ehg137
      DOUBLE PRECISION ehg128
      external ehg128
      DOUBLE PRECISION DDOT
      external DDOT
      save execnt
      data execnt /0/
c     V -> e
c     Identity -> identi
c     X -> b
      execnt=execnt+1
c     l2fit with trace(L)
      if(.not.(k.le.nf-1))then
         call ehg182(104)
      end if
      if(.not.(k.le.15))then
         call ehg182(105)
      end if
      if(trl.ne.0)then
         do 3 i5=1,n
            diagl(i5)=0
    3    continue
         do 4 i6=1,nv
            do 5 i5=0,d
               vval2(i5,i6)=0
    5       continue
    4    continue
      end if
      do 6 identi=1,n
         psi(identi)=identi
    6 continue
      do 7 l=1,nv
         do 8 i5=1,d
            q(i5)=v(l,i5)
    8    continue
         call ehg127(q,n,d,nf,f,x,psi,y,rw,kernel,k,dist,eta,b,od,w,rcon
     +d,sing,sigma,u,e,DGAMMA,qraux,work,tol,dd,tdeg,cdeg,s(0,l))
         if(trl.ne.0)then
c           invert $psi$
            do 9 i5=1,n
               phi(i5)=0
    9       continue
            do 10 i=1,nf
               phi(psi(i))=i
   10       continue
            do 11 i5=1,d
               z(i5)=v(l,i5)
   11       continue
            call ehg137(z,vhit(l),leaf,nleaf,d,nv,nvmax,ncmax,vc,a,xi,lo
     +,hi,c,v)
            do 12 ileaf=1,nleaf
               do 13 ii=lo(leaf(ileaf)),hi(leaf(ileaf))
                  i=phi(pi(ii))
                  if(i.ne.0)then
                     if(.not.(psi(i).eq.pi(ii)))then
                        call ehg182(194)
                     end if
                     do 14 i5=1,nf
                        eta(i5)=0
   14                continue
                     eta(i)=w(i)
c                    $eta = Q sup T W e sub i$
                     call DQRSL(b,nf,nf,k,qraux,eta,work,eta,eta,work,wo
     +rk,1000,info)
                     do 15 j=1,k
                        if(tol.lt.sigma(j))then
                           i4=DDOT(k,u(1,j),1,eta,1)/sigma(j)
                        else
                           i4=0.D0
                        end if
                       DGAMMA(j)=i4
   15                continue
                     do 16 j=1,d+1
c                       bug fix 2006-07-15 for k=1, od>1.   (thanks btyner@gmail.com)
                        if(j.le.k)then
                           vval2(j-1,l)=DDOT(k,e(j,1),15,DGAMMA,1)
                        else
                           vval2(j-1,l)=0
                        end if
   16                continue
                     do 17 i5=1,d
                        z(i5)=x(pi(ii),i5)
   17                continue
                     term=ehg128(z,d,ncmax,vc,a,xi,lo,hi,c,v,nvmax,vval2
     +)
                     diagl(pi(ii))=diagl(pi(ii))+term
                     do 18 i5=0,d
                        vval2(i5,l)=0
   18                continue
                  end if
   13          continue
   12       continue
         end if
         if(setlf)then
c           $Lf sub {:,l,:} = V SIGMA sup {+} U sup T Q sup T W$
            if(.not.(k.ge.d+1))then
               call ehg182(196)
            end if
            do 19 i5=1,nf
               lq(l,i5)=psi(i5)
   19       continue
            do 20 i6=1,nf
               do 21 i5=0,d
                  lf(i5,l,i6)=0
   21          continue
   20       continue
            do 22 j=1,k
               do 23 i5=1,nf
                  eta(i5)=0
   23          continue
               do 24 i5=1,k
                  eta(i5)=u(i5,j)
   24          continue
               call DQRSL(b,nf,nf,k,qraux,eta,eta,work,work,work,work,10
     +000,info)
               if(tol.lt.sigma(j))then
                  scale=1.D0/sigma(j)
               else
                  scale=0.D0
               end if
               do 25 i5=1,nf
                  eta(i5)=eta(i5)*(scale*w(i5))
   25          continue
               do 26 i=1,nf
                  i7=eta(i)
                  do 27 i5=0,d
                     if(i5.lt.k)then
                        lf(i5,l,i)=lf(i5,l,i)+e(1+i5,j)*i7
                     else
                        lf(i5,l,i)=0
                     end if
   27             continue
   26          continue
   22       continue
         end if
    7 continue
      if(trl.ne.0)then
         if(n.le.0)then
            trl=0.D0
         else
            i3=n
            i1=diagl(i3)
            do 28 i2=i3-1,1,-1
               i1=diagl(i2)+i1
   28       continue
            trl=i1
         end if
      end if
      return
      end


C----------------------------------------------------------------------C
      subroutine dqrdc(x,ldx,n,p,qraux,jpvt,work,job)
      integer ldx,n,p,job
      integer jpvt(1)
      double precision x(ldx,1),qraux(1),work(1)
c
c     dqrdc uses householder transformations to compute the qr
c     factorization of an n by p matrix x.  column pivoting
c     based on the 2-norms of the reduced columns may be
c     performed at the users option.
c
c     on entry
c
c        x       double precision(ldx,p), where ldx .ge. n.
c                x contains the matrix whose decomposition is to be
c                computed.
c
c        ldx     integer.
c                ldx is the leading dimension of the array x.
c
c        n       integer.
c                n is the number of rows of the matrix x.
c
c        p       integer.
c                p is the number of columns of the matrix x.
c
c        jpvt    integer(p).
c                jpvt contains integers that control the selection
c                of the pivot columns.  the k-th column x(k) of x
c                is placed in one of three classes according to the
c                value of jpvt(k).
c
c                   if jpvt(k) .gt. 0, then x(k) is an initial
c                                      column.
c
c                   if jpvt(k) .eq. 0, then x(k) is a free column.
c
c                   if jpvt(k) .lt. 0, then x(k) is a final column.
c
c                before the decomposition is computed, initial columns
c                are moved to the beginning of the array x and final
c                columns to the end.  both initial and final columns
c                are frozen in place during the computation and only
c                free columns are moved.  at the k-th stage of the
c                reduction, if x(k) is occupied by a free column
c                it is interchanged with the free column of largest
c                reduced norm.  jpvt is not referenced if
c                job .eq. 0.
c
c        work    double precision(p).
c                work is a work array.  work is not referenced if
c                job .eq. 0.
c
c        job     integer.
c                job is an integer that initiates column pivoting.
c                if job .eq. 0, no pivoting is done.
c                if job .ne. 0, pivoting is done.
c
c     on return
c
c        x       x contains in its upper triangle the upper
c                triangular matrix r of the qr factorization.
c                below its diagonal x contains information from
c                which the orthogonal part of the decomposition
c                can be recovered.  note that if pivoting has
c                been requested, the decomposition is not that
c                of the original matrix x but that of x
c                with its columns permuted as described by jpvt.
c
c        qraux   double precision(p).
c                qraux contains further information required to recover
c                the orthogonal part of the decomposition.
c
c        jpvt    jpvt(k) contains the index of the column of the
c                original matrix that has been interchanged into
c                the k-th column, if pivoting was requested.
c
c     linpack. this version dated 08/14/78 .
c     g.w. stewart, university of maryland, argonne national lab.
c
c     dqrdc uses the following functions and subprograms.
c
c     blas daxpy,ddot,dscal,dswap,dnrm2
c     fortran dabs,dmax1,min0,dsqrt
c
c     internal variables
c
      integer j,jp,l,lp1,lup,maxj,pl,pu
      double precision maxnrm,dnrm2,tt
      double precision ddot,nrmxl,t
      logical negj,swapj
c
c
      pl = 1
      pu = 0
      if (job .eq. 0) go to 60
c
c        pivoting has been requested.  rearrange the columns
c        according to jpvt.
c
         do 20 j = 1, p
            swapj = jpvt(j) .gt. 0
            negj = jpvt(j) .lt. 0
            jpvt(j) = j
            if (negj) jpvt(j) = -j
            if (.not.swapj) go to 10
               if (j .ne. pl) call dswap(n,x(1,pl),1,x(1,j),1)
               jpvt(j) = jpvt(pl)
               jpvt(pl) = j
               pl = pl + 1
   10       continue
   20    continue
         pu = p
         do 50 jj = 1, p
            j = p - jj + 1
            if (jpvt(j) .ge. 0) go to 40
               jpvt(j) = -jpvt(j)
               if (j .eq. pu) go to 30
                  call dswap(n,x(1,pu),1,x(1,j),1)
                  jp = jpvt(pu)
                  jpvt(pu) = jpvt(j)
                  jpvt(j) = jp
   30          continue
               pu = pu - 1
   40       continue
   50    continue
   60 continue
c
c     compute the norms of the free columns.
c
      if (pu .lt. pl) go to 80
      do 70 j = pl, pu
         qraux(j) = dnrm2(n,x(1,j),1)
         work(j) = qraux(j)
   70 continue
   80 continue
c
c     perform the householder reduction of x.
c
      lup = min0(n,p)
      do 200 l = 1, lup
         if (l .lt. pl .or. l .ge. pu) go to 120
c
c           locate the column of largest norm and bring it
c           into the pivot position.
c
            maxnrm = 0.0d0
            maxj = l
            do 100 j = l, pu
               if (qraux(j) .le. maxnrm) go to 90
                  maxnrm = qraux(j)
                  maxj = j
   90          continue
  100       continue
            if (maxj .eq. l) go to 110
               call dswap(n,x(1,l),1,x(1,maxj),1)
               qraux(maxj) = qraux(l)
               work(maxj) = work(l)
               jp = jpvt(maxj)
               jpvt(maxj) = jpvt(l)
               jpvt(l) = jp
  110       continue
  120    continue
         qraux(l) = 0.0d0
         if (l .eq. n) go to 190
c
c           compute the householder transformation for column l.
c
            nrmxl = dnrm2(n-l+1,x(l,l),1)
            if (nrmxl .eq. 0.0d0) go to 180
               if (x(l,l) .ne. 0.0d0) nrmxl = dsign(nrmxl,x(l,l))
               call dscal(n-l+1,1.0d0/nrmxl,x(l,l),1)
               x(l,l) = 1.0d0 + x(l,l)
c
c              apply the transformation to the remaining columns,
c              updating the norms.
c
               lp1 = l + 1
               if (p .lt. lp1) go to 170
               do 160 j = lp1, p
                  t = -ddot(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
                  call daxpy(n-l+1,t,x(l,l),1,x(l,j),1)
                  if (j .lt. pl .or. j .gt. pu) go to 150
                  if (qraux(j) .eq. 0.0d0) go to 150
                     tt = 1.0d0 - (dabs(x(l,j))/qraux(j))**2
                     tt = dmax1(tt,0.0d0)
                     t = tt
                     tt = 1.0d0 + 0.05d0*tt*(qraux(j)/work(j))**2
                     if (tt .eq. 1.0d0) go to 130
                        qraux(j) = qraux(j)*dsqrt(t)
                     go to 140
  130                continue
                        qraux(j) = dnrm2(n-l,x(l+1,j),1)
                        work(j) = qraux(j)
  140                continue
  150             continue
  160          continue
  170          continue
c
c              save the transformation.
c
               qraux(l) = x(l,l)
               x(l,l) = -nrmxl
  180       continue
  190    continue
  200 continue
      return
      end

C---------------------------------------------------------------------C
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n .lt. 1 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end

C----------------------------------------------------------------------C
C     build kd tree
      subroutine lowesb(xx,yy,ww,diagl,infl,iv,liv,lv,wv)
      logical infl,setlf
      integer execnt
      integer iv(*)
      DOUBLE PRECISION trl
      DOUBLE PRECISION diagl(*),wv(*),ww(*),xx(*),yy(*)
      external ehg131,ehg182,ehg183
      integer ifloor
      external ifloor
      save execnt
      data execnt /0/
      execnt=execnt+1
      if(.not.(iv(28).ne.173))then
         call ehg182(174)
      end if
      if(iv(28).ne.172)then
         if(.not.(iv(28).eq.171))then
            call ehg182(171)
         end if
      end if
      iv(28)=173
      if(infl)then
         trl=1.D0
      else
         trl=0.D0
      end if
      setlf=(iv(27).ne.iv(25))
      call ehg131(xx,yy,ww,trl,diagl,iv(20),iv(29),iv(3),iv(2),iv(5),iv(
     +17),iv(4),iv(6),iv(14),iv(19),wv(1),iv(iv(7)),iv(iv(8)),iv(iv(9)),
     +iv(iv(10)),iv(iv(22)),iv(iv(27)),wv(iv(11)),iv(iv(23)),wv(iv(13)),
     +wv(iv(12)),wv(iv(15)),wv(iv(16)),wv(iv(18)),ifloor(iv(3)*wv(2)),wv
     +(3),wv(iv(26)),wv(iv(24)),wv(4),iv(30),iv(33),iv(32),iv(41),iv(iv(
     +25)),wv(iv(34)),setlf)
      if(iv(14).lt.iv(6)+DFLOAT(iv(4))/2.D0)then
         call ehg183('Warning. k-d tree limited by memory; nvmax=',iv(14
     +),1,1)
      else
         if(iv(17).lt.iv(5)+2)then
            call ehg183('Warning. k-d tree limited by memory. ncmax=',iv
     +(17),1,1)
         end if
      end if
      return
      end

C----------------------------------------------------------------------C
C     setup workspace
      subroutine lowesd(versio,iv,liv,lv,v,d,n,f,ideg,nvmax,setlf)
      logical setlf
      integer bound,d,execnt,i,i1,i2,ideg,j,liv,lv,n,ncmax,nf,nvmax,vc,v
     +ersio
      integer iv(liv)
      double precision f
      double precision v(lv)
      external ehg182
      integer ifloor
      external ifloor
      save execnt
      data execnt /0/
c     version -> versio
      execnt=execnt+1
      if(.not.(versio.eq.106))then
         call ehg182(100)
      end if
      iv(28)=171
      iv(2)=d
      iv(3)=n
      vc=2**d
      iv(4)=vc
      if(.not.(0.lt.f))then
         call ehg182(120)
      end if
      nf=min(n,ifloor(n*f))
      iv(19)=nf
      iv(20)=1
      if(ideg.eq.0)then
         i1=1
      else
         if(ideg.eq.1)then
            i1=d+1
         else
            if(ideg.eq.2)then
               i1=dfloat((d+2)*(d+1))/2.d0
            end if
         end if
      end if
      iv(29)=i1
      iv(21)=1
      iv(14)=nvmax
      ncmax=nvmax
      iv(17)=ncmax
      iv(30)=0
      iv(32)=ideg
      if(.not.(ideg.ge.0))then
         call ehg182(195)
      end if
      if(.not.(ideg.le.2))then
         call ehg182(195)
      end if
      iv(33)=d
      do 3 i2=41,49
         iv(i2)=ideg
    3 continue
      iv(7)=50
      iv(8)=iv(7)+ncmax
      iv(9)=iv(8)+vc*ncmax
      iv(10)=iv(9)+ncmax
      iv(22)=iv(10)+ncmax
c     initialize permutation
      j=iv(22)-1
      do 4 i=1,n
         iv(j+i)=i
    4 continue
      iv(23)=iv(22)+n
      iv(25)=iv(23)+nvmax
      if(setlf)then
         iv(27)=iv(25)+nvmax*nf
      else
         iv(27)=iv(25)
      end if
      bound=iv(27)+n
      if(.not.(bound-1.le.liv))then
         call ehg182(102)
      end if
      iv(11)=50
      iv(13)=iv(11)+nvmax*d
      iv(12)=iv(13)+(d+1)*nvmax
      iv(15)=iv(12)+ncmax
      iv(16)=iv(15)+n
      iv(18)=iv(16)+nf
      iv(24)=iv(18)+iv(29)*nf
      iv(34)=iv(24)+(d+1)*nvmax
      if(setlf)then
         iv(26)=iv(34)+(d+1)*nvmax*nf
      else
         iv(26)=iv(34)
      end if
      bound=iv(26)+nf
      if(.not.(bound-1.le.lv))then
         call ehg182(103)
      end if
      v(1)=f
      v(2)=0.05d0
      v(3)=0.d0
      v(4)=1.d0
      return
      end

C----------------------------------------------------------------------C
C     evaluate smooth at z
      subroutine lowese(iv,liv,lv,wv,m,z,s)
      integer execnt,m
      integer iv(*)
      double precision s(m),wv(*),z(m,1)
      external ehg133,ehg182
      save execnt
      data execnt /0/
      execnt=execnt+1
      if(.not.(iv(28).ne.172))then
         call ehg182(172)
      end if
      if(.not.(iv(28).eq.173))then
         call ehg182(173)
      end if
      call ehg133(iv(3),iv(2),iv(4),iv(14),iv(5),iv(17),iv(iv(7)),iv(iv(
     +8)),iv(iv(9)),iv(iv(10)),wv(iv(11)),wv(iv(13)),wv(iv(12)),m,z,s)
      return
      end

C----------------------------------------------------------------------C
C     slow smooth at z
      subroutine lowesf(xx,yy,ww,iv,liv,lv,wv,m,z,l,ihat,s)
      logical i1
      integer execnt,ihat,m,n
      integer iv(*)
      double precision l(m,*),s(m),wv(*),ww(*),xx(*),yy(*),z(m,1)
      external ehg182,ehg136
      save execnt
      data execnt /0/
      execnt=execnt+1
      if(171.le.iv(28))then
         i1=(iv(28).le.174)
      else
         i1=.false.
      end if
      if(.not.i1)then
         call ehg182(171)
      end if
      iv(28)=172
      if(.not.(iv(14).ge.iv(19)))then
         call ehg182(186)
      end if
      call ehg136(z,m,m,iv(3),iv(2),iv(19),wv(1),xx,iv(iv(22)),yy,ww,iv(
     +20),iv(29),wv(iv(15)),wv(iv(16)),wv(iv(18)),0,l,ihat,wv(iv(26)),wv
     +(4),iv(30),iv(33),iv(32),iv(41),s)
      return
      end

C----------------------------------------------------------------------C
C     explicit hat matrix mapping y to z
      subroutine lowesl(iv,liv,lv,wv,m,z,l)
      integer execnt,m,n
      integer iv(*)
      double precision l(m,*),wv(*),z(m,1)
      external ehg182,ehg191
      save execnt
      data execnt /0/
      execnt=execnt+1
      if(.not.(iv(28).ne.172))then
         call ehg182(172)
      end if
      if(.not.(iv(28).eq.173))then
         call ehg182(173)
      end if
      if(.not.(iv(26).ne.iv(34)))then
         call ehg182(175)
      end if
      call ehg191(m,z,l,iv(2),iv(3),iv(19),iv(6),iv(17),iv(4),iv(iv(7)),
     +wv(iv(12)),iv(iv(10)),iv(iv(9)),iv(iv(8)),wv(iv(11)),iv(14),wv(iv(
     +24)),wv(iv(34)),iv(iv(25)))
      return
      end

C----------------------------------------------------------------------C
C     rebuild with new data values (does not change y)
      subroutine lowesr(yy,iv,liv,lv,wv)
      integer execnt
      integer iv(*)
      DOUBLE PRECISION wv(*),yy(*)
      external ehg182,ehg192
      save execnt
      data execnt /0/
      execnt=execnt+1
      if(.not.(iv(28).ne.172))then
         call ehg182(172)
      end if
      if(.not.(iv(28).eq.173))then
         call ehg182(173)
      end if
      call ehg192(yy,iv(2),iv(3),iv(19),iv(6),iv(14),wv(iv(13)),wv(iv(34
     +)),iv(iv(25)))
      return
      end

C----------------------------------------------------------------------C
C     robustness weights
      subroutine lowesw(res,n,rw,pi)
      integer identi,execnt,i,i1,n,nh
      integer pi(n)
      double precision cmad,rsmall
      double precision res(n),rw(n)
      external ehg106
      integer ifloor
      external ifloor
      double precision d1mach
      external d1mach
      save execnt
      data execnt /0/
c     Identity -> identi
      execnt=execnt+1
c     tranliterated from Devlin's ratfor
c     find median of absolute residuals
      do 3 i1=1,n
         rw(i1)=dabs(res(i1))
    3 continue
      do 4 identi=1,n
         pi(identi)=identi
    4 continue
      nh=ifloor(dfloat(n)/2.d0)+1
c     partial sort to find 6*mad
      call ehg106(1,n,nh,1,rw,pi,n)
      if((n-nh)+1.lt.nh)then
         call ehg106(1,nh-1,nh-1,1,rw,pi,n)
         cmad=3*(rw(pi(nh))+rw(pi(nh-1)))
      else
         cmad=6*rw(pi(nh))
      end if
      rsmall=d1mach(1)
      if(cmad.lt.rsmall)then
         do 5 i1=1,n
            rw(i1)=1
    5    continue
      else
         do 6 i=1,n
            if(cmad*0.999d0.lt.rw(i))then
               rw(i)=0
            else
               if(cmad*0.001d0.lt.rw(i))then
                  rw(i)=(1-(rw(i)/cmad)**2)**2
               else
                  rw(i)=1
               end if
            end if
    6    continue
      end if
      return
      end

C----------------------------------------------------------------------C
C     pseudovalues
      subroutine lowesp(n,y,yhat,pwgts,rwgts,pi,ytilde)
      integer identi,execnt,i2,i3,i5,m,n
      integer pi(n)
      double precision c,i1,i4,mad
      double precision pwgts(n),rwgts(n),y(n),yhat(n),ytilde(n)
      external ehg106
      integer ifloor
      external ifloor
      save execnt
      data execnt /0/
c     Identity -> identi
      execnt=execnt+1
c     median absolute deviation
      do 3 i5=1,n
         ytilde(i5)=dabs(y(i5)-yhat(i5))*dsqrt(pwgts(i5))
    3 continue
      do 4 identi=1,n
         pi(identi)=identi
    4 continue
      m=ifloor(dfloat(n)/2.d0)+1
      call ehg106(1,n,m,1,ytilde,pi,n)
      if((n-m)+1.lt.m)then
         call ehg106(1,m-1,m-1,1,ytilde,pi,n)
         mad=(ytilde(pi(m-1))+ytilde(pi(m)))/2
      else
         mad=ytilde(pi(m))
      end if
c     magic constant
      c=(6*mad)**2/5
      do 5 i5=1,n
         ytilde(i5)=1-((y(i5)-yhat(i5))**2*pwgts(i5))/c
    5 continue
      do 6 i5=1,n
         ytilde(i5)=ytilde(i5)*dsqrt(rwgts(i5))
    6 continue
      if(n.le.0)then
         i4=0.d0
      else
         i3=n
         i1=ytilde(i3)
         do 7 i2=i3-1,1,-1
            i1=ytilde(i2)+i1
    7    continue
         i4=i1
      end if
      c=n/i4
c     pseudovalues
      do 8 i5=1,n
         ytilde(i5)=yhat(i5)+(c*rwgts(i5))*(y(i5)-yhat(i5))
    8 continue
      return
      end

C----------------------------------------------------------------------C
C     rbuild
      subroutine ehg124(ll,uu,d,n,nv,nc,ncmax,vc,x,pi,a,xi,lo,hi,c,v,vhi
     +t,nvmax,fc,fd,dd)
      logical i1,i2,i3,leaf
      integer d,dd,execnt,fc,i4,inorm2,k,l,ll,m,n,nc,ncmax,nv,nvmax,p,u,
     +uu,vc,lower,upper,check,offset
      integer a(ncmax),c(vc,ncmax),hi(ncmax),lo(ncmax),pi(n),vhit(nvmax)
      DOUBLE PRECISION diam,fd
      DOUBLE PRECISION diag(8),sigma(8),v(nvmax,d),x(n,d),xi(ncmax)
      external ehg125,ehg106,ehg129
      integer IDAMAX
      external IDAMAX
      save execnt
      data execnt /0/
      execnt=execnt+1
      p=1
      l=ll
      u=uu
      lo(p)=l
      hi(p)=u
c     top of while loop
    3 if(.not.(p.le.nc))goto 4
         do 5 i4=1,dd
            diag(i4)=v(c(vc,p),i4)-v(c(1,p),i4)
    5    continue
         diam=0
         do 6 inorm2=1,dd
            diam=diam+diag(inorm2)**2
    6    continue
         diam=DSQRT(diam)
         if((u-l)+1.le.fc)then
            i1=.true.
         else
            i1=(diam.le.fd)
         end if
         if(i1)then
            leaf=.true.
         else
            if(ncmax.lt.nc+2)then
               i2=.true.
            else
               i2=(nvmax.lt.nv+DFLOAT(vc)/2.D0)
            end if
            leaf=i2
         end if
         if(.not.leaf)then
            call ehg129(l,u,dd,x,pi,n,sigma)
            k=IDAMAX(dd,sigma,1)
            m=DFLOAT(l+u)/2.D0
            call ehg106(l,u,m,1,x(1,k),pi,n)

c bug fix from btyner@gmail.com 2006-07-20
      offset = 0
    7 if(((m+offset).ge.u).or.((m+offset).lt.l))then
         goto 8
      else
        if(offset.lt.0)then
          lower = l
          check = m + offset
          upper = check
        else
          lower = m + offset + 1
          check = lower
          upper = u
        end if
        call ehg106(lower,upper,check,1,x(1,k),pi,n)
        if(x(pi(m + offset),k).eq.x(pi(m+offset+1),k))then
          offset = (-offset)
          if(offset.ge.0)then
            offset = offset + 1
          end if
      goto 7
        else
          m = m + offset
          goto 8
        end if
      end if

    8       if(v(c(1,p),k).eq.x(pi(m),k))then
               leaf=.true.
            else
               leaf=(v(c(vc,p),k).eq.x(pi(m),k))
            end if
         end if
         if(leaf)then
            a(p)=0
         else
            a(p)=k
            xi(p)=x(pi(m),k)
c           left son
            nc=nc+1
            lo(p)=nc
            lo(nc)=l
            hi(nc)=m
c           right son
            nc=nc+1
            hi(p)=nc
            lo(nc)=m+1
            hi(nc)=u
            call ehg125(p,nv,v,vhit,nvmax,d,k,xi(p),2**(k-1),2**(d-k),c(
     +1,p),c(1,lo(p)),c(1,hi(p)))
         end if
         p=p+1
         l=lo(p)
         u=hi(p)
         goto 3
c     bottom of while loop
    4 return
      end

C----------------------------------------------------------------------C
C     spread
      subroutine ehg129(l,u,d,x,pi,n,sigma)
      integer d,execnt,i,k,l,n,u
      integer pi(n)
      DOUBLE PRECISION machin,alpha,beta,t
      DOUBLE PRECISION sigma(d),x(n,d)
      DOUBLE PRECISION D1MACH
      external D1MACH
      save machin,execnt
      data execnt /0/
c     MachInf -> machin
      execnt=execnt+1
      if(execnt.eq.1)then
         machin=D1MACH(2)
      end if
      do 3 k=1,d
         alpha=machin
         beta=-machin
         do 4 i=l,u
            t=x(pi(i),k)
            alpha=min(alpha,x(pi(i),k))
            beta=max(beta,t)
    4    continue
         sigma(k)=beta-alpha
    3 continue
      return
      end

C----------------------------------------------------------------------C
C     vleaf
      subroutine ehg137(z,kappa,leaf,nleaf,d,nv,nvmax,ncmax,vc,a,xi,lo,h
     +i,c,v)
      integer d,execnt,nc,ncmax,nleaf,p,stackt
      integer a(ncmax),hi(ncmax),leaf(256),lo(ncmax),pstack(20)
      DOUBLE PRECISION xi(ncmax),z(d)
      external ehg182
      save execnt
      data execnt /0/
c     stacktop -> stackt
      execnt=execnt+1
c     find leaf cells affected by $z$
      stackt=0
      p=1
      nleaf=0
c     top of while loop
    3 if(.not.(0.lt.p))goto 4
         if(a(p).eq.0)then
c           leaf
            nleaf=nleaf+1
            leaf(nleaf)=p
c           Pop
            if(stackt.ge.1)then
               p=pstack(stackt)
            else
               p=0
            end if
            stackt=max(0,stackt-1)
         else
            if(z(a(p)).eq.xi(p))then
c              Push
               stackt=stackt+1
               if(.not.(stackt.le.20))then
                  call ehg182(187)
               end if
               pstack(stackt)=hi(p)
               p=lo(p)
            else
               if(z(a(p)).le.xi(p))then
                  p=lo(p)
               else
                  p=hi(p)
               end if
            end if
         end if
         goto 3
c     bottom of while loop
    4 if(.not.(nleaf.le.256))then
         call ehg182(185)
      end if
      return
      end
