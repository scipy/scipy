      subroutine mc21b(n,icn,licn,ip,lenr,iperm,numnz,pr,arp,cv,out)
      integer ip(n)
c   pr(i) is the previous row to i in the depth first search.
c it is used as a work array in the sorting algorithm.
c   elements (iperm(i),i) i=1, ... n  are non-zero at the end of the
c algorithm unless n assignments have not been made.  in which case
c (iperm(i),i) will be zero for n-numnz entries.
c   cv(i) is the most recent row extension at which column i
c was visited.
c   arp(i) is one less than the number of non-zeros in row i
c which have not been scanned when looking for a cheap assignment.
c   out(i) is one less than the number of non-zeros in row i
c which have not been scanned during one pass through the main loop.
      integer icn(licn),lenr(n),iperm(n),pr(n),cv(n),
     1arp(n),out(n)
c
c   initialization of arrays.
      do 10 i=1,n
      arp(i)=lenr(i)-1
      cv(i)=0
   10 iperm(i)=0
      numnz=0
c
c
c   main loop.
c   each pass round this loop either results in a new assignment
c or gives a row with no assignment.
      do 130 jord=1,n
      j=jord
      pr(j)=-1
      do 100 k=1,jord
c look for a cheap assignment
      in1=arp(j)
      if (in1.lt.0) go to 60
      in2=ip(j)+lenr(j)-1
      in1=in2-in1
      do 50 ii=in1,in2
      i=icn(ii)
      if (iperm(i).eq.0) go to 110
   50 continue
c   no cheap assignment in row.
      arp(j)=-1
c   begin looking for assignment chain starting with row j.
   60 out(j)=lenr(j)-1
c inner loop.  extends chain by one or backtracks.
      do 90 kk=1,jord
      in1=out(j)
      if (in1.lt.0) go to 80
      in2=ip(j)+lenr(j)-1
      in1=in2-in1
c forward scan.
      do 70 ii=in1,in2
      i=icn(ii)
      if (cv(i).eq.jord) go to 70
c   column i has not yet been accessed during this pass.
      j1=j
      j=iperm(i)
      cv(i)=jord
      pr(j)=j1
      out(j1)=in2-ii-1
      go to 100
   70 continue
c
c   backtracking step.
   80 j=pr(j)
      if (j.eq.-1) go to 130
   90 continue
c
  100 continue
c
c   new assignment is made.
  110 iperm(i)=j
      arp(j)=in2-ii-1
      numnz=numnz+1
      do 120 k=1,jord
      j=pr(j)
      if (j.eq.-1) go to 130
      ii=ip(j)+lenr(j)-out(j)-2
      i=icn(ii)
      iperm(i)=j
  120 continue
c
  130 continue
c
c   if matrix is structurally singular, we now complete the
c permutation iperm.
      if (numnz.eq.n) return
      do 140 i=1,n
  140 arp(i)=0
      k=0
      do 160 i=1,n
      if (iperm(i).ne.0) go to 150
      k=k+1
      out(k)=i
      go to 160
  150 j=iperm(i)
      arp(j)=i
  160 continue
      k=0
      do 170 i=1,n
      if (arp(i).ne.0) go to 170
      k=k+1
      ioutk=out(k)
      iperm(ioutk)=i
  170 continue
      return
      end
