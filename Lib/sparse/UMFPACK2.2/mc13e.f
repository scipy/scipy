      subroutine mc13e(n,icn,licn,ip,lenr,arp,ib,num,lowl,numb,prev)
      integer stp,dummy
      integer ip(n)
c
c arp(i) is one less than the number of unsearched edges leaving
c     node i.  at the end of the algorithm it is set to a
c     permutation which puts the matrix in block lower
c     triangular form.
c ib(i) is the position in the ordering of the start of the ith
c     block.  ib(n+1-i) holds the node number of the ith node
c     on the stack.
c lowl(i) is the smallest stack position of any node to which a path
c     from node i has been found.  it is set to n+1 when node i
c     is removed from the stack.
c numb(i) is the position of node i in the stack if it is on
c     it, is the permuted order of node i for those nodes
c     whose final position has been found and is otherwise zero.
c prev(i) is the node at the end of the path when node i was
c     placed on the stack.
      integer icn(licn),lenr(n),arp(n),ib(n),lowl(n),numb(n),
     1prev(n)
c
c
c   icnt is the number of nodes whose positions in final ordering have
c     been found.
      icnt=0
c num is the number of blocks that have been found.
      num=0
      nnm1=n+n-1
c
c initialization of arrays.
      do 20 j=1,n
      numb(j)=0
      arp(j)=lenr(j)-1
   20 continue
c
c
      do 120 isn=1,n
c look for a starting node
      if (numb(isn).ne.0) go to 120
      iv=isn
c ist is the number of nodes on the stack ... it is the stack pointer.
      ist=1
c put node iv at beginning of stack.
      lowl(iv)=1
      numb(iv)=1
      ib(n)=iv
c
c the body of this loop puts a new node on the stack or backtracks.
      do 110 dummy=1,nnm1
      i1=arp(iv)
c have all edges leaving node iv been searched.
      if (i1.lt.0) go to 60
      i2=ip(iv)+lenr(iv)-1
      i1=i2-i1
c
c look at edges leaving node iv until one enters a new node or
c     all edges are exhausted.
      do 50 ii=i1,i2
      iw=icn(ii)
c has node iw been on stack already.
      if (numb(iw).eq.0) go to 100
c update value of lowl(iv) if necessary.
  50  lowl(iv)=min0(lowl(iv),lowl(iw))
c
c there are no more edges leaving node iv.
      arp(iv)=-1
c is node iv the root of a block.
   60 if (lowl(iv).lt.numb(iv)) go to 90
c
c order nodes in a block.
      num=num+1
      ist1=n+1-ist
      lcnt=icnt+1
c peel block off the top of the stack starting at the top and
c     working down to the root of the block.
      do 70 stp=ist1,n
      iw=ib(stp)
      lowl(iw)=n+1
      icnt=icnt+1
      numb(iw)=icnt
      if (iw.eq.iv) go to 80
   70 continue
   80 ist=n-stp
      ib(num)=lcnt
c are there any nodes left on the stack.
      if (ist.ne.0) go to 90
c have all the nodes been ordered.
      if (icnt.lt.n) go to 120
      go to 130
c
c backtrack to previous node on path.
   90 iw=iv
      iv=prev(iv)
c update value of lowl(iv) if necessary.
      lowl(iv)=min0(lowl(iv),lowl(iw))
      go to 110
c
c put new node on the stack.
 100  arp(iv)=i2-ii-1
      prev(iw)=iv
      iv=iw
      ist=ist+1
      lowl(iv)=ist
      numb(iv)=ist
      k=n+1-ist
      ib(k)=iv
  110 continue
c
  120 continue
c
c
c put permutation in the required form.
  130 do 140 i=1,n
      ii=numb(i)
 140  arp(ii)=i
      return
      end
