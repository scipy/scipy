      recursive subroutine fpadno(maxtr,up,left,right,info,count,
     *   merk,jbind,n1,ier)
      implicit none
c  subroutine fpadno adds a branch of length n1 to the triply linked
c  tree,the information of which is kept in the arrays up,left,right
c  and info. the information field of the nodes of this new branch is
c  given in the array jbind. in linking the new branch fpadno takes
c  account of the property of the tree that
c    info(k) < info(right(k)) ; info(k) < info(left(k))
c  if necessary the subroutine calls subroutine fpfrno to collect the
c  free nodes of the tree. if no computer words are available at that
c  moment, the error parameter ier is set to 1.
c  ..
c  ..scalar arguments..
      integer maxtr,count,merk,n1,ier
c  ..array arguments..
      integer up(maxtr),left(maxtr),right(maxtr),info(maxtr),jbind(n1)
c  ..local scalars..
      integer k,niveau,point
      logical bool
c  ..subroutine references..
c    fpfrno
c  ..
      point = 1
      niveau = 1
  10  k = left(point)
      bool = .true.
  20  if(k.eq.0) go to 50
      if (info(k)-jbind(niveau).lt.0) go to 30
      if (info(k)-jbind(niveau).eq.0) go to 40
      go to 50
  30  point = k
      k = right(point)
      bool = .false.
      go to 20
  40  point = k
      niveau = niveau+1
      go to 10
  50  if(niveau.gt.n1) go to 90
      count = count+1
      if(count.le.maxtr) go to 60
      call fpfrno(maxtr,up,left,right,info,point,merk,n1,count,ier)
      if(ier.ne.0) go to 100
  60  info(count) = jbind(niveau)
      left(count) = 0
      right(count) = k
      if(bool) go to 70
      bool = .true.
      right(point) = count
      up(count) = up(point)
      go to 80
  70  up(count) = point
      left(point) = count
  80  point = count
      niveau = niveau+1
      k = 0
      go to 50
  90  ier = 0
 100  return
      end
