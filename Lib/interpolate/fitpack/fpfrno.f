      subroutine fpfrno(maxtr,up,left,right,info,point,merk,n1,
     * count,ier)
c  subroutine fpfrno collects the free nodes (up field zero) of the
c  triply linked tree the information of which is kept in the arrays
c  up,left,right and info. the maximal length of the branches of the
c  tree is given by n1. if no free nodes are found, the error flag
c  ier is set to 1.
c  ..
c  ..scalar arguments..
      integer maxtr,point,merk,n1,count,ier
c  ..array arguments..
      integer up(maxtr),left(maxtr),right(maxtr),info(maxtr)
c  ..local scalars
      integer i,j,k,l,n,niveau
c  ..
      ier = 1
      if(n1.eq.2) go to 140
      niveau = 1
      count = 2
  10  j = 0
      i = 1
  20  if(j.eq.niveau) go to 30
      k = 0
      l = left(i)
      if(l.eq.0) go to 110
      i = l
      j = j+1
      go to 20
  30  if(i-count) 110,100,40
  40  if(up(count).eq.0) go to 50
      count = count+1
      go to 30
  50  up(count) = up(i)
      left(count) = left(i)
      right(count) = right(i)
      info(count) = info(i)
      if(merk.eq.i) merk = count
      if(point.eq.i) point = count
      if(k.eq.0) go to 60
      right(k) = count
      go to 70
  60  n = up(i)
      left(n) = count
  70  l = left(i)
  80  if(l.eq.0) go to 90
      up(l) = count
      l = right(l)
      go to 80
  90  up(i) = 0
      i = count
 100  count = count+1
 110  l = right(i)
      k = i
      if(l.eq.0) go to 120
      i = l
      go to 20
 120  l = up(i)
      j = j-1
      if(j.eq.0) go to 130
      i = l
      go to 110
 130  niveau = niveau+1
      if(niveau.le.n1) go to 10
      if(count.gt.maxtr) go to 140
      ier = 0
 140  return
      end
