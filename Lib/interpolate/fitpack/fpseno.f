      subroutine fpseno(maxtr,up,left,right,info,merk,ibind,nbind)
c  subroutine fpseno fetches a branch of a triply linked tree the
c  information of which is kept in the arrays up,left,right and info.
c  the branch has a specified length nbind and is determined by the
c  parameter merk which points to its terminal node. the information
c  field of the nodes of this branch is stored in the array ibind. on
c  exit merk points to a new branch of length nbind or takes the value
c  1 if no such branch was found.
c  ..
c  ..scalar arguments..
      integer maxtr,merk,nbind
c  ..array arguments..
      integer up(maxtr),left(maxtr),right(maxtr),info(maxtr),
     * ibind(nbind)
c  ..scalar arguments..
      integer i,j,k
c  ..
      k = merk
      j = nbind
      do 10 i=1,nbind
        ibind(j) = info(k)
        k = up(k)
        j = j-1
  10  continue
  20  k = right(merk)
      if(k.ne.0) go to 30
      merk = up(merk)
      if(merk-1) 40,40,20
  30  merk = k
      k = left(merk)
      if(k.ne.0) go to 30
  40  return
      end
