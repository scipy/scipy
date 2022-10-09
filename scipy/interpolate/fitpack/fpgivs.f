      recursive subroutine fpgivs(piv,ww,cos,sin)
      implicit none
c  subroutine fpgivs calculates the parameters of a givens
c  transformation .
c  ..
c  ..scalar arguments..
      real*8 piv,ww,cos,sin
c  ..local scalars..
      real*8 dd,one,store
c  ..function references..
      real*8 abs,sqrt
c  ..
      one = 0.1e+01
      store = abs(piv)
      if(store.ge.ww) dd = store*sqrt(one+(ww/piv)**2)
      if(store.lt.ww) dd = ww*sqrt(one+(piv/ww)**2)
      cos = ww/dd
      sin = piv/dd
      ww = dd
      return
      end
