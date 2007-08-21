      subroutine xsetf (mflag)
c
c this routine resets the print control flag mflag.
c
      integer mflag, mesflg, lunit
      common /eh0001/ mesflg, lunit
c
      if (mflag .eq. 0 .or. mflag .eq. 1) mesflg = mflag
      return
c----------------------- end of subroutine xsetf -----------------------
      end
