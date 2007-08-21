      subroutine dluno
     +   (lun, fn)

      integer lun
      character*(*) fn

      open(unit=lun, file=fn, status='new')

      return

      end

      subroutine dlunc
     +   (lun)

      integer lun

      close(unit=lun)

      return

      end
