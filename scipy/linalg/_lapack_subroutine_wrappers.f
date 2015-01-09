      subroutine dlamchwrapper(ret, cmach)
        external dlamch
        double precision dlamch
        double precision ret(1)
        character cmach
        ret(1) = dlamch(cmach)
      end

      subroutine slamchwrapper(ret, cmach)
        external wslamch
        real wslamch
        real ret(1)
        character cmach
        ret(1) = wslamch(cmach)
      end
