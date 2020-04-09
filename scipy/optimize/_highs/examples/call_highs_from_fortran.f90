program fortrantest
  use, intrinsic :: iso_c_binding
  use highs_lp_solver
  implicit none

  integer ( c_int ), parameter :: n = 2
  integer ( c_int ), parameter :: m = 2
  integer ( c_int ), parameter :: nz = 4

  real ( c_double ) cc(n)
  real ( c_double ) cl(n)
  real ( c_double ) cu(n)
  real ( c_double ) rl(m)
  real ( c_double ) ru(m)
  integer ( c_int ) as(n+1)
  integer ( c_int ) ai(nz)
  real ( c_double ) av(nz)

  real ( c_double ) cv(n)
  real ( c_double ) cd(n)
  real ( c_double ) rv(m)
  real ( c_double ) rd(m)
  integer ( c_int ) cbs(n)
  integer ( c_int ) rbs(m)
  integer ( c_int ) ms

  cc(1) = 1
  cc(2) = -2
  cl(1) = 0
  cl(2) = 0
  cu(1) = 1000
  cu(2) = 1000
  rl(1) = 0.0
  rl(2) = 0.0
  ru(1) = 10.0
  ru(2) = 10.0
  as(1) = 0
  as(2) = 2
  as(3) = 4
  ai(1) = 0
  ai(2) = 1
  ai(3) = 0
  ai(4) = 1
  av(1) = 1
  av(2) = -1
  av(3) = 3
  av(4) = 0.2

  call Highs_call( n, m, nz, cc, cl, cu, rl, ru, as, ai, av, cv, cd, rv, rd, cbs, rbs, ms)
      

end program fortrantest