
!A = numpy.array([[12.0, 28.0, 76.0, 220.0], [16.0, 32.0, 80.0, 224.0], [24.0, 40.0, 88.0, 232.0], [40.0, 56.0, 104.0, 248.0]], dtype='float64')
! B = numpy.array([[2.0, 4.0, 10.0, 28.0], [3.0, 5.0, 11.0, 29.0], [5.0, 7.0, 13.0, 31.0], [9.0, 11.0, 17.0, 35.0]], dtype='float64')
! D, V = scipy.linalg.eig(A, B); D 
implicit none
integer, parameter :: n = 4

integer :: lda, ldb, ldvr, ldvl, lwork, info
character*1 :: jobvl, jobvr
real*8 :: alphar(n)
real*8 :: alphai(n)
real*8 :: beta(n)
real*8 :: vl(n, n)
real*8 :: vr(n, n)
real*8, allocatable :: work(:)


real*8 :: a(n, n)
real*8 :: b(n, n)

a(1, :) = (/12.0, 28.0, 76.0, 220.0/)
a(2, :) = (/16.0, 32.0, 80.0, 224.0/)
a(3, :) = (/24.0, 40.0, 88.0, 232.0/)
a(4, :) = (/40.0, 56.0, 104.0, 248.0/)

b(1, :) = (/2.0, 4.0, 10.0, 28.0/)
b(2, :) = (/3.0, 5.0, 11.0, 29.0/)
b(3, :) = (/5.0, 7.0, 13.0, 31.0/)
b(4, :) = (/9.0, 11.0, 17.0, 35.0/)

lda = n
ldb = n
ldvr = n
ldvl = n
jobvr = 'V'
jobvl = 'V'

! workspace query

allocate(work(1:8*n)) ! min value
lwork = -1

print*, 'workspace query: lwork = ', lwork

call dggev(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info)

print*, 'info = ', info
lwork = int(work(1))
print*, 'opt lwork =', lwork

! do the work

deallocate(work)
allocate(work(1:lwork))

call dggev(jobvl, jobvr, n, a, lda, b, ldb, alphar, alphai, beta, vl, ldvl, vr, ldvr, work, lwork, info)

print*
print*, 'info = ', info
print*, 'alphar = ', alphar
print*, 'alphai = ', alphai
print*, 'beta = ', beta
print*
print*, 'Re(eigv) = ', alphar / beta
print*, 'Im(eigv) = ', alphai / beta

deallocate(work)
end
