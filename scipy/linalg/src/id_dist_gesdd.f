c
c ILP64 LAPACK wrappers for id_dist
c
      subroutine dgesdd_id_dist(jobz, m, n, a, lda, s, u, ldu,
     *     vt, ldvt, work, lwork, iwork, info)
      implicit none
      character*1 jobz
      integer m, n, lda, ldu, ldvt, lwork, iwork(*), info
      double precision a, s, u, vt, work

      integer*8 m2, n2, lda2, ldu2, ldvt2, lwork2, info2
      integer*8 iwork2(8*min(m,n))

      m2 = m
      n2 = n
      lda2 = lda
      ldu2 = ldu
      ldvt2 = ldvt
      lwork2 = lwork
      call dgesdd(jobz, m2, n2, a, lda2, s, u, ldu2,
     *     vt, ldvt2, work, lwork2, iwork2, info2)
      info = int(info2)
      return
      end

      subroutine zgesdd_id_dist(jobz, m, n, a, lda, s, u, ldu,
     *     vt, ldvt, work, lwork, rwork, iwork, info)
      implicit none
      character*1 jobz
      integer m, n, lda, ldu, ldvt, lwork, iwork(*), info
      double complex a, s, u, vt, work
      double precision rwork

      integer*8 m2, n2, lda2, ldu2, ldvt2, lwork2, info2
      integer*8 iwork2(8*min(m,n))

      m2 = m
      n2 = n
      lda2 = lda
      ldu2 = ldu
      ldvt2 = ldvt
      lwork2 = lwork
      call zgesdd(jobz, m2, n2, a, lda2, s, u, ldu2,
     *     vt, ldvt2, work, lwork2, rwork, iwork2, info2)
      info = int(info2)
      return
      end
