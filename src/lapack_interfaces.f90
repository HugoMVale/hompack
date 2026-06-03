module lapack_interfaces
  !! Double precision interfaces for the LAPACK procedures used by `hompack`.

   use hompack_kinds, only: dp
   implicit none

   interface

      subroutine dgeqpf(m, n, a, lda, jpvt, tau, work, info)
      !! QR factorization with column pivoting: `A*P = Q*R`.
         import :: dp
         integer, intent(in) :: m
         integer, intent(in) :: n
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(in) :: lda
         integer, intent(inout) :: jpvt(*)
         real(dp), intent(out) :: tau(*)
         real(dp), intent(out) :: work(*)
         integer, intent(out) :: info
      end subroutine dgeqpf

      subroutine dgeqr2(m, n, a, lda, tau, work, info)
      !! Unblocked QR factorization: `A = Q*R`.
         import :: dp
         integer, intent(in) :: m
         integer, intent(in) :: n
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(in) :: lda
         real(dp), intent(out) :: tau(*)
         real(dp), intent(out) :: work(*)
         integer, intent(out) :: info
      end subroutine dgeqr2

      subroutine dgeqrf(m, n, a, lda, tau, work, lwork, info)
      !! Blocked QR factorization: `A = Q*R`.
         import :: dp
         integer, intent(in) :: m
         integer, intent(in) :: n
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(in) :: lda
         real(dp), intent(out) :: tau(*)
         real(dp), intent(out) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(out) :: info
      end subroutine dgeqrf

      pure logical function disnan(din)
      !! Test if a double precision value is NaN.
         import :: dp
         real(dp), intent(in) :: din
      end function disnan

      pure subroutine dlacpy(uplo, m, n, a, lda, b, ldb)
      !! Copy all or part of a matrix: `B = A`.
         import :: dp
         character, intent(in) :: uplo
         integer, intent(in) :: m
         integer, intent(in) :: n
         real(dp), intent(in) :: a(lda, *)
         integer, intent(in) :: lda
         real(dp), intent(out) :: b(ldb, *)
         integer, intent(in) :: ldb
      end subroutine dlacpy

      pure subroutine dlaic1(job, j, x, sest, w, gamma, sestpr, s, c)
      !! Apply one step of incremental condition estimation.
         import :: dp
         integer, intent(in) :: job
         integer, intent(in) :: j
         real(dp), intent(in) :: x(j)
         real(dp), intent(in) :: sest
         real(dp), intent(in) :: w(j)
         real(dp), intent(in) :: gamma
         real(dp), intent(out) :: sestpr
         real(dp), intent(out) :: s
         real(dp), intent(out) :: c
      end subroutine dlaic1

      pure logical function dlaisnan(din1, din2)
      !! Test if two double precision values differ (NaN helper).
         import :: dp
         real(dp), intent(in) :: din1
         real(dp), intent(in) :: din2
      end function dlaisnan

      pure real(dp) function dlamch(cmach)
      !! Determine double precision machine parameters.
         import :: dp
         character, intent(in) :: cmach
      end function dlamch

      pure real(dp) function dlapy2(x, y)
      !! Compute `sqrt(x**2 + y**2)` without undue overflow or underflow.
         import :: dp
         real(dp), intent(in) :: x
         real(dp), intent(in) :: y
      end function dlapy2

      pure subroutine dlarf(side, m, n, v, incv, tau, c, ldc, work)
      !! Apply a Householder reflector `H = I - tau*v*v**T` to a matrix.
         import :: dp
         character, intent(in) :: side
         integer, intent(in) :: m
         integer, intent(in) :: n
         real(dp), intent(in) :: v(*)
         integer, intent(in) :: incv
         real(dp), intent(in) :: tau
         real(dp), intent(inout) :: c(ldc, *)
         integer, intent(in) :: ldc
         real(dp), intent(out) :: work(*)
      end subroutine dlarf

      pure subroutine dlarf1f(side, m, n, v, incv, tau, c, ldc, work)
      !! Apply a Householder reflector with first element 1.
         import :: dp
         character, intent(in) :: side
         integer, intent(in) :: m
         integer, intent(in) :: n
         real(dp), intent(in) :: v(*)
         integer, intent(in) :: incv
         real(dp), intent(in) :: tau
         real(dp), intent(inout) :: c(ldc, *)
         integer, intent(in) :: ldc
         real(dp), intent(out) :: work(*)
      end subroutine dlarf1f

      pure subroutine dlarfb(side, trans, direct, storev, m, n, k, v, ldv, t, ldt, c, ldc, work, ldwork)
      !! Apply a block Householder reflector `H = I - V*T*V**T` to a matrix.
         import :: dp
         character, intent(in) :: side
         character, intent(in) :: trans
         character, intent(in) :: direct
         character, intent(in) :: storev
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: k
         real(dp), intent(in) :: v(ldv, *)
         integer, intent(in) :: ldv
         real(dp), intent(in) :: t(ldt, *)
         integer, intent(in) :: ldt
         real(dp), intent(inout) :: c(ldc, *)
         integer, intent(in) :: ldc
         real(dp), intent(out) :: work(ldwork, *)
         integer, intent(in) :: ldwork
      end subroutine dlarfb

      pure subroutine dlarfg(n, alpha, x, incx, tau)
      !! Generate a Householder reflector `H` such that `H*[alpha; x] = [beta; 0]`.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(inout) :: alpha
         real(dp), intent(inout) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(out) :: tau
      end subroutine dlarfg

      pure subroutine dlarft(direct, storev, n, k, v, ldv, tau, t, ldt)
      !! Form the triangular factor `T` of a block Householder reflector.
         import :: dp
         character, intent(in) :: direct
         character, intent(in) :: storev
         integer, intent(in) :: n
         integer, intent(in) :: k
         real(dp), intent(in) :: v(ldv, *)
         integer, intent(in) :: ldv
         real(dp), intent(in) :: tau(*)
         real(dp), intent(out) :: t(ldt, *)
         integer, intent(in) :: ldt
      end subroutine dlarft

      subroutine dorg2r(m, n, k, a, lda, tau, work, info)
      !! Unblocked generation of Q from QR factorization.
         import :: dp
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: k
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(in) :: lda
         real(dp), intent(in) :: tau(*)
         real(dp), intent(out) :: work(*)
         integer, intent(out) :: info
      end subroutine dorg2r

      subroutine dorgqr(m, n, k, a, lda, tau, work, lwork, info)
      !! Blocked generation of Q from QR factorization.
         import :: dp
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: k
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(in) :: lda
         real(dp), intent(in) :: tau(*)
         real(dp), intent(out) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(out) :: info
      end subroutine dorgqr

      subroutine dorm2r(side, trans, m, n, k, a, lda, tau, c, ldc, work, info)
      !! Unblocked application of Q from QR factorization to a matrix.
         import :: dp
         character, intent(in) :: side
         character, intent(in) :: trans
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: k
         real(dp), intent(in) :: a(lda, *)
         integer, intent(in) :: lda
         real(dp), intent(in) :: tau(*)
         real(dp), intent(inout) :: c(ldc, *)
         integer, intent(in) :: ldc
         real(dp), intent(out) :: work(*)
         integer, intent(out) :: info
      end subroutine dorm2r

      subroutine dormqr(side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
      !! Blocked application of Q from QR factorization to a matrix.
         import :: dp
         character, intent(in) :: side
         character, intent(in) :: trans
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: k
         real(dp), intent(in) :: a(lda, *)
         integer, intent(in) :: lda
         real(dp), intent(in) :: tau(*)
         real(dp), intent(inout) :: c(ldc, *)
         integer, intent(in) :: ldc
         real(dp), intent(out) :: work(*)
         integer, intent(in) :: lwork
         integer, intent(out) :: info
      end subroutine dormqr

      pure integer function ieeeck(ispec, zero, one)
      !! Check IEEE arithmetic support.
         integer, intent(in) :: ispec
         real, intent(in) :: zero
         real, intent(in) :: one
      end function ieeeck

      pure integer function iladlc(m, n, a, lda)
      !! Find the last non-zero column of a matrix.
         import :: dp
         integer, intent(in) :: m
         integer, intent(in) :: n
         real(dp), intent(in) :: a(lda, *)
         integer, intent(in) :: lda
      end function iladlc

      pure integer function iladlr(m, n, a, lda)
      !! Find the last non-zero row of a matrix.
         import :: dp
         integer, intent(in) :: m
         integer, intent(in) :: n
         real(dp), intent(in) :: a(lda, *)
         integer, intent(in) :: lda
      end function iladlr

      pure integer function ilaenv(ispec, name, opts, n1, n2, n3, n4)
      !! Return LAPACK environment parameters (e.g. optimal block size).
         integer, intent(in) :: ispec
         character(len=*), intent(in) :: name
         character(len=*), intent(in) :: opts
         integer, intent(in) :: n1
         integer, intent(in) :: n2
         integer, intent(in) :: n3
         integer, intent(in) :: n4
      end function ilaenv

      pure integer function iparmq(ispec, name, opts, n, ilo, ihi, lwork)
      !! Return shift parameters for the double-shift QR iteration.
         integer, intent(in) :: ispec
         character(len=*), intent(in) :: name
         character(len=*), intent(in) :: opts
         integer, intent(in) :: n
         integer, intent(in) :: ilo
         integer, intent(in) :: ihi
         integer, intent(in) :: lwork
      end function iparmq

      subroutine xerbla(srname, info)
      !! Error handler called by LAPACK/BLAS routines.
         character(len=*), intent(in) :: srname
         integer, intent(in) :: info
      end subroutine xerbla

   end interface

end module lapack_interfaces
