module blas_interfaces
!! Double precision interfaces for the BLAS procedures used by `hompack`.

   use hompack_kinds, only: dp
   implicit none

   interface

      pure subroutine daxpy(n, a, x, incx, y, incy)
      !! Computation `Y = A*X + Y`.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: a
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: incy
      end subroutine daxpy

      pure subroutine dcopy(n, x, incx, y, incy)
      !! Vector copy `Y = X`.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: incy
      end subroutine dcopy

      pure real(dp) function ddot(n, x, incx, y, incy)
      !! Inner product of vectors.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(in) :: y(*)
         integer, intent(in) :: incy
      end function ddot

      pure subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
      !! Matrix-matrix operation `C = alpha*op(A)*op(B) + beta*C`.
         import :: dp
         character, intent(in) :: transa
         character, intent(in) :: transb
         integer, intent(in) :: m
         integer, intent(in) :: n
         integer, intent(in) :: k
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: a(lda, *)
         integer, intent(in) :: lda
         real(dp), intent(in) :: b(ldb, *)
         integer, intent(in) :: ldb
         real(dp), intent(in) :: beta
         real(dp), intent(inout) :: c(ldc, *)
         integer, intent(in) :: ldc
      end subroutine dgemm

      pure real(dp) function dasum(n, x, incx)
      !! Sum of magnitudes of vector components.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
      end function dasum

      pure subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
      !! Matrix-vector operation `y = alpha*A*x + beta*y`.
         import :: dp
         character, intent(in) :: trans
         integer, intent(in) :: m
         integer, intent(in) :: n
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: a(lda, *)
         integer, intent(in) :: lda
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(in) :: beta
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: incy
      end subroutine dgemv

      pure subroutine dger(m, n, alpha, x, incx, y, incy, a, lda)
      !! Rank-1 update `A = alpha*x*y**T + A`.
         import :: dp
         integer, intent(in) :: m
         integer, intent(in) :: n
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(in) :: y(*)
         integer, intent(in) :: incy
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(in) :: lda
      end subroutine dger

      pure real(dp) function dnrm2(n, x, incx)
      !! Euclidean length (L2 Norm) of vector.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
      end function dnrm2

      pure subroutine dscal(n, a, x, incx)
      !! Vector scale `X = A*X`.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: a
         real(dp), intent(inout) :: x(*)
         integer, intent(in) :: incx
      end subroutine dscal

      pure subroutine dswap(n, x, incx, y, incy)
      !! Interchange vectors.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
         real(dp), intent(inout) :: y(*)
         integer, intent(in) :: incy
      end subroutine dswap

      pure subroutine dtpmv(uplo, trans, diag, n, ap, x, incx)
      !! Triangular packed matrix-vector operation `x = A*x` or `x = A**T*x`.
         import :: dp
         character, intent(in) :: uplo
         character, intent(in) :: trans
         character, intent(in) :: diag
         integer, intent(in) :: n
         real(dp), intent(in) :: ap(*)
         real(dp), intent(inout) :: x(*)
         integer, intent(in) :: incx
      end subroutine dtpmv

      pure subroutine dtpsv(uplo, trans, diag, n, ap, x, incx)
      !! Solve triangular packed system `A*x = b` or `A**T*x = b`.
         import :: dp
         character, intent(in) :: uplo
         character, intent(in) :: trans
         character, intent(in) :: diag
         integer, intent(in) :: n
         real(dp), intent(in) :: ap(*)
         real(dp), intent(inout) :: x(*)
         integer, intent(in) :: incx
      end subroutine dtpsv

      pure subroutine dtrmm(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)
      !! Triangular matrix-matrix operation `B = alpha*op(A)*B` or `B = alpha*B*op(A)`.
         import :: dp
         character, intent(in) :: side
         character, intent(in) :: uplo
         character, intent(in) :: transa
         character, intent(in) :: diag
         integer, intent(in) :: m
         integer, intent(in) :: n
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: a(lda, *)
         integer, intent(in) :: lda
         real(dp), intent(inout) :: b(ldb, *)
         integer, intent(in) :: ldb
      end subroutine dtrmm

      pure subroutine dtrsv(uplo, trans, diag, n, a, lda, x, incx)
      !! Solve triangular system `A*x = b` or `A**T*x = b`.
         import :: dp
         character, intent(in) :: uplo
         character, intent(in) :: trans
         character, intent(in) :: diag
         integer, intent(in) :: n
         real(dp), intent(in) :: a(lda, *)
         integer, intent(in) :: lda
         real(dp), intent(inout) :: x(*)
         integer, intent(in) :: incx
      end subroutine dtrsv

      pure integer function idamax(n, x, incx)
      !! Find largest component of vector.
         import :: dp
         integer, intent(in) :: n
         real(dp), intent(in) :: x(*)
         integer, intent(in) :: incx
      end function idamax

   end interface

end module blas_interfaces
