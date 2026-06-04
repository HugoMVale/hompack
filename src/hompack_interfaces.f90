module hompack_interfaces
! Interfaces for user written subroutines.

   use hompack_kinds, only: dp
   implicit none

   interface

      subroutine f(x, v)
         ! interface for subroutine that evaluates f(x) and returns it in the vector v.
         import :: dp
         real(dp), dimension(:), intent(in) :: x
         real(dp), dimension(:), intent(out) :: v
      end subroutine f

      subroutine fjac(x, v, k)
         ! interface for subroutine that returns in v the k-th column of the jacobian
         ! matrix of f(x) evaluated at x.
         import :: dp
         real(dp), dimension(:), intent(in) :: x
         real(dp), dimension(:), intent(out) :: v
         integer, intent(in) :: k
      end subroutine fjac

      subroutine rho(a, lambda, x, v)
         ! interface for subroutine that evaluates rho(a,lambda,x) and returns it
         ! in the vector v.
         import :: dp
         real(dp), intent(in) :: a(:), x(:)
         real(dp), intent(in out) :: lambda
         real(dp), intent(out) :: v(:)
      end subroutine rho

      ! the following code is specifically for the polynomial system driver
      ! polsys1h, and should be used verbatim with polsys1h in the external
      ! subroutine rho.
      !     use hompack90_global, only: ipar, par  ! for polsys1h only.
      !     interface
      !       subroutine hfunp(n,a,lambda,x)
      !       use real_precision
      !       integer, intent(in) :: n
      !       real(dp), intent(in) :: a(2*n), lambda, x(2*n)
      !       end subroutine hfunp
      !     end interface
      !     integer :: j, npol
      ! force predicted point to have  lambda .ge. 0  .
      !     if (lambda .lt. 0.0) lambda=0.0
      !     npol=ipar(1)
      !     call hfunp(npol,a,lambda,x)
      !     do j=1,2*npol
      !       v(j)=par(ipar(3 + (4-1)) + (j-1))
      !     end do
      !     return
      ! if calling fixp?? or step?? directly, supply appropriate replacement
      ! code in the external subroutine rho.

      subroutine rhoa(a, lambda, x)
         ! interface for subroutine that calculates and returns in a the vector
         ! z such that rho(z,lambda,x) = 0 .
         import :: dp
         real(dp), dimension(:), intent(out) :: a
         real(dp), intent(in) :: lambda, x(:)
      end subroutine rhoa

      subroutine rhojac(a, lambda, x, v, k)
         ! interface for subroutine that returns in the vector v the kth column
         ! of the jacobian matrix [d rho/d lambda, d rho/dx] evaluated at the
         ! point (a, lambda, x).
         import :: dp
         real(dp), intent(in) :: a(:), x(:)
         real(dp), intent(in out) :: lambda
         real(dp), intent(out) :: v(:)
         integer, intent(in) :: k
      end subroutine rhojac

      ! the following code is specifically for the polynomial system driver
      ! polsys1h, and should be used verbatim with polsys1h in the external
      ! subroutine rhojac.
      !     use hompack90_global, only: ipar, par  ! for polsys1h only.
      !     interface
      !       subroutine hfunp(n,a,lambda,x)
      !       use real_precision
      !       integer, intent(in) :: n
      !       real(dp), intent(in) :: a(2*n), lambda, x(2*n)
      !       end subroutine hfunp
      !     end interface
      !     integer :: j, npol, n2
      !     npol=ipar(1)
      !     n2=2*npol
      !     if (k .eq. 1) then
      ! force predicted point to have  lambda .ge. 0  .
      !       if (lambda .lt. 0.0) lambda=0.0
      !       call hfunp(npol,a,lambda,x)
      !       do j=1,n2
      !         v(j)=par(ipar(3 + (6-1)) + (j-1))
      !       end do
      !       return
      !     else
      !       do j=1,n2
      !         v(j)=par(ipar(3 + (5-1)) + (j-1) + n2*(k-2))
      !       end do
      !     endif
      !
      !     return
      ! if calling fixp?? or step?? directly, supply appropriate replacement
      ! code in the external subroutine rhojac.

      subroutine fjacs(x)
         ! interface for subroutine that evaluates a sparse jacobian matrix of
         ! f(x) at x, and operates as follows:
         !
         ! if mode = 1,
         ! evaluate the n x n symmetric jacobian matrix of f(x) at x, and return
         ! the result in packed skyline storage format in qrsparse.  lenqr is the
         ! length of qrsparse, and rowpos contains the indices of the diagonal
         ! elements of the jacobian matrix within qrsparse.  rowpos(n+1) and
         ! rowpos(n+2) are set by subroutine fodeds.  the allocatable array colpos
         ! is not used by this storage format.
         !
         ! if mode = 2,
         ! evaluate the n x n jacobian matrix of f(x) at x, and return the result
         ! in sparse row storage format in qrsparse.  lenqr is the length of
         ! qrsparse, rowpos contains the indices of where each row begins within
         ! qrsparse, and colpos (of length lenqr) contains the column indices of
         ! the corresponding elements in qrsparse.  even if zero, the diagonal
         ! elements of the jacobian matrix must be stored in qrsparse.
         use hompack_global, only: qrsparse, rowpos, colpos
         import :: dp
         real(dp), dimension(:), intent(in) :: x
      end subroutine fjacs

      subroutine rhojs(a, lambda, x)
         ! interface for subroutine that evaluates a sparse jacobian matrix of
         ! rho(a,x,lambda) at (a,x,lambda), and operates as follows:
         !
         ! if mode = 1,
         ! evaluate the n x n symmetric jacobian matrix [d rho/dx] at
         ! (a,x,lambda), and return the result in packed skyline storage format in
         ! qrsparse.  lenqr is the length of qrsparse, and rowpos contains the
         ! indices of the diagonal elements of [d rho/dx] within qrsparse.  pp
         ! contains -[d rho/d lambda] evaluated at (a,x,lambda).  note the minus
         ! sign in the definition of pp.  the allocatable array colpos is not used
         ! in this storage format.
         !
         ! if mode = 2,
         ! evaluate the n x (n+1) jacobian matrix [d rho/dx, d rho/dlambda] at
         ! (a,x,lambda), and return the result in sparse row storage format in
         ! qrsparse.  lenqr is the length of qrsparse, rowpos contains the indices
         ! of where each row begins within qrsparse, and colpos (of length lenqr)
         ! contains the column indices of the corresponding elements in qrsparse.
         ! even if zero, the diagonal elements of the jacobian matrix must be
         ! stored in qrsparse.  the allocatable array pp is not used in this
         ! storage format.
         use hompack_global, only: qrsparse, rowpos, colpos
         import :: dp
         real(dp), intent(in) :: a(:), lambda, x(:)
      end subroutine rhojs

   end interface

end module hompack_interfaces
