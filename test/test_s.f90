module switch
!  ROWSET  IS USED TO INITIALIZE SPARSE MATRIX DATA STRUCTURES ONLY
!  ONCE, AFTER THEY ARE ALLOCATED.

   logical::  rowset
end module switch

program test_s
!! Main program to test `fixpqs`, `fixpns`, and `fixpds.
!!
!! The output from this routine should be as follows, with the execution times
!! corresponding to a DEC AXP 3000/600.
!!
!! ```
!!       TESTING FIXPQS WITH STORAGE MODE = 1
!!
!! LAMBDA = 1.00000000  FLAG = 1      33 JACOBIAN EVALUATIONS
!! EXECUTION TIME(SECS) =     0.119    ARC LENGTH =     1.274
!!   4.00864019E-01  2.65454893E-01  8.40421103E-02  4.83042527E-01
!!   3.01797132E-01  2.32508994E-01  4.96639853E-01  3.00908894E-01
!!
!!       TESTING FIXPNS WITH STORAGE MODE = 1
!!
!! LAMBDA = 1.00000000  FLAG = 1      20 JACOBIAN EVALUATIONS
!! EXECUTION TIME(SECS) =     0.013    ARC LENGTH =     1.275
!!   4.00864019E-01  2.65454893E-01  8.40421103E-02  4.83042527E-01
!!   3.01797132E-01  2.32508994E-01  4.96639853E-01  3.00908894E-01
!!
!!       TESTING FIXPDS WITH STORAGE MODE = 1
!!
!! LAMBDA = 1.00000000  FLAG = 1      70 JACOBIAN EVALUATIONS
!! EXECUTION TIME(SECS) =     0.031    ARC LENGTH =     1.281
!!   4.00864019E-01  2.65454893E-01  8.40421103E-02  4.83042527E-01
!!   3.01797132E-01  2.32508994E-01  4.96639853E-01  3.00908894E-01
!!
!!       TESTING FIXPQS WITH STORAGE MODE = 2
!!
!! LAMBDA = 1.00000000  FLAG = 1      33 JACOBIAN EVALUATIONS
!! EXECUTION TIME(SECS) =     0.015    ARC LENGTH =     1.274
!!   4.00864019E-01  2.65454893E-01  8.40421103E-02  4.83042527E-01
!!   3.01797132E-01  2.32508994E-01  4.96639853E-01  3.00908894E-01
!!
!!       TESTING FIXPNS WITH STORAGE MODE = 2
!!
!! LAMBDA = 1.00000000  FLAG = 1      20 JACOBIAN EVALUATIONS
!! EXECUTION TIME(SECS) =     0.011    ARC LENGTH =     1.275
!!   4.00864019E-01  2.65454893E-01  8.40421103E-02  4.83042527E-01
!!   3.01797132E-01  2.32508994E-01  4.96639853E-01  3.00908894E-01
!!
!!       TESTING FIXPDS WITH STORAGE MODE = 2
!!
!! LAMBDA = 1.00000000  FLAG = 1      70 JACOBIAN EVALUATIONS
!! EXECUTION TIME(SECS) =     0.022    ARC LENGTH =     1.281
!!   4.00864019E-01  2.65454893E-01  8.40421103E-02  4.83042527E-01
!!   3.01797132E-01  2.32508994E-01  4.96639853E-01  3.00908894E-01
!! ```

   use switch
   use hompack_kinds, only: dp, zero, one
   use hompack, only: fixpds, fixpns, fixpqs
   implicit none

   integer, parameter :: n = 8, ndima = 8
   real(dp) :: a(n), ansae, ansre, arcae, arcre, arclen, dtime, sspar(8), y(n + 1)
   integer :: iflag, ii, j, lenqr, mode, nfe, np1, timenew(8), timeold(8), trace
   character(len=6) :: name

   ! TEST EACH OF THE THREE ALGORITHMS WITH BOTH STORAGE MODES
   do mode = 1, 2

      select case (mode)
      case (1)
         lenqr = 18
      case (2)
         lenqr = 36
      end select

      do ii = 1, 3

         ! DEFINE ARGUMENTS FOR CALL TO HOMPACK PROCEDURE
         np1 = n + 1
         arcre = 0.5d-4
         arcae = 0.5d-4
         ansre = 1.0d-12
         ansae = 1.0d-12
         trace = 0
         sspar = zero
         iflag = -mode
         y(1:n) = 0.5_dp
         if (iflag .eq. -2) a = y(1:n)
         rowset = .false.

         ! GET CURRENT DATE AND TIME
         call date_and_time(values=timeold)

         ! CALL TO HOMPACK ROUTINE
         if (ii .eq. 1) then
            name = 'FIXPQS'
            call fixpqs(n, y, iflag, arcre, arcae, ansre, ansae, trace, a, &
                        nfe, arclen, mode, lenqr, sspar)
         else if (ii .eq. 2) then
            name = 'FIXPNS'
            call fixpns(n, y, iflag, arcre, arcae, ansre, ansae, trace, a, &
                        nfe, arclen, mode, lenqr, sspar)
         else
            name = 'FIXPDS'
            call fixpds(n, y, iflag, arcre, ansre, trace, a, ndima, nfe, &
                        arclen, mode, lenqr)
         end if

         ! CALCULATE EXECUTION TIME
         call date_and_time(values=timenew)

         if (timenew(8) .lt. timeold(8)) then
            timenew(8) = timenew(8) + 1000
            timenew(7) = timenew(7) - 1
         end if

         if (timenew(7) .lt. timeold(7)) then
            timenew(7) = timenew(7) + 60
            timenew(6) = timenew(6) - 1
         end if

         if (timenew(6) .lt. timeold(6)) then
            timenew(6) = timenew(6) + 60
            timenew(5) = timenew(5) - 1
         end if

         if (timenew(5) .lt. timeold(5)) timenew(5) = timenew(5) + 24
         dtime = dot_product(timenew(5:8) - timeold(5:8), (/3600000, 60000, 1000, 1/))

         write (6, 45) name, mode
45       format(//, 7x, 'TESTING', 1x, a6, ' WITH STORAGE MODE =', i2)
         write (6, 50) y(np1), iflag, nfe, dtime, arclen, (y(j), j=1, n)
50       format(/' LAMBDA =', f11.8, '  FLAG =', i2, i8, ' JACOBIAN ', &
                 'EVALUATIONS', /, 1x, 'EXECUTION TIME(MILIS) =', e10.2, 4x, &
                 'ARC LENGTH =', f10.3/(1x, 4es16.8))
      end do
   end do

end program test_s

subroutine f(x, v)
!! Evaluate `f(x)` and return in the vector `v`.

   use hompack_kinds, only: dp
   implicit none

   real(dp), intent(in) :: x(:)
   real(dp), intent(out) :: v(:)

   v(1) = x(1)**3 + 6*x(2)*x(3) - 1d0 + 2*x(1)
   v(2) = 6*x(1)*x(3) + x(2)**4*x(5) - 1d0 + 3*x(2)
   v(3) = 6*x(1)*x(2) + x(3)*x(5) - 1d0 + 4*x(3)
   v(4) = x(4)**3*x(8) - 1d0 + 2*x(4)
   v(5) = (x(2)**5)/5 + (x(3)**2)/2 + x(8)*x(5) - 1d0 + 3*x(5)
   v(6) = x(6)*x(8) - 1d0 + 4*x(6)
   v(7) = x(7)**2*x(8)**3 - 1d0 + 2*x(7)
   v(8) = (x(4)**4)/4 + (x(5)**2)/2 + (x(6)**2)/2 + x(7)**3* &
          x(8)**2 - 1d0 + 3*x(8)

end subroutine f

subroutine fjacs(x)
!! Compute the jacobian matrix of `f` at the point `x`, returning the jacobian matrix in
!! packed skyline form (`mode=1`) in the arrays `qrsparse` and `rowpos`.

   use hompack_kinds, only: dp
   use hompack_global, only: qrsparse, rowpos
   use switch
   implicit none

   real(dp), intent(in) :: x(:)

   integer:: n

   ! If MODE = 1,
   ! evaluate the N x N symmetric Jacobian matrix of F(X) at X, and return
   ! the result in packed skyline storage format in QRSPARSE.  LENQR is the
   ! length of QRSPARSE, and ROWPOS contains the indices of the diagonal
   ! elements of the Jacobian matrix within QRSPARSE.  ROWPOS(N+1) and
   ! ROWPOS(N+2) are set by subroutine FODEDS. The allocatable array COLPOS
   ! is not used by this storage format.
   !
   ! If MODE = 2,
   ! evaluate the N x N Jacobian matrix of F(X) at X, and return the result
   ! in sparse row storage format in QRSPARSE.  LENQR is the length of
   ! QRSPARSE, ROWPOS contains the indices of where each row begins within
   ! QRSPARSE, and COLPOS (of length LENQR) contains the column indices of
   ! the corresponding elements in QRSPARSE.  Even if zero, the diagonal
   ! elements of the Jacobian matrix must be stored in QRSPARSE.

   n = size(x)
   if (.not. rowset) then
      rowset = .true.
      rowpos(1:n + 1) = (/1, 2, 4, 7, 8, 12, 13, 14, 19/)
   end if

   qrsparse(1) = 3*x(1)**2 + 2d0
   qrsparse(2) = 4*x(2)**3*x(5) + 3d0
   qrsparse(3) = 6*x(3)
   qrsparse(4) = x(5) + 4d0
   qrsparse(5) = 6*x(1)
   qrsparse(6) = 6*x(2)
   qrsparse(7) = 3*x(4)**2*x(8) + 2d0
   qrsparse(8) = x(8) + 3d0
   qrsparse(9) = 0d0
   qrsparse(10) = x(3)
   qrsparse(11) = x(2)**4
   qrsparse(12) = x(8) + 4d0
   qrsparse(13) = 2*x(7)*x(8)**3 + 2d0
   qrsparse(14) = 2*x(7)**3*x(8) + 3d0
   qrsparse(15) = 3*x(7)**2*x(8)**2
   qrsparse(16) = x(6)
   qrsparse(17) = x(5)
   qrsparse(18) = x(4)**3

end subroutine fjacs

subroutine rho(a, lambda, x, v)
!! Evaluate `rho(a,lambda,x)` and return in the vector `v`.
   use hompack_kinds, only: dp
   use hompack_interfaces, only: f
   implicit none

   real(dp), intent(in):: a(:), x(:)
   real(dp), intent(in out):: lambda
   real(dp), intent(out):: v(:)

   integer:: n
   n = size(x)
   call f(x(1:n), v(1:n))
   v(1:n) = lambda*v(1:n) + (1d0 - lambda)*(x(1:n) - a(1:n))

end subroutine rho

subroutine rhoa(a, lambda, x)
!! Calculate and return in `a` the vector `z` such that `rho(z,lambda,x) = 0`.

   use hompack_kinds, only: dp
   use hompack_interfaces, only: f
   implicit none

   real(dp), intent(out) :: a(:)
   real(dp), intent(in) :: lambda, x(:)

   ! N=NDIMA FOR THIS TEST PROBLEM.
   integer:: n
   n = size(x)
   call f(x(1:n), a(1:n))
   a(1:n) = lambda*a(1:n)/(1d0 - lambda) + x(1:n)

end subroutine rhoa

subroutine rhojs(a, lambda, x)
!! Computes the Jacobian matrix of `rho(a,x,lambda) = lambda*f(x) + (1 - lambda)*(x - a)`
!! at the point `(a,x,lambda)`, returning the Jacobian matrix in sparse row storage format
!! (`mode=2`) in the arrays `qrsparse`, `rowpos`, and `colpos`.

   use hompack_kinds, only: dp, zero
   use hompack_interfaces, only: f
   use hompack_global, only: qrsparse, rowpos, colpos
   use switch
   implicit none

   real(dp), intent(in) :: a(:), lambda, x(:)

   integer, parameter :: n = 8
   integer :: j, jpos, elem(n) = (/4, 9, 14, 18, 24, 27, 30, 36/)
   real(dp) :: drhodl(n)

   ! If MODE = 1,
   ! evaluate the N x N symmetric Jacobian matrix of F(X) at X, and return
   ! the result in packed skyline storage format in QRSPARSE.  LENQR is the
   ! length of QRSPARSE, and ROWPOS contains the indices of the diagonal
   ! elements of the Jacobian matrix within QRSPARSE.  ROWPOS(N+1) and
   ! ROWPOS(N+2) are set by subroutine FODEDS.  The allocatable array COLPO
   ! is not used by this storage format.
   !
   ! If MODE = 2,
   ! evaluate the N x N Jacobian matrix of F(X) at X, and return the result
   ! in sparse row storage format in QRSPARSE.  LENQR is the length of
   ! QRSPARSE, ROWPOS contains the indices of where each row begins within
   ! QRSPARSE, and COLPOS (of length LENQR) contains the column indices of
   ! the corresponding elements in QRSPARSE.  Even if zero, the diagonal
   ! elements of the Jacobian matrix must be stored in QRSPARSE.
   !
   !---------------------------------------------------------------------
   !    [QRSPARSE] = [ LAMBDA*DF(X) + (1 - LAMBDA)*I | F(X) - X + A ] .
   !---------------------------------------------------------------------

   if (.not. rowset) then
      rowset = .true.
      rowpos(1:n + 1) = (/1, 5, 10, 15, 19, 25, 28, 31, 37/)
      colpos(1:36) = (/1, 2, 3, 9, 1, 2, 3, 5, 9, 1, 2, 3, 5, 9, 4, 5, 8, 9, &
                       2, 3, 4, 5, 8, 9, 6, 8, 9, 7, 8, 9, 4, 5, 6, 7, 8, 9/)
   end if

   qrsparse = zero
   ! ROW 1.
   qrsparse(1) = 3*x(1)**2 + 2d0
   qrsparse(2) = 6*x(3)
   qrsparse(3) = 6*x(2)
   ! ROW 2.
   qrsparse(5) = 6*x(3)
   qrsparse(6) = 4*x(2)**3*x(5) + 3d0
   qrsparse(7) = 6*x(1)
   qrsparse(8) = x(2)**4
   ! ROW 3.
   qrsparse(10) = 6*x(2)
   qrsparse(11) = 6*x(1)
   qrsparse(12) = x(5) + 4d0
   qrsparse(13) = x(3)
   ! ROW 4.
   qrsparse(15) = 3*x(4)**2*x(8) + 2d0
   qrsparse(16) = zero
   qrsparse(17) = x(4)**3
   ! ROW 5.
   qrsparse(19) = x(2)**4
   qrsparse(20) = x(3)
   qrsparse(21) = zero
   qrsparse(22) = x(8) + 3d0
   qrsparse(23) = x(5)
   ! ROW 6.
   qrsparse(25) = x(8) + 4d0
   qrsparse(26) = x(6)
   colpos(25) = 6
   colpos(26) = 8
   colpos(27) = 9
   ! ROW 7.
   qrsparse(28) = 2*x(7)*x(8)**3 + 2d0
   qrsparse(29) = 3*x(7)**2*x(8)**2
   ! ROW 8.
   qrsparse(31) = x(4)**3
   qrsparse(32) = x(5)
   qrsparse(33) = x(6)
   qrsparse(34) = 3*x(7)**2*x(8)**2
   qrsparse(35) = 2*x(7)**3*x(8) + 3d0

   qrsparse = lambda*qrsparse

   ! FIND INDEX JPOS OF DIAGONAL ELEMENT IN JTH ROW OF QR
   do j = 1, n
      jpos = rowpos(j)
      do
         if (colpos(jpos) .eq. j) exit
         jpos = jpos + 1
      end do
      qrsparse(jpos) = qrsparse(jpos) + 1d0 - lambda
   end do

   ! INITIALIZE (N+1)ST COLUMN.
   call f(x(1:n), drhodl(1:n))
   drhodl = drhodl - x(1:n) + a(1:n)
   qrsparse(elem) = drhodl(1:8)

end subroutine rhojs

! THE REST OF THESE SUBROUTINES ARE NOT USED BY PROGRAM TEST_S, AND ARE
! INCLUDED HERE SIMPLY FOR COMPLETENESS AND AS TEMPLATES FOR THEIR USE.

subroutine fjac(x, v, k)
!! Return in `v` the `k`-th column of the Jacobian matrix of `f(x)` evaluated at `x`.

   use hompack_kinds, only: dp
   implicit none

   real(dp), intent(in) :: x(:)
   real(dp), intent(out) :: v(:)
   integer, intent(in) :: k

   v(1) = x(1)

end subroutine fjac

subroutine rhojac(a, lambda, x, v, k)
!! Return in the vector `v` the `k`-th column of the jacobian matrix
!! `[d rho / d lambda, d rho / d x]` evaluated at the point `(a, lambda, x)`.

   use hompack_kinds, only: dp, zero
   use hompack_core, only: hfunp
   use hompack_global, only: par, ipar
   implicit none

   real(dp), intent(in) :: a(:), x(:)
   real(dp), intent(inout) :: lambda
   real(dp), intent(out) :: v(:)
   integer, intent(in) :: k

   integer:: j, npol, n2

   ! THE FOLLOWING CODE IS SPECIFICALLY FOR THE POLYNOMIAL SYSTEM DRIVER
   !  POLSYS1H , AND SHOULD BE USED VERBATIM WITH  POLSYS1H .  IF THE USER
   ! CALLING  FIXP??  OR   STEP??  DIRECTLY, HE MUST SUPPLY APPROPRIATE
   ! REPLACEMENT CODE HERE.

   npol = ipar(1)
   n2 = 2*npol
   if (k .eq. 1) then
      ! FORCE PREDICTED POINT TO HAVE  LAMBDA .GE. 0  .
      if (lambda .lt. zero) lambda = zero
      ! CALL HFUNP(NPOL,A,LAMBDA,X)
      do j = 1, n2
         v(j) = par(ipar(3 + (6 - 1)) + (j - 1))
      end do
      return
   else
      do j = 1, n2
         v(j) = par(ipar(3 + (5 - 1)) + (j - 1) + n2*(k - 2))
      end do
   end if

end subroutine rhojac
