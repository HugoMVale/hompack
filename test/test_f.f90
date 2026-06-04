program test_f
!! Main program to test `fixpqf`, `fixpnf`, `fixpdf`, `stepnx`, and `rootnx`.
!! brown's function, zero finding.
!!
!! This program tests the hompack routines `fixpqf`, `fixpnf`, `fixpdf`,
!! `stepnx`, and `rootnx`.
!!
!! The output from this routine should be as follows, with the execution times
!! corresponding to a DEC AXP 3000/600.
!!
!!       TESTING FIXPQF
!!
!! LAMBDA = 1.00000000  FLAG = 1       6 JACOBIAN EVALUATIONS
!! EXECUTION TIME(SECS) =     0.106    ARCLEN =     2.693
!!   1.00000000E+00  1.00000000E+00  1.00000000E+00  1.00000000E+00
!!   1.00000000E+00
!!
!!
!!       TESTING FIXPNF
!!
!! LAMBDA = 1.00000000  FLAG = 1      19 JACOBIAN EVALUATIONS
!! EXECUTION TIME(SECS) =     0.005    ARCLEN =     2.676
!!   1.00000000E+00  1.00000000E+00  1.00000000E+00  1.00000000E+00
!!   1.00000000E+00
!!
!!
!!       TESTING FIXPDF
!!
!! LAMBDA = 1.00000000  FLAG = 1      71 JACOBIAN EVALUATIONS
!! EXECUTION TIME(SECS) =     0.016    ARCLEN =     2.712
!!   1.00000000E+00  1.00000000E+00  1.00000000E+00  1.00000000E+00
!!   1.00000000E+00
!!
!!
!!       TESTING STEPNX AND ROOTNX
!!
!! LAMBDA = 1.00000000  FLAG = -1      80 JACOBIAN EVALUATIONS
!! EXECUTION TIME(SECS) =     0.020    ARCLEN =     2.711
!!   1.00000000E+00  1.00000000E+00  1.00000000E+00  1.00000000E+00
!!   1.00000000E+00

   use hompack_kinds, only: dp
   use hompack, only: fixpdf, fixpnf, fixpqf
   implicit none

   integer, parameter :: n = 5, ndima = 5
   real(dp) :: a(n), ansae, ansre, arcae, arcre, &
               arclen, dtime, sspar(8), y(n + 1)
   integer :: iflag, ii, j, nfe, np1, timenew(8), timeold(8), trace
   character(len=6) :: name

   interface
      subroutine mainx
      end subroutine mainx
   end interface

   ! TEST EACH OF THE THREE ALGORITHMS
   do ii = 1, 3

      ! DEFINE ARGUMENTS FOR CALL TO HOMPACK PROCEDURE
      np1 = n + 1
      arcre = 0.5d-4
      arcae = 0.5d-4
      ansre = 1.0d-10
      ansae = 1.0d-10
      trace = 0
      sspar = 0d0
      iflag = -1
      y(2:np1) = 0d0

      ! GET CURRENT DATE AND TIME
      call date_and_time(values=timeold)

      ! CALL TO HOMPACK ROUTINE
      if (ii .eq. 1) then
         name = 'FIXPQF'
         call fixpqf(n, y, iflag, arcre, arcae, ansre, ansae, trace, a, &
                     sspar, nfe, arclen)
      else if (ii .eq. 2) then
         name = 'FIXPNF'
         call fixpnf(n, y, iflag, arcre, arcae, ansre, ansae, trace, a, &
                     sspar, nfe, arclen)
      else
         name = 'FIXPDF'
         call fixpdf(n, y, iflag, arcre, ansre, trace, a, ndima, nfe, arclen)
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

      write (6, 45) name
45    format(//, 7x, 'TESTING', 1x, a6)
      write (6, 50) y(1), iflag, nfe, dtime, arclen, (y(j), j=2, np1)
50    format(/' LAMBDA =', f11.8, '  FLAG =', i2, i8, ' JACOBIAN ', &
              'EVALUATIONS', /, 1x, 'EXECUTION TIME(MILIS) =', e10.2, 4x, &
              'ARCLEN =', f10.3/(1x, 4es16.8))
   end do

   ! TEST REVERSE CALL SUBROUTINES  STEPNX  AND  ROOTNX  ON THE SAME PROBLEM
   call mainx

end program test_f

! SAMPLE USER WRITTEN HOMOTOPY SUBROUTINES FOR TESTING FIXP*F.

subroutine f(x, v)
!********************************************************************
!
!      SUBROUTINE F(X,V) -- EVALUATES BROWN'S FUNCTION AT THE POINT
!         X, AND RETURNS THE VALUE IN V.
!
!********************************************************************
   use hompack_kinds, only: dp
   implicit none

   real(dp), intent(in) :: x(:)
   real(dp), intent(out) :: v(:)
   integer:: n

   n = size(x)
   v(1) = product(x) - 1d0
   v(2:n) = sum(x) - real(n + 1, kind=dp) + x(2:n)

end subroutine f

subroutine fjac(x, v, k)
!********************************************************************
!
!      SUBROUTINE FJAC(X,V,K)  --  EVALUATES THE K-TH COLUMN OF
!         THE JACOBIAN MATRIX FOR BROWN'S FUNCTION EVALUATED AT
!         THE POINT X, RETURNING THE VALUE IN V.
!
!********************************************************************
   use hompack_kinds, only: dp
   implicit none

   real(dp), intent(in) :: x(:)
   real(dp), intent(out) :: v(:)
   integer, intent(in) :: k
   integer :: j, n
   real(dp) :: prod

   n = size(x)
   prod = 1d0
   do j = 1, k - 1
      prod = prod*x(j)
   end do
   do j = k + 1, n
      prod = prod*x(j)
   end do
   v(1) = prod
   v(2:n) = 1d0
   if (k .gt. 1) v(k) = v(k) + 1d0

end subroutine fjac

! **********************************************************************
! THE REST OF THESE SUBROUTINES ARE NOT USED BY PROGRAM TEST_F, AND ARE
! INCLUDED HERE SIMPLY FOR COMPLETENESS AND AS TEMPLATES FOR THEIR USE.
! *********************************************************************

subroutine rho(a, lambda, x, v)
!! Evaluate `rho(a,lambda,x)` and return in the vector `v`.

   use hompack_kinds, only: dp, zero
   use hompack_core, only: hfunp
   use hompack_global, only: par, ipar
   implicit none

   real(dp), intent(in) :: a(:), x(:)
   real(dp), intent(inout) :: lambda
   real(dp), intent(out) :: v(:)

   integer:: j, npol

   ! THE FOLLOWING CODE IS SPECIFICALLY FOR THE POLYNOMIAL SYSTEM DRIVER
   !  POLSYS1H , AND SHOULD BE USED VERBATIM WITH  POLSYS1H .  IF THE USER
   ! CALLING  FIXP??  OR   STEP??  DIRECTLY, HE MUST SUPPLY APPROPRIATE
   ! REPLACEMENT CODE HERE.

   ! FORCE PREDICTED POINT TO HAVE  LAMBDA .GE. 0
   if (lambda .lt. zero) lambda = zero
   npol = ipar(1)
   ! CALL HFUNP(NPOL,A,LAMBDA,X)
   do j = 1, 2*npol
      v(j) = par(ipar(3 + (4 - 1)) + (j - 1))
   end do

end subroutine rho

subroutine rhoa(a, lambda, x)
!! Calculate and return in `a` the vector `z` such that `rho(z,lambda,x) = 0`.

   use hompack_kinds, only: dp
   implicit none

   real(dp), intent(out) :: a(:)
   real(dp), intent(in) :: lambda, x(:)

   a(1) = lambda

end subroutine rhoa

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

subroutine fjacs(x)
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

   use hompack_kinds, only: dp
   use hompack_global
   implicit none

   real(dp), intent(in) :: x(:)

end subroutine fjacs

subroutine rhojs(a, lambda, x)
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

   use hompack_kinds, only: dp
   use hompack_global
   implicit none

   real(dp), intent(in) :: a(:), lambda, x(:)

end subroutine rhojs

subroutine mainx
!!  Subroutine to test the reverse call subroutines `stepnx` and `rootnx`. The test
!!  problem is Brown's function, zero finding. The output is similar to that from the
!!  test of  `fixpnf`, except with more Jacobian evaluations since the undefined function
!!  option of `stepnx`  is used to force smaller steps.

   use hompack_kinds, only: dp, zero, one
   use hompack_core, only: rootnx, stepnx, tangnf
   implicit none

   integer, parameter :: n = 5, ndima = 5
   real(dp) :: a(ndima), abserr, alpha(3*n + 3), &
               ansae, ansre, arcae, arcre, arclen, dtime, gofw, h, hold, &
               qr(n, n + 2), relerr, rholen, s, sspar(8), tz(n + 1), w(n + 1), &
               wp(n + 1), y(n + 1), yold(n + 1), yolds(n + 1), yp(n + 1), ypold(n + 1)
   integer :: iflag, iter = 0, j, nfe, nfec = 0, np1, pivot(n + 1), &
              timenew(8), timeold(8), trace
   logical :: crash, start

   ! DEFINE ARGUMENTS FOR CALL TO HOMPACK PROCEDURE
   np1 = n + 1
   nfe = 0
   arcre = 0.5d-4
   arcae = 0.5d-4
   ansre = 1.0d-10
   ansae = 1.0d-10
   abserr = arcae; relerr = arcre
   trace = 0
   sspar = zero
   iflag = -1
   a = zero
   y(1:np1) = zero
   yp(1) = one; yp(2:np1) = zero
   yold = y; ypold = yp
   start = .true.
   crash = .false.
   hold = one
   h = 0.1_dp
   s = zero

   ! GET CURRENT DATE AND TIME
   call date_and_time(values=timeold)

   ! TRACK CURVE TILL LAMBDA > 1.0
   track: do while (y(1) < 1.0_dp)
      call stepnx(n, nfe, iflag, start, crash, hold, h, relerr, &
                  abserr, s, y, yp, yold, ypold, a, tz, w, wp, rholen, sspar)
      if (crash) cycle track
      select case (iflag)
      case (-2:0)
         if (trace .gt. 0) then
            iter = iter + 1
            write (trace, 11) iter, nfe, s, y(1), y(2:np1)
11          format(/' STEP', i5, 3x, 'NFE =', i5, 3x, 'ARC LENGTH =', &
                    f9.4, 3x, 'LAMBDA =', f7.4, 5x, 'X VECTOR:'/(1x, 6es12.4))
         end if
         cycle track
      case (-12:-10)   ! TANGENT VECTOR
         if (h > .1_dp) then
            iflag = iflag - 100
            cycle track
         end if
         rholen = 0d0
         call tangnf(rholen, w, wp, ypold, a, qr, alpha, tz, pivot, &
                     nfec, n, iflag)
      case (-32:-20)   ! TANGENT VECTOR AND NEWTON STEP
         if (h > .1_dp) then
            iflag = iflag - 100
            cycle track
         end if
         rholen = -1d0
         call tangnf(rholen, w, wp, ypold, a, qr, alpha, tz, pivot, &
                     nfec, n, iflag)
      case (4, 6, 7)
         write (6, 13) iflag
13       format(/' FATAL ERROR OCCURRED DURING TRACKING WITH', &
                 ' FLAG =', i2, //)
         stop
      end select
   end do track

   ! CLEAN UP WORKING STORAGE
   iflag = iflag - 40
   call stepnx(n, nfe, iflag, start, crash, hold, h, relerr, &
               abserr, s, y, yp, yold, ypold, a, tz, w, wp, rholen, sspar)

   ! SAVE  YOLD  FOR ARC LENGTH CALCULATION LATER
   yolds = yold

   ! FIND POINT ON HOMOTOPY ZERO CURVE SATISFYING G(Y(S)) = LAMBDA(S) - 1 = 0
   abserr = ansae
   relerr = ansre
   end_game: do
      call rootnx(n, nfe, iflag, relerr, abserr, y, yp, yold, &
                  ypold, a, gofw, tz, w, wp)
      select case (iflag)
      case (-42:-10)   ! G(W)
         gofw = w(1) - 1d0
      case (-52:-50)   ! TANGENT VECTOR AND NEWTON STEP
         rholen = -1d0
         call tangnf(rholen, w, wp, ypold, a, qr, alpha, tz, pivot, &
                     nfec, n, iflag)
      case (-2:0, 4, 6, 7)
         exit end_game
      end select
   end do end_game

   ! CALCULATE FINAL ARC LENGTH
   w = y - yolds
   arclen = s - hold + sqrt(dot_product(w, w))

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
!
   write (6, 45)
45 format(//, 7x, 'TESTING STEPNX AND ROOTNX')
   write (6, 50) y(1), iflag, nfe, dtime, arclen, (y(j), j=2, np1)
50 format(/' LAMBDA =', f11.8, '  FLAG =', i3, i8, ' JACOBIAN ', &
           'EVALUATIONS', /, 1x, 'EXECUTION TIME(MILIS) =', e10.2, 4x, &
           'ARCLEN =', f10.3/(1x, 4es16.8))

end subroutine mainx
