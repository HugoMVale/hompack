module test_f_mod
!! Module to hold the user-supplied function and Jacobian evaluation subroutines

   use iso_c_binding, only: c_ptr
   use hompack_kinds, only: dp
   implicit none

contains

   subroutine f2(x, v, data)

      use hompack_core_legacy, only: f ! # TEMPORARY
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: v(:)
      type(c_ptr), value :: data

      call f(x, v)

   end subroutine f2

   subroutine fjac2(x, v, k, data)

      use hompack_core_legacy, only: fjac ! # TEMPORARY
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: v(:)
      integer, intent(in) :: k
      type(c_ptr), value :: data

      call fjac(x, v, k)

   end subroutine fjac2

end module test_f_mod

program test_f
!! Main program to test `fixpqf`, `fixpnf`, `fixpdf`, `stepnx`, `stepnx`, and `rootnx`.
!!
!! Brown's function, zero finding.
!!
!! The output from this routine should be as follows, with the execution times
!! corresponding to a DEC AXP 3000/600.
!!
!! ```
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
!! ```

   use hompack_kinds, only: dp
   use hompack_df, only: fixpdf
   use hompack_qf, only: fixpqf
   use hompack_nf, only: fixpnf, hompack_callbacks, fixnpf_state
   use hompack_core_legacy, only: f, fjac
   use test_f_mod, only: f2, fjac2
   implicit none

   integer, parameter :: n = 5, ndima = 5
   real(dp) :: a(n), ansae, ansre, arcae, arcre, arclen, dtime, sspar(8), y(n + 1), x(n)
   integer :: iflag, ii, j, nfe, np1, timenew(8), timeold(8), trace
   character(len=6) :: name
   type(hompack_callbacks) :: callbacks
   type(fixnpf_state) :: state

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
      x = y(2:np1)

      ! GET CURRENT DATE AND TIME
      call date_and_time(values=timeold)

      ! CALL TO HOMPACK ROUTINE
      if (ii .eq. 1) then
         name = 'FIXPQF'
         call fixpqf(n, y, iflag, arcre, arcae, ansre, ansae, trace, a, &
                     sspar, nfe, arclen)
      else if (ii .eq. 2) then
         name = 'FIXPNF'
         callbacks%f => f2
         callbacks%fjac => fjac2
         call fixpnf(state, callbacks, n, x, iflag, arcre, arcae, ansre, ansae, &
                     sspar=sspar, a=a, lunit=trace)
         nfe = state%nfe
         arclen = state%s
         y = state%y
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

! *************************************************************
! THE REST OF THESE SUBROUTINES ARE NOT USED BY PROGRAM TEST_F
! *************************************************************

subroutine rho(a, lambda, x, v)

   use hompack_kinds, only: dp
   implicit none

   real(dp), intent(in) :: a(:), x(:)
   real(dp), intent(inout) :: lambda
   real(dp), intent(out) :: v(:)

end subroutine rho

subroutine rhoa(a, lambda, x)

   use hompack_kinds, only: dp
   implicit none

   real(dp), intent(out) :: a(:)
   real(dp), intent(in) :: lambda, x(:)

   a(1) = lambda

end subroutine rhoa

subroutine rhojac(a, lambda, x, v, k)

   use hompack_kinds, only: dp
   implicit none

   real(dp), intent(in) :: a(:), x(:)
   real(dp), intent(inout) :: lambda
   real(dp), intent(out) :: v(:)
   integer, intent(in) :: k

end subroutine rhojac

subroutine fjacs(x)

   use hompack_kinds, only: dp
   use hompack_global_legacy
   implicit none

   real(dp), intent(in) :: x(:)

end subroutine fjacs

subroutine rhojs(a, lambda, x)

   use hompack_kinds, only: dp
   use hompack_global_legacy
   implicit none

   real(dp), intent(in) :: a(:), lambda, x(:)

end subroutine rhojs

subroutine mainx
!!  Subroutine to test the reverse call subroutines `stepnx` and `rootnx`. The test
!!  problem is Brown's function, zero finding. The output is similar to that from the
!!  test of  `fixpnf`, except with more Jacobian evaluations since the undefined function
!!  option of `stepnx`  is used to force smaller steps.

   use hompack_kinds, only: dp, zero, one
   use hompack_core_legacy, only: rootnx, stepnx
   use hompack_core_legacy, only: f, fjac
   use hompack_nf, only: hompack_callbacks, tangnf
   use test_f_mod, only: f2, fjac2
   implicit none

   integer, parameter :: n = 5, ndima = 5
   real(dp) :: a(ndima), abserr, alpha(3*n + 3), &
               ansae, ansre, arcae, arcre, arclen, dtime, gofw, h, hold, &
               qr(n, n + 2), relerr, rholen, s, sspar(8), tz(n + 1), w(n + 1), &
               wp(n + 1), y(n + 1), yold(n + 1), yolds(n + 1), yp(n + 1), ypold(n + 1)
   integer :: iflag, iter = 0, j, nfe, nfec = 0, np1, pivot(n + 1), &
              timenew(8), timeold(8), trace
   logical :: crash, start
   type(hompack_callbacks) :: callbacks

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

   ! CALLBACKS
   callbacks%f => f2
   callbacks%fjac => fjac2

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
         call tangnf(callbacks, rholen, n, w, wp, ypold, a, nfec, iflag, &
                     qr, alpha, tz, pivot)
      case (-32:-20)   ! TANGENT VECTOR AND NEWTON STEP
         if (h > .1_dp) then
            iflag = iflag - 100
            cycle track
         end if
         rholen = -1d0
         call tangnf(callbacks, rholen, n, w, wp, ypold, a, nfec, iflag, &
                     qr, alpha, tz, pivot)
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
         call tangnf(callbacks, rholen, n, w, wp, ypold, a, nfec, iflag, &
                     qr, alpha, tz, pivot)
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
