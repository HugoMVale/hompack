program test_p
!! Main program to test `polsys1h`.
!!
!! This routine requires one input file, `innhp.dat`.
!!
!! A sample input file and associated output are given in the comments that follow. This
!! sample problem is cited in the hompack report.
!!
!! *Sample input data (read from the file 'innhp.dat')*
!!
!! ```
!! ' TWO QUADRICS, NO SOLUTIONS AT INFINITY, TWO REAL SOLUTIONS.'
!! &PROBLEM
!!       IFLGHM = 1
!!       IFLGSC = 1
!!       TOTDG  = 4
!!       MAXT = 6
!!       EPSBIG = 1.D-04
!!       EPSSML = 1.D-14
!!       SSPAR(5) = 1.D+00
!!       NUMRR = 10
!!       N = 2 /
!! 00006                     NUMTRM(1)
!! 00002                     DEG(1,1,1)
!! 00000                     DEG(1,2,1)
!!            -.00098D+00
!! 00000                     DEG(1,1,2)
!! 00002                     DEG(1,2,2)
!!            978000.D+00
!! 00001                     DEG(1,1,3)
!! 00001                     DEG(1,2,3)
!!               -9.8D+00
!! 00001                     DEG(1,1,4)
!! 00000                     DEG(1,2,4)
!!             -235.0D+00
!! 00000                     DEG(1,1,5)
!! 00001                     DEG(1,2,5)
!!            88900.0D+00
!! 00000                     DEG(1,1,6)
!! 00000                     DEG(1,2,6)
!!             -1.000D+00
!! 00006                     NUMTRM(2)
!! 00002                     DEG(2,1,1)
!! 00000                     DEG(2,2,1)
!!             -.0100D+00
!! 00000                     DEG(2,1,2)
!! 00002                     DEG(2,2,2)
!!             -.9840D+00
!! 00001                     DEG(2,1,3)
!! 00001                     DEG(2,2,3)
!!             -29.70D+00
!! 00001                     DEG(2,1,4)
!! 00000                     DEG(2,2,4)
!!             .00987D+00
!! 00000                     DEG(2,1,5)
!! 00001                     DEG(2,2,5)
!!             -.1240D+00
!! 00000                     DEG(2,1,6)
!! 00000                     DEG(2,2,6)
!!             -.2500D+00
!! ```
!!
!! *Associated sample output (written to the file 'outhp.dat')*
!!
!! ```
!!     POLSYS1H TEST ROUTINE 7/7/95
!!
!! TWO QUADRICS, NO SOLUTIONS AT INFINITY, TWO REAL SOLUTIONS.
!!
!! IF IFLGHM=1, HOMOGENEOUS; IF IFLGHM=0, INHOMOGENEOUS; IFLGHM= 1
!!
!! IF IFLGSC=1, SCLGNP USED; IF IFLGSC=0, NO SCALING; IFLGSC=    1
!!
!! TOTDG=    4          MAXT=    6
!!
!! EPSBIG, EPSSML =  1.000000000000000E-04  1.000000000000000E-14
!!
!! SSPAR(5) =  1.000000000000000E+00
!!
!! NUMBER OF EQUATIONS =    2
!!
!! NUMBER OF RECALLS WHEN IFLAG=3:   10
!!
!!
!!  ****** COEFFICIENT TABLEAU ******
!!
!!  NUMT( 1) =    6
!!  KDEG( 1, 1, 1) =    2
!!  KDEG( 1, 2, 1) =    0
!!  COEF( 1, 1) =  -9.800000000000000E-04
!!  KDEG( 1, 1, 2) =    0
!!  KDEG( 1, 2, 2) =    2
!!  COEF( 1, 2) =   9.780000000000000E+05
!!  KDEG( 1, 1, 3) =    1
!!  KDEG( 1, 2, 3) =    1
!!  COEF( 1, 3) =  -9.800000000000001E+00
!!  KDEG( 1, 1, 4) =    1
!!  KDEG( 1, 2, 4) =    0
!!  COEF( 1, 4) =  -2.350000000000000E+02
!!  KDEG( 1, 1, 5) =    0
!!  KDEG( 1, 2, 5) =    1
!!  COEF( 1, 5) =   8.890000000000000E+04
!!  KDEG( 1, 1, 6) =    0
!!  KDEG( 1, 2, 6) =    0
!!  COEF( 1, 6) =  -1.000000000000000E+00
!!
!!  NUMT( 2) =    6
!!  KDEG( 2, 1, 1) =    2
!!  KDEG( 2, 2, 1) =    0
!!  COEF( 2, 1) =  -1.000000000000000E-02
!!  KDEG( 2, 1, 2) =    0
!!  KDEG( 2, 2, 2) =    2
!!  COEF( 2, 2) =  -9.840000000000000E-01
!!  KDEG( 2, 1, 3) =    1
!!  KDEG( 2, 2, 3) =    1
!!  COEF( 2, 3) =  -2.970000000000000E+01
!!  KDEG( 2, 1, 4) =    1
!!  KDEG( 2, 2, 4) =    0
!!  COEF( 2, 4) =   9.870000000000000E-03
!!  KDEG( 2, 1, 5) =    0
!!  KDEG( 2, 2, 5) =    1
!!  COEF( 2, 5) =  -1.240000000000000E-01
!!  KDEG( 2, 1, 6) =    0
!!  KDEG( 2, 2, 6) =    0
!!  COEF( 2, 6) =  -2.500000000000000E-01
!!
!!
!!
!!  IFLG1 =   11
!!
!!  PATH NUMBER =    1
!!
!!  FINAL VALUES FOR PATH
!!
!!  ARCLEN =   1.005533190562901E+01
!!  NFE =   53
!!  IFLG2 =  1
!! REAL, FINITE SOLUTION
!!  LAMBDA =  1.000000000000003E+00
!! X( 1) =  2.342338519591276E+03  8.841149143431121E-13
!! X( 2) = -7.883448240941412E-01 -9.356862757018485E-16
!!
!! X( 3) = -9.493594594086552E-03 -1.064475509002627E-03
!!
!!
!!  PATH NUMBER =    2
!!
!!  FINAL VALUES FOR PATH
!!
!!  ARCLEN =   1.721129286057142E+00
!!  NFE =   37
!!  IFLG2 =  1
!! COMPLEX, FINITE SOLUTION
!!  LAMBDA =  1.000000000000006E+00
!! X( 1) =  1.614785792344189E-02  1.684969554988811E+00
!! X( 2) =  2.679947396144760E-04  4.428029939736605E-03
!!
!! X( 3) = -3.819489729424030E-01  3.720689434572830E-01
!!
!!
!!  PATH NUMBER =    3
!!
!!  FINAL VALUES FOR PATH
!!
!!  ARCLEN =   2.023295279367267E+00
!!  NFE =   35
!!  IFLG2 =  1
!! COMPLEX, FINITE SOLUTION
!!  LAMBDA =  1.000000000000000E+00
!! X( 1) =  1.614785792343521E-02 -1.684969554988812E+00
!! X( 2) =  2.679947396144598E-04 -4.428029939736611E-03
!!
!! X( 3) = -3.293704938476598E-01  5.566197755230126E-01
!!
!!
!!  PATH NUMBER =    4
!!
!!  FINAL VALUES FOR PATH
!!
!!  ARCLEN =   4.163266156958467E+00
!!  NFE =   46
!!  IFLG2 =  1
!! REAL, FINITE SOLUTION
!!  LAMBDA =  9.999999999999998E-01
!! X( 1) =  9.089212296153869E-02  1.153793567884107E-16
!! X( 2) = -9.114970981974997E-02  1.887399041592030E-17
!!
!! X( 3) = -5.736733957279616E-02  1.362436637092185E-01
!!
!!
!! TOTAL NFE OVER ALL PATHS =        171
!!
!! ```
!!
!! ```
!!  PROGRAM DESCRIPTION:  1. READS IN DATA.
!!                        2. GENERATES POLSYS1H INPUT.
!!                        3. CALLS POLSYS1H.
!!                        4. WRITES POLSYS1H OUTPUT.
!!
!! DIMENSIONS SHOULD BE SET AS FOLLOWS:
!!
!!     DIMENSION NUMT(NN),COEF(NN,MMAXT),KDEG(NN,NN+1,MMAXT)
!!     DIMENSION IFLG2(TTOTDG)
!!     DIMENSION LAMBDA(TTOTDG),ROOTS(2,NN+1,TTOTDG),ARCLEN(TTOTDG),
!!    & NFE(TTOTDG)
!!
!! WHERE:
!!    N   IS THE NUMBER OF EQUATIONS.  NN .GE. N.
!!    MAXT  IS THE MAXIMUM NUMBER OF TERMS IN ANY ONE EQUATION.
!!       MMAXT  .GE.  MAXT.
!!    TOTDG  IS THE TOTAL DEGREE OF THE SYSTEM.  TTOTDG .GE. TOTDG.
!!
!! THIS TEST CODE HAS DIMENSIONS SET AS FOLLOWS:
!!
!! NN=10, MMAXT=30, TTOTDG=1024
!! ```

   use hompack_kinds, only: dp, zero
   use hompack, only: polsys1h
   implicit none

   integer, parameter :: nn = 10, mmaxt = 30, ttotdg = 1024
   integer :: iflg1, iflg2(ttotdg), iflghm, iflgsc, itotit, j, k, &
              kdeg(nn, nn + 1, mmaxt), l, m, maxt, n, nfe(ttotdg), np1, nt, &
              numrr, numt(nn), totdg
   real(dp) :: arclen(ttotdg), coef(nn, mmaxt), epsbig, epssml, &
               lambda(ttotdg), roots(2, nn + 1, ttotdg), sspar(8)
   character(len=72) :: title

   namelist /problem/ iflghm, iflgsc, totdg, maxt, n, epsbig, epssml, sspar, numrr

   open (unit=7, file='./test/innhp.dat', action='READ', &
         position='REWIND', delim='APOSTROPHE', status='OLD')
   open (unit=6, file='./test/outhp.dat', action='WRITE', &
         delim='APOSTROPHE', status='REPLACE')

   sspar(1:8) = zero
   read (7, *) title
   write (6, 10) title
10 format(5x, 'POLSYS1H TEST ROUTINE 7/7/95', //, a72)

   read (7, nml=problem)

   write (6, 100) iflghm
100 format(/ &
      ' IF IFLGHM=1, HOMOGENEOUS; IF IFLGHM=0, INHOMOGENEOUS; IFLGHM=' &
      , i2)
   write (6, 102) iflgsc
102 format(/ &
      ' IF IFLGSC=1, SCLGNP USED; IF IFLGSC=0, NO SCALING; IFLGSC=', i5)
   write (6, 104) totdg, maxt
104 format(/, ' TOTDG=', i5, 10x, 'MAXT=', i5)
   write (6, 106) epsbig, epssml, sspar(5), n, numrr
106 format(/, ' EPSBIG, EPSSML =', 2es22.14, &
           //, ' SSPAR(5) =', es22.14, &
           //, ' NUMBER OF EQUATIONS =', i5, &
           //, ' NUMBER OF RECALLS WHEN IFLAG=3:', i5)

   np1 = n + 1

   ! NOTE THAT THE DEGREES OF VARIABLES IN EACH TERM OF EACH EQUATION
   ! ARE DEFINED BY THE FOLLOWING INDEXING SCHEME:
   !
   !     KDEG(J,  L,  K)
   !
   !          ^   ^   ^
   !
   !          E   V   T
   !          Q   A   E
   !          U   R   R
   !          A   I   M
   !          T   A
   !          I   B
   !          O   L
   !          N   E

   write (6, 200)
200 format(//, '  ****** COEFFICIENT TABLEAU ******')
   kdeg = 0  !SET UNUSED DEGREES TO ZERO
   eqn: do j = 1, n
      read (7, 1000) numt(j)
      write (6, 210) j, numt(j)
210   format(/, '  NUMT(', i2, ') =', i5)
      nt = numt(j)
      terms: do k = 1, nt
         vars: do l = 1, n
            read (7, 1000) kdeg(j, l, k)
            write (6, 220) j, l, k, kdeg(j, l, k)
220         format('  KDEG(', i2, ',', i2, ',', i2, ') =', i5)
         end do vars
         read (7, 2000) coef(j, k)
         write (6, 230) j, k, coef(j, k)
230      format('  COEF(', i2, ',', i2, ') =', es22.14)
      end do terms
   end do eqn
   write (6, fmt="(//)")

   iflg1 = 10*iflghm + iflgsc
   do m = 1, totdg
      iflg2(m) = -2
   end do
   call polsys1h(n, numt(1:n), coef(1:n, 1:maxt), &
                 kdeg(1:n, 1:n + 1, 1:maxt), iflg1, iflg2(1:totdg), epsbig, epssml, &
                 sspar, numrr, lambda(1:totdg), roots(1:2, 1:n + 1, 1:totdg), &
                 arclen(1:totdg), nfe(1:totdg))

   write (6, 240) iflg1
240 format('  IFLG1 =', i5,/)
   itotit = sum(nfe(1:totdg))
   do m = 1, totdg
      write (6, 260) m
260   format('  PATH NUMBER =', i5, //'  FINAL VALUES FOR PATH'/)
      write (6, 280) arclen(m)
280   format('  ARCLEN =', es22.14)
      write (6, 290) nfe(m)
290   format('  NFE =', i5)
      write (6, 300) iflg2(m)
300   format('  IFLG2 =', i3)

      ! DESIGNATE SOLUTIONS "REAL" OR "COMPLEX"
      if (any(abs(roots(2, 1:n, m)) .ge. 1.0d-4)) then
         write (6, 779, advance='NO')
779      format(' COMPLEX, ')
      else
         write (6, 780, advance='NO')
780      format(' REAL, ')
      end if

      ! DESIGNATE SOLUTION "FINITE" OR "INFINITE"
      if (sum(abs(roots(1:2, np1, m))) .lt. 1.0d-6) then
         write (6, 781)
781      format('INFINITE SOLUTION')
      else
         write (6, 782)
782      format('FINITE SOLUTION')
      end if
!
      write (6, 320) lambda(m), (j, (roots(l, j, m), l=1, 2), j=1, n)
320   format('  LAMBDA =', es22.14, /, (' X(', i2, ') =', 2es22.14))
      write (6, 330) np1, roots(1:2, np1, m)
330   format(/, ' X(', i2, ') =', 2es22.14, //)
   end do
   write (6, 400) itotit
400 format(' TOTAL NFE OVER ALL PATHS = ', i10)
   stop
1000 format(i5)
2000 format(es22.14)

end program test_p

! HOMOTOPY subroutines for the polynomial system driver POLSYS1H.
! These subroutines should be used verbatim with POLSYS1H for solving
! polynomial systems of equations. The polynomial coefficients, defined
! as input to POLSYS1H, are accessed by the routines here via the global
! arrays in HOMPACK_GLOBAL.
!
! ###################################################################
! ONLY THE SUBROUTINES RHO AND RHOJAC ARE USED BY THE POLYNOMIAL
! SYSTEM DRIVER POLSYS1H.  ALL THE OTHER ROUTINES HERE ARE PROVIDED
! SIMPLY AS TEMPLATES.
! ###################################################################

subroutine f(x, v)
!! Evaluate `f(x)` and return in the vector `v`.

   use hompack_kinds, only: dp
   implicit none

   real(dp), intent(in) :: x(:)
   real(dp), intent(out) :: v(:)

   v(1) = x(1)

end subroutine f

subroutine fjac(x, v, k)
!! Return in `v` the `k`-th column of the Jacobian matrix of `f(x)` evaluated at `x`.

   use hompack_kinds, only: dp
   implicit none

   real(dp), intent(in) :: x(:)
   real(dp), intent(out) :: v(:)
   integer, intent(in) :: k

   v(1) = x(1)

end subroutine fjac

subroutine rho(a, lambda, x, v)
!! Evaluate `rho(a,lambda,x)` and return in the vector `v`.

   use hompack_kinds, only: dp, zero
   use hompack_core, only: hfunp
   use hompack_global, only: ipar, par
   implicit none

   real(dp), intent(in) :: a(:), x(:)
   real(dp), intent(inout) :: lambda
   real(dp), intent(out) :: v(:)

   integer:: j, npol

   ! THE FOLLOWING CODE IS SPECIFICALLY FOR THE POLYNOMIAL SYSTEM DRIVER
   ! POLSYS1H , AND SHOULD BE USED VERBATIM WITH  POLSYS1H .  IF THE USER
   ! CALLING  FIXP??  OR   STEP??  DIRECTLY, HE MUST SUPPLY APPROPRIATE
   ! REPLACEMENT CODE HERE.

   ! FORCE PREDICTED POINT TO HAVE LAMBDA .GE. 0
   if (lambda .lt. zero) lambda = zero
   npol = ipar(1)
   call hfunp(npol, a, lambda, x)
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
   use hompack_global, only: ipar, par
   implicit none

   real(dp), intent(in) :: a(:), x(:)
   real(dp), intent(inout) :: lambda
   real(dp), intent(out) :: v(:)
   integer, intent(in) :: k

   ! THE FOLLOWING CODE IS SPECIFICALLY FOR THE POLYNOMIAL SYSTEM DRIVER
   !  POLSYS1H , AND SHOULD BE USED VERBATIM WITH  POLSYS1H .  IF THE USER
   ! CALLING  FIXP??  OR   STEP??  DIRECTLY, HE MUST SUPPLY APPROPRIATE
   ! REPLACEMENT CODE HERE.

   integer:: j, npol, n2
   npol = ipar(1)
   n2 = 2*npol
   if (k .eq. 1) then
      ! FORCE PREDICTED POINT TO HAVE  LAMBDA .GE. 0
      if (lambda .lt. zero) lambda = zero
      call hfunp(npol, a, lambda, x)
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
   ! use hompack_global, only: qrsparse, rowpos, colpos
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
   ! use hompack_global, only: qrsparse, rowpos, colpos
   implicit none

   real(dp), intent(in) :: a(:), lambda, x(:)

end subroutine rhojs
