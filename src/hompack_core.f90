module hompack_core

   use hompack_kinds, only: dp
   implicit none

   type root_state
      real(dp) :: t
      real(dp) :: ft
      real(dp) :: b
      real(dp) :: c
      real(dp) :: relerr
      real(dp) :: abserr
      integer :: iflag
   end type root_state

contains

   subroutine root(t, ft, b, c, relerr, abserr, iflag)
!
!  ROOT COMPUTES A ROOT OF THE NONLINEAR EQUATION F(X)=0
!  WHERE F(X) IS A CONTINOUS REAL FUNCTION OF A SINGLE REAL
!  VARIABLE X.  THE METHOD USED IS A COMBINATION OF BISECTION
!  AND THE SECANT RULE.
!
!  NORMAL INPUT CONSISTS OF A CONTINUOS FUNCTION F AND AN
!  INTERVAL (B,C) SUCH THAT F(B)*F(C).LE.0.0.  EACH ITERATION
!  FINDS NEW VALUES OF B AND C SUCH THAT THE INTERVAL(B,C) IS
!  SHRUNK AND F(B)*F(C).LE.0.0.  THE STOPPING CRITERION IS
!
!          DABS(B-C).LE.2.0*(RELERR*DABS(B)+ABSERR)
!
!  WHERE RELERR=RELATIVE ERROR AND ABSERR=ABSOLUTE ERROR ARE
!  INPUT QUANTITIES.  SET THE FLAG, IFLAG, POSITIVE TO INITIALIZE
!  THE COMPUTATION.  AS B,C AND IFLAG ARE USED FOR BOTH INPUT AND
!  OUTPUT, THEY MUST BE VARIABLES IN THE CALLING PROGRAM.
!
!  IF 0 IS A POSSIBLE ROOT, ONE SHOULD NOT CHOOSE ABSERR=0.0.
!
!  THE OUTPUT VALUE OF B IS THE BETTER APPROXIMATION TO A ROOT
!  AS B AND C ARE ALWAYS REDEFINED SO THAT DABS(F(B)).LE.DABS(F(C)).
!
!  TO SOLVE THE EQUATION, ROOT MUST EVALUATE F(X) REPEATEDLY. THIS
!  IS DONE IN THE CALLING PROGRAM.  WHEN AN EVALUATION OF F IS
!  NEEDED AT T, ROOT RETURNS WITH IFLAG NEGATIVE.  EVALUATE FT=F(T)
!  AND CALL ROOT AGAIN.  DO NOT ALTER IFLAG.
!
!  WHEN THE COMPUTATION IS COMPLETE, ROOT RETURNS TO THE CALLING
!  PROGRAM WITH IFLAG POSITIVE=
!
!     IFLAG=1  IF F(B)*F(C).LT.0 AND THE STOPPING CRITERION IS MET.
!
!          =2  IF A VALUE B IS FOUND SUCH THAT THE COMPUTED VALUE
!              F(B) IS EXACTLY ZERO.  THE INTERVAL (B,C) MAY NOT
!              SATISFY THE STOPPING CRITERION.
!
!          =3  IF DABS(F(B)) EXCEEDS THE INPUT VALUES DABS(F(B)),
!              DABS(F(C)).  IN THIS CASE IT IS LIKELY THAT B IS CLOSE
!              TO A POLE OF F.
!
!          =4  IF NO ODD ORDER ROOT WAS FOUND IN THE INTERVAL.  A
!              LOCAL MINIMUM MAY HAVE BEEN OBTAINED.
!
!          =5  IF TOO MANY FUNCTION EVALUATIONS WERE MADE.
!              (AS PROGRAMMED, 500 ARE ALLOWED.)
!
!  THIS CODE IS A MODIFICATION OF THE CODE ZEROIN WHICH IS COMPLETELY
!  EXPLAINED AND DOCUMENTED IN THE TEXT  NUMERICAL COMPUTING:  AN
!  INTRODUCTION,  BY L. F. SHAMPINE AND R. C. ALLEN.
!
      use hompack_kinds, only: one, zero
      implicit none
!
      real(dp):: a, abserr, acbs, acmb, ae, b, c, cmb, fa, fb, &
                 fc, ft, fx, p, q, re, relerr, t, tol, u
      integer ic, iflag, kount
      save

      if (iflag .ge. 0) go to 100
      iflag = abs(iflag)
      go to(200, 300, 400), iflag

100   u = epsilon(one)
      re = max(relerr, u)
      ae = max(abserr, zero)
      ic = 0
      acbs = abs(b - c)
      a = c
      t = a
      iflag = -1
      return
200   fa = ft
      t = b
      iflag = -2
      return
300   fb = ft
      fc = fa
      kount = 2
      fx = max(abs(fb), abs(fc))
1     if (abs(fc) .ge. abs(fb)) go to 2

      ! INTERCHANGE B AND C SO THAT ABS(F(B)).LE.ABS(F(C)).
      a = b
      fa = fb
      b = c
      fb = fc
      c = a
      fc = fa
2     cmb = (c - b)/2
      acmb = abs(cmb)
      tol = re*abs(b) + ae

      ! TEST STOPPING CRITERION AND FUNCTION COUNT
      if (acmb .le. tol) go to 8
      if (kount .ge. 500) go to 12

      ! CALCULATE NEW ITERATE EXPLICITLY AS B+P/Q
      ! WHERE WE ARRANGE P.GE.0.  THE IMPLICIT
      ! FORM IS USED TO PREVENT OVERFLOW.
      p = (b - a)*fb
      q = fa - fb
      if (p .ge. zero) go to 3
      p = -p
      q = -q

      ! UPDATE A, CHECK IF REDUCTION IN THE SIZE OF BRACKETING
      ! INTERVAL IS SATISFACTORY. IF NOT BISECT UNTIL IT IS.
3     a = b
      fa = fb
      ic = ic + 1
      if (ic .lt. 4) go to 4
      if (8*acmb .ge. acbs) go to 6
      ic = 0
      acbs = acmb

      ! TEST FOR TOO SMALL A CHANGE
4     if (p .gt. abs(q)*tol) go to 5

      ! INCREMENT BY TOLERANCE
      b = b + sign(tol, cmb)
      go to 7

      !  ROOT OUGHT TO BE BETWEEN B AND (C+B)/2
5     if (p .ge. cmb*q) go to 6

      ! USE SECANT RULE
      b = b + p/q
      go to 7

      ! USE BISECTION
6     b = (c + b)/2

      ! HAVE COMPLETED COMPUTATION FOR NEW ITERATE B
7     t = b
      iflag = -3
      return
400   fb = ft
      if (fb .eq. zero) go to 9
      kount = kount + 1
      if (sign(one, fb) .ne. sign(one, fc)) go to 1
      c = a
      fc = fa
      go to 1

      ! FINISHED. SET IFLAG.
8     if (sign(one, fb) .eq. sign(one, fc)) go to 11
      if (abs(fb) .gt. fx) go to 10
      iflag = 1
      return
9     iflag = 2
      return
10    iflag = 3
      return
11    iflag = 4
      return
12    iflag = 5
      return

   end subroutine root

end module hompack_core
