module hompack_core

   use hompack_kinds, only: dp
   implicit none

   type root_state
      real(dp) :: a
      real(dp) :: acbs
      real(dp) :: acmb
      real(dp) :: ae
      real(dp) :: cmb
      real(dp) :: fa
      real(dp) :: fb
      real(dp) :: fc
      real(dp) :: fx
      real(dp) :: p
      real(dp) :: q
      real(dp) :: re
      real(dp) :: tol
      real(dp) :: u
      integer  :: ic
      integer  :: kount
   end type root_state

contains

   subroutine root(t, ft, b, c, relerr, abserr, iflag)
   !! This subroutine computes a root of the nonlinear equation `f(x) = 0` where `f(x)`
   !! is a continuous real function of a single real variable `x`. The method used is a
   !! combination of bisection and the secant rule.
   !!
   !! Normal input consists of a continuous function `f` and an interval `(b, c)` such
   !! that `f(b) * f(c) <= 0.0` . Each iteration finds new values of `b` and `c` such that
   !! the interval `(b, c)` is shrunk and `f(b) * f(c) <= 0.0` . The stopping criterion is
   !!
   !!      `abs(b - c) <= 2.0 * (relerr * dabs(b) + abserr)`
   !!
   !! where `relerr` = relative error and `abserr` = absolute error are input quantities.
   !! Set the flag, `iflag`, positive to initialize the computation. As `b`, `c` and
   !! `iflag` are used for both input and output, they must be variables in the calling
   !! program.
   !!
   !! If 0 is a possible root, one should not choose `abserr = 0.0` .
   !!
   !! The output value of `b` is the better approximation to a root as `b` and `c` are
   !! always redefined so that `abs(f(b)) <= abs(f(c))` .
   !!
   !! To solve the equation, `root` must evaluate `f(x)` repeatedly. This is done in the
   !! calling program. When an evaluation of `f` is needed at `t`, `root` returns with
   !! `iflag` negative. Evaluate `ft = f(t)` and call `root` again. Do not alter `iflag`.
   !! When the computation is complete, `root` returns to the calling program with `iflag`
   !! positive.
   !!
   !! This code is a modification of the code `zeroin` which is completely explained and
   !! documented in the text: "Numerical Computing: An Introduction", by L. F. Shampine and
   !! R. C. Allen.

      use hompack_kinds, only: one, zero
      implicit none

      real(dp), intent(inout) :: t
         !! Point at which the function is to be evaluated.
         !! On output with `iflag < 0`, `t` contains the next point where the caller must
         !! evaluate the function and set `ft = f(t)` before calling `root` again.
         !! When the iteration terminates successfully, `t` contains the final evaluation
         !! point associated with the computed root approximation.
      real(dp), intent(inout) :: ft
         !! Function value at `t`.
         !! When `iflag < 0`, the caller must compute `ft = f(t)` and call `root` again
         !! without modifying `iflag`.
      real(dp), intent(inout) :: b
         !! Lower endpoint of the current bracketing interval.
         !! On input, together with `c`, defines an interval containing a root, typically
         !! satisfying `f(b)*f(c) <= 0`.
         !! On output, contains the best approximation to the root. The algorithm
         !! maintains `abs(f(b)) <= abs(f(c))`.
      real(dp), intent(inout) :: c
         !! Upper endpoint of the current bracketing interval.
         !! On output, contains the second endpoint of the final bracketing interval.
      real(dp), intent(in) :: relerr
         !! Relative error tolerance.
         !! Convergence is declared when `abs(b-c) <= 2*(relerr*abs(b) + abserr)`.
      real(dp), intent(in) :: abserr
         !! Absolute error tolerance.
         !! Used together with `relerr` in the convergence test.
         !! A nonzero value is recommended when a root near zero is possible.
      integer, intent(inout) :: iflag
         !! Reverse-communication control and status flag.
         !!
         !! On input, set to a positive value to initialize the computation.
         !! During iteration, negative values indicate that the caller must evaluate the
         !! function at the point returned in `t` and then call `root` again without
         !! changing `iflag`.
         !!
         !! On successful or terminating return:
         !!
         !! * `1` : root bracketed and convergence criterion satisfied.
         !! * `2` : a point was found for which the computed function value is exactly
         !!         zero.
         !! * `3` : `abs(f(b))` increased relative to the initial values, suggesting
         !!         proximity to a pole or singularity.
         !! * `4` : no odd-order root was detected in the interval; a local minimum
         !!         may have been encountered.
         !! * `5` : maximum number of function evaluations (500) exceeded.

      real(dp):: a, acbs, acmb, ae, cmb, fa, fb, fc, fx, p, q, re, tol, u
      integer ic, kount
      save

      if (iflag >= 0) go to 100
      iflag = abs(iflag)
      if (iflag == 1) go to 200
      if (iflag == 2) go to 300
      if (iflag == 3) go to 400

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
1     if (abs(fc) >= abs(fb)) go to 2

      ! INTERCHANGE B AND C SO THAT ABS(F(B))<=ABS(F(C)).
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
      if (acmb <= tol) go to 8
      if (kount >= 500) go to 12

      ! CALCULATE NEW ITERATE EXPLICITLY AS B+P/Q
      ! WHERE WE ARRANGE P>=0.  THE IMPLICIT
      ! FORM IS USED TO PREVENT OVERFLOW.
      p = (b - a)*fb
      q = fa - fb
      if (p >= zero) go to 3
      p = -p
      q = -q

      ! UPDATE A, CHECK IF REDUCTION IN THE SIZE OF BRACKETING
      ! INTERVAL IS SATISFACTORY. IF NOT BISECT UNTIL IT IS.
3     a = b
      fa = fb
      ic = ic + 1
      if (ic < 4) go to 4
      if (8*acmb >= acbs) go to 6
      ic = 0
      acbs = acmb

      ! TEST FOR TOO SMALL A CHANGE
4     if (p > abs(q)*tol) go to 5

      ! INCREMENT BY TOLERANCE
      b = b + sign(tol, cmb)
      go to 7

      !  ROOT OUGHT TO BE BETWEEN B AND (C+B)/2
5     if (p >= cmb*q) go to 6

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
      if (fb == zero) go to 9
      kount = kount + 1
      if (sign(one, fb) .ne. sign(one, fc)) go to 1
      c = a
      fc = fa
      go to 1

      ! FINISHED. SET IFLAG.
8     if (sign(one, fb) == sign(one, fc)) go to 11
      if (abs(fb) > fx) go to 10
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
