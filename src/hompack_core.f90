module hompack_core

   use hompack_kinds, only: dp
   implicit none

   type root_state
   !! State variables for [[root]].
      real(dp) :: a
         !! Previous iterate point used for secant step calculation.
      real(dp) :: acbs
         !! Width of the bracketing interval at the last forced bisection checkpoint.
      real(dp) :: acmb
         !! Absolute value of half the current interval length, `|(c - b)/2|`.
      real(dp) :: ae
         !! Internal absolute error tolerance bound, `max(abserr, 0.0)`.
      real(dp) :: cmb
         !! Half-interval vector representation, `(c - b)/2`.
      real(dp) :: fa
         !! Function value evaluated at point `a`, `f(a)`.
      real(dp) :: fb
         !! Function value evaluated at point `b`, `f(b)`.
      real(dp) :: fc
         !! Function value evaluated at point `c`, `f(c)`.
      real(dp) :: fx
         !! Maximum initial absolute function value, used to detect poles or divergence.
      real(dp) :: p
         !! Numerator of the proposed secant step update fraction `(p/q)`.
      real(dp) :: q
         !! Denominator of the proposed secant step update fraction `(p/q)`.
      real(dp) :: re
         !! Internal relative error tolerance bound, `max(relerr, epsilon)`.
      real(dp) :: tol
         !! Dynamic convergence tolerance threshold for the current iteration step.
      integer  :: ic
         !! Counter of sequential secant steps used to force bisection if convergence stalls.
      integer  :: fcount
         !! Total number of function evaluations completed so far.
   end type root_state

contains

   subroutine root(t, ft, b, c, relerr, abserr, iflag, state)
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
      type(root_state), intent(inout) :: state
         !! Internal state of the root-finding iteration. The caller should not modify it.
         !! Used internally by `root` to maintain the state of the iteration across calls.

      real(dp), parameter :: u = epsilon(one)
      integer, parameter :: max_fcount = 500

      associate (a => state%a, acbs => state%acbs, acmb => state%acmb, &
                 ae => state%ae, cmb => state%cmb, fa => state%fa, &
                 fb => state%fb, fc => state%fc, fx => state%fx, &
                 p => state%p, q => state%q, re => state%re, &
                 tol => state%tol, ic => state%ic, fcount => state%fcount)

         if (iflag >= 0) go to 100
         iflag = abs(iflag)
         if (iflag == 1) go to 200
         if (iflag == 2) go to 300
         if (iflag == 3) go to 400

100      re = max(relerr, u)
         ae = max(abserr, zero)
         ic = 0
         acbs = abs(b - c)
         a = c
         t = a
         iflag = -1
         return
200      fa = ft
         t = b
         iflag = -2
         return
300      fb = ft
         fc = fa
         fcount = 2
         fx = max(abs(fb), abs(fc))
1        if (abs(fc) >= abs(fb)) go to 2

         ! Interchange 'b' and 'c' so that 'abs(f(b))<=abs(f(c))'
         a = b
         fa = fb
         b = c
         fb = fc
         c = a
         fc = fa
2        cmb = (c - b)/2
         acmb = abs(cmb)
         tol = re*abs(b) + ae

         ! Test stopping criterion and function count
         if (acmb <= tol) go to 8
         if (fcount >= max_fcount) go to 12

         ! Calculate new iterate explicitly as 'b+p/q' where we arrange 'p>=0'.
         ! The implicit form is used to prevent overflow.
         p = (b - a)*fb
         q = fa - fb
         if (p >= zero) go to 3
         p = -p
         q = -q

         ! Update 'a', check if reduction in the size of bracketing interval is
         ! satisfactory. If not bisect until it is.
3        a = b
         fa = fb
         ic = ic + 1
         if (ic < 4) go to 4
         if (8*acmb >= acbs) go to 6
         ic = 0
         acbs = acmb

         ! Test for too small a change
4        if (p > abs(q)*tol) go to 5

         ! Increment by tolerance
         b = b + sign(tol, cmb)
         go to 7

         ! Root ought to be between 'b' and '(c+b)/2'
5        if (p >= cmb*q) go to 6

         ! Use secant rule
         b = b + p/q
         go to 7

         ! Use bisection
6        b = (c + b)/2

         ! Have completed computation for new iterate 'b'
7        t = b
         iflag = -3
         return
400      fb = ft
         if (fb == zero) go to 9
         fcount = fcount + 1
         if (sign(one, fb) .ne. sign(one, fc)) go to 1
         c = a
         fc = fa
         go to 1

         ! Finished. Set 'iflag'.
8        if (sign(one, fb) == sign(one, fc)) go to 11
         if (abs(fb) > fx) go to 10
         iflag = 1
         return
9        iflag = 2
         return
10       iflag = 3
         return
11       iflag = 4
         return
12       iflag = 5
         return

      end associate

   end subroutine root

end module hompack_core
