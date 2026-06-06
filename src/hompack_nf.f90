module hompack_nf
   use hompack_kinds, only: dp
   use iso_c_binding, only: c_ptr, c_null_ptr
   implicit none

   abstract interface
      subroutine f_t(x, v, data)
         import :: dp, c_ptr
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: v(:)
         type(c_ptr), value :: data
      end subroutine f_t

      subroutine fjac_t(x, v, k, data)
         import :: dp, c_ptr
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: v(:)
         integer, intent(in) :: k
         type(c_ptr), value :: data
      end subroutine fjac_t

      subroutine rho_t(a, lambda, x, v, data)
         import :: dp, c_ptr
         real(dp), intent(in) :: a(:)
         real(dp), intent(in) :: lambda
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: v(:)
         type(c_ptr), value :: data
      end subroutine rho_t

      subroutine rhojac_t(a, lambda, x, v, k, data)
         import :: dp, c_ptr
         real(dp), intent(in) :: a(:)
         real(dp), intent(in) :: lambda
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: v(:)
         integer, intent(in) :: k
         type(c_ptr), value :: data
      end subroutine rhojac_t
   end interface

   type hompack_callbacks
   !! Type to hold user-supplied function and Jacobian evaluation subroutines for
   !! f-suffix (dense Jacobian) methods.
      procedure(f_t), nopass, pointer :: f => null()
      procedure(fjac_t), nopass, pointer :: fjac => null()
      procedure(rho_t), nopass, pointer :: rho => null()
      procedure(rhojac_t), nopass, pointer :: rhojac => null()
      type(c_ptr) :: data = c_null_ptr
   end type hompack_callbacks

contains

   impure subroutine fixpnf( &
      callbacks, &
      n, y, iflag, arcre, arcae, ansre, ansae, trace, a, sspar, nfe, arclen, poly_switch)
   !! This subroutine finds a fixed point or zero of the N-dimensional vector function
   !! \( F(x) \), or tracks a zero curve of a general homotopy map \( \rho(a,\lambda,x) \).
   !!
   !! For the fixed-point problem, \( F(x) \) is assumed to be a \(C^2\) map of some ball
   !! into itself. The equation \( x = F(x) \) is solved by following the zero curve of
   !! the homotopy map
   !!
   !! $$ \rho_a(\lambda, x) = \lambda (x - F(x)) + (1 - \lambda) (x - a) $$
   !!
   !! starting from \( \lambda=0, x=a \). The curve is parameterized by arc length \(s\),
   !! and is followed by solving the ordinary differential equation \( d \rho/ d s = 0 \)
   !! for \( y(s) = (\lambda(s), x(s)) \) using a Hermite cubic predictor and a corrector
   !! which returns to the zero curve along the flow normal to the Davidenko flow (which
   !! consists of the integral curves of \( d \rho/ d s \)).
   !!
   !! For the zero-finding problem, \( F(x) \) is assumed to be a \(C^2\) map such that
   !! for some \(r > 0\),  \(x F(x) \ge 0\) whenever \( ||x|| = r \). The equation
   !! \( F(x) = 0 \) is solved by following the zero curve of the homotopy map
   !!
   !! $$ \rho_a(\lambda, x) = \lambda F(x) + (1 - \lambda) (x - a) $$
   !!
   !! emanating from \( \lambda = 0 \), \(x = a \). Parameter \(a\)  must be an interior
   !! point of the above mentioned balls.
   !!
   !! For the curve tracking problem, \( \rho(a,\lambda,x) \) is assumed to be a \(C^2\)
   !! map from \( E^m \times [0,1) \times E^n \) into \( E^n \), which for almost all
   !! parameter vectors \( a \) in some nonempty open subset of \( E^m \) satisfies
   !!
   !!  $$\operatorname{rank} \left[ \frac{\partial \rho(a,\lambda,x)}{\partial \lambda}, \frac{\partial \rho(a,\lambda,x)}{\partial x} \right] = N$$
   !!
   !! for all points \( (\lambda,x) \) such that \( \rho(a,\lambda,x)=0 \). It is further
   !! assumed that
   !!
   !! $$\operatorname{rank} \left[ \frac{\partial \rho(a,0,x_0)}{\partial x} \right] = N $$
   !!
   !! With \(a\) fixed, the zero curve of \( \rho(a,\lambda,x) \) emanating from
   !! \( \lambda=0 \), \( x=x_0 \) is tracked until \( \lambda=1 \) by solving the
   !! ordinary differential equation \( d\rho(a,\lambda(s),x(s))/ds = 0 \) for
   !! \( y(s) = (\lambda(s), x(s)) \), where \(s\) is the arc length along the zero curve.
   !! Also the homotopy map \( \rho(a,\lambda,x) \) is assumed to be constructed such that
   !! \( d\lambda(0) / ds > 0 \).
   !!
   !! For the fixed point and zero finding problems, the user must supply a subroutine
   !! `f(x,v)` which evaluates \(F(x)\) at \(x\) and returns the vector \(F(x)\) in `v`,
   !! and a subroutine `fjac(x,v,k)` which returns in `v` the `k`th column of the Jacobian
   !! matrix of \(F(x)\) evaluated at \(x\).
   !!
   !! For the curve tracking problem, the user must supply a subroutine `rho(a,lambda,x,v)`
   !! which evaluates the homotopy map \(\rho\) at \( (a,\lambda,x) \) and returns the
   !! corresponding vector in `v`, and a subroutine `rhojac(a,lambda,x,v,k)` which
   !! returns in `v` the `k`th column of the \( N \times (N+1) \) Jacobian matrix
   !! \( [\partial \rho/\partial \lambda, \partial \rho/\partial x] \) evaluated
   !! at \( (a,\lambda,x) \).

      ! ON INPUT:
      !
      ! N  is the dimension of X, F(X), and RHO(A,LAMBDA,X).
      !
      ! Y(:)  is an array of length  N + 1.  (Y(2),...,Y(N+1)) = A  is the
      !    starting point for the zero curve for the fixed point and
      !    zero finding problems.  (Y(2),...,Y(N+1)) = X0  for the curve
      !    tracking problem.
      !
      ! IFLAG  can be -2, -1, 0, 2, or 3.  IFLAG  should be 0 on the
      !    first call to  FIXPNF  for the problem  X=F(X), -1 for the
      !    problem  F(X)=0, and -2 for the problem  RHO(A,LAMBDA,X)=0.
      !    In certain situations  IFLAG  is set to 2 or 3 by  FIXPNF,
      !    and  FIXPNF  can be called again without changing  IFLAG.
      !
      ! ARCRE , ARCAE  are the relative and absolute errors, respectively,
      !    allowed the normal flow iteration along the zero curve.  If
      !    ARC?E .LE. 0.0  on input it is reset to  .5*SQRT(ANS?E) .
      !    Normally  ARC?E should be considerably larger than  ANS?E .
      !
      ! ANSRE , ANSAE  are the relative and absolute error values used for
      !    the answer at LAMBDA = 1.  The accepted answer  Y = (LAMBDA, X)
      !    satisfies
      !
      !       |Y(1) - 1|  .LE.  ANSRE + ANSAE           .AND.
      !
      !       ||Z||  .LE.  ANSRE*||X|| + ANSAE          where
      !
      !    (.,Z) is the Newton step to Y.
      !
      ! TRACE  is an integer specifying the logical I/O unit for
      !    intermediate output.  If  TRACE .GT. 0  the points computed on
      !    the zero curve are written to I/O unit  TRACE .
      !
      ! A(:)  contains the parameter vector  A .  For the fixed point
      !    and zero finding problems, A  need not be initialized by the
      !    user, and is assumed to have length  N.  For the curve
      !    tracking problem, A  must be initialized by the user.
      !
      ! SSPAR(1:8) = (LIDEAL, RIDEAL, DIDEAL, HMIN, HMAX, BMIN, BMAX, P)  is
      !    a vector of parameters used for the optimal step size estimation.
      !    If  SSPAR(J) .LE. 0.0  on input, it is reset to a default value
      !    by  FIXPNF .  Otherwise the input value of  SSPAR(J)  is used.
      !    See the comments below and in  STEPNF  for more information about
      !    these constants.
      !
      ! POLY_SWITCH  is an optional logical variable used only by the driver
      !    POLSYS1H  for polynomial systems.
      !
      !
      ! ON OUTPUT:
      !
      ! N , TRACE , A  are unchanged.
      !
      ! Y(1) = LAMBDA, (Y(2),...,Y(N+1)) = X, and Y is an approximate
      !    zero of the homotopy map.  Normally LAMBDA = 1 and X is a
      !    fixed point(zero) of F(X).  In abnormal situations LAMBDA
      !    may only be near 1 and X is near a fixed point(zero).
      !
      ! IFLAG =
      !  -2   causes  FIXPNF  to initialize everything for the problem
      !       RHO(A,LAMBDA,X) = 0 (use on first call).
      !
      !  -1   causes  FIXPNF  to initialize everything for the problem
      !       F(X) = 0 (use on first call).
      !
      !   0   causes  FIXPNF  to initialize everything for the problem
      !       X = F(X) (use on first call).
      !
      !   1   Normal return.
      !
      !   2   Specified error tolerance cannot be met.  Some or all of
      !       ARCRE , ARCAE , ANSRE , ANSAE  have been increased to
      !       suitable values.  To continue, just call  FIXPNF  again
      !       without changing any parameters.
      !
      !   3   STEPNF  has been called 1000 times.  To continue, call
      !       FIXPNF  again without changing any parameters.
      !
      !   4   Jacobian matrix does not have full rank.  The algorithm
      !       has failed (the zero curve of the homotopy map cannot be
      !       followed any further).
      !
      !   5   The tracking algorithm has lost the zero curve of the
      !       homotopy map and is not making progress.  The error tolerances
      !       ARC?E  and  ANS?E  were too lenient.  The problem should be
      !       restarted by calling  FIXPNF  with smaller error tolerances
      !       and  IFLAG = 0 (-1, -2).
      !
      !   6   The normal flow Newton iteration in  STEPNF  or  ROOTNF
      !       failed to converge.  The error tolerances  ANS?E  may be too
      !       stringent.
      !
      !   7   Illegal input parameters, a fatal error.
      !
      !   8   Memory allocation error, fatal.
      !
      ! ARCRE , ARCAE , ANSRE , ANSAE  are unchanged after a normal return
      !    (IFLAG = 1).  They are increased to appropriate values on the
      !    return  IFLAG = 2 .
      !
      ! NFE  is the number of function evaluations (= number of
      !    Jacobian matrix evaluations).
      !
      ! ARCLEN  is the length of the path followed.
      !
      ! Allocatable and automatic work arrays:
      !
      ! YP(1:N+1)  is a work array containing the tangent vector to
      !    the zero curve at the current point  Y .
      !
      ! YOLD(1:N+1)  is a work array containing the previous point found
      !    on the zero curve.
      !
      ! YPOLD(1:N+1)  is a work array containing the tangent vector to
      !    the zero curve at  YOLD .
      !
      ! QR(1:N,1:N+2), ALPHA(1:3*N+3), TZ(1:N+1), PIVOT(1:N+1), W(1:N+1),
      !    WP(1:N+1), Z0(1:N+1), Z1(1:N+1)  are all work arrays used by
      !    STEPNF  to calculate the tangent vectors and Newton steps.

      use hompack_kinds, only: zero, one
      use blas_interfaces, only: dnrm2
      implicit none

      type(hompack_callbacks) :: callbacks
         !! User-supplied function and Jacobian evaluation subroutines.
      integer, intent(in) :: n
         !! Dimension of `x`, `f(x)`, and `rho(a,lambda,x)`.
      real(dp), intent(inout) :: y(:)
         !! Homotopy solution vector. `Shape: (n+1)`.
         !! On input, `y(2:n+1)` contains the starting point:
         !! for fixed-point and zero-finding problems this is the initial point `a`,
         !! and for curve tracking it is the initial solution `x0`.
         !! On output, `y(1)=lambda` and `y(2:n+1)=x`, an approximate zero of the
         !! homotopy map. Normally `lambda=1` and `x` is a fixed point or zero of
         !! the target function.
      integer, intent(inout) :: iflag
         !! Problem type and status flag.
         !!
         !! On input:
         !! * `0`  : solve `x = f(x)` (first call).
         !! * `-1` : solve `f(x) = 0` (first call).
         !! * `-2` : track a zero curve of `rho(a,lambda,x) = 0` (first call).
         !!
         !! On output:
         !! * `1` : normal return.
         !! * `2` : requested tolerances cannot be achieved; tolerances have been
         !!         increased and the routine may be called again.
         !! * `3` : iteration limit reached; call again to continue.
         !! * `4` : Jacobian matrix lost full rank.
         !! * `5` : tracking algorithm lost the zero curve.
         !! * `6` : Newton iteration failed to converge.
         !! * `7` : illegal input parameters.
         !! * `8` : memory allocation failure.
      real(dp), intent(inout) :: arcre
         !! Relative error tolerance for the normal-flow iteration used while
         !! tracking the zero curve. If nonpositive on input, a default value is
         !! chosen. May be increased when `iflag=2`.
      real(dp), intent(inout) :: arcae
         !! Absolute error tolerance for the normal-flow iteration used while
         !! tracking the zero curve. If nonpositive on input, a default value is
         !! chosen. May be increased when `iflag=2`.
      real(dp), intent(inout) :: ansre
         !! Relative error tolerance required of the final solution at `lambda=1`.
         !! May be increased when `iflag=2`.
      real(dp), intent(inout) :: ansae
         !! Absolute error tolerance required of the final solution at `lambda = 1`.
         !! May be increased when `iflag=2`.
      integer, intent(in) :: trace
         !! Logical I/O unit for intermediate output. If `trace > 0`, points
         !! computed along the zero curve are written to this unit.
      real(dp), intent(inout) :: a(:)
         !! Parameter vector `a`.
         !! For fixed-point and zero-finding problems, the array is assumed to have
         !! length `n` and need not be initialized by the user.
         !! For curve-tracking problems, it must be initialized on input.
         !! Unchanged on output.
      real(dp), intent(inout) :: sspar(8)
         !! Step-size control parameters:
         !! `(LIDEAL, RIDEAL, DIDEAL, HMIN, HMAX, BMIN, BMAX, P)`.
         !! Parameters used by the optimal step-size estimation algorithm.
         !! Elements that are nonpositive on input are replaced by default values.
      integer, intent(out) :: nfe
         !! Number of homotopy/Jacobian evaluations performed.
         !! This counter is incremented once for each call to `tangnf`, where the
         !! homotopy map and its Jacobian are assembled and used to compute a tangent
         !! vector and Newton correction.
      real(dp), intent(out) :: arclen
         !! Total arc length of the solution path followed by the algorithm.
      logical, intent(in), optional :: poly_switch
         !! Optional flag used only by the polynomial-system driver `POLSYS1H`.

      real(dp), save :: abserr, curtol, h, hold, relerr, s
      integer, save :: iflagc, iter, jw, limit, nc, nfec, np1
      logical, save :: crash, polsys, start
      real(dp), dimension(:), allocatable, save :: yold, yp, ypold

      real(dp) :: alpha(3*n + 3), qr(n, n + 2), tz(n + 1), &
                  w(n + 1), wp(n + 1), z0(n + 1), z1(n + 1)
      integer :: pivot(n + 1)

      ! Upper bound on the number of steps
      integer, parameter :: limitd = 1000

      ! Switch from the tolerance arc?e to the (finer) tolerance ans?e if the curvature
      ! of any component of y exceeds cursw
      real(dp), parameter :: cursw = 10.0_dp

      ! Check callbacks are present
      if (iflag == 0 .or. iflag == -1) then
         if (.not. associated(callbacks%f) .or. .not. associated(callbacks%fjac)) then
            iflag = 7
            return
         end if
      else if (iflag == -2) then
         if (.not. associated(callbacks%rho) .or. .not. associated(callbacks%rhojac)) then
            iflag = 7
            return
         end if
      end if

      ! Test logical switch to reflect intended usage of 'fixpnf'
      if (present(poly_switch)) then
         polsys = poly_switch
      else
         polsys = .false.
      end if

      if (n .le. 0 .or. ansre .le. zero .or. ansae .lt. zero &
          .or. (n + 1) .ne. size(y) .or. &
          ((iflag .eq. -1 .or. iflag .eq. 0) .and. n .ne. size(a))) &
         iflag = 7
      if (iflag .ge. -2 .and. iflag .le. 0) go to 20
      if (iflag .eq. 2) go to 120
      if (iflag .eq. 3) go to 90

      ! Only valid input for 'iflag' is -2, -1, 0, 2, 3.
      iflag = 7
      return

      ! Initialization block
20    arclen = zero
      if (arcre .le. zero) arcre = sqrt(ansre)/2
      if (arcae .le. zero) arcae = sqrt(ansae)/2

      nc = n
      nfec = 0
      iflagc = iflag
      np1 = n + 1

      call cleanup
      allocate (yp(np1), yold(np1), ypold(np1), stat=jw)

      if (jw /= 0) then
         iflag = 8
         return
      end if

      ! Set initial conditions for first call to 'stepnf'
      start = .true.
      crash = .false.
      hold = one
      h = 0.1_dp
      s = zero
      ypold(1) = one
      yp(1) = one
      y(1) = zero
      ypold(2:np1) = zero
      yp(2:np1) = zero

      ! Set optimal step size estimation parameters
      ! Let Z[K] denote the Newton iterates along the flow normal to the Davidenko flow
      ! and Y their limit
      ! Ideal contraction factor: ||Z[2] - Z[1]|| / ||Z[1] - Z[0]||
      if (sspar(1) .le. zero) sspar(1) = 0.5_dp
      ! Ideal residual factor:  ||RHO(A, Z[1])|| / ||RHO(A, Z[0])||
      if (sspar(2) .le. zero) sspar(2) = 0.01_dp
      ! Ideal distance factor:  ||Z[1] - Y|| / ||Z[0] - Y||
      if (sspar(3) .le. zero) sspar(3) = 0.5_dp
      ! Minimum step size HMIN
      if (sspar(4) .le. zero) sspar(4) = (sqrt(real(n + 1, dp)) + 4.0_dp)*epsilon(one)
      ! Maximum step size HMAX
      if (sspar(5) .le. zero) sspar(5) = one
      ! Minimum step size reduction factor BMIN
      if (sspar(6) .le. zero) sspar(6) = 0.1_dp
      ! Maximum step size expansion factor BMAX
      if (sspar(7) .le. zero) sspar(7) = 3.0_dp
      ! Assumed operating order P
      if (sspar(8) .le. zero) sspar(8) = 2.0_dp

      ! Load 'A' for the fixed point and zero finding problems
      if (iflagc .ge. -1) then
         a = y(2:np1)
      end if

90    limit = limitd

      ! Main loop
120   main_loop: do iter = 1, limit

         if (y(1) .lt. zero) then
            arclen = s
            iflag = 5
            call cleanup
            return
         end if

         ! Set different error tolerance if the trajectory Y(S) has any high
         ! curvature components
         curtol = cursw*hold
         relerr = arcre
         abserr = arcae
         if (any(abs(yp - ypold) .gt. curtol)) then
            relerr = ansre
            abserr = ansae
         end if

         ! Take a step along the curve
         call stepnf(callbacks, nc, nfec, iflagc, start, crash, hold, h, relerr, abserr, &
                     s, y, yp, yold, ypold, a, qr, alpha, tz, pivot, w, wp, z0, z1, sspar)

         ! Print latest point on curve if requested
         if (trace .gt. 0) then
            write (trace, 217) iter, nfec, s, y(1), (y(jw), jw=2, np1)
217         format(/' STEP', i5, 3x, 'NFE =', i5, 3x, 'ARC LENGTH =', f9.4, 3x, &
                    'LAMBDA =', f7.4, 5x, 'X VECTOR:'/(1x, 6es12.4))
         end if
         nfe = nfec

         ! Check if the step was successful
         if (iflagc .gt. 0) then
            arclen = s
            iflag = iflagc
            call cleanup
            return
         end if

         if (crash) then
            ! Return code for error tolerance too small
            iflag = 2
            ! Change error tolerances
            if (arcre .lt. relerr) arcre = relerr
            if (ansre .lt. relerr) ansre = relerr
            if (arcae .lt. abserr) arcae = abserr
            if (ansae .lt. abserr) ansae = abserr
            ! Change limit on number of iterations
            limit = limit - iter
            return
         end if

         ! Use Hermite cubic interpolation and Newton iteration to get the
         ! answer at lambda = 1.0
         if (y(1) .ge. one) then

            ! Save 'yold' for arc length calculation later
            z0 = yold
            call rootnf(callbacks, nc, nfec, iflagc, ansre, ansae, y, yp, yold, ypold, &
                        a, qr, alpha, tz, pivot, w, wp)
            nfe = nfec
            iflag = 1

            ! Set error flag if 'rootnf' could not get the point on the zero
            ! curve at lambda = 1.0
            if (iflagc .gt. 0) iflag = iflagc

            ! Calculate final arc length
            w = y - z0
            arclen = s - hold + dnrm2(np1, w, 1)
            call cleanup
            return

         end if

         ! For polynomial systems and the POLSYS1H homotopy map, D LAMBDA/DS >= 0
         ! necessarily; this condition is forced here if the 'poly_switch' variable is
         ! present
         if (polsys) then
            if (yp(1) .lt. zero) then
               ! Reverse tangent direction so D LAMBDA/DS = YP(1) > 0
               yp = -yp
               ypold = yp
               ! Force 'stepnf' to use the linear predictor for the next step only
               start = .true.
            end if
         end if

      end do main_loop

      ! Lambda has not reached 1 in 'limitd' steps
      iflag = 3
      arclen = s

   contains

      subroutine cleanup
         if (allocated(yp)) deallocate (yp)
         if (allocated(yold)) deallocate (yold)
         if (allocated(ypold)) deallocate (ypold)
      end subroutine cleanup

   end subroutine fixpnf

   subroutine rootnf( &
      callbacks, &
      n, nfe, iflag, relerr, abserr, y, yp, yold, &
      ypold, a, qr, alpha, tz, pivot, w, wp)
   !! This subroutine finds the point `ybar = (1, xbar)` on the zero curve of the homotopy
   !! map. It starts with two points `yold = (lambdaold, xold)` and `y = (lambda, x)` such
   !! that `lambdaold < 1 <= lambda` , and alternates between secant estimates of `ybar`
   !! and Newton iteration until convergence.

      use hompack_kinds, only: zero, one
      use hompack_core, only: root
      use blas_interfaces, only: dnrm2
      implicit none

      type(hompack_callbacks) :: callbacks
         !! User-supplied function and Jacobian evaluation subroutines.
      integer, intent(in) :: n
         !! Problem dimension.
      integer, intent(inout) :: nfe
         !! Number of homotopy/Jacobian evaluations performed.
         !! On input, contains the current count.
         !! Updated on output to include evaluations performed by `rootnf`.
      integer, intent(inout) :: iflag
         !! Problem type and return status flag.
         !!
         !! On input:
         !! * `0`  : fixed-point problem, `F(x) = x`.
         !! * `-1` : zero-finding problem, `F(x) = 0`.
         !! * `-2` : general homotopy curve-tracking problem.
         !!
         !! On output:
         !! * unchanged (`0`, `-1`, or `-2`) on normal return.
         !! * `4` : Jacobian matrix lost full rank (`rank < n`); iteration not completed.
         !! * `6` : iteration failed to converge.
      real(dp), intent(in) :: relerr
         !! Relative convergence tolerance.
         !! Iteration is considered converged when `|y(1)-1| <= relerr + abserr` and the
         !! Newton correction satisfies `||z|| <= relerr*||x|| + abserr`.
      real(dp), intent(in) :: abserr
         !! Absolute convergence tolerance.
         !! Used together with `relerr` in the convergence criteria.
      real(dp), intent(inout) :: y(:)
         !! Current point on the zero curve. `Shape: (n+1)`.
         !! Contains `(lambda, x)` on input.
         !! On successful return, contains the point on the zero curve of the homotopy map
         !! at `lambda = 1`.
      real(dp), intent(inout) :: yp(:)
         !! Unit tangent vector to the zero curve at `y`. `Shape: (n+1)`.
      real(dp), intent(inout) :: yold(:)
         !! Previous point on the zero curve distinct from `y`. `Shape: (n+1)`.
         !! If convergence fails (`iflag=6`), `y` and `yold` contain the last
         !! two points computed on the zero curve.
      real(dp), intent(inout) :: ypold(:)
         !! Unit tangent vector to the zero curve at `yold`. `Shape: (n+1)`.
      real(dp), intent(in) :: a(:)
         !! Parameter vector used in the homotopy map.
      real(dp), intent(inout) :: qr(:, :)
         !! Workspace for QR factorizations used in Newton-step computations.
         !! `Shape: (n, n+2)`.
      real(dp), intent(inout) :: alpha(:)
         !! Workspace used during interpolation and Newton iteration.
         !! `Shape: (3*n+3)`.
      real(dp), intent(inout) :: tz(:)
         !! Workspace array used in QR-factorization and Newton-step calculations.
         !! `Shape: (n+1)`.
      integer, intent(inout) :: pivot(:)
         !! Pivot indices used by the QR factorization.
         !! `Shape: (n+1)`.
      real(dp), intent(inout) :: w(:)
         !! Workspace array used for interpolation and Newton-step calculations.
         !! `Shape: (n+1)`.
      real(dp), intent(inout) :: wp(:)
         !! Workspace array used for interpolation and Newton-step calculations.
         !! `Shape: (n+1)`.

      real(dp) :: dd001, dd0011, dd01, dd011, dels, f0, f1, fp0, fp1, qofs, qsout, aerr, &
                  rerr, s, sa, sb, sout, u
      integer :: judy, jw, lcode, limit, np1
      logical :: bracket

      ! Definition of hermite cubic interpolant via divided differences
      ! Note: these weird things are fortran statement functions
      dd01(f0, f1, dels) = &
         (f1 - f0)/dels
      dd001(f0, fp0, f1, dels) = &
         (dd01(f0, f1, dels) - fp0)/dels
      dd011(f0, f1, fp1, dels) = &
         (fp1 - dd01(f0, f1, dels))/dels
      dd0011(f0, fp0, f1, fp1, dels) = &
         (dd011(f0, f1, fp1, dels) - dd001(f0, fp0, f1, dels))/dels
      qofs(f0, fp0, f1, fp1, dels, s) = &
         ((dd0011(f0, fp0, f1, fp1, dels)*(s - dels) + &
           dd001(f0, fp0, f1, dels))*s + fp0)*s + f0

      u = epsilon(one)
      rerr = max(relerr, u)
      aerr = max(abserr, zero)
      np1 = n + 1

      ! The limit on the number of iterations allowed may be changed by changing the
      ! following parameter statement
      limit = 2*(int(abs(log10(aerr + rerr))) + 1)

      tz = y - yold
      dels = dnrm2(np1, tz, 1)

      ! Using two points and tangents on the homotopy zero curve, construct the Hermite
      ! cubic interpolant q(s). Then use 'root' to find the 's' corresponding to
      ! 'lambda=1'. The two points on the zero curve are always chosen to bracket
      ! 'lambda=1', with the bracketing interval always being [0, dels].
      sa = zero
      sb = dels
      lcode = 1
130   call root(sout, qsout, sa, sb, rerr, aerr, lcode)
      if (lcode .gt. 0) go to 140
      qsout = qofs(yold(1), ypold(1), y(1), yp(1), dels, sout) - one
      go to 130

      ! If lambda=1 were bracketed, root cannot fail
140   if (lcode .gt. 2) then
         iflag = 6
         return
      end if

      ! Calculate q(sa) as the initial point for a Newton iteration
      do jw = 1, np1
         w(jw) = qofs(yold(jw), ypold(jw), y(jw), yp(jw), dels, sa)
      end do

      ! Tangent information 'yp' is no longer needed. Hereafter, 'yp' represents the most
      ! recent point which is on the opposite side of the hyperplane 'lambda=1' from 'y'

      ! Prepare for main loop
      yp = yold

      ! Initialize bracket to indicate that the points 'y' and 'yold' bracket 'lambda=1',
      ! thus 'yold = yp'
      bracket = .true.

      ! Main loop
      do judy = 1, limit

         ! Calculate Newton step at current estimate 'w'
         call tangnf(callbacks, sa, w, wp, ypold, a, qr, alpha, tz, pivot, nfe, n, iflag)
         if (iflag .gt. 0) return

         ! Next point = current point + Newton step
         w = w + tz

         ! Check for convergence
         if ((abs(w(1) - one) .le. rerr + aerr) .and. &
             (dnrm2(np1, tz, 1) .le. rerr*dnrm2(n, w(2:np1), 1) + aerr)) then
            y = w
            return
         end if

         ! Prepare for next iteration
         if (abs(w(1) - one) .le. rerr + aerr) then
            ypold = wp
            cycle
         end if

         ! Update 'y' and 'yold'
         yold = y
         y = w

         ! Update 'yp' such that 'yp' is the most recent point opposite of 'lambda=1'
         ! from 'y'. Set bracket=.true. iff 'y' and 'yold' bracket 'lambda=1' so that
         ! yp = yold .
         if ((y(1) - one)*(yold(1) - one) .gt. 0) then
            bracket = .false.
         else
            bracket = .true.
            yp = yold
         end if

         ! Compute dels=||y-yp||
         tz = y - yp
         dels = dnrm2(np1, tz, 1)

         ! Compute tz for the linear predictor w = y + tz, where tz = sa*(yold-y).
         sa = (one - y(1))/(yold(1) - y(1))
         tz = sa*(yold - y)

         ! To insure stability, the linear prediction must be no farther from y than yp is.
         ! This is guaranteed if bracket=true. If linear prediction is too far away, use
         ! bracketing points to compute linear prediction.
         if (.not. bracket) then
            if (dnrm2(np1, tz, 1) .gt. dels) then
               ! Compute tz = sa*(yp-y)
               sa = (one - y(1))/(yp(1) - y(1))
               tz = sa*(yp - y)
            end if
         end if

         ! Compute estimate  w = y + tz  and save old tangent vector.
         w = w + tz
         ypold = wp

      end do

      ! The alternating secant estimation and Newton iteration has not converged in limit
      ! steps. Error return.
      iflag = 6

   end subroutine rootnf

   subroutine stepnf( &
      callbacks, &
      n, nfe, iflag, start, crash, hold, h, relerr, &
      abserr, s, y, yp, yold, ypold, a, qr, alpha, tz, pivot, w, wp, &
      z0, z1, sspar)
   !! This subroutine takes one step along the zero curve of the homotopy map using a
   !! predictor-corrector algorithm. The predictor uses a Hermite cubic interpolant, and
   !! the corrector returns to the zero curve along the flow normal to the Davidenko flow.
   !! `stepnf` also estimates a step size `h` for the next step along the zero curve.
   !! Normally, `stepnf` is used indirectly through `fixpnf`, and should be called
   !! directly only if it is necessary to modify the stepping algorithm's parameters.

      use hompack_kinds, only: one, zero
      use blas_interfaces, only: dnrm2
      implicit none

      type(hompack_callbacks), intent(in) :: callbacks
         !! User-supplied function and Jacobian evaluation subroutines.
      integer, intent(in) :: n
         !! Problem dimension.
      integer, intent(inout) :: nfe
         !! Number of homotopy/Jacobian evaluations performed.
         !! On input, contains the current count.
         !! Updated on output to include evaluations performed by `stepnf`.
      integer, intent(inout) :: iflag
         !! Problem type and return status flag.
         !!
         !! On input:
         !! * `0`  : fixed-point problem, `F(x) = x`.
         !! * `-1` : zero-finding problem, `F(x) = 0`.
         !! * `-2` : general homotopy curve-tracking problem.
         !!
         !! On output:
         !! * unchanged (`0`, `-1`, or `-2`) on normal return.
         !! * `4` : Jacobian matrix lost full rank (`rank < n`); iteration not completed.
         !! * `6` : iteration failed to converge.
      logical, intent(inout) :: start
         !! Indicates whether this is the first call to `stepnf`.
         !! Set to `.true.` on the initial call and changed to `.false.` after successful
         !! initialization.
      logical, intent(out) :: crash
         !! Status flag indicating that continuation cannot proceed with the
         !! current parameters.
         !!
         !! * `.false.` on a normal return.
         !! * `.true.` if the step size `h` is too small and has been increased to an
         !!    acceptable value.
         !! * `.true.` if `relerr` and/or `abserr` are too stringent and have been
         !!    relaxed to acceptable values.
      real(dp), intent(inout) :: hold
         !! Distance between the two most recent points on the zero curve, `||y - yold||`.
         !! Should not be modified by the user.
      real(dp), intent(inout) :: h
         !! Upper bound on the step length to be attempted. Must be initialized to a
         !! positive value on the first call. Thereafter, `stepnf` computes an estimated
         !! optimal value for the next continuation step.
      real(dp), intent(inout) :: relerr
         !! Relative convergence tolerance for the Newton corrector.
         !! On normal return the value is unchanged.
         !! If `crash=.true.`, the value may be increased to a more appropriate level.
      real(dp), intent(inout) :: abserr
         !! Absolute convergence tolerance for the Newton corrector.
         !! On normal return the value is unchanged.
         !! If `crash=.true.`, the value may be increased to a more appropriate level.
      real(dp), intent(inout) :: s
         !! Approximate arc length along the homotopy zero curve.
         !! On output, corresponds to the latest point returned in `y`.
      real(dp), intent(inout) :: y(:)
         !! Current point on the zero curve. `Shape: (n+1)`.
         !! Contains `(lambda, x)` on input.
         !! On output, updated to the latest point found by the continuation algorithm.
      real(dp), intent(inout) :: yp(:)
         !! Unit tangent vector to the zero curve at `y`. `Shape: (n+1)`.
      real(dp), intent(inout) :: yold(:)
         !! Previous point on the zero curve preceding `y`. `Shape: (n+1)`.
      real(dp), intent(inout) :: ypold(:)
         !! Unit tangent vector to the zero curve at `yold`. `Shape: (n+1)`.
      real(dp), intent(in) :: a(:)
         !! Parameter vector used in the homotopy map.
      real(dp), intent(inout) :: qr(:, :)
         !! Workspace for QR factorizations used in Newton-step calculations.
         !! `Shape: (n, n+2)`.
      real(dp), intent(inout) :: alpha(:)
         !! Workspace array used during interpolation and Newton-step calculations.
         !! `Shape: (3*n+3)`.
      real(dp), intent(inout) :: tz(:)
         !! Workspace array used in QR-factorization and Newton-step computations.
         !! `Shape: (n+1)`.
      integer, intent(inout) :: pivot(:)
         !! Pivot indices used by the QR factorization. `Shape: (n+1)`.
      real(dp), intent(inout) :: w(:)
         !! Workspace array used during interpolation and Newton-step calculations.
         !! `Shape: (n+1)`.
      real(dp), intent(inout) :: wp(:)
         !! Workspace array used during interpolation and Newton-step calculations.
         !! `Shape: (n+1)`.
      real(dp), intent(inout) :: z0(:)
         !! Workspace array used for estimating the optimal next step size. `Shape: (n+1)`.
      real(dp), intent(inout) :: z1(:)
         !! Workspace array used for estimating the optimal next step size. `Shape: (n+1)`.
      real(dp), intent(in) :: sspar(8)
         !! Step-size estimation parameters:
         !! `(lideal, rideal, dideal, hmin, hmax, bmin, bmax, p)`.
         !! Controls the adaptive continuation step-size strategy.

      real(dp) :: dcalc, dd001, dd0011, dd01, dd011, dels, f0, f1, fouru, fp0, fp1, &
                  hfail, ht, lcalc, qofs, rcalc, rholen, temp, twou
      integer :: itnum, j, judy, np1
      logical :: fail

      ! The limit on the number of Newton iterations allowed before reducing the step
      ! size 'h' may be changed by changing the following parameter
      integer, parameter:: litfh = 4

      ! Definition of Hermite cubic interpolant via divided differences
      dd01(f0, f1, dels) = &
         (f1 - f0)/dels
      dd001(f0, fp0, f1, dels) = &
         (dd01(f0, f1, dels) - fp0)/dels
      dd011(f0, f1, fp1, dels) = &
         (fp1 - dd01(f0, f1, dels))/dels
      dd0011(f0, fp0, f1, fp1, dels) = &
         (dd011(f0, f1, fp1, dels) - dd001(f0, fp0, f1, dels))/dels
      qofs(f0, fp0, f1, fp1, dels, s) = &
         ((dd0011(f0, fp0, f1, fp1, dels)*(s - dels) + &
           dd001(f0, fp0, f1, dels))*s + fp0)*s + f0

      twou = 2*epsilon(one)
      fouru = 2*twou
      np1 = n + 1
      crash = .true.

      ! The arclength 's' must be nonnegative
      if (s .lt. zero) return

      ! If step size is too small, determine an acceptable one
      if (h .lt. fouru*(one + s)) then
         h = fouru*(one + s)
         return
      end if

      ! If error tolerances are too small, increase them to acceptable values
      temp = dnrm2(np1, y, 1) + one
      if (0.5_dp*(relerr*temp + abserr) .lt. twou*temp) then
         if (relerr .ne. zero) then
            relerr = fouru*(one + fouru)
            abserr = max(abserr, zero)
         else
            abserr = fouru*temp
         end if
         return
      end if

      ! STARTUP SECTION (FIRST STEP ALONG ZERO CURVE)
      crash = .false.
      startup: if (start) then
         fail = .false.
         start = .false.

         ! Determine suitable initial step size
         h = min(h, 0.1_dp, sqrt(sqrt(relerr*temp + abserr)))

         ! Use linear predictor along tangent direction to start Newton iteration
         ypold(1) = one
         ypold(2:np1) = zero
         call tangnf(callbacks, s, y, yp, ypold, a, qr, alpha, tz, pivot, nfe, n, iflag)
         if (iflag .gt. 0) return

         lp: do
            w = y + h*yp
            z0 = w
            do judy = 1, litfh
               rholen = -one

               ! Calculate the Newton step 'tz' at the current point 'w'
               call tangnf(callbacks, &
                           rholen, w, wp, ypold, a, qr, alpha, tz, pivot, nfe, n, iflag)
               if (iflag .gt. 0) return

               ! Take Newton step and check convergence
               w = w + tz
               itnum = judy

               ! Compute quantities used for optimal step size estimation
               if (judy .eq. 1) then
                  lcalc = dnrm2(np1, tz, 1)
                  rcalc = rholen
                  z1 = w
               else if (judy .eq. 2) then
                  lcalc = dnrm2(np1, tz, 1)/lcalc
                  rcalc = rholen/rcalc
               end if

               ! Go to mop-up section after convergence
               if (dnrm2(np1, tz, 1) .le. relerr*dnrm2(np1, w, 1) + abserr) go to 600

            end do

            ! No convergence in litfh iterations. Reduce h and try again.
            if (h .le. fouru*(one + s)) then
               iflag = 6
               return
            end if
            h = h/2

         end do lp
      end if startup

      ! PREDICTOR SECTION
      fail = .false.
      hp: do

         ! Compute point predicted by Hermite interpolant. Use step size 'h' computed on
         ! last call to 'stepnf'.
         do j = 1, np1
            w(j) = qofs(yold(j), ypold(j), y(j), yp(j), hold, hold + h)
         end do
         z0 = w

         ! CORRECTOR SECTION
         corrector: do judy = 1, litfh

            ! Calculate the Newton step 'tz' at the current point 'w'
            rholen = -one
            call tangnf(callbacks, &
                        rholen, w, wp, yp, a, qr, alpha, tz, pivot, nfe, n, iflag)
            if (iflag .gt. 0) return

            ! Take Newton step and check convergence
            w = w + tz
            itnum = judy

            ! Compute quantities used for optimal step size estimation.
            if (judy .eq. 1) then
               lcalc = dnrm2(np1, tz, 1)
               rcalc = rholen
               z1 = w
            else if (judy .eq. 2) then
               lcalc = dnrm2(np1, tz, 1)/lcalc
               rcalc = rholen/rcalc
            end if

            ! Go to mop-up section after convergence.
            if (dnrm2(np1, tz, 1) .le. relerr*dnrm2(np1, w, 1) + abserr) go to 600

         end do corrector

         ! No convergence in 'litfh' iterations. Record failure at calculated 'h'
         ! Save this step size, reduce 'h' and try again.
         fail = .true.
         hfail = h
         if (h .le. fouru*(one + s)) then
            iflag = 6
            return
         end if
         h = h/2

      end do hp

      ! MOP-UP SECTION

      ! 'yold'  and  'y'  always contain the last two points found on the zero curve of
      ! the homotopy map. 'ypold' and 'yp' contain the tangent vectors to the zero curve
      ! at  'yold'  and  'y' , respectively.
600   ypold = yp
      yold = y
      y = w
      yp = wp
      w = y - yold

      ! Update arc length
      hold = dnrm2(np1, w, 1)
      s = s + hold

      ! OPTIMAL STEP SIZE ESTIMATION SECTION

      ! Calculate the distance factor 'dcalc'
      tz = z0 - y
      w = z1 - y
      dcalc = dnrm2(np1, tz, 1)
      if (dcalc .ne. zero) dcalc = dnrm2(np1, w, 1)/dcalc

      ! The optimal step size hbar is defined by
      !
      !   ht = hold * [min(lideal/lcalc, rideal/rcalc, dideal/dcalc)]**(1/p)
      !
      !   hbar = min [ max(ht, bmin*hold, hmin), bmax*hold, hmax ]

      ! If convergence had occurred after 1 iteration, set the contraction factor 'lcalc'
      ! to zero.
      if (itnum .eq. 1) lcalc = zero

      ! Formula for optimal step size
      if (lcalc + rcalc + dcalc .eq. zero) then
         ht = sspar(7)*hold
      else
         ht = (one/max(lcalc/sspar(1), rcalc/sspar(2), dcalc/sspar(3))) &
              **(one/sspar(8))*hold
      end if

      ! 'ht' contains the estimated optimal step size. Now put it within reasonable bounds
      h = min(max(ht, sspar(6)*hold, sspar(4)), sspar(7)*hold, sspar(5))

      if (itnum .eq. 1) then
         ! If convergence had occurred after 1 iteration, don't decrease 'h'.
         h = max(h, hold)
      else if (itnum .eq. litfh) then
         ! If convergence required the maximum 'litfh' iterations, don't increase 'h'.
         h = min(h, hold)
      end if

      ! If convergence did not occur in 'litfh' iterations for a particular 'h = hfail',
      ! don't choose the new step size larger than 'hfail'.
      if (fail) h = min(h, hfail)

   end subroutine stepnf

   subroutine tangnf( &
      callbacks, &
      rholen, y, yp, ypold, a, qr, alpha, tz, pivot, nfe, n, iflag)
   !! This subroutine builds the Jacobian matrix of the homotopy map, computes a QR
   !! decomposition of that matrix, and then calculates the (unit) tangent vector and the
   !! Newton step.

      use hompack_kinds, only: zero, one
      use blas_interfaces, only: dnrm2
      use lapack_interfaces, only: dgeqpf, dormqr
      implicit none

      type(hompack_callbacks) :: callbacks
         !! User-supplied function and Jacobian evaluation subroutines.
      real(dp), intent(inout) :: rholen
         !! Controls computation of the homotopy residual norm.
         !!
         !! On input:
         !! * `rholen < 0` requests evaluation of the homotopy map norm.
         !! * `rholen >= 0` suppresses norm evaluation.
         !!
         !! On output, if `rholen < 0` on entry:
         !! `rholen = ||rho(a, lambda, x)||`. Otherwise the value is unchanged.
      real(dp), intent(in) :: y(:)
         !! Current point on the homotopy zero curve. `Shape: (n+1)`.
         !! Contains `(lambda, x)`.
      real(dp), intent(out) :: yp(:)
         !! Unit tangent vector to the zero curve at `y`. `Shape: (n+1)`.
         !! Represents `dy/ds`, where `s` is arc length along the zero curve.
      real(dp), intent(in) :: ypold(:)
         !! Unit tangent vector at the previous point on the zero curve. `Shape: (n+1)`.
         !! Used to determine a consistent orientation for the newly computed tangent
         !! vector.
      real(dp), intent(in) :: a(:)
         !! Parameter vector used in the homotopy map.
      real(dp), intent(inout) :: qr(:, :)
         !! Workspace containing the Jacobian matrix and its QR factorization.
         !! `Shape: (n, n+2)`.
      real(dp), intent(inout) :: alpha(:)
         !! Workspace array used during QR factorization and related linear algebra
         !! operations. `Shape: (3*n+3)`.
      real(dp), intent(out) :: tz(:)
         !! Newton correction vector. `Shape: (n+1)`.
         !! Equal to the negative pseudoinverse of the homotopy Jacobian applied to the
         !! homotopy residual.
      integer, intent(inout) :: pivot(:)
         !! Pivot indices produced by the QR factorization. `Shape: (n+1)`.
      integer, intent(inout) :: nfe
         !! Number of homotopy/Jacobian evaluations performed.
         !! Incremented by one on every successful call.
      integer, intent(in) :: n
         !! Problem dimension.
      integer, intent(inout) :: iflag
         !! Problem type and return status flag.
         !!
         !! On input:
         !! * `0`  : fixed-point problem, `F(x) = x`.
         !! * `-1` : zero-finding problem, `F(x) = 0`.
         !! * `-2` : general homotopy curve-tracking problem.
         !!
         !! On output:
         !! * unchanged (`0`, `-1`, or `-2`) on normal return.
         !! * `4` : Jacobian matrix lost full rank (`rank < n`); iteration not completed.

      real(dp) :: lambda, sigma, ypnorm
      integer :: i, j, k, kp1, np1, np2

      lambda = y(1)
      np1 = n + 1
      np2 = n + 2
      nfe = nfe + 1

      ! Compute the jacobian matrix, store it and homotopy map in QR
      if (iflag .eq. -2) then

         ! QR = [ d rho(a,lambda,x)/d lambda , d rho(a,lambda,x)/dx , rho(a,lambda,x) ]
         do k = 1, np1
            call callbacks%rhojac(a, lambda, y(2:np1), qr(:, k), k, callbacks%data)
         end do
         call callbacks%rho(a, lambda, y(2:np1), qr(:, np2), callbacks%data)

      else

         call callbacks%f(y(2:np1), tz(1:n), callbacks%data)

         if (iflag .eq. 0) then

            ! QR = [ a - f(x), i - lambdadf(x), x - a + lambda(a - f(x)) ]
            qr(:, 1) = a - tz(1:n)
            qr(:, np2) = y(2:np1) - a + lambda*qr(:, 1)
            do k = 1, n
               call callbacks%fjac(y(2:np1), tz(1:n), k, callbacks%data)
               kp1 = k + 1
               qr(:, kp1) = -lambda*tz(1:n)
               qr(k, kp1) = one + qr(k, kp1)
            end do

         else

            ! QR = [ f(x) - x + a, lambda*df(x) + (1 - lambda)i , x - a + lambda(f(x) - x + a) ]
            qr(:, 1) = tz(1:n) - y(2:np1) + a
            qr(:, np2) = y(2:np1) - a + lambda*qr(:, 1)
            do k = 1, n
               call callbacks%fjac(y(2:np1), tz(1:n), k, callbacks%data)
               kp1 = k + 1
               qr(:, kp1) = lambda*tz(1:n)
               qr(k, kp1) = one - lambda + qr(k, kp1)
            end do

         end if

      end if

      ! Compute the norm of the homotopy map if it was requested
      if (rholen .lt. zero) rholen = dnrm2(n, qr(:, np2), 1)

      ! Reduce the Jacobian matrix to upper triangular form
      pivot = 0
      call dgeqpf(n, np1, qr, n, pivot, yp, alpha, k)

      if (abs(qr(n, n)) .le. abs(qr(1, 1))*epsilon(one)) then
         iflag = 4
         return
      end if

      call dormqr('L', 'T', n, 1, n, qr, n, yp, qr(:, np2), n, alpha, 3*n + 3, k)

      do i = 1, n
         alpha(i) = qr(i, i)
      end do

      ! Compute kernel of Jacobian, which specifies yp=dy/ds.
      tz(np1) = one
      do i = n, 1, -1
         j = i + 1
         tz(i) = -dot_product(qr(i, j:np1), tz(j:np1))/alpha(i)
      end do
      ypnorm = dnrm2(np1, tz, 1)
      yp(pivot) = tz/ypnorm
      if (dot_product(yp, ypold) .lt. zero) yp = -yp

      ! 'yp' is the unit tangent vector in the correct direction
      !
      ! Compute the minimum norm solution of [d rho(y(s))] v = -rho(y(s)).
      ! 'v' is given by 'p - (p,q)q', where 'p' is any solution of [d rho] v = -rho
      ! and 'q' is a unit vector in the kernel of [d rho].
      alpha(np1) = one
      do i = n, 1, -1
         j = i + 1
         alpha(i) = -(dot_product(qr(i, j:np1), alpha(j:np1)) + qr(i, np2))/alpha(i)
      end do
      tz(pivot) = alpha(1:np1)

      ! 'tz' now contains a particular solution 'p', and 'yp' contains a vector 'q'
      ! in the kernel (the tangent)
      sigma = dot_product(tz, yp)

      ! 'tz' is the Newton step from the current point y(s) = (lambda(s), x(s)).
      tz = tz - sigma*yp

   end subroutine tangnf

end module hompack_nf
