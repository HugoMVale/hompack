module hompack_nf
!! Specific routines for the [[fixpnf]] solver.

   use iso_c_binding, only: c_ptr, c_null_ptr
   use iso_fortran_env, only: output_unit
   use hompack_kinds, only: dp
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

   type fixnpf_workspace
   !! Linear-algebra workspace for [[fixpnf]].
      real(dp), allocatable :: alpha(:)
         !! Array used during interpolation and Newton-step calculations.
      real(dp), allocatable :: qr(:, :)
         !! Matrix for QR factorizations used in Newton-step calculations.
      real(dp), allocatable :: tz(:)
         !! Array used in QR-factorization and Newton-step computations.
      real(dp), allocatable :: w(:)
         !! Array used during interpolation and Newton-step calculations.
      real(dp), allocatable :: wp(:)
         !! Array used during interpolation and Newton-step calculations.
      real(dp), allocatable :: z0(:)
         !! Array used for estimating the optimal next step size.
      real(dp), allocatable :: z1(:)
         !! Array used for estimating the optimal next step size.
      integer, allocatable :: pivot(:)
         !! Pivot indices used by the QR factorization.
   contains
      procedure :: init => init_workspace
   end type fixnpf_workspace

   type fixnpf_state
   !! State variables for [[fixpnf]].
      real(dp) :: abserr
         !! Absolute error tolerance.
      real(dp) :: relerr
         !! Relative error tolerance.
      real(dp) :: curtol
         !! Curvature-based error tolerance.
      real(dp) :: h
         !! Optimal step size for the next step to be attempted by [[stepnf]].
      real(dp) :: hold
         !! ||yp - ypold|| at the previous step.
      real(dp) :: s
         !! Total arc length of the solution path followed by the algorithm.
      integer :: iflag
         !! Problem type and status flag.
      integer :: limit
         !! Limit on the number of steps allowed in the main loop of [[fixpnf]].
      integer :: n
         !! Problem dimension.
      integer :: nfe
         !! Number of homotopy/Jacobian evaluations performed.
         !! This counter is incremented once for each call to `tangnf`, where the
         !! homotopy map and its Jacobian are assembled and used to compute a tangent
         !! vector and Newton correction.
      integer :: lunit
         !! Logical I/O unit for intermediate output.
      logical :: start
         !! Flag to indicate the first call to [[stepnf]].
      logical :: crash
         !! Flag to indicate that the return state of [[stepnf]] is a crash.
      logical :: ispoly
         !! Flag to indicate polynomial mode.
      character(:), allocatable :: message
         !! Error message.
      real(dp) :: sspar(8)
         !! Step-size control parameters.
         !! `(lideal, rideal, dideal, hmin, hmax, bmin, bmax, p)`
      real(dp), allocatable :: a(:)
      real(dp), allocatable :: y(:)
         !! Current point on the zero curve.
      real(dp), allocatable :: yold(:)
         !! Previous point found on the zero curve.
      real(dp), allocatable :: yp(:)
         !! Unit tangent vector to the zero curve at `y`.
      real(dp), allocatable :: ypold(:)
         !! Unit tangent vector to the zero curve at `yold`.
      type(fixnpf_workspace) :: workspace
         !! Linear-algebra workspace.
   contains
      procedure :: init => init_state
   end type fixnpf_state

contains

   subroutine fixpnf( &
      state, callbacks, n, x, iflag, &
      arcre, arcae, ansre, ansae, sspar, a, lunit, ispoly)
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
      !    ARC?E <= 0.0  on input it is reset to  .5*SQRT(ANS?E) .
      !    Normally  ARC?E should be considerably larger than  ANS?E .
      !
      ! ANSRE , ANSAE  are the relative and absolute error values used for
      !    the answer at LAMBDA = 1.  The accepted answer  Y = (LAMBDA, X)
      !    satisfies
      !
      !       |Y(1) - 1|  <=  ANSRE + ANSAE           .AND.
      !
      !       ||Z||  <=  ANSRE*||X|| + ANSAE          where
      !
      !    (.,Z) is the Newton step to Y.
      !
      ! TRACE  is an integer specifying the logical I/O unit for
      !    intermediate output.  If  TRACE > 0  the points computed on
      !    the zero curve are written to I/O unit  TRACE .
      !
      ! A(:)  contains the parameter vector  A .  For the fixed point
      !    and zero finding problems, A  need not be initialized by the
      !    user, and is assumed to have length  N.  For the curve
      !    tracking problem, A  must be initialized by the user.
      !
      ! SSPAR(1:8) = (LIDEAL, RIDEAL, DIDEAL, HMIN, HMAX, BMIN, BMAX, P)  is
      !    a vector of parameters used for the optimal step size estimation.
      !    If  SSPAR(J) <= 0.0  on input, it is reset to a default value
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

      use hompack_kinds, only: zero, one, eps64
      implicit none

      type(fixnpf_state), intent(inout) :: state
         !! State variables for [[fixpnf]]. Initialized on the first call, and updated on
         !! subsequent calls when `iflag=2` or `iflag=3`.
      type(hompack_callbacks), intent(in) :: callbacks
         !! User-supplied function and Jacobian evaluation subroutines.
      integer, intent(in) :: n
         !! Problem dimension.
      real(dp), intent(inout) :: x(:)
         !! On input, `x` is the initial point `a` for the fixed-point and zero-finding
         !! problems, and the initial solution `x0` for the curve-tracking problem.
         !! On output, `x` is the approximate solution of the problem at `lambda = 1`.
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
      real(dp), intent(inout), optional :: sspar(8)
         !! Step-size control parameters:
         !! `(lideal, rideal, dideal, hmin, hmax, bmin, bmax, p)`.
         !! Parameters used by the optimal step-size estimation algorithm.
         !! Elements that are nonpositive on input are replaced by default values.
      real(dp), intent(in), optional :: a(:)
         !! Parameter vector `a`.
         !! For fixed-point and zero-finding problems, the array is assumed to have
         !! length `n` and need not be initialized by the user.
         !! For curve-tracking problems, it must be initialized on input.
         !! Unchanged on output.
      integer, intent(in), optional :: lunit
         !! Logical I/O unit for intermediate output.
         !! * `0` : No output is printed (default).
         !! * `6` : Standard output (mapped internally to `output_unit`)
         !! * otherwise : Output to the specified I/O unit.
      logical, intent(in), optional :: ispoly
         !! Optional flag used only by the polynomial-system driver [[POLSYS1H]].

      ! Switch from the tolerance arc?e to the (finer) tolerance ans?e if the curvature
      ! of any component of y exceeds cursw
      real(dp), parameter :: cursw = 10.0_dp

      ! Upper bound on the number of steps
      integer, parameter :: limitd = 1000

      integer :: ierr, iter, np1, dima

      np1 = n + 1
      state%message = ""

      if (n <= 0) then
         state%message = "Illegal input: `n` must be greater or equal than one."
         iflag = 7
         return
      end if

      if (size(x) /= n) then
         state%message = "Illegal input: length of `x` must be equal to `n`."
         iflag = 7
         return
      end if

      if (ansre <= zero) then
         state%message = "Illegal input: `ansre` must be greater than zero."
         iflag = 7
         return
      end if

      if (ansae < zero) then
         state%message = "Illegal input: `ansae` must be greater or equal than zero."
         iflag = 7
         return
      end if

      if (iflag == 0 .or. iflag == -1) then
         ! First run for fixed-point or zero-finding problem
         if (.not. associated(callbacks%f)) then
            state%message = "Illegal input: callback `f` must be provided for fixed-point and zero-finding problems."
            iflag = 7
            return
         end if
         if (.not. associated(callbacks%fjac)) then
            state%message = "Illegal input: callback `fjac` must be provided for fixed-point and zero-finding problems."
            iflag = 7
            return
         end if
         dima = n
         go to 20
      else if (iflag == -2) then
         ! First run for curve-tracking problem
         if (.not. associated(callbacks%rho)) then
            state%message = "Illegal input: callback `rho` must be provided for curve-tracking problems."
            iflag = 7
            return
         end if
         if (.not. associated(callbacks%rhojac)) then
            state%message = "Illegal input: callback `rhojac` must be provided for curve-tracking problems."
            iflag = 7
            return
         end if
         if (.not. present(a)) then
            state%message = "Illegal input: parameter vector `a` must be provided for curve-tracking problems."
            iflag = 7
            return
         end if
         dima = size(a)
         if (dima <= 0) then
            state%message = "Illegal input: length of parameter vector `a` must be greater than zero."
            iflag = 7
            return
         end if
         go to 20
      else if (iflag == 2) then
         ! Restart after error tolerances were increased
         if (state%n /= n .or. state%nfe < 1) then
            state%message = "Restart error: possible corrupted state or invalid `iflag` value."
            iflag = 7
            return
         end if
         go to 120
      else if (iflag == 3) then
         ! Restart after iteration limit was reached
         if (state%n /= n .or. state%nfe < 1) then
            state%message = "Restart error: possible corrupted state or invalid `iflag` value."
            iflag = 7
            return
         end if
         go to 90
      else
         state%message = "Illegal input: `iflag` must be -2, -1, 0, 2, or 3."
         iflag = 7
         return
      end if

      ! INITIALIZATION (FIRST CALL ONLY)

      ! Initialize state
20    call state%init(n, dima, ierr)

      if (ierr /= 0) then
         iflag = 8
         return
      end if

      state%iflag = iflag

      state%y(1) = zero
      state%y(2:np1) = x
      state%yp(1) = one
      state%ypold(1) = one

      ! Load 'a'
      if (state%iflag == 0 .or. state%iflag == -1) then
         state%a = x
      else
         state%a = a
      end if

      ! Default arc tolerances
      if (arcre <= zero) arcre = sqrt(ansre)/2
      if (arcae <= zero) arcae = sqrt(ansae)/2

      ! Set optimal step size estimation parameters:
      ! `(lideal, rideal, dideal, hmin, hmax, bmin, bmax, p)`.
      ! Let 'z[k]' denote the Newton iterates along the flow normal to the Davidenko flow
      ! and 'y' their limit
      if (present(sspar)) state%sspar = sspar
      ! Ideal contraction factor: ||z[2] - z[1]|| / ||z[1] - z[0]||
      if (state%sspar(1) <= zero) state%sspar(1) = 0.5_dp
      ! Ideal residual factor: ||rho(A, z[1])|| / ||rho(a, z[0])||
      if (state%sspar(2) <= zero) state%sspar(2) = 0.01_dp
      ! Ideal distance factor: ||z[1] - y|| / ||z[0] - y||
      if (state%sspar(3) <= zero) state%sspar(3) = 0.5_dp
      ! Minimum step size 'hmin'
      if (state%sspar(4) <= zero) state%sspar(4) = (sqrt(real(n + 1, dp)) + 4.0_dp)*eps64
      ! Maximum step size 'hmax'
      if (state%sspar(5) <= zero) state%sspar(5) = one
      ! Minimum step size reduction factor 'bmin'
      if (state%sspar(6) <= zero) state%sspar(6) = 0.1_dp
      ! Maximum step size expansion factor 'bmax'
      if (state%sspar(7) <= zero) state%sspar(7) = 3.0_dp
      ! Assumed operating order 'p'
      if (state%sspar(8) <= zero) state%sspar(8) = 2.0_dp

      ! Set the trace output unit
      if (present(lunit)) then
         if (lunit == 6) then
            state%lunit = output_unit
         else
            state%lunit = lunit
         end if
      end if

      ! Set special mode for 'polsys1h'
      if (present(ispoly)) then
         state%ispoly = ispoly
      else
         state%ispoly = .false.
      end if

      ! COMMON PART FOR FIRST CALL AND RESTARTS

      ! Set default iteration limit
90    state%limit = limitd

      ! Main loop
120   do iter = 1, state%limit

         ! Tracking algorithm lost the zero curve
         if (state%y(1) < zero) then
            iflag = 5
            return
         end if

         ! Set different error tolerance if the trajectory y(s) has any high curvature
         ! components
         state%curtol = cursw*state%hold
         state%relerr = arcre
         state%abserr = arcae
         if (any(abs(state%yp - state%ypold) > state%curtol)) then
            state%relerr = ansre
            state%abserr = ansae
         end if

         ! Take a step along the curve
         call stepnf(state, callbacks)

         ! Print latest point on curve if requested
         if (state%lunit /= 0) then
            write (lunit, 217) iter, state%nfe, state%s, state%y(1), state%y(2:np1)
217         format(/' STEP', i5, 3x, 'NFE =', i5, 3x, 'ARC LENGTH =', f9.4, 3x, &
                    'LAMBDA =', f7.4, 5x, 'X VECTOR:'/(1x, 6es12.4))
         end if

         ! Check if the step was successful
         if (state%iflag > 0) then
            iflag = state%iflag
            return
         end if

         if (state%crash) then
            ! Return code for error tolerance too small
            iflag = 2
            ! Change error tolerances
            if (arcre < state%relerr) arcre = state%relerr
            if (ansre < state%relerr) ansre = state%relerr
            if (arcae < state%abserr) arcae = state%abserr
            if (ansae < state%abserr) ansae = state%abserr
            ! Change limit on number of iterations
            state%limit = state%limit - iter
            return
         end if

         ! Use Hermite cubic interpolation and Newton iteration to get the answer at lambda=1
         if (state%y(1) >= one) then

            associate (ws => state%workspace)

               ! Save 'yold' for arc length calculation later
               ws%z0 = state%yold
               call rootnf(state, callbacks, ansre, ansae)
               iflag = 1

               ! Set error flag if 'rootnf' could not get the point on the zero curve at
               ! lambda=1
               if (state%iflag > 0) iflag = state%iflag

               ! Calculate final arc length
               ws%w = state%y - ws%z0
               state%s = state%s - state%hold + norm2(ws%w)
               return

            end associate

         end if

         ! For polynomial systems and the 'polsys1h' homotopy map, dlambda/ds>= 0
         ! necessarily; this condition is enforced here
         if (state%ispoly) then
            if (state%yp(1) < zero) then
               ! Reverse tangent direction so D LAMBDA/DS = YP(1) > 0
               state%yp = -state%yp
               state%ypold = state%yp
               ! Force 'stepnf' to use the linear predictor for the next step only
               state%start = .true.
            end if
         end if

      end do

      ! Lambda has not reached 1 in 'limitd' steps
      iflag = 3

   end subroutine fixpnf

   subroutine rootnf(state, callbacks, relerr, abserr)
   !! This subroutine finds the point `ybar = (1, xbar)` on the zero curve of the homotopy
   !! map. It starts with two points `yold = (lambdaold, xold)` and `y = (lambda, x)` such
   !! that `lambdaold < 1 <= lambda` , and alternates between secant estimates of `ybar`
   !! and Newton iteration until convergence.

      use hompack_kinds, only: zero, one, eps64
      use hompack_core, only: root, root_state, qofs
      implicit none

      type(fixnpf_state), intent(inout) :: state
         !! State variables for [[fixpnf]].
      type(hompack_callbacks) :: callbacks
         !! User-supplied function and Jacobian evaluation subroutines.
      real(dp), intent(in) :: relerr
         !! Relative convergence tolerance.
         !! Iteration is considered converged when `|y(1)-1| <= relerr + abserr` and the
         !! Newton correction satisfies `||z|| <= relerr*||x|| + abserr`.
      real(dp), intent(in) :: abserr
         !! Absolute convergence tolerance.
         !! Used together with `relerr` in the convergence criteria.

      real(dp) :: dels, qsout, aerr, rerr, sa, sb, sout
      integer :: judy, jw, lcode, limit, np1
      logical :: bracket
      type(root_state) :: state_root

      rerr = max(relerr, eps64)
      aerr = max(abserr, zero)
      np1 = state%n + 1

      ! The limit on the number of iterations allowed may be changed by changing the
      ! following parameter statement
      limit = 2*(int(abs(log10(aerr + rerr))) + 1)

      associate (ws => state%workspace)

         ws%tz = state%y - state%yold
         dels = norm2(ws%tz)

         ! Using two points and tangents on the homotopy zero curve, construct the Hermite
         ! cubic interpolant q(s). Then use 'root' to find the 's' corresponding to
         ! 'lambda=1'. The two points on the zero curve are always chosen to bracket
         ! 'lambda=1', with the bracketing interval always being [0, dels].
         sa = zero
         sb = dels
         lcode = 1 ! forces initialization of 'root'
         do
            call root(sout, qsout, sa, sb, rerr, aerr, lcode, state_root)
            if (lcode > 0) exit
            qsout = qofs(state%yold(1), state%ypold(1), state%y(1), state%yp(1), dels, sout) - one
         end do

         ! If lambda=1 were bracketed, root cannot fail
         if (lcode > 2) then
            state%iflag = 6
            return
         end if

         ! Calculate 'q(sa)' as the initial point for a Newton iteration
         do jw = 1, np1
            ws%w(jw) = qofs(state%yold(jw), state%ypold(jw), state%y(jw), state%yp(jw), dels, sa)
         end do

         ! Tangent information 'yp' is no longer needed. Hereafter, 'yp' represents the most
         ! recent point which is on the opposite side of the hyperplane 'lambda=1' from 'y'

         ! Prepare for main loop
         state%yp = state%yold

         ! Initialize bracket to indicate that the points 'y' and 'yold' bracket 'lambda=1',
         ! thus 'yold = yp'
         bracket = .true.

         ! Main loop
         do judy = 1, limit

            ! Calculate Newton step at current estimate 'w'
            call tangnf(callbacks, sa, &
                        state%n, ws%w, ws%wp, state%ypold, state%a, state%nfe, state%iflag, &
                        ws%qr, ws%alpha, ws%tz, ws%pivot)
            if (state%iflag > 0) return

            ! Next point = current point + Newton step
            ws%w = ws%w + ws%tz

            ! Check for convergence
            if ((abs(ws%w(1) - one) <= rerr + aerr) .and. &
                (norm2(ws%tz) <= rerr*norm2(ws%w(2:np1)) + aerr)) then
               state%y = ws%w
               return
            end if

            ! Prepare for next iteration
            if (abs(ws%w(1) - one) <= rerr + aerr) then
               state%ypold = ws%wp
               cycle
            end if

            ! Update 'y' and 'yold'
            state%yold = state%y
            state%y = ws%w

            ! Update 'yp' such that 'yp' is the most recent point opposite of 'lambda=1'
            ! from 'y'. Set bracket=.true. iff 'y' and 'yold' bracket 'lambda=1' so that
            ! yp = yold .
            if ((state%y(1) - one)*(state%yold(1) - one) > 0) then
               bracket = .false.
            else
               bracket = .true.
               state%yp = state%yold
            end if

            ! Compute dels=||y-yp||
            ws%tz = state%y - state%yp
            dels = norm2(ws%tz)

            ! Compute tz for the linear predictor w = y + tz, where tz = sa*(yold-y).
            sa = (one - state%y(1))/(state%yold(1) - state%y(1))
            ws%tz = sa*(state%yold - state%y)

            ! To insure stability, the linear prediction must be no farther from y than
            ! yp is. This is guaranteed if bracket=true. If linear prediction is too far
            ! away, use bracketing points to compute linear prediction.
            if (.not. bracket) then
               if (norm2(ws%tz) > dels) then
                  ! Compute tz = sa*(yp-y)
                  sa = (one - state%y(1))/(state%yp(1) - state%y(1))
                  ws%tz = sa*(state%yp - state%y)
               end if
            end if

            ! Compute estimate w = y + tz  and save old tangent vector.
            ws%w = ws%w + ws%tz
            state%ypold = ws%wp

         end do

         ! The alternating secant estimation and Newton iteration has not converged in
         ! 'limit' iterations
         state%iflag = 6

      end associate

   end subroutine rootnf

   subroutine stepnf(state, callbacks)
   !! This subroutine takes one step along the zero curve of the homotopy map using a
   !! predictor-corrector algorithm. The predictor uses a Hermite cubic interpolant, and
   !! the corrector returns to the zero curve along the flow normal to the Davidenko flow.
   !! [[stepnf]] also estimates a step size `h` for the next step along the zero curve.
   !! Normally, [[stepnf]] is used indirectly through [[fixpnf]], and should be called
   !! directly only if it is necessary to modify the stepping algorithm's parameters.

      use hompack_kinds, only: one, zero, eps64
      use hompack_core, only: qofs
      implicit none

      type(fixnpf_state), intent(inout) :: state
         !! State variables for [[fixpnf]].
      type(hompack_callbacks), intent(in) :: callbacks
         !! User-supplied function and Jacobian evaluation subroutines.

      real(dp), parameter :: twou = 2*eps64, fouru = 4*eps64
      real(dp) :: dcalc, hfail, ht, lcalc, rcalc, rholen, temp
      integer :: itnum, j, judy, np1
      logical :: fail

      ! The limit on the number of Newton iterations allowed before reducing the step
      ! size 'h' may be changed by changing the following parameter
      integer, parameter:: litfh = 4

      np1 = state%n + 1
      state%crash = .true.

      associate (ws => state%workspace, sspar => state%sspar)

         ! The arclength 's' must be nonnegative
         if (state%s < zero) return

         ! If step size is too small, determine an acceptable one
         if (state%h < fouru*(one + state%s)) then
            state%h = fouru*(one + state%s)
            return
         end if

         ! If error tolerances are too small, increase them to acceptable values
         temp = norm2(state%y) + one
         if (0.5_dp*(state%relerr*temp + state%abserr) < twou*temp) then
            if (state%relerr /= zero) then
               state%relerr = fouru*(one + fouru)
               state%abserr = max(state%abserr, zero)
            else
               state%abserr = fouru*temp
            end if
            return
         end if

         ! STARTUP SECTION (FIRST STEP ALONG ZERO CURVE)

         state%crash = .false.
         if (state%start) then
            fail = .false.
            state%start = .false.

            ! Determine suitable initial step size
            state%h = min(state%h, 0.1_dp, sqrt(sqrt(state%relerr*temp + state%abserr)))

            ! Use linear predictor along tangent direction to start Newton iteration
            state%ypold(1) = one
            state%ypold(2:np1) = zero
            call tangnf(callbacks, state%s, state%n, state%y, state%yp, state%ypold, state%a, &
                        state%nfe, state%iflag, &
                        ws%qr, ws%alpha, ws%tz, ws%pivot)

            if (state%iflag > 0) return

            do

               ws%w = state%y + state%h*state%yp
               ws%z0 = ws%w

               do judy = 1, litfh

                  ! Calculate the Newton step 'tz' at the current point 'w'
                  rholen = -one
                  call tangnf(callbacks, rholen, &
                              state%n, ws%w, ws%wp, state%ypold, state%a, &
                              state%nfe, state%iflag, &
                              ws%qr, ws%alpha, ws%tz, ws%pivot)
                  if (state%iflag > 0) return

                  ! Take Newton step and check convergence
                  ws%w = ws%w + ws%tz
                  itnum = judy

                  ! Compute quantities used for optimal step size estimation
                  if (judy == 1) then
                     lcalc = norm2(ws%tz)
                     rcalc = rholen
                     ws%z1 = ws%w
                  else if (judy == 2) then
                     lcalc = norm2(ws%tz)/lcalc
                     rcalc = rholen/rcalc
                  end if

                  ! Go to mop-up section after convergence
                  if (norm2(ws%tz) <= state%relerr*norm2(ws%w) + state%abserr) then
                     go to 600
                  end if

               end do

               ! No convergence in litfh iterations. Reduce h and try again.
               if (state%h <= fouru*(one + state%s)) then
                  state%iflag = 6
                  return
               end if
               state%h = state%h/2

            end do

         end if

         ! PREDICTOR SECTION

         fail = .false.

         do

            ! Compute point predicted by Hermite interpolant. Use step size 'h' computed on
            ! last call to 'stepnf'.
            do j = 1, np1
               ws%w(j) = qofs(state%yold(j), state%ypold(j), state%y(j), state%yp(j), &
                              state%hold, state%hold + state%h)
            end do

            ws%z0 = ws%w

            ! CORRECTOR SECTION
            do judy = 1, litfh

               ! Calculate the Newton step 'tz' at the current point 'w'
               rholen = -one
               call tangnf(callbacks, rholen, &
                           state%n, ws%w, ws%wp, state%yp, state%a, state%nfe, &
                           state%iflag, ws%qr, ws%alpha, ws%tz, ws%pivot)
               if (state%iflag > 0) return

               ! Take Newton step and check convergence
               ws%w = ws%w + ws%tz
               itnum = judy

               ! Compute quantities used for optimal step size estimation.
               if (judy == 1) then
                  lcalc = norm2(ws%tz)
                  rcalc = rholen
                  ws%z1 = ws%w
               else if (judy == 2) then
                  lcalc = norm2(ws%tz)/lcalc
                  rcalc = rholen/rcalc
               end if

               ! Go to mop-up section after convergence
               if (norm2(ws%tz) <= state%relerr*norm2(ws%w) + state%abserr) then
                  go to 600
               end if

            end do

            ! No convergence in 'litfh' iterations. Record failure at calculated 'h'
            ! Save this step size, reduce 'h' and try again.

            fail = .true.
            hfail = state%h

            if (state%h <= fouru*(one + state%s)) then
               state%iflag = 6
               return
            end if

            state%h = state%h/2

         end do

         ! MOP-UP SECTION

         ! 'yold'  and  'y'  always contain the last two points found on the zero curve of
         ! the homotopy map. 'ypold' and 'yp' contain the tangent vectors to the zero curve
         ! at  'yold'  and  'y' , respectively.
600      state%ypold = state%yp
         state%yold = state%y
         state%y = ws%w
         state%yp = ws%wp
         ws%w = state%y - state%yold

         ! Update arc length
         state%hold = norm2(ws%w)
         state%s = state%s + state%hold

         ! OPTIMAL STEP SIZE ESTIMATION SECTION

         ! Calculate the distance factor 'dcalc'
         ws%tz = ws%z0 - state%y
         ws%w = ws%z1 - state%y
         dcalc = norm2(ws%tz)
         if (dcalc /= zero) dcalc = norm2(ws%w)/dcalc

         ! The optimal step size hbar is defined by
         !
         !   ht = hold * [min(lideal/lcalc, rideal/rcalc, dideal/dcalc)]**(1/p)
         !
         !   hbar = min [ max(ht, bmin*hold, hmin), bmax*hold, hmax ]

         ! If convergence had occurred after 1 iteration, set the contraction factor 'lcalc'
         ! to zero.
         if (itnum == 1) lcalc = zero

         ! Formula for optimal step size
         if (lcalc + rcalc + dcalc == zero) then
            ht = sspar(7)*state%hold
         else
            ht = (one/max(lcalc/sspar(1), rcalc/sspar(2), &
                          dcalc/sspar(3)))**(one/sspar(8))*state%hold
         end if

         ! 'ht' contains the estimated optimal step size. Put it within reasonable bounds
         state%h = min( &
                   max(ht, sspar(6)*state%hold, sspar(4)), &
                   sspar(7)*state%hold, &
                   sspar(5) &
                   )

         if (itnum == 1) then
            ! If convergence had occurred after 1 iteration, don't decrease 'h'.
            state%h = max(state%h, state%hold)
         else if (itnum == litfh) then
            ! If convergence required the maximum 'litfh' iterations, don't increase 'h'.
            state%h = min(state%h, state%hold)
         end if

         ! If convergence did not occur in 'litfh' iterations for a particular 'h = hfail',
         ! don't choose the new step size larger than 'hfail'.
         if (fail) state%h = min(state%h, hfail)

      end associate

   end subroutine stepnf

   subroutine tangnf( &
      callbacks, rholen, n, y, yp, ypold, a, nfe, iflag, qr, alpha, tz, pivot)
   !! This subroutine builds the Jacobian matrix of the homotopy map, computes a QR
   !! decomposition of that matrix, and then calculates the (unit) tangent vector and the
   !! Newton step.

      use hompack_kinds, only: zero, one, eps64
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
      integer, intent(in) :: n
         !! Problem dimension.
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
      integer, intent(inout) :: nfe
         !! Number of homotopy/Jacobian evaluations performed.
         !! Incremented by one on every successful call.
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

      real(dp) :: lambda, sigma, ypnorm
      integer :: i, j, k, kp1, np1, np2, info

      lambda = y(1)
      np1 = n + 1
      np2 = n + 2
      nfe = nfe + 1

      ! Compute the jacobian matrix, store it and homotopy map in QR
      if (iflag == -2) then

         ! QR = [ d rho(a,lambda,x)/d lambda , d rho(a,lambda,x)/dx , rho(a,lambda,x) ]
         do k = 1, np1
            call callbacks%rhojac(a, lambda, y(2:np1), qr(:, k), k, callbacks%data)
         end do

         call callbacks%rho(a, lambda, y(2:np1), qr(:, np2), callbacks%data)

      else

         call callbacks%f(y(2:np1), tz(1:n), callbacks%data)

         if (iflag == 0) then

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
      if (rholen < zero) rholen = norm2(qr(:, np2))

      ! Reduce the Jacobian matrix to upper triangular form
      pivot = 0
      call dgeqpf(n, np1, qr, n, pivot, yp, alpha, info)

      if (abs(qr(n, n)) <= abs(qr(1, 1))*eps64) then
         iflag = 4
         return
      end if

      call dormqr('L', 'T', n, 1, n, qr, n, yp, qr(:, np2), n, alpha, 3*n + 3, info)

      do i = 1, n
         alpha(i) = qr(i, i)
      end do

      ! Compute kernel of Jacobian, which specifies yp=dy/ds.
      tz(np1) = one
      do i = n, 1, -1
         j = i + 1
         tz(i) = -dot_product(qr(i, j:np1), tz(j:np1))/alpha(i)
      end do
      ypnorm = norm2(tz)
      yp(pivot) = tz/ypnorm
      if (dot_product(yp, ypold) < zero) yp = -yp

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

   pure subroutine init_workspace(self, n, stat)
   !! Initializes [[fixnpf_workspace]].

      use hompack_kinds, only: zero

      class(fixnpf_workspace), intent(inout) :: self
         !! Workspace.
      integer, intent(in) :: n
         !! Problem dimension.
      integer, intent(out), optional :: stat
         !! Error status of the allocation.

      integer :: ierr(8)

      if (present(stat)) stat = 0

      ! Deallocate any previously allocated arrays
      if (allocated(self%alpha)) deallocate (self%alpha)
      if (allocated(self%qr)) deallocate (self%qr)
      if (allocated(self%tz)) deallocate (self%tz)
      if (allocated(self%w)) deallocate (self%w)
      if (allocated(self%wp)) deallocate (self%wp)
      if (allocated(self%z0)) deallocate (self%z0)
      if (allocated(self%z1)) deallocate (self%z1)
      if (allocated(self%pivot)) deallocate (self%pivot)

      ! Allocate/initialize workspace arrays
      allocate (self%alpha(3*n + 3), source=zero, stat=ierr(1))
      allocate (self%qr(n, n + 2), source=zero, stat=ierr(2))
      allocate (self%tz(n + 1), source=zero, stat=ierr(3))
      allocate (self%w(n + 1), source=zero, stat=ierr(4))
      allocate (self%wp(n + 1), source=zero, stat=ierr(5))
      allocate (self%z0(n + 1), source=zero, stat=ierr(6))
      allocate (self%z1(n + 1), source=zero, stat=ierr(7))
      allocate (self%pivot(n + 1), source=0, stat=ierr(8))

      if (any(ierr /= 0)) then
         if (present(stat)) then
            stat = ierr(findloc(ierr /= 0, .true., dim=1))
         else
            error stop "Error: Allocation failed in fixnpf_workspace%init()."
         end if
      end if

   end subroutine init_workspace

   pure subroutine init_state(self, n, dima, stat)
   !! Initializes [[fixnpf_state]].

      use hompack_kinds, only: zero, one

      class(fixnpf_state), intent(inout) :: self
         !! State.
      integer, intent(in)  :: n
         !! Problem dimension.
      integer, intent(in)  :: dima
         !! Dimension of the parameter vector `a`.
      integer, intent(out), optional   :: stat
         !! Error status of the allocation.

      integer :: ierr(6)

      if (present(stat)) stat = 0

      ! Initialize scalar state variables
      self%abserr = zero
      self%curtol = zero
      self%h = 0.1_dp
      self%hold = zero
      self%s = zero
      self%iflag = 0
      self%limit = 0
      self%n = n
      self%nfe = 0
      self%lunit = 0
      self%ispoly = .false.
      self%start = .true.
      self%crash = .false.
      self%sspar = -one
      self%message = ""

      ! Deallocate any previously allocated arrays
      if (allocated(self%y)) deallocate (self%y)
      if (allocated(self%yold)) deallocate (self%yold)
      if (allocated(self%yp)) deallocate (self%yp)
      if (allocated(self%ypold)) deallocate (self%ypold)
      if (allocated(self%a)) deallocate (self%a)

      ! Allocate/initialize state arrays
      allocate (self%y(n + 1), source=zero, stat=ierr(1))
      allocate (self%yold(n + 1), source=zero, stat=ierr(2))
      allocate (self%yp(n + 1), source=zero, stat=ierr(3))
      allocate (self%ypold(n + 1), source=zero, stat=ierr(4))
      allocate (self%a(dima), source=zero, stat=ierr(5))

      ! Initialize workspace
      call self%workspace%init(n, stat=ierr(6))

      if (any(ierr /= 0)) then
         if (present(stat)) then
            stat = ierr(findloc(ierr /= 0, .true., dim=1))
         else
            error stop "Error: Allocation failed in fixnpf_state%init()."
         end if
      end if

   end subroutine init_state

end module hompack_nf
