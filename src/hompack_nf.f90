module hompack_nf
!! Specific routines for the [[fixpnf]] solver.

   use iso_c_binding, only: c_ptr, c_null_ptr
   use iso_fortran_env, only: output_unit
   use hompack_kinds, only: dp, zero, one, eps64
   use hompack_core
   implicit none

   integer, parameter :: fixpnf_success = hompack_success
   integer, parameter :: fixpnf_tol_increased = 2
   integer, parameter :: fixpnf_iter_limit = 3
   integer, parameter :: fixpnf_rank_loss = 4
   integer, parameter :: fixpnf_lost_curve = 5
   integer, parameter :: fixpnf_newton_failed = 6
   integer, parameter :: fixpnf_input_error = 7
   integer, parameter :: fixpnf_memory_error = 8

   type :: nf_workspace
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
      procedure :: alloc => alloc_workspace
   end type

   type :: nf_state
   !! State variables for [[fixpnf]].
      real(dp) :: abserr = zero
         !! Absolute error tolerance.
      real(dp) :: relerr = zero
         !! Relative error tolerance.
      real(dp) :: curtol = zero
         !! Curvature-based error tolerance.
      real(dp) :: h = zero
         !! Optimal step size for the next step to be attempted by [[stepnf]].
      real(dp) :: hold = zero
         !! ||yp - ypold|| at the previous step.
      real(dp) :: s = zero
         !! Total arc length of the solution path followed by the algorithm.
      integer :: iflag = -100
         !! Problem type and status flag.
      integer :: limit = -100
         !! Limit on the number of steps allowed in the main loop of [[fixpnf]].
      integer :: n = 0
         !! Problem dimension.
      integer :: nfe = 0
         !! Number of homotopy/Jacobian evaluations performed.
         !! This counter is incremented once for each call to `tangnf`, where the
         !! homotopy map and its Jacobian are assembled and used to compute a tangent
         !! vector and Newton correction.
      logical :: start = .true.
         !! Flag to indicate the first call to [[stepnf]].
      logical :: crash = .false.
         !! Flag to indicate that the return state of [[stepnf]] is a crash.
      real(dp), allocatable :: a(:)
         !! Parameter vector \(a\) for curve-tracking problems, or initial point for the
         !! fixed-point and zero-finding problems.
      real(dp), allocatable :: y(:)
         !! Current point on the zero curve.
      real(dp), allocatable :: yold(:)
         !! Previous point found on the zero curve.
      real(dp), allocatable :: yp(:)
         !! Unit tangent vector to the zero curve at `y`.
      real(dp), allocatable :: ypold(:)
         !! Unit tangent vector to the zero curve at `yold`.
      type(nf_workspace) :: workspace
         !! Linear-algebra workspace.
   contains
      procedure :: alloc => alloc_state
   end type

   type :: nf_config
   !! Container to hold the configuration for a call to [[fixpnf]].
      real(dp) :: arcre = zero
         !! Relative error tolerance for the normal-flow iteration used while tracking the
         !! zero curve. If nonpositive, the default value `arcre=0.5*sqrt(ansre)` is used.
      real(dp) :: arcae = zero
         !! Absolute error tolerance for the normal-flow iteration used while tracking the
         !! zero curve. If nonpositive, the default value `arcae=0.5*sqrt(ansae)` is used.
      real(dp) :: ansre = 1e-10_dp
         !! Relative error tolerance required of the final solution at `lambda=1`.
      real(dp) :: ansae = 1e-10_dp
         !! Absolute error tolerance required of the final solution at `lambda=1`.
      real(dp) :: cursw = 10.0_dp
         !! Curvature switch. If the curvature of any component of `y` exceeds this value, the
         !! tolerances for the normal-flow iteration are switched from `arcre` and `arcae` to
         !! the finer tolerances `ansre` and `ansae`, respectively.
      real(dp), dimension(8) :: sspar = zero
         !! Step-size control parameters:
         !! `(lideal, rideal, dideal, hmin, hmax, bmin, bmax, p)`.
         !! Elements that are nonpositive on input are replaced by default values.
      integer :: max_steps = 1000
         !! Maximum number of steps allowed in the main loop.
      integer :: lunit = 0
         !! Logical I/O unit for intermediate output.
         !! * `0` : No output is printed (default).
         !! * `6` : Standard output (mapped internally to `output_unit`)
         !! * otherwise : Output to the specified I/O unit.
      logical :: polsys = .false.
         !! Optional flag used only by the polynomial-system driver [[POLSYS1H]].
   end type

   type :: nf_solver
   !! Container to hold the state variables and user-supplied callbacks for [[fixpnf]].
      integer :: problem_type = -100
      logical :: initialized = .false.
      type(nf_state) :: state
      type(nf_config) :: config
      type(hompack_f_callbacks) :: callbacks
   contains
      procedure :: initialize => nf_solver_initialize
      procedure :: solve => nf_solver_solve
      procedure :: restart => nf_solver_restart
   end type

   type :: nf_result
   !! Container to hold the results of a call to [[fixpnf]].
      real(dp) :: arc_length = zero
         !! Arc length of the zero curve until the solution was found.
      real(dp) :: lambda = zero
         !! Final value of the homotopy parameter \(\lambda\) at the solution.
      real(dp), allocatable :: x(:)
         !! Final value of the independent variable \(x\) at the solution.
      integer :: nfe = 0
         !! Number of homotopy/Jacobian evaluations.
      type(hompack_status) :: status
         !! Status.
   end type

contains

   type(hompack_status) function nf_solver_initialize( &
      self, problem_type, callbacks, n, dima, config) result(status)
   !! This function initializes a [[nf_solver]] object for a specified problem type,
   !! user-supplied callbacks, and configuration parameters.
   !! It validates the inputs, allocates necessary state variables (deallocation included),
   !! and binds the callbacks to the solver object.

      class(nf_solver), intent(inout) :: self
         !! Solver object.
      integer, intent(in) :: problem_type
         !! Problem type. Must be one of the following:
         !! * `fix_point` : solve \( F(x) = x \).
         !! * `zero_find` : solve \( F(x) = 0 \).
         !! * `curve_track` : track a zero curve of \( \rho(a,\lambda,x) = 0 \).
      type(hompack_f_callbacks), intent(in) :: callbacks
         !! User-supplied function and Jacobian evaluation subroutines. Required callbacks
         !! depend on `problem_type` as follows:
         !! * For `problem_type = fix_point` or `zero_find`, callbacks `f` and `fjac` must be
         !!   provided.
         !! * For `problem_type = curve_track`, callbacks `rho` and `rhojac` must be provided.
      integer, intent(in) :: n
         !! Problem dimension, i.e., the dimension of the independent variable \(x\).
      integer, intent(in), optional :: dima
         !! Dimension of the parameter vector \(a\) for curve-tracking problems. Required if
         !! `problem_type` is `curve_track`.
      type(nf_config), intent(in), optional :: config
         !! Configuration parameters. If not provided, default values are used.

      integer :: dima_, ierr

      ! Validate 'problem_type'
      if (problem_type /= problem_fix_point .and. problem_type /= problem_zero_find .and. &
          problem_type /= problem_curve_track) then
         status%message = "Illegal input: `prob_type` must be `fix_point`, `zero_find`, or `curve_track`."
         status%code = fixpnf_input_error
         return
      end if

      ! Validate 'n'
      if (n <= 0) then
         status%message = "Illegal input: `n` must be greater or equal than one."
         status%code = fixpnf_input_error
         return
      end if

      ! Validate callbacks and parameter vector 'a'
      if (problem_type == problem_fix_point .or. problem_type == problem_zero_find) then
         if (.not. associated(callbacks%f)) then
            status%message = "Illegal input: callback `f` must be provided for fixed-point and zero-finding problems."
            status%code = fixpnf_input_error
            return
         end if
         if (.not. associated(callbacks%fjac)) then
            status%message = "Illegal input: callback `fjac` must be provided for fixed-point and zero-finding problems."
            status%code = fixpnf_input_error
            return
         end if
         dima_ = n
      else if (problem_type == problem_curve_track) then
         if (.not. associated(callbacks%rho)) then
            status%message = "Illegal input: callback `rho` must be provided for curve-tracking problems."
            status%code = fixpnf_input_error
            return
         end if
         if (.not. associated(callbacks%rhojac)) then
            status%message = "Illegal input: callback `rhojac` must be provided for curve-tracking problems."
            status%code = fixpnf_input_error
            return
         end if
         if (.not. present(dima)) then
            status%message = "Illegal input: parameter vector dimension `dima` must be provided for curve-tracking problems."
            status%code = fixpnf_input_error
            return
         end if
         if (dima <= 0) then
            status%message = "Illegal input: parameter vector dimension `dima` must be greater than zero."
            status%code = fixpnf_input_error
            return
         end if
         dima_ = dima
      end if

      ! Validate configuration parameters
      if (present(config)) then
         if (config%max_steps < 1) then
            status%message = "Illegal input: `config.max_steps` must be greater than zero."
            status%code = fixpnf_input_error
            return
         end if
         if (config%ansre <= zero) then
            status%message = "Illegal input: `config.ansre` must be greater than zero."
            status%code = fixpnf_input_error
            return
         end if
         if (config%ansae < zero) then
            status%message = "Illegal input: `config.ansae` must be greater or equal than zero."
            status%code = fixpnf_input_error
            return
         end if
      end if

      ! Allocate state variables
      call self%state%alloc(n, dima_, stat=ierr)
      if (ierr /= 0) then
         status%message = "Memory allocation failure during initialization."
         status%code = fixpnf_memory_error
         return
      end if

      ! Update object fields now that initialization is successful
      self%problem_type = problem_type
      self%callbacks = callbacks ! copy is intentional
      if (present(config)) self%config = config
      self%initialized = .true.

      ! Finalize
      status%message = "Initialization successful."
      status%code = fixpnf_success

   end function nf_solver_initialize

   type(nf_result) function nf_solver_solve(self, x0, a) result(result)

      class(nf_solver), intent(inout) :: self
         !! Solver object.
      real(dp), intent(in) :: x0(:)
         !! Initial point for the solver. For fixed-point and zero-finding problems, this is
         !! the initial guess \(a\). For curve-tracking problems, this is the initial solution
         !! \(x_0\) at \(\lambda=0\).
      real(dp), intent(in), optional :: a(:)
         !! Parameter vector \(a\) for curve-tracking problems. Required if the solver was
         !! initialized with `problem_type = curve_track`.

      integer :: iflag, n
      real(dp), allocatable :: a_(:)

      n = self%state%n

      associate (message => result%status%message, code => result%status%code)

         ! Validate that the solver has been initialized
         if (.not. self%initialized) then
            message = "Solver has not been initialized. Call `initialize` before `solve`."
            code = fixpnf_input_error
            return
         end if

         ! Validate the initial point 'x0'
         if (size(x0) /= n) then
            message = "Illegal input: length of `x0` must be equal to `n` specified during initialization."
            code = fixpnf_input_error
            return
         end if

         ! Validate the parameter vector 'a' for curve-tracking problems
         if (self%problem_type == problem_curve_track) then
            if (.not. present(a)) then
               message = "Illegal input: parameter vector `a` must be provided for curve-tracking problems."
               code = fixpnf_input_error
               return
            end if
            if (size(a) /= size(self%state%a)) then
               message = "Illegal input: length of parameter vector `a` must match the dimension specified during initialization."
               code = fixpnf_input_error
               return
            end if
            a_ = a
         else
            a_ = [zero]
         end if

      end associate

      ! Determine 'iflag'
      if (self%problem_type == problem_fix_point) then
         iflag = 0
      else if (self%problem_type == problem_zero_find) then
         iflag = -1
      else if (self%problem_type == problem_curve_track) then
         iflag = -2
      end if

      ! Call the solver
      call fixpnf(self%state, self%callbacks, self%state%n, x0, iflag, self%config, a_)

      ! Populate the result object
      result%arc_length = self%state%s
      result%lambda = self%state%y(1)
      result%x = self%state%y(2:n + 1)
      result%nfe = self%state%nfe
      result%status%code = iflag

   end function nf_solver_solve

   type(nf_result) function nf_solver_restart(self) result(result)

      class(nf_solver), intent(inout) :: self
         !! Solver object.

      integer :: iflag, n

      ! n = self%state%n

      ! associate (message => result%status%message, code => result%status%code)

      !    ! Validate that the solver has been initialized
      !    if (.not. self%initialized) then
      !       message = "Solver has not been initialized. Call `initialize` before `solve`."
      !       code = fixpnf_input_error
      !       return
      !    end if

      !    ! Validate the initial point 'x0'
      !    if (size(x0) /= n) then
      !       message = "Illegal input: length of `x0` must be equal to `n` specified during initialization."
      !       code = fixpnf_input_error
      !       return
      !    end if

      !    ! Validate the parameter vector 'a' for curve-tracking problems
      !    if (self%problem_type == curve_track) then
      !       if (.not. present(a)) then
      !          message = "Illegal input: parameter vector `a` must be provided for curve-tracking problems."
      !          code = fixpnf_input_error
      !          return
      !       end if
      !       if (size(a) /= size(self%state%a)) then
      !          message = "Illegal input: length of parameter vector `a` must match the dimension specified during initialization."
      !          code = fixpnf_input_error
      !          return
      !       end if
      !       a_ = a
      !    else
      !       a_ = zero
      !    end if

      ! end associate

      ! ! Determine 'iflag'
      ! if (self%problem_type == fix_point) then
      !    iflag = 0
      ! else if (self%problem_type == zero_find) then
      !    iflag = -1
      ! else if (self%problem_type == curve_track) then
      !    iflag = -2
      ! end if

      ! ! Call the solver
      ! call fixpnf(self%state, self%callbacks, self%state%n, x0, iflag, self%config, a_)

      ! ! Populate the result object
      ! result%arc_length = self%state%s
      ! result%lambda = self%state%y(1)
      ! result%x = self%state%y(2:n + 1)
      ! result%nfe = self%state%nfe
      ! result%status%code = iflag

   end function nf_solver_restart

   subroutine fixpnf(state, callbacks, n, x0, iflag, config, a)
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

      implicit none

      type(nf_state), intent(inout) :: state
         !! Preallocated state for [[fixpnf]]. Set to start conditions on the first call, and
         !! updated on subsequent calls when `iflag=2` or `iflag=3`.
      type(hompack_f_callbacks), intent(in) :: callbacks
         !! User-supplied function and Jacobian evaluation subroutines.
      integer, intent(in) :: n
         !! Problem dimension.
      real(dp), intent(in) :: x0(:)
         !! Initial point for the solver. For fixed-point and zero-finding problems, this is
         !! identical to the initial guess \(a\). For curve-tracking problems, this is the
         !! initial solution at \(\lambda=0\).
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
         !! * `5` : tracking algorithm lost the zero curve and is not making progress; the
         !!         tolerances `arc?e` and `ans?e` were too lenient. Retry with smaller
         !!         tolerances.
         !! * `6` : normal flow Newton iteration failed to converge; tolerances `ans?e` may be
         !!         too stringent.
         !! * `7` : illegal input parameters.
         !! * `8` : memory allocation failure.
      type(nf_config), intent(inout) :: config
         !! Configuration parameters for the solver. Optional settings will be set to default
         !! values. Tolerances may be updated (increased) on output if `iflag=2`.
      real(dp), intent(in), optional :: a(:)
         !! Parameter vector `a` for curve-tracking problems. Ignored for fixed-point and
         !! zero-finding problems.

      integer :: iter, np1, dima

      np1 = n + 1

      ! INPUT CHECKS (intentionally light; extensive in the callers)

      ! Validate problem dimension
      if (n < 1 .or. size(x0) /= n .or. state%n /= n) then
         iflag = fixpnf_input_error
         return
      end if

      ! Validate 'iflag' and launch job
      if (iflag == 0 .or. iflag == -1) then
         dima = n
         goto 20
      else if (iflag == -2) then
         if (.not. present(a)) then
            iflag = fixpnf_input_error
            return
         end if
         dima = size(a)
         goto 20
      else if (iflag == 2) then
         goto 120
      else if (iflag == 3) then
         goto 90
      else
         iflag = fixpnf_input_error
         return
      end if

      ! FIRST CALL

      ! Set state variables to starting values
20    state%abserr = zero
      state%relerr = zero
      state%curtol = zero
      state%h = 0.1_dp
      state%hold = zero
      state%s = zero
      state%nfe = 0
      state%start = .true.
      state%crash = .false.

      state%iflag = iflag

      state%y = zero
      state%yold = zero
      state%yp = zero
      state%ypold = zero

      state%y(2:np1) = x0
      state%yp(1) = one
      state%ypold(1) = one

      ! Load 'a'
      if (state%iflag == 0 .or. state%iflag == -1) then
         state%a = x0
      else
         state%a = a
      end if

      ! Default arc tolerances
      if (config%arcre <= zero) config%arcre = sqrt(config%ansre)/2
      if (config%arcae <= zero) config%arcae = sqrt(config%ansae)/2

      ! Set optimal step size estimation parameters
      ! 'z[k]' denote the Newton iterates along the flow normal to the Davidenko flow
      ! and 'y' their limit

      associate (sspar => config%sspar)
         ! Ideal contraction factor: ||z[2] - z[1]|| / ||z[1] - z[0]||
         if (sspar(1) <= zero) sspar(1) = 0.5_dp
         ! Ideal residual factor: ||rho(A, z[1])|| / ||rho(a, z[0])||
         if (sspar(2) <= zero) sspar(2) = 0.01_dp
         ! Ideal distance factor: ||z[1] - y|| / ||z[0] - y||
         if (sspar(3) <= zero) sspar(3) = 0.5_dp
         ! Minimum step size 'hmin'
         if (sspar(4) <= zero) sspar(4) = (sqrt(real(n + 1, dp)) + 4.0_dp)*eps64
         ! Maximum step size 'hmax'
         if (sspar(5) <= zero) sspar(5) = one
         ! Minimum step size reduction factor 'bmin'
         if (sspar(6) <= zero) sspar(6) = 0.1_dp
         ! Maximum step size expansion factor 'bmax'
         if (sspar(7) <= zero) sspar(7) = 3.0_dp
         ! Assumed operating order 'p'
         if (sspar(8) <= zero) sspar(8) = 2.0_dp
      end associate

      ! COMMON PART FOR FIRST CALL AND RESTARTS

      ! Set default iteration limit
90    state%limit = config%max_steps

      ! Main loop
120   do iter = 1, state%limit

         ! Tracking algorithm lost the zero curve
         if (state%y(1) < zero) then
            iflag = 5
            return
         end if

         ! Set different error tolerance if the trajectory y(s) has any high curvature
         ! components
         state%curtol = config%cursw*state%hold
         state%relerr = config%arcre
         state%abserr = config%arcae
         if (any(abs(state%yp - state%ypold) > state%curtol)) then
            state%relerr = config%ansre
            state%abserr = config%ansae
         end if

         ! Take a step along the curve
         call stepnf(state, callbacks, config%sspar)

         ! Print latest point on curve if requested
         if (config%lunit /= 0) then
            write (config%lunit, 217) iter, state%nfe, state%s, state%y(1), state%y(2:np1)
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
            if (config%arcre < state%relerr) config%arcre = state%relerr
            if (config%ansre < state%relerr) config%ansre = state%relerr
            if (config%arcae < state%abserr) config%arcae = state%abserr
            if (config%ansae < state%abserr) config%ansae = state%abserr
            ! Change limit on number of iterations
            state%limit = state%limit - iter
            return
         end if

         ! Use Hermite cubic interpolation and Newton iteration to get the answer at lambda=1
         if (state%y(1) >= one) then

            associate (ws => state%workspace)

               ! Save 'yold' for arc length calculation later
               ws%z0 = state%yold
               call rootnf(state, callbacks, config%ansre, config%ansae)
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
         if (config%polsys) then
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

      use hompack_core, only: root, root_state, qofs
      implicit none

      type(nf_state), intent(inout) :: state
         !! State variables for [[fixpnf]].
      type(hompack_f_callbacks) :: callbacks
         !! User-supplied function and Jacobian evaluation subroutines.
      real(dp), intent(in) :: relerr
         !! Relative convergence tolerance.
         !! Iteration is considered converged when `|y(1) - 1| <= relerr + abserr` and the
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

   subroutine stepnf(state, callbacks, sspar)
   !! This subroutine takes one step along the zero curve of the homotopy map using a
   !! predictor-corrector algorithm. The predictor uses a Hermite cubic interpolant, and
   !! the corrector returns to the zero curve along the flow normal to the Davidenko flow.
   !! [[stepnf]] also estimates a step size `h` for the next step along the zero curve.
   !! Normally, [[stepnf]] is used indirectly through [[fixpnf]], and should be called
   !! directly only if it is necessary to modify the stepping algorithm's parameters.

      use hompack_core, only: qofs
      implicit none

      type(nf_state), intent(inout) :: state
         !! State variables for [[fixpnf]].
      type(hompack_f_callbacks), intent(in) :: callbacks
         !! User-supplied function and Jacobian evaluation subroutines.
      real(dp), intent(in) :: sspar(:)
         !! Step-size control parameters:
         !! `(lideal, rideal, dideal, hmin, hmax, bmin, bmax, p)`.

      real(dp), parameter :: twou = 2*eps64, fouru = 4*eps64
      real(dp) :: dcalc, hfail, ht, lcalc, rcalc, rholen, temp
      integer :: itnum, j, judy, np1
      logical :: fail

      ! Newton iterations allowed before reducing the step size 'h'
      integer, parameter:: max_newton = 4

      np1 = state%n + 1
      state%crash = .true.

      associate (ws => state%workspace)

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

               do judy = 1, max_newton

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

               ! No convergence in 'max_newton' iterations. Reduce h and try again.
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
            do judy = 1, max_newton

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

            ! No convergence in 'max_newton' iterations. Record failure at calculated 'h'
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
         else if (itnum == max_newton) then
            ! If convergence required the maximum 'max_newton' iterations, don't increase 'h'.
            state%h = min(state%h, state%hold)
         end if

         ! If convergence did not occur in 'max_newton' iterations for a particular
         ! 'h = hfail', don't choose the new step size larger than 'hfail'.
         if (fail) state%h = min(state%h, hfail)

      end associate

   end subroutine stepnf

   subroutine tangnf( &
      callbacks, rholen, n, y, yp, ypold, a, nfe, iflag, qr, alpha, tz, pivot)
   !! This subroutine builds the Jacobian matrix of the homotopy map, computes a QR
   !! decomposition of that matrix, and then calculates the (unit) tangent vector and the
   !! Newton step.

      use lapack_interfaces, only: dgeqpf, dormqr
      implicit none

      type(hompack_f_callbacks) :: callbacks
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

   pure subroutine alloc_workspace(self, n, stat)
   !! Initializes [[nf_workspace]], i.e., (re)allocates all allocatable arrays and sets them
   !! to zero.

      class(nf_workspace), intent(inout) :: self
         !! Workspace object.
      integer, intent(in) :: n
         !! Problem dimension.
      integer, intent(out), optional :: stat
         !! Error status of the allocation.

      integer :: ierr(8)

      if (n <= 0) then
         error stop "Error: 'n' must be positive in hompack_nf_state%alloc()."
      end if

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
            error stop "Error: Allocation failed in hompack_nf_workspace%alloc()."
         end if
      end if

   end subroutine alloc_workspace

   pure subroutine alloc_state(self, n, dima, stat)
   !! Initializes [[nf_state]], i.e., (re)allocates all allocatable arrays and sets them to
   !! zero.

      class(nf_state), intent(inout) :: self
         !! State object.
      integer, intent(in)  :: n
         !! Problem dimension.
      integer, intent(in)  :: dima
         !! Dimension of the parameter vector `a`.
      integer, intent(out), optional   :: stat
         !! Error status of the allocation.

      integer :: ierr(6)

      if (n <= 0) then
         error stop "Error: 'n' must be positive in hompack_nf_state%alloc()."
      end if

      if (dima <= 0) then
         error stop "Error: 'dima' must be positive in hompack_nf_state%alloc()."
      end if

      if (present(stat)) stat = 0

      ! Deallocate any previously allocated arrays
      self%n = 0
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
      call self%workspace%alloc(n, stat=ierr(6))

      if (any(ierr /= 0)) then
         if (present(stat)) then
            stat = ierr(findloc(ierr /= 0, .true., dim=1))
         else
            error stop "Error: Allocation failed in nf_state%alloc()."
         end if
      end if

      ! Set problem dimension now that allocation is successful
      self%n = n

   end subroutine alloc_state

end module hompack_nf
