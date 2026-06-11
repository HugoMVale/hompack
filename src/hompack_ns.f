      MODULE HOMPACK_NS
C
      USE HOMPACK_KINDS, ONLY: DP
      IMPLICIT NONE
C
      CONTAINS
C
      SUBROUTINE FIXPNS(N,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,A,
     &   NFE,ARCLEN,MODE,LENQR,SSPAR)
C
C Subroutine  FIXPNS  finds a fixed point or zero of the
C N-dimensional vector function F(X), or tracks a zero curve
C of a general homotopy map RHO(A,X,LAMBDA).  For the fixed 
C point problem F(X) is assumed to be a C2 map of some ball 
C into itself.  The equation  X = F(X)  is solved by
C following the zero curve of the homotopy map
C
C  LAMBDA*(X - F(X)) + (1 - LAMBDA)*(X - A)  ,
C
C starting from LAMBDA = 0, X = A.  The curve is parameterized
C by arc length S, and is followed by solving the ordinary
C differential equation  D(HOMOTOPY MAP)/DS = 0  for
C Y(S) = (X(S), LAMBDA(S)) using a Hermite cubic predictor and a
C corrector which returns to the zero curve along the flow normal
C to the Davidenko flow (which consists of the integral curves of
C D(HOMOTOPY MAP)/DS ).
C
C For the zero finding problem F(X) is assumed to be a C2 map
C such that for some R > 0,  X*F(X) >= 0  whenever NORM(X) = R.
C The equation  F(X) = 0  is solved by following the zero curve
C of the homotopy map
C
C   LAMBDA*F(X) + (1 - LAMBDA)*(X - A)
C
C emanating from LAMBDA = 0, X = A.
C
C  A  must be an interior point of the above mentioned balls.
C
C For the curve tracking problem RHO(A,X,LAMBDA) is assumed to
C be a C2 map from E**M X E**N X [0,1) into E**N, which for
C almost all parameter vectors A in some nonempty open subset
C of E**M satisfies
C
C  rank [D RHO(A,X,LAMBDA)/DX , D RHO(A,X,LAMBDA)/D LAMBDA] = N
C
C for all points (X,LAMBDA) such that RHO(A,X,LAMBDA)=0.  It is
C further assumed that
C
C           rank [ D RHO(A,X0,0)/DX ] = N  .
C
C With A fixed, the zero curve of RHO(A,X,LAMBDA) emanating
C from  LAMBDA = 0, X = X0  is tracked until  LAMBDA = 1  by
C solving the ordinary differential equation
C D RHO(A,X(S),LAMBDA(S))/DS = 0  for  Y(S) = (X(S), LAMBDA(S)),
C where S is arc length along the zero curve.  Also the homotopy
C map RHO(A,X,LAMBDA) is assumed to be constructed such that
C
C              D LAMBDA(0)/DS > 0  .
C
C
C For the fixed point and zero finding problems, the user
C must supply a subroutine  F(X,V)  which evaluates F(X) at X
C and returns the vector F(X) in V, and a subroutine
C  FJACS(X)  which evaluates, if
C MODE = 1,
C   the (symmetric) Jacobian matrix of F(X) at X, and returns the
C   symmetric Jacobian matrix in packed skyline storage format in
C   QR, or if
C MODE = 2,
C   returns the (nonsymmetric) Jacobian matrix in sparse row format
C   in QR.  The MODE 1 format is defined by QR, LENQR, ROWPOS; the
C   MODE 2 format is defined by QR, LENQR, ROWPOS, COLPOS.
C
C For the curve tracking problem, the user must supply a subroutine
C  RHO(A,LAMBDA,X,V)  which evaluates the homotopy map RHO 
C at (A,X,LAMBDA) and returns the vector RHO(A,X,LAMBDA) in V, and 
C a subroutine  RHOJS(A,LAMBDA,X)  which, if
C MODE = 1,
C   returns in QR the symmetric N X N Jacobian matrix [D RHO/DX] 
C   evaluated at (A,X,LAMBDA) and stored in packed skyline format, 
C   and returns in PP the vector -(D RHO/D LAMBDA) evaluated at 
C   (A,X,LAMBDA).  This data structure is described by QR, LENQR,
C   ROWPOS, PP.  *** Note the minus sign in the definition of PP. ***  If
C MODE = 2,
C   the (nonsymmetric) N X (N+1) Jacobian matrix [D RHO/DX, D RHO/DLAMBDA]
C   evaluated at (A,X,LAMBDA) is returned in a data structure described
C   by QR, LENQR, ROWPOS, COLPOS.
C
C Whichever of the routines  F,  FJACS,  RHO,  RHOJS  are required
C should be supplied as external subroutines, conforming with the
C interfaces in the module  HOMOTOPY.
C
C
C Functions and subroutines directly or indirectly called by FIXPNS:
C F (or  RHO ), FJACS (or  RHOJS ), GMFADS , GMRES , GMRILUDS ,
C ILUFDS , ILUSOLVDS , MULTDS , MULT2DS , PCGDS , ROOT , ROOTNS ,
C SOLVDS , STEPNS , TANGNS , and the BLAS functions  DDOT , DLAIC1 ,
C DLAMCH , DNRM2 .  The module  REAL_PRECISION  specifies 64-bit
C real arithmetic, which the user may want to change.
C 
C
C ON INPUT:
C
C N  is the dimension of X, F(X), and RHO(A,X,LAMBDA).
C
C Y  is an array of length  N + 1.  (Y(1),...,Y(N)) = A  is the
C    starting point for the zero curve for the fixed point and 
C    zero finding problems.  (Y(1),...,Y(N)) = X0  for the curve
C    tracking problem.
C
C IFLAG  can be -2, -1, 0, 2, or 3.  IFLAG  should be 0 on the 
C    first call to  FIXPNS  for the problem  X=F(X), -1 for the
C    problem  F(X)=0, and -2 for the problem  RHO(A,X,LAMBDA)=0.
C    In certain situations  IFLAG  is set to 2 or 3 by  FIXPNS,
C    and  FIXPNS  can be called again without changing  IFLAG.
C
C ARCRE , ARCAE  are the relative and absolute errors, respectively,
C    allowed the normal flow iteration along the zero curve.  If
C    ARC?E .LE. 0.0  on input it is reset to  .5*SQRT(ANS?E) .
C    Normally  ARC?E should be considerably larger than  ANS?E .
C
C ANSRE , ANSAE  are the relative and absolute error values used for
C    the answer at LAMBDA = 1.  The accepted answer  Y = (X, LAMBDA)
C    satisfies
C
C       |Y(NP1) - 1|  .LE.  ANSRE + ANSAE           .AND.
C
C       ||Z||  .LE.  ANSRE*||X|| + ANSAE          where
C
C    (Z,.) is the Newton step to Y.
C
C TRACE  is an integer specifying the logical I/O unit for
C    intermediate output.  If  TRACE .GT. 0  the points computed on
C    the zero curve are written to I/O unit  TRACE .
C
C A(:)  contains the parameter vector  A.  For the fixed point
C    and zero finding problems, A  need not be initialized by the
C    user, and is assumed to have length  N.  For the curve
C    tracking problem, A  must be initialized by the user.
C
C MODE = 1 if the Jacobian matrix is symmetric and stored in a packed
C          skyline format;
C      = 2 if the Jacobian matrix is stored in a sparse row format.
C
C LENQR  is the number of nonzero entries in the sparse Jacobian
C    matrices, used to determine the sparse matrix data structures.
C
C SSPAR(1:8) = (LIDEAL, RIDEAL, DIDEAL, HMIN, HMAX, BMIN, BMAX, P)  is
C    a vector of parameters used for the optimal step size estimation.
C    If  SSPAR(J) .LE. 0.0  on input, it is reset to a default value
C    by  FIXPNS .  Otherwise the input value of  SSPAR(J)  is used.
C    See the comments below and in  STEPNS  for more information about
C    these constants.
C
C
C ON OUTPUT:
C
C N , TRACE , A  are unchanged.
C
C (Y(1),...,Y(N)) = X, Y(NP1) = LAMBDA, and Y is an approximate
C    zero of the homotopy map.  Normally LAMBDA = 1 and X is a
C    fixed point(zero) of F(X).  In abnormal situations LAMBDA
C    may only be near 1 and X is near a fixed point(zero).
C
C IFLAG =
C  -2   causes  FIXPNS  to initialize everything for the problem
C       RHO(A,X,LAMBDA) = 0 (use on first call).
C
C  -1   causes  FIXPNS  to initialize everything for the problem
C       F(X) = 0 (use on first call).
C
C   0   causes  FIXPNS  to initialize everything for the problem
C       X = F(X) (use on first call).
C
C   1   Normal return.
C
C   2   Specified error tolerance cannot be met.  Some or all of
C       ARCRE , ARCAE , ANSRE , ANSAE  have been increased to 
C       suitable values.  To continue, just call  FIXPNS  again 
C       without changing any parameters.
C
C   3   STEPNS  has been called 1000 times.  To continue, call
C       FIXPNS  again without changing any parameters.
C
C   4   The preconditioned conjugate gradient iteration failed to
C       converge, or the Jacobian matrix does not have full rank
C       or has a zero on the diagonal.  The algorithm has failed
C       (the zero curve of the homotopy map cannot be followed any
C       further).
C
C   5   The tracking algorithm has lost the zero curve of the
C       homotopy map and is not making progress.  The error tolerances
C       ARC?E  and  ANS?E  were too lenient.  The problem should be
C       restarted by calling  FIXPNS  with smaller error tolerances
C       and  IFLAG = 0 (-1, -2).
C
C   6   The normal flow Newton iteration in  STEPNS  or  ROOTNS
c       failed to converge.  The error tolerances  ANS?E  may be too
C       stringent.
C
C   7   Illegal input parameters, a fatal error.
C
C ARCRE , ARCAE , ANSRE , ANSAE  are unchanged after a normal return 
C    (IFLAG = 1).  They are increased to appropriate values on the 
C    return  IFLAG = 2 .
C
C NFE  is the number of function evaluations (= number of
C    Jacobian evaluations).
C
C ARCLEN  is the length of the path followed.
C
C
C Allocatable and automatic work arrays:
C
C YP(1:N+1)  is a work array containing the tangent vector to 
C    the zero curve at the current point  Y .
C
C YOLD(1:N+1)  is a work array containing the previous point found
C    on the zero curve.
C
C YPOLD(1:N+1)  is a work array containing the tangent vector to 
C    the zero curve at  YOLD .
C
C QR(1:LENQR), PP(1:N), ROWPOS(1:N+2), COLPOS(1:LENQR) are all work
C    arrays used to define the sparse Jacobian matrices, allocated
C    here, and distributed via the module  HOMOTOPY.
C
C
      USE HOMPACK_GLOBAL_LEGACY, ONLY: QR => QRSPARSE, ROWPOS, COLPOS,
     & PP, PAR, IPAR
      USE HOMPACK_CORE_LEGACY, ONLY: ROOTNS
      USE BLAS_INTERFACES, ONLY: DNRM2
      IMPLICIT NONE
C
      INTEGER, INTENT(IN)::LENQR,MODE,N,TRACE
      REAL (dp), DIMENSION(:), INTENT(IN OUT)::A,Y
      INTEGER, INTENT(IN OUT)::IFLAG
      REAL (dp), INTENT(IN OUT)::ANSAE,ANSRE,ARCAE,ARCRE,SSPAR(8)
      INTEGER, INTENT(OUT)::NFE
      REAL (dp), INTENT(OUT)::ARCLEN
C
C *****  LOCAL VARIABLES.  *****
C
      REAL (dp), SAVE:: ABSERR,CURTOL,H,HOLD,RELERR,S
      INTEGER, SAVE:: IFLAGC,ITER,JW,LIMIT,NC,NFEC,NP1
      LOGICAL, SAVE:: CRASH,START
C ***** WORK ARRAYS. *****
      REAL (dp), ALLOCATABLE, DIMENSION(:), SAVE:: YP,YOLD,YPOLD
      REAL (dp):: TZ(N+1),W(N+1),WP(N+1),Z0(N+1),Z1(N+1)
C
C LIMITD  IS AN UPPER BOUND ON THE NUMBER OF STEPS.  IT MAY BE
C CHANGED BY CHANGING THE FOLLOWING PARAMETER STATEMENT:
      INTEGER, PARAMETER:: LIMITD=1000
C
C SWITCH FROM THE TOLERANCE  ARC?E  TO THE (FINER) TOLERANCE  ANS?E  IF
C THE CURVATURE OF ANY COMPONENT OF  Y  EXCEEDS  CURSW.
      REAL (dp), PARAMETER:: CURSW=10.0
C
C ***** END OF SPECIFICATION INFORMATION. *****
C
C :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :
      IF (N .LE. 0  .OR.  ANSRE .LE. 0.0  .OR.  ANSAE .LT. 0.0
     &  .OR.  N+1 .NE. SIZE(Y) .OR.
     &  ((IFLAG .EQ. -1  .OR.  IFLAG .EQ. 0) .AND.  N .NE. SIZE(A))
     &  .OR.  MODE .LE. 0  .OR.  MODE .GE. 3)
     &                                                     IFLAG=7
      IF (IFLAG .GE. -2  .AND.  IFLAG .LE. 0) GO TO 20
      IF (IFLAG .EQ. 2) GO TO 120
      IF (IFLAG .EQ. 3) GO TO 90
C ONLY VALID INPUT FOR  IFLAG  IS -2, -1, 0, 2, 3.
      IFLAG=7
      RETURN
C
C *****  INITIALIZATION BLOCK.  *****
C
20    ARCLEN=0.0
      IF (ARCRE .LE. 0.0) ARCRE=.5*SQRT(ANSRE)
      IF (ARCAE .LE. 0.0) ARCAE=.5*SQRT(ANSAE)
      NC=N
      NFEC=0
      IFLAGC=IFLAG
      NP1=N+1
C ALLOCATE SAVED WORK ARRAYS.
      IF (ALLOCATED(YOLD)) DEALLOCATE(YOLD)
      IF (ALLOCATED(YP)) DEALLOCATE(YP)
      IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
      ALLOCATE(YOLD(NP1),YP(NP1),YPOLD(NP1))
C SET INITIAL CONDITIONS FOR FIRST CALL TO  STEPNS .
      START=.TRUE.
      CRASH=.FALSE.
      HOLD=1.0
      H=.1
      S=0.0
      YPOLD(NP1)=1.0
      YP(NP1)=1.0
      Y(NP1)=0.0
      YPOLD(1:N)=0.0
      YP(1:N)=0.0
C SET OPTIMAL STEP SIZE ESTIMATION PARAMETERS.
C LET Z[K] DENOTE THE NEWTON ITERATES ALONG THE FLOW NORMAL TO THE
C DAVIDENKO FLOW AND Y THEIR LIMIT.
C IDEAL CONTRACTION FACTOR:  ||Z[2] - Z[1]|| / ||Z[1] - Z[0]||
      IF (SSPAR(1) .LE. 0.0) SSPAR(1)= .5_dp
C IDEAL RESIDUAL FACTOR:  ||RHO(A, Z[1])|| / ||RHO(A, Z[0])||
      IF (SSPAR(2) .LE. 0.0) SSPAR(2)= .01_dp
C IDEAL DISTANCE FACTOR:  ||Z[1] - Y|| / ||Z[0] - Y||
      IF (SSPAR(3) .LE. 0.0) SSPAR(3)= .5_dp
C MINIMUM STEP SIZE  HMIN .
      IF (SSPAR(4) .LE. 0.0) SSPAR(4)=(SQRT(N+1.0)+4.0)*EPSILON(1.0_dp)
C MAXIMUM STEP SIZE  HMAX .
      IF (SSPAR(5) .LE. 0.0) SSPAR(5)= 1.0
C MINIMUM STEP SIZE REDUCTION FACTOR  BMIN .
      IF (SSPAR(6) .LE. 0.0) SSPAR(6)= .1_dp
C MAXIMUM STEP SIZE EXPANSION FACTOR  BMAX .
      IF (SSPAR(7) .LE. 0.0) SSPAR(7)= 3.0
C ASSUMED OPERATING ORDER  P .
      IF (SSPAR(8) .LE. 0.0) SSPAR(8)= 2.0
C
C LOAD  A  FOR THE FIXED POINT AND ZERO FINDING PROBLEMS.
      IF (IFLAGC .GE. -1) A(1:N) = Y(1:N)
90    LIMIT=LIMITD
C ALLOCATE ARRAYS FOR SPARSE JACOBIAN MATRIX DATA STRUCTURE.
120   SELECT CASE (MODE)
        CASE (1)
          IF (.NOT. ALLOCATED(QR)) ALLOCATE(QR(LENQR))
          IF (.NOT. ALLOCATED(ROWPOS)) ALLOCATE(ROWPOS(N+2))
          IF (.NOT. ALLOCATED(PP)) ALLOCATE(PP(N))
        CASE (2)
          IF (.NOT. ALLOCATED(QR)) ALLOCATE(QR(LENQR))
          IF (.NOT. ALLOCATED(ROWPOS)) ALLOCATE(ROWPOS(N+2))
          IF (.NOT. ALLOCATED(COLPOS)) ALLOCATE(COLPOS(LENQR))
          IF ((.NOT. ALLOCATED(PP)) .AND. (IFLAGC .GE. -1))
     &      ALLOCATE(PP(N))
      END SELECT
C
C *****  END OF INITIALIZATION BLOCK.  *****
C
      MAIN_LOOP: DO ITER=1,LIMIT  ! *****  MAIN LOOP.  *****
      IF (Y(NP1) .LT. 0.0) THEN
        ARCLEN=S
        IFLAG=5
        CALL CLEANUPALL
        RETURN
      ENDIF
C
C SET DIFFERENT ERROR TOLERANCE IF THE TRAJECTORY Y(S) HAS ANY HIGH 
C CURVATURE COMPONENTS.
      CURTOL=CURSW*HOLD
      RELERR=ARCRE
      ABSERR=ARCAE
        IF (ANY(ABS(YP-YPOLD) .GT. CURTOL)) THEN
          RELERR=ANSRE
          ABSERR=ANSAE
        ENDIF
C
C TAKE A STEP ALONG THE CURVE.
      CALL STEPNS(NC,NFEC,IFLAGC,START,CRASH,HOLD,H,RELERR,
     &     ABSERR,S,Y,YP,YOLD,YPOLD,A,MODE,LENQR,SSPAR,TZ,W,WP,Z0,Z1)
C PRINT LATEST POINT ON CURVE IF REQUESTED.
      IF (TRACE .GT. 0) THEN
        WRITE (TRACE,217) ITER,NFEC,S,Y(NP1),(Y(JW),JW=1,NC)
217     FORMAT(/' STEP',I5,3X,'NFE =',I5,3X,'ARC LENGTH =',F9.4,3X,
     &  'LAMBDA =',F7.4,5X,'X vector:'/(1X,6ES12.4))
      ENDIF
      NFE=NFEC
C CHECK IF THE STEP WAS SUCCESSFUL.
      IF (IFLAGC .GT. 0) THEN
        ARCLEN=S
        IFLAG=IFLAGC
        CALL CLEANUPALL
        RETURN
      ENDIF
      IF (CRASH) THEN
C RETURN CODE FOR ERROR TOLERANCE TOO SMALL.
        IFLAG=2
C CHANGE ERROR TOLERANCES.
        IF (ARCRE .LT. RELERR) ARCRE=RELERR
        IF (ANSRE .LT. RELERR) ANSRE=RELERR
        IF (ARCAE .LT. ABSERR) ARCAE=ABSERR
        IF (ANSAE .LT. ABSERR) ANSAE=ABSERR
C CHANGE LIMIT ON NUMBER OF ITERATIONS.
        LIMIT=LIMIT-ITER
        CALL CLEANUP
        RETURN
      ENDIF
C
      IF (Y(NP1) .GE. 1.0) THEN
C
C USE HERMITE CUBIC INTERPOLATION AND NEWTON ITERATION TO GET THE 
C ANSWER AT LAMBDA = 1.0 .
C
C SAVE  YOLD  FOR ARC LENGTH CALCULATION LATER.
        W=YOLD
C
        CALL ROOTNS(NC,NFEC,IFLAGC,ANSRE,ANSAE,Y,YP,YOLD,YPOLD,
     &              A,MODE,LENQR)
C
        NFE=NFEC
        IFLAG=1
C SET ERROR FLAG IF  ROOTNS  COULD NOT GET THE POINT ON THE ZERO
C CURVE AT  LAMBDA = 1.0  .
        IF (IFLAGC .GT. 0) IFLAG=IFLAGC
C CALCULATE FINAL ARC LENGTH.
        W = Y - W
        ARCLEN = S - HOLD + DNRM2(NP1,W,1)
        CALL CLEANUPALL
        RETURN
      ENDIF
C
      END DO MAIN_LOOP  !  *****  END OF MAIN LOOP.  *****
C
C LAMBDA HAS NOT REACHED 1 IN 1000 STEPS.
      IFLAG=3
      ARCLEN=S
      CALL CLEANUP
      RETURN
C
      CONTAINS
      SUBROUTINE CLEANUPALL
      IF (ALLOCATED(YOLD)) DEALLOCATE(YOLD)
      IF (ALLOCATED(YP)) DEALLOCATE(YP)
      IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
      CALL CLEANUP
      RETURN
      END SUBROUTINE CLEANUPALL
      SUBROUTINE CLEANUP
      IF (ALLOCATED(QR)) DEALLOCATE(QR)
      IF (ALLOCATED(ROWPOS)) DEALLOCATE(ROWPOS)
      IF (ALLOCATED(COLPOS)) DEALLOCATE(COLPOS)
      IF (ALLOCATED(PP)) DEALLOCATE(PP)
      IF (ALLOCATED(PAR)) DEALLOCATE(PAR)
      IF (ALLOCATED(IPAR)) DEALLOCATE(IPAR)
      RETURN
      END SUBROUTINE CLEANUP
      END SUBROUTINE FIXPNS
C
C
      SUBROUTINE STEPNS(N,NFE,IFLAG,START,CRASH,HOLD,H,RELERR,
     &   ABSERR,S,Y,YP,YOLD,YPOLD,A,MODE,LENQR,SSPAR,TZ,W,WP,Z0,Z1)
C
C  STEPNS  TAKES ONE STEP ALONG THE ZERO CURVE OF THE HOMOTOPY MAP
C USING A PREDICTOR-CORRECTOR ALGORITHM.  THE PREDICTOR USES A HERMITE
C CUBIC INTERPOLANT, AND THE CORRECTOR RETURNS TO THE ZERO CURVE ALONG
C THE FLOW NORMAL TO THE DAVIDENKO FLOW.  STEPNS  ALSO ESTIMATES A
C STEP SIZE H FOR THE NEXT STEP ALONG THE ZERO CURVE.  NORMALLY
C  STEPNS  IS USED INDIRECTLY THROUGH  FIXPNS , AND SHOULD BE CALLED
C DIRECTLY ONLY IF IT IS NECESSARY TO MODIFY THE STEPPING ALGORITHM'S
C PARAMETERS.  SEE ALSO THE REVERSE CALL ROUTINE  STEPNX .
C
C THE CALLING PROGRAM MUST INCLUDE THE FOLLOWING INTERFACE BLOCK:
C
C     INTERFACE
C       SUBROUTINE STEPNS(N,NFE,IFLAG,START,CRASH,HOLD,H,RELERR,
C    &    ABSERR,S,Y,YP,YOLD,YPOLD,A,MODE,LENQR,SSPAR,TZ,W,WP,Z0,Z1)
C       USE REAL_PRECISION
C       INTEGER, INTENT(IN):: LENQR,MODE,N
C       INTEGER, INTENT(IN OUT):: IFLAG,NFE
C       LOGICAL, INTENT(IN OUT):: CRASH,START
C       REAL(DP), INTENT(IN):: A(:),SSPAR(8)
C       REAL(DP), INTENT(IN OUT):: ABSERR,H,HOLD,RELERR,S,
C    &    Y(:),YOLD(:),YP(:),YPOLD(:)
C       REAL(DP), INTENT(OUT), DIMENSION(:):: TZ,W,WP,Z0,Z1
C       END SUBROUTINE STEPNS
C     END INTERFACE
C
C
C ON INPUT:
C
C N = DIMENSION OF X AND THE HOMOTOPY MAP.
C
C NFE = NUMBER OF JACOBIAN MATRIX EVALUATIONS.
C
C IFLAG = -2, -1, OR 0, INDICATING THE PROBLEM TYPE.
C
C START = .TRUE. ON FIRST CALL TO  STEPNS , .FALSE. OTHERWISE.
C
C HOLD = ||Y - YOLD||; SHOULD NOT BE MODIFIED BY THE USER.
C
C H = UPPER LIMIT ON LENGTH OF STEP THAT WILL BE ATTEMPTED.  H  MUST BE
C    SET TO A POSITIVE NUMBER ON THE FIRST CALL TO  STEPNS .
C    THEREAFTER  STEPNS  CALCULATES AN OPTIMAL VALUE FOR  H , AND  H
C    SHOULD NOT BE MODIFIED BY THE USER.
C
C RELERR, ABSERR = RELATIVE AND ABSOLUTE ERROR VALUES.  THE ITERATION IS
C    CONSIDERED TO HAVE CONVERGED WHEN A POINT W=(X,LAMBDA) IS FOUND 
C    SUCH THAT
C
C    ||Z|| <= RELERR*||W|| + ABSERR  ,          WHERE
C
C    Z IS THE NEWTON STEP TO W=(X,LAMBDA).
C
C S = (APPROXIMATE) ARC LENGTH ALONG THE HOMOTOPY ZERO CURVE UP TO
C    Y(S) = (X(S), LAMBDA(S)).
C
C Y(1:N+1) = PREVIOUS POINT (X(S), LAMBDA(S)) FOUND ON THE ZERO CURVE OF 
C    THE HOMOTOPY MAP.
C
C YP(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY MAP
C    AT  Y .
C
C YOLD(1:N+1) = A POINT BEFORE  Y  ON THE ZERO CURVE OF THE HOMOTOPY MAP.
C
C YPOLD(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY
C    MAP AT  YOLD .
C
C A(:) = PARAMETER VECTOR IN THE HOMOTOPY MAP.
C
C MODE = 1 IF THE JACOBIAN MATRIX IS SYMMETRIC AND STORED IN A PACKED
C          SKYLINE FORMAT;
C      = 2 IF THE JACOBIAN MATRIX IS STORED IN A SPARSE ROW FORMAT.
C
C LENQR  IS THE NUMBER OF NONZERO ENTRIES IN THE SPARSE JACOBIAN
C    MATRICES, USED TO DETERMINE THE SPARSE MATRIX DATA STRUCTURES.
C
C SSPAR(1:8) = (LIDEAL, RIDEAL, DIDEAL, HMIN, HMAX, BMIN, BMAX, P)  IS
C    A VECTOR OF PARAMETERS USED FOR THE OPTIMAL STEP SIZE ESTIMATION.
C
C TZ(1:N+1), W(1:N+1), WP(1:N+1), Z0(1:N+1), AND  Z1(1:N+1)  ARE WORK
C    ARRAYS USED FOR THE CALCULATION OF THE JACOBIAN MATRIX KERNEL, THE
C    NEWTON STEP, INTERPOLATION, AND THE ESTIMATION OF THE NEXT STEP
C    SIZE  H .
C
C
C ON OUTPUT:
C
C N , A , SSPAR  ARE UNCHANGED.
C
C NFE  HAS BEEN UPDATED.
C
C IFLAG  
C    = -2, -1, OR 0 (UNCHANGED) ON A NORMAL RETURN.
C
C    = 4 IF THE CONJUGATE GRADIENT ITERATION FAILED TO CONVERGE
C        (MOST LIKELY DUE TO A JACOBIAN MATRIX WITH RANK < N).  THE
C        ITERATION WAS NOT COMPLETED.
C
C    = 6 IF THE NEWTON ITERATION FAILED TO CONVERGE.  W  CONTAINS 
C        THE LAST NEWTON ITERATE.
C
C START = .FALSE. ON A NORMAL RETURN.
C
C CRASH 
C    = .FALSE. ON A NORMAL RETURN.
C
C    = .TRUE. IF THE STEP SIZE  H  WAS TOO SMALL.  H  HAS BEEN
C      INCREASED TO AN ACCEPTABLE VALUE, WITH WHICH  STEPNS  MAY BE
C      CALLED AGAIN.
C
C    = .TRUE. IF  RELERR  AND/OR  ABSERR  WERE TOO SMALL.  THEY HAVE
C      BEEN INCREASED TO ACCEPTABLE VALUES, WITH WHICH  STEPNS  MAY
C      BE CALLED AGAIN.
C
C HOLD = ||Y - YOLD||.
C
C H = OPTIMAL VALUE FOR NEXT STEP TO BE ATTEMPTED.  NORMALLY  H  SHOULD
C    NOT BE MODIFIED BY THE USER.
C
C RELERR, ABSERR  ARE UNCHANGED ON A NORMAL RETURN.
C
C S = (APPROXIMATE) ARC LENGTH ALONG THE ZERO CURVE OF THE HOMOTOPY MAP 
C    UP TO THE LATEST POINT FOUND, WHICH IS RETURNED IN  Y .
C
C Y, YP, YOLD, YPOLD  CONTAIN THE TWO MOST RECENT POINTS AND TANGENT
C    VECTORS FOUND ON THE ZERO CURVE OF THE HOMOTOPY MAP.
C
C
C CALLS  DNRM2 , TANGNS .
      USE HOMPACK_KINDS, ONLY: ZERO, ONE
      USE HOMPACK_CORE_LEGACY, ONLY: TANGNS
      USE BLAS_INTERFACES, ONLY: DNRM2
      IMPLICIT NONE
C
      INTEGER, INTENT(IN):: LENQR,MODE,N
      INTEGER, INTENT(IN OUT):: IFLAG,NFE
      LOGICAL, INTENT(IN OUT):: CRASH,START
      REAL(DP), INTENT(IN):: A(:),SSPAR(8)
      REAL(DP), INTENT(IN OUT):: ABSERR,H,HOLD,RELERR,S,
     &  Y(:),YOLD(:),YP(:),YPOLD(:)
      REAL(DP), INTENT(OUT), DIMENSION(:):: TZ,W,WP,Z0,Z1
C
C *****  LOCAL VARIABLES.  *****
C
      REAL(DP):: DCALC,DD001,DD0011,DD01,DD011,DELS,F0,F1,
     &   FOURU,FP0,FP1,HFAIL,HT,LCALC,QOFS,RCALC,RHOLEN,TEMP,TWOU
      INTEGER:: ITNUM,J,JUDY,NP1
      LOGICAL:: FAIL
C
C THE LIMIT ON THE NUMBER OF NEWTON ITERATIONS ALLOWED BEFORE REDUCING
C THE STEP SIZE  H  MAY BE CHANGED BY CHANGING THE FOLLOWING PARAMETER 
C STATEMENT:
      INTEGER, PARAMETER:: LITFH=4
C
C DEFINITION OF HERMITE CUBIC INTERPOLANT VIA DIVIDED DIFFERENCES.
C
      DD01(F0,F1,DELS)=(F1-F0)/DELS
      DD001(F0,FP0,F1,DELS)=(DD01(F0,F1,DELS)-FP0)/DELS
      DD011(F0,F1,FP1,DELS)=(FP1-DD01(F0,F1,DELS))/DELS
      DD0011(F0,FP0,F1,FP1,DELS)=(DD011(F0,F1,FP1,DELS) - 
     &                            DD001(F0,FP0,F1,DELS))/DELS
      QOFS(F0,FP0,F1,FP1,DELS,S)=((DD0011(F0,FP0,F1,FP1,DELS)*(S-DELS) +
     &   DD001(F0,FP0,F1,DELS))*S + FP0)*S + F0
C
C ***** END OF SPECIFICATION INFORMATION. *****
C
C
      TWOU=2*EPSILON(ONE)
      FOURU=TWOU+TWOU
      NP1=N+1
      CRASH=.TRUE.
C THE ARCLENGTH  S  MUST BE NONNEGATIVE.
      IF (S .LT. ZERO) RETURN
C IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE.
      IF (H .LT. FOURU*(ONE+S)) THEN
        H=FOURU*(ONE+S)
        RETURN
      ENDIF
C IF ERROR TOLERANCES ARE TOO SMALL, INCREASE THEM TO ACCEPTABLE VALUES.
      TEMP=DNRM2(NP1,Y,1)+ONE
      IF (0.5_DP*(RELERR*TEMP+ABSERR) .LT. TWOU*TEMP) THEN
        IF (RELERR .NE. ZERO) THEN
          RELERR=FOURU*(ONE+FOURU)
          ABSERR=MAX(ABSERR,ZERO)
        ELSE
          ABSERR=FOURU*TEMP
        ENDIF
        RETURN
      ENDIF
      CRASH=.FALSE.
      STARTUP: IF (START) THEN
C
C *****  STARTUP SECTION (FIRST STEP ALONG ZERO CURVE).  *****
C
      FAIL=.FALSE.
      START=.FALSE.
C DETERMINE SUITABLE INITIAL STEP SIZE.
      H=MIN(H, .10_dp, SQRT(SQRT(RELERR*TEMP+ABSERR)))
C USE LINEAR PREDICTOR ALONG TANGENT DIRECTION TO START NEWTON ITERATION.
      YPOLD(NP1)=ONE
      YPOLD(1:N)=ZERO
      CALL TANGNS(S,Y,YP,TZ,YPOLD,A,MODE,LENQR,NFE,N,IFLAG)
      IF (IFLAG .GT. 0) RETURN
      LP: DO
      W=Y + H*YP
      Z0=W
      DO JUDY=1,LITFH
        RHOLEN = -ONE
C CALCULATE THE NEWTON STEP  TZ  AT THE CURRENT POINT  W .
        CALL TANGNS(RHOLEN,W,WP,TZ,YPOLD,A,MODE,LENQR,NFE,N,IFLAG)
        IF (IFLAG .GT. 0) RETURN
C
C TAKE NEWTON STEP AND CHECK CONVERGENCE.
        W=W + TZ
        ITNUM=JUDY
C COMPUTE QUANTITIES USED FOR OPTIMAL STEP SIZE ESTIMATION.
        IF (JUDY .EQ. 1) THEN
          LCALC=DNRM2(NP1,TZ,1)
          RCALC=RHOLEN
          Z1=W
        ELSE IF (JUDY .EQ. 2) THEN
          LCALC=DNRM2(NP1,TZ,1)/LCALC
          RCALC=RHOLEN/RCALC
        ENDIF
C GO TO MOP-UP SECTION AFTER CONVERGENCE.
        IF (DNRM2(NP1,TZ,1) .LE. RELERR*DNRM2(NP1,W,1)+ABSERR)
     &                                                 GO TO 600
C
      END DO
C
C NO CONVERGENCE IN  LITFH  ITERATIONS.  REDUCE  H  AND TRY AGAIN.
      IF (H .LE. FOURU*(ONE + S)) THEN
        IFLAG=6
        RETURN
      ENDIF
      H = H/2
      END DO LP
      END IF STARTUP
C
C ***** END OF STARTUP SECTION. *****
C
C ***** PREDICTOR SECTION. *****
C
      FAIL=.FALSE.
C COMPUTE POINT PREDICTED BY HERMITE INTERPOLANT.  USE STEP SIZE  H
C COMPUTED ON LAST CALL TO  STEPNF .
      HP: DO
      DO J=1,NP1
        W(J)=QOFS(YOLD(J),YPOLD(J),Y(J),YP(J),HOLD,HOLD+H)
      END DO
      Z0=W 
C
C ***** END OF PREDICTOR SECTION. *****
C
C ***** CORRECTOR SECTION. *****
C
      CORRECTOR: DO JUDY=1,LITFH
        RHOLEN = -ONE
C CALCULATE THE NEWTON STEP  TZ  AT THE CURRENT POINT  W .
        CALL TANGNS(RHOLEN,W,WP,TZ,YP,A,MODE,LENQR,NFE,N,IFLAG)
        IF (IFLAG .GT. 0) RETURN
C
C TAKE NEWTON STEP AND CHECK CONVERGENCE.
        W=W + TZ
        ITNUM=JUDY
C COMPUTE QUANTITIES USED FOR OPTIMAL STEP SIZE ESTIMATION.
        IF (JUDY .EQ. 1) THEN
          LCALC=DNRM2(NP1,TZ,1)
          RCALC=RHOLEN
          Z1=W
        ELSE IF (JUDY .EQ. 2) THEN
          LCALC=DNRM2(NP1,TZ,1)/LCALC
          RCALC=RHOLEN/RCALC
        ENDIF
C GO TO MOP-UP SECTION AFTER CONVERGENCE.
        IF (DNRM2(NP1,TZ,1) .LE. RELERR*DNRM2(NP1,W,1)+ABSERR)
     &                                                 GO TO 600
C
      END DO CORRECTOR
C
C NO CONVERGENCE IN  LITFH  ITERATIONS.  RECORD FAILURE AT CALCULATED  H , 
C SAVE THIS STEP SIZE, REDUCE  H  AND TRY AGAIN.
      FAIL=.TRUE.
      HFAIL=H
      IF (H .LE. FOURU*(ONE + S)) THEN
        IFLAG=6
        RETURN
      ENDIF
      H = H/2
      END DO HP
C
C ***** END OF CORRECTOR SECTION. *****
C
C ***** MOP-UP SECTION. *****
C
C YOLD  AND  Y  ALWAYS CONTAIN THE LAST TWO POINTS FOUND ON THE ZERO
C CURVE OF THE HOMOTOPY MAP.  YPOLD  AND  YP  CONTAIN THE TANGENT
C VECTORS TO THE ZERO CURVE AT  YOLD  AND  Y , RESPECTIVELY.
C
600   YPOLD=YP
      YOLD=Y
      Y=W
      YP=WP
      W=Y - YOLD
C UPDATE ARC LENGTH.
      HOLD=DNRM2(NP1,W,1)
      S=S+HOLD
C
C ***** END OF MOP-UP SECTION. *****
C
C ***** OPTIMAL STEP SIZE ESTIMATION SECTION. *****
C
C CALCULATE THE DISTANCE FACTOR  DCALC .
      TZ=Z0 - Y
      W=Z1 - Y
      DCALC=DNRM2(NP1,TZ,1)
      IF (DCALC .NE. ZERO) DCALC=DNRM2(NP1,W,1)/DCALC
C
C THE OPTIMAL STEP SIZE HBAR IS DEFINED BY
C
C   HT=HOLD * [MIN(LIDEAL/LCALC, RIDEAL/RCALC, DIDEAL/DCALC)]**(1/P)
C
C     HBAR = MIN [ MAX(HT, BMIN*HOLD, HMIN), BMAX*HOLD, HMAX ]
C
C IF CONVERGENCE HAD OCCURRED AFTER 1 ITERATION, SET THE CONTRACTION
C FACTOR  LCALC  TO ZERO.
      IF (ITNUM .EQ. 1) LCALC = ZERO
C FORMULA FOR OPTIMAL STEP SIZE.
      IF (LCALC+RCALC+DCALC .EQ. ZERO) THEN
        HT = SSPAR(7) * HOLD
      ELSE 
        HT = (ONE/MAX(LCALC/SSPAR(1), RCALC/SSPAR(2), DCALC/SSPAR(3)))
     &       **(ONE/SSPAR(8)) * HOLD
      ENDIF
C  HT  CONTAINS THE ESTIMATED OPTIMAL STEP SIZE.  NOW PUT IT WITHIN
C REASONABLE BOUNDS.
      H=MIN(MAX(HT,SSPAR(6)*HOLD,SSPAR(4)), SSPAR(7)*HOLD, SSPAR(5))
      IF (ITNUM .EQ. 1) THEN
C IF CONVERGENCE HAD OCCURRED AFTER 1 ITERATION, DON'T DECREASE  H .
        H=MAX(H,HOLD)
      ELSE IF (ITNUM .EQ. LITFH) THEN
C IF CONVERGENCE REQUIRED THE MAXIMUM  LITFH  ITERATIONS, DON'T
C INCREASE  H .
        H=MIN(H,HOLD)
      ENDIF
C IF CONVERGENCE DID NOT OCCUR IN  LITFH  ITERATIONS FOR A PARTICULAR
C H = HFAIL , DON'T CHOOSE THE NEW STEP SIZE LARGER THAN  HFAIL .
      IF (FAIL) H=MIN(H,HFAIL)
C
      END SUBROUTINE STEPNS
C
      END MODULE HOMPACK_NS