      MODULE HOMPACK_QS
C
      USE HOMPACK_KINDS, ONLY: DP
      IMPLICIT NONE
C
      CONTAINS
C
      SUBROUTINE FIXPQS(N,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,A,
     &     NFE,ARCLEN,MODE,LENQR,SSPAR)
C
C Subroutine  FIXPQS  finds a fixed point or zero of the 
C N-dimensional vector function  F(X), or tracks a zero curve of a 
C general homotopy map  RHO(A,X,LAMBDA).  For the fixed point problem
C F(X) is assumed to be a C2 map of some ball into itself.  The 
C equation  X=F(X)  is solved by following the zero curve of the 
C homotopy map
C
C  LAMBDA*(X - F(X)) + (1 - LAMBDA)*(X - A),
C
C starting from  LAMBDA = 0, X = A.   The curve is parameterized
C by arc length  S, and is followed by solving the ordinary 
C differential equation  D(HOMOTOPY MAP)/DS = 0  for  
C Y(S) = (X(S),LAMBDA(S)).  This is done by using a Hermite cubic 
C predictor and a corrector which returns to the zero curve in a 
C hyperplane perpendicular to the tangent to the zero curve at the 
C most recent point.
C
C For the zero finding problem  F(X)  is assumed to be a C2 map such
C that for some  R > 0,  X*F(X) >= 0  whenever  NORM(X) = R.
C The equation  F(X) = 0  is solved by following the zero curve of
C the homotopy map
C
C  LAMBDA*F(X) + (1 - LAMBDA)*(X - A)
C
C emanating from  LAMBDA = 0, X = A.
C
C A  must be an interior point of the above mentioned balls.
C
C For the curve tracking problem RHO(A,X,LAMBDA) is assumed to 
C be a C2 map from  E**M X [0,1) X E**N  into  E**N, which for 
C almost all parameter vectors  A  in some nonempty open subset
C of E**M satisfies
C
C  rank [D RHO(A,X,LAMBDA)/D LAMBDA, D RHO(A,X,LAMBDA)/DX] = N
C
C for all points  (X,LAMBDA)  such that  RHO(A,X,LAMBDA) = 0.  It is
C further assumed that
C
C         rank [ D RHO(A,X0,0)/DX ] = N.
C
C With  A  fixed, the zero curve of  RHO(A,X,LAMBDA)  emanating from
C LAMBDA = 0, X = X0  is tracked until  LAMBDA = 1  by solving the 
C ordinary differential equation  D RHO(A,X(S),LAMBDA(S))/DS = 0
C for  Y(S) = (X(S),LAMBDA(S)), where  S  is arc length along the
C zero curve.  Also the homotopy map  RHO(A,X,LAMBDA)  is assumed to
C be constructed such that
C
C         D LAMBDA(0)/DS > 0.
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
C FIXPQS directly or indirectly uses the subroutines  
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
C Y(1:N+1)  contains the starting point for tracking the homotopy map.
C    (Y(1),...,Y(N)) = A  for the fixed point and zero finding 
C    problems.  (Y(1),...,Y(N)) = X0  for the curve tracking problem.
C    Y(N+1)  need not be defined by the user.
C
C IFLAG  can be -2, -1, 0, 2, or 3.  IFLAG  should be 0 on the first
C    call to  FIXPQS  for the problem  X=F(X), -1 for the problem
C    F(X)=0, and -2 for the problem  RHO(A,X,LAMBDA)=0.   In certain
C    situations  IFLAG  is set to 2 or 3 by  FIXPQS, and  FIXPQS  can
C    be called again without changing  IFLAG.
C
C ARCRE, ARCAE  are the relative and absolute errors, respectively,
C    allowed the iteration along the zero curve.  If
C    ARC?E .LE. 0.0  on input, it is reset to  .5*SQRT(ANS?E).
C    Normally  ARC?E  should be considerably larger than  ANS?E.
C
C ANSRE, ANSAE  are the relative and absolute error values used for 
C    the answer at  LAMBDA = 1.  The accepted answer  Y = (X,LAMBDA)
C    satisfies
C
C      |Y(1) - 1| .LE. ANSRE + ANSAE      .AND.
C  
C      ||DZ|| .LE. ANSRE*||Y|| + ANSAE      where
C
C      DZ is the Newton step to Y.
C
C TRACE  is an integer specifying the logical I/O unit for
C    intermediate output.  If  TRACE .GT. 0  the points computed on
C    the zero curve are written to I/O unit  TRACE .
C
C A(:)  contains the parameter vector  A.  For the fixed point
C    and zero finding problems,  A  need not be initialized by the 
C    user, and is assumed to have length  N.  For the curve
C    tracking problem,  A  must be initialized by the user.
C
C MODE = 1 if the Jacobian matrix is symmetric and stored in a packed
C          skyline format;
C      = 2 if the Jacobian matrix is stored in a sparse row format.
C
C LENQR  is the number of nonzero entries in the sparse Jacobian
C    matrices, used to determine the sparse matrix data structures.
C
C SSPAR(1:4) =  (HMIN, HMAX, BMIN, BMAX)  is a vector of parameters 
C    used for the optimal step size estimation.  A default value
C    can be specified for any of these four parameters by setting it
C    .LE. 0.0  on input.  See the comments in  STEPQS  for more
C    information about these parameters.
C
C
C ON OUTPUT:
C
C N , TRACE , A , LENQR  are unchanged.
C
C Y(N+1) = LAMBDA, (Y(1),...,Y(N)) = X, and  Y  is an approximate
C    zero of the homotopy map.  Normally  LAMBDA = 1  and  X  is a
C    fixed point or zero of  F(X).   In abnormal situations,  LAMBDA
C    may only be near 1 and  X  near a fixed point or zero.
C
C IFLAG =
C
C   1   Normal return.
C
C   2   Specified error tolerance cannot be met.  Some or all of
C       ARCRE, ARCAE, ANSRE, ANSAE  have been increased to 
C       suitable values.  To continue, just call  FIXPQS  again
C       without changing any parameters.
C
C   3   STEPQS  has been called 1000 times.  To continue, call
C       FIXPQS  again without changing any parameters.
C
C   4   Jacobian matrix does not have full rank.  The algorithm
C       has failed (the zero curve of the homotopy map cannot be
C       followed any further).
C
C   5   The tracking algorithm has lost the zero curve of the 
C       homotopy map and is not making progress.  The error 
C       tolerances  ARC?E  and  ANS?E  were too lenient.  The problem 
C       should be restrarted by calling  FIXPQS  with smaller error 
C       tolerances and  IFLAG = 0 (-1, -2).
C
C   6   The Newton iteration in  STEPQS  or  ROOTNS  failed to
C       converge.  The error tolerances  ANS?E  may be too stringent.
C
C   7   Illegal input parameters, a fatal error.
C
C ARCRE, ARCAE, ANSRE, ANSAE  are unchanged after a normal return
C    (IFLAG = 1).  They are increased to appropriate values on the
C    return  IFLAG = 2.
C
C NFE  is the number of Jacobian evaluations.
C
C ARCLEN  is the approximate length of the zero curve.  
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
C    here, and distributed via the module  HOMOTOPY .
C
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
      REAL (dp), INTENT(IN OUT)::ANSAE,ANSRE,ARCAE,ARCRE,SSPAR(4)
      INTEGER, INTENT(OUT)::NFE
      REAL (dp), INTENT(OUT)::ARCLEN
C
C     LOCAL VARIABLES 
C
      REAL (dp), SAVE:: ABSERR, H, HOLD, RELERR, S, WK 
      INTEGER, SAVE:: IFLAGC, ITER, JW, LIMIT, NP1
      LOGICAL, SAVE:: CRASH, START       
C
C     WORK ARRAYS 
C
      REAL (dp), ALLOCATABLE, DIMENSION(:), SAVE:: YP,YOLD,YPOLD
      REAL (dp):: DZ(N+1),T(N+1),Z0(N+1) 
C
C LIMITD  IS AN UPPER BOUND ON THE NUMBER OF STEPS.  IT MAY BE 
C CHANGED BY CHANGING THE FOLLOWING PARAMETER STATEMENT:
      INTEGER, PARAMETER:: LIMITD = 1000
C
C ***** FIRST EXECUTABLE STATEMENT *****
C
C CHECK IFLAG
C
      IF (N .LE. 0  .OR.  ANSRE .LE. 0.0  .OR.  ANSAE .LT. 0.0
     &  .OR.  N+1 .NE. SIZE(Y) .OR.
     &  ((IFLAG .EQ. -1  .OR.  IFLAG .EQ. 0) .AND.  N .NE. SIZE(A))
     &  .OR.  MODE .LE. 0  .OR.  MODE .GE. 3)
     &                                                     IFLAG=7
      IF (IFLAG .GE. -2 .AND. IFLAG .LE. 0) GO TO 10
      IF (IFLAG .EQ. 2) GO TO 50
      IF (IFLAG .EQ. 3) GO TO 40
C
C ONLY VALID INPUT FOR IFLAG IS -2, -1, 0, 2, 3.
C
      IFLAG = 7
      RETURN
C
C ***** INITIALIZATION BLOCK  *****
C
 10   ARCLEN = 0.0
      IF (ARCRE .LE. 0.0) ARCRE = .5*SQRT(ANSRE)
      IF (ARCAE .LE. 0.0) ARCAE = .5*SQRT(ANSAE)
      NFE=0
      IFLAGC = IFLAG
      NP1=N+1
C ALLOCATE SAVED WORK ARRAYS.
      IF (ALLOCATED(YOLD)) DEALLOCATE(YOLD)
      IF (ALLOCATED(YP)) DEALLOCATE(YP)
      IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
      ALLOCATE(YOLD(NP1),YP(NP1),YPOLD(NP1))
C 
C SET INITIAL CONDITIONS FOR FIRST CALL TO STEPQS.
C
        START=.TRUE.
        CRASH=.FALSE.
        RELERR = ARCRE
        ABSERR = ARCAE
        HOLD=1.0
        H=0.1
        S=0.0
        YPOLD(NP1) = 1.0
        Y(NP1) = 0.0
        YPOLD(1:N)=0.0
C
C SET OPTIMAL STEP SIZE ESTIMATION PARAMETERS.
C
C     MINIMUM STEP SIZE HMIN
      IF (SSPAR(1) .LE. 0.0) SSPAR(1)=(SQRT(N+1.0)+4.0)*EPSILON(1.0_dp)
C     MAXIMUM STEP SIZE HMAX
      IF (SSPAR(2) .LE. 0.0) SSPAR(2)= 1.0
C     MINIMUM STEP REDUCTION FACTOR BMIN
      IF (SSPAR(3) .LE. 0.0) SSPAR(3)= 0.1_dp
C     MAXIMUM STEP EXPANSION FACTOR BMAX
      IF (SSPAR(4) .LE. 0.0) SSPAR(4)= 7.0
C
C LOAD  A  FOR THE FIXED POINT AND ZERO FINDING PROBLEMS.
C
      IF (IFLAGC .GE. -1) A(1:N) = Y(1:N)
C
 40   LIMIT=LIMITD
C ALLOCATE ARRAYS FOR SPARSE JACOBIAN MATRIX DATA STRUCTURE.
 50   SELECT CASE (MODE)
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
C ***** END OF INITIALIZATION BLOCK. *****
C
      MAIN_LOOP: DO ITER=1,LIMIT ! ***** MAIN LOOP. *****
        IF (Y(NP1) .LT. 0.0) THEN
          ARCLEN = S
          IFLAG = 5
          CALL CLEANUPALL
          RETURN
        END IF
C
C TAKE A STEP ALONG THE CURVE.
C
        CALL STEPQS(N,NFE,IFLAGC,MODE,LENQR,START,CRASH,HOLD,H,
     &    WK,RELERR,ABSERR,S,Y,YP,YOLD,YPOLD,A,Z0,DZ,T,SSPAR)
C
C PRINT LATEST POINT ON CURVE IF REQUESTED.
C
        IF (TRACE .GT. 0) THEN
          WRITE (TRACE,217) ITER,NFE,S,Y(NP1),(Y(JW),JW=1,N)
217       FORMAT(/' STEP',I5,3X,'NFE =',I5,3X,'ARC LENGTH =',F9.4,3X,
     &    'LAMBDA =',F7.4,5X,'X vector:'/(1X,6ES12.4))
        ENDIF
C
C CHECK IF THE STEP WAS SUCCESSFUL.
C
        IF (IFLAGC .GT. 0) THEN
          ARCLEN=S
          IFLAG=IFLAGC
          CALL CLEANUPALL
          RETURN
        END IF
C
        IF (CRASH) THEN
C
C         RETURN CODE FOR ERROR TOLERANCE TOO SMALL.
C      
          IFLAG=2
C
C         CHANGE ERROR TOLERANCES.
C
          IF (ARCRE .LT. RELERR) THEN
            ARCRE=RELERR
            ANSRE=RELERR
          ENDIF
          IF (ARCAE .LT. ABSERR) ARCAE=ABSERR
C
C         CHANGE LIMIT ON NUMBER OF ITERATIONS.
C
          LIMIT = LIMIT - ITER
          CALL CLEANUP
          RETURN
        END IF
C
C IF  LAMBDA >= 1.0,  USE  ROOTNS  TO FIND SOLUTION.
C
        IF (Y(NP1) .GE. 1.0) GO TO 500
C
      END DO MAIN_LOOP   ! ***** END OF MAIN LOOP *****
C
C DID NOT CONVERGE IN  LIMIT  ITERATIONS, SET  IFLAG  AND RETURN.
C
      ARCLEN = S
      IFLAG = 3
      CALL CLEANUP
      RETURN
C
C ***** FINAL STEP -- FIND SOLUTION AT LAMBDA=1 *****
C
C SAVE  YOLD  FOR ARC LENGTH CALCULATION LATER.
C
 500  T = YOLD
C
C FIND SOLUTION.
C
      CALL ROOTNS(N,NFE,IFLAGC,ANSRE,ANSAE,Y,YP,
     &    YOLD,YPOLD,A,MODE,LENQR)
C
C CHECK IF SOLUTION WAS FOUND AND SET  IFLAG  ACCORDINGLY.
C
      IFLAG=1
C
C     SET ERROR FLAG IF ROOTNS COULD NOT GET THE POINT ON THE ZERO
C     CURVE AT  LAMBDA = 1.0 .
C
      IF (IFLAGC .GT. 0) IFLAG=IFLAGC
C
C CALCULATE FINAL ARC LENGTH.
C
      DZ = Y - T
      ARCLEN = S - HOLD + DNRM2(NP1,DZ,1)
C
C ***** END OF FINAL STEP *****
C
      CALL CLEANUPALL
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
      END SUBROUTINE FIXPQS
C
      SUBROUTINE STEPQS(N,NFE,IFLAG,MODE,LENQR,START,CRASH,HOLD,H,
     &  WK,RELERR,ABSERR,S,Y,YP,YOLD,YPOLD,A,Z0,DZ,T,SSPAR)
C
C SUBROUTINE  STEPQS  TAKES ONE STEP ALONG THE ZERO CURVE OF THE 
C HOMOTOPY MAP  RHO(X,LAMBDA)  USING A PREDICTOR-CORRECTOR ALGORITHM.
C THE PREDICTOR USES A HERMITE CUBIC INTERPOLANT, AND THE CORRECTOR 
C RETURNS TO THE ZERO CURVE USING A NEWTON ITERATION, REMAINING
C IN A HYPERPLANE PERPENDICULAR TO THE MOST RECENT TANGENT VECTOR.
C  STEPQS  ALSO ESTIMATES A STEP SIZE  H  FOR THE NEXT STEP ALONG THE 
C ZERO CURVE.  SEE ALSO THE REVERSE CALL ROUTINE  STEPNX .
C 
C THE CALLING PROGRAM MUST CONTAIN THE FOLLOWING INTERFACE BLOCK:
C
C     INTERFACE
C       SUBROUTINE STEPQS(N,NFE,IFLAG,MODE,LENQR,START,CRASH,HOLD,H,
C    &    WK,RELERR,ABSERR,S,Y,YP,YOLD,YPOLD,A,Z0,DZ,T,SSPAR)
C       USE REAL_PRECISION
C       INTEGER, INTENT(IN):: LENQR,MODE,N
C       INTEGER, INTENT(IN OUT):: IFLAG,NFE
C       LOGICAL, INTENT(IN OUT):: CRASH,START
C       REAL(DP), INTENT(IN):: A(:),SSPAR(4)
C       REAL(DP), INTENT(IN OUT):: ABSERR,H,HOLD,RELERR,S,WK,
C    &    Y(:),YOLD(:),YP(:),YPOLD(:)
C       REAL(DP), INTENT(OUT), DIMENSION(:):: DZ,T,Z0
C       END SUBROUTINE STEPQS
C     END INTERFACE
C
C
C ON INPUT:
C 
C N = DIMENSION OF  X. 
C
C NFE = NUMBER OF JACOBIAN MATRIX EVALUATIONS.
C
C IFLAG = -2, -1, OR 0, INDICATING THE PROBLEM TYPE.
C
C MODE = 1 IF THE JACOBIAN MATRIX IS SYMMETRIC AND STORED IN A PACKED
C          SKYLINE FORMAT;
C      = 2 IF THE JACOBIAN MATRIX IS STORED IN A SPARSE ROW FORMAT.
C
C LENQR  IS THE NUMBER OF NONZERO ENTRIES IN THE SPARSE JACOBIAN
C    MATRICES, USED TO DETERMINE THE SPARSE MATRIX DATA STRUCTURES.
C
C START = .TRUE. ON FIRST CALL TO  STEPQS, .FALSE. OTHERWISE.
C         SHOULD NOT BE MODIFIED BY THE USER AFTER THE FIRST CALL.
C
C HOLD = ||Y - YOLD|| ; SHOULD NOT BE MODIFIED BY THE USER.
C
C H = UPPER LIMIT ON LENGTH OF STEP THAT WILL BE ATTEMPTED.  H  MUST
C    BE SET TO A POSITIVE NUMBER ON THE FIRST CALL TO  STEPQS.
C    THEREAFTER,  STEPQS  CALCULATES AN OPTIMAL VALUE FOR  H, AND  H
C    SHOULD NOT BE MODIFIED BY THE USER.
C
C WK = APPROXIMATE CURVATURE FOR THE LAST STEP (COMPUTED BY PREVIOUS
C    CALL TO  STEPQS).  UNDEFINED ON FIRST CALL.  SHOULD NOT BE
C    MODIFIED BY THE USER.
C  
C RELERR, ABSERR = RELATIVE AND ABSOLUTE ERROR VALUES.  THE ITERATION
C    IS CONSIDERED TO HAVE CONVERGED WHEN A POINT  Z=(X,LAMBDA)  IS 
C    FOUND SUCH THAT
C       ||DZ|| .LE. RELERR*||Z|| + ABSERR,
C    WHERE  DZ  IS THE LAST NEWTON STEP.
C
C S  = (APPROXIMATE) ARC LENGTH ALONG THE HOMOTOPY ZERO CURVE UP TO
C    Y(S) = (X(S),LAMBDA(S)).
C
C Y(1:N+1) = PREVIOUS POINT (X(S),LAMBDA(S)) FOUND ON THE ZERO CURVE
C    OF THE HOMOTOPY MAP.
C
C YP(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE HOMOTOPY
C    MAP AT  Y.  INPUT IN THIS VECTOR IS NOT USED ON THE FIRST CALL
C    TO  STEPQS.
C
C YOLD(1:N+1) = A POINT BEFORE  Y  ON THE ZERO CURVE OF THE HOMOTOPY
C    MAP.  INPUT IN THIS VECTOR IS NOT USED ON THE FIRST CALL TO 
C    STEPQS.
C
C YPOLD(1:N+1) = UNIT TANGENT VECTOR TO THE ZERO CURVE OF THE 
C    HOMOTOPY MAP AT  YOLD.
C
C A(:) = PARAMETER VECTOR IN THE HOMOTOPY MAP.
C
C QR(1:LENQR), PP(1:N), ROWPOS(1:N+2), COLPOS(1:LENQR) ARE ALL WORK
C    ARRAYS USED TO DEFINE THE SPARSE JACOBIAN MATRICES, ALLOCATED
C    IN FIXPQS, AND DISTRIBUTED VIA THE MODULE  HOMOTOPY .
C
C Z0(1:N+1), DZ(1:N+1), T(1:N+1)  ARE ALL WORK ARRAYS USED TO
C    CALCULATE THE TANGENT VECTORS AND NEWTON STEPS.
C    
C SSPAR(1:4) = PARAMETERS USED FOR COMPUTATION OF THE OPTIMAL STEP SIZE.
C    SSPAR(1) = HMIN, SSPAR(2) = HMAX, SSPAR(3) = BMIN, SSPAR(4) = BMAX.
C    THE OPTIMAL STEP  H  IS RESTRICTED SUCH THAT 
C       HMIN .LE. H .LE. HMAX, AND  BMIN*HOLD .LE. H .LE. BMAX*HOLD.
C
C
C ON OUTPUT:
C
C N, LENQR, A  ARE UNCHANGED.
C
C NFE HAS BEEN UPDATED.
C
C IFLAG
C
C    = -2, -1, OR 0 (UNCHANGED) ON A NORMAL RETURN.
C
C    = 4 IF A JACOBIAN MATRIX WITH RANK <  N  HAS OCCURRED.  THE
C        ITERATION WAS NOT COMPLETED.
C
C    = 6 IF THE ITERATION FAILED TO CONVERGE. 
C
C START = .FALSE. ON A NORMAL RETURN.
C
C CRASH 
C
C    = .FALSE. ON A NORMAL RETURN.
C
C    = .TRUE. IF THE STEP SIZE  H  WAS TOO SMALL.  H  HAS BEEN
C      INCREASED TO AN ACCEPTABLE VALUE, WITH WHICH  STEPQS  MAY BE
C      CALLED AGAIN.
C
C    = .TRUE. IF  RELERR  AND/OR  ABSERR  WERE TOO SMALL.  THEY HAVE
C      BEEN INCREASED TO ACCEPTABLE VALUES, WITH WHICH  STEPQS  MAY
C      BE CALLED AGAIN.
C
C HOLD = ||Y-YOLD||.
C
C H = OPTIMAL VALUE FOR NEXT STEP TO BE ATTEMPTED.  NORMALLY  H  SHOULD
C     NOT BE MODIFIED BY THE USER.
C
C WK = APPROXIMATE CURVATURE FOR THE STEP TAKEN BY  STEPQS.
C
C S = (APPROXIMATE) ARC LENGTH ALONG THE ZERO CURVE OF THE HOMOTOPY 
C     MAP UP TO THE LATEST POINT FOUND, WHICH IS RETURNED IN  Y.
C
C RELERR, ABSERR  ARE UNCHANGED ON A NORMAL RETURN.  THEY ARE POSSIBLY
C     CHANGED IF  CRASH  = .TRUE. (SEE DESCRIPTION OF  CRASH  ABOVE).
C
C Y, YP, YOLD, YPOLD  CONTAIN THE TWO MOST RECENT POINTS AND TANGENT
C     VECTORS FOUND ON THE ZERO CURVE OF THE HOMOTOPY MAP.
C
C
C CALLS  DNRM2, TANGNS.
C
      USE HOMPACK_KINDS, ONLY: ZERO, ONE
      USE HOMPACK_CORE_LEGACY, ONLY: TANGNS
      USE BLAS_INTERFACES, ONLY: DNRM2 
      IMPLICIT NONE
C
      INTEGER, INTENT(IN):: LENQR,MODE,N
      INTEGER, INTENT(IN OUT):: IFLAG,NFE
      LOGICAL, INTENT(IN OUT):: CRASH,START
      REAL(DP), INTENT(IN):: A(:),SSPAR(4)
      REAL(DP), INTENT(IN OUT):: ABSERR,H,HOLD,RELERR,S,WK,
     &    Y(:),YOLD(:),YP(:),YPOLD(:)
      REAL(DP), INTENT(OUT), DIMENSION(:):: DZ,T,Z0
C
      REAL(DP):: DD001,DD0011,DD01,DD011,QOFS
C
C     LOCAL VARIABLES.
C
      REAL(DP), SAVE:: ACOF(12), ALPHA, CORDIS, DELS, FOURU,
     &  GAMMA, HFAIL, HTEMP, IDLERR, OMEGA, P0, P1, PP0, PP1, 
     &  SIGMA, TEMP, THETA, TWOU, WKOLD, WRGE(8), XSTEP
      INTEGER:: I, ITCNT, LK, LST, NP1
      LOGICAL:: FAILED
      DATA WRGE  /
     &   .8735115E+00_dp, .1531947E+00_dp, .3191815E-01_dp,
     &   .3339946E-10_dp, .4677788E+00_dp, .6970123E-03_dp,
     &   .1980863E-05_dp, .1122789E-08_dp/
      DATA ACOF  /
     &   .9043128E+00_dp, -.7075675E+00_dp, -.4667383E+01_dp,
     &  -.3677482E+01_dp,  .8516099E+00_dp, -.1953119E+00_dp,
     &  -.4830636E+01_dp, -.9770528E+00_dp,  .1040061E+01_dp,
     &   .3793395E-01_dp,  .1042177E+01_dp,  .4450706E-01_dp/
C
C THE LIMIT ON THE NUMBER OF NEWTON ITERATIONS ALLOWED BEFORE REDUCING
C THE STEP SIZE  H  MAY BE CHANGED BY CHANGING THE FOLLOWING PARAMETER 
C STATEMENT:
      INTEGER, PARAMETER:: LITFH = 10
C
C DEFINITION OF HERMITE CUBIC INTERPOLANT VIA DIVIDED DIFFERENCES.
C
      DD01(P0,P1,DELS) = (P1-P0)/DELS
      DD001(P0,PP0,P1,DELS) = (DD01(P0,P1,DELS)-PP0)/DELS
      DD011(P0,P1,PP1,DELS) = (PP1-DD01(P0,P1,DELS))/DELS
      DD0011(P0,PP0,P1,PP1,DELS) = (DD011(P0,P1,PP1,DELS) -
     &  DD001(P0,PP0,P1,DELS))/DELS
      QOFS(P0,PP0,P1,PP1,DELS,S) = ((DD0011(P0,PP0,P1,PP1,DELS)*
     &  (S-DELS) + DD001(P0,PP0,P1,DELS))*S + PP0)*S + P0
C
C ***** END OF SPECIFICATION SECTION. *****
C
C ***** INITIALIZATION. *****
C
      TWOU = 2*EPSILON(ONE)
      FOURU = TWOU + TWOU
      NP1 = N+1
      FAILED = .FALSE.
      CRASH = .TRUE.
C 
C CHECK THAT ALL INPUT PARAMETERS ARE CORRECT.
C
C     THE ARCLENGTH  S  MUST BE NONNEGATIVE.
C
      IF (S .LT. ZERO) RETURN
C
C     IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE.
C   
      IF (H .LT. FOURU*(ONE+S)) THEN
          H=FOURU*(ONE + S)
          RETURN
      END IF
C
C     IF ERROR TOLERANCES ARE TOO SMALL, INCREASE THEM TO ACCEPTABLE 
C     VALUES.
C
      TEMP=DNRM2(NP1,Y,1) + ONE
      IF (0.5_DP*(RELERR*TEMP+ABSERR) .LT. TWOU*TEMP) THEN
          IF (RELERR .NE. ZERO) THEN
            RELERR = FOURU*(ONE+FOURU)
            TEMP = ZERO
            ABSERR = MAX(ABSERR,TEMP)
          ELSE
            ABSERR=FOURU*TEMP
          END IF
          RETURN
      END IF
C
C     INPUT PARAMETERS WERE ALL ACCEPTABLE.
C
      CRASH = .FALSE.
C
C COMPUTE  YP  ON FIRST CALL.
C
      IF (START) THEN
C
C         INITIALIZE THE IDEAL ERROR USED FOR STEP SIZE ESTIMATION.
C
          IDLERR=SQRT(SQRT(ABSERR))
C
          CALL TANGNS(S,Y,YP,DZ,YPOLD,A,MODE,LENQR,NFE,N,IFLAG)
          IF (IFLAG .GT. 0) RETURN
      END IF
C
      CONV: DO
C
C ***** COMPUTE PREDICTOR POINT Z0. *****
C
        IF (START) THEN
C           
C         COMPUTE Z0 WITH LINEAR PREDICTOR USING Y, YP --
C         
          Z0 = Y + H*YP
        ELSE
C
C         COMPUTE Z0 WITH CUBIC PREDICTOR.
C
          DO I=1,NP1
            Z0(I) = QOFS(YOLD(I),YPOLD(I),Y(I),YP(I),HOLD,HOLD+H) 
          END DO
        END IF
C
C ***** END OF PREDICTOR SECTION. *****
C
        NEWTON: DO ITCNT = 1,LITFH   ! ***** NEWTON ITERATION. *****
C
C COMPUTE TANGENT  T  AND MINIMUM NORM NEWTON STEP  DZ  AT THE
C CURRENT POINT  Z0 .
C
          TEMP = -ONE
          CALL TANGNS(TEMP,Z0,T,DZ,YP,A,MODE,LENQR,NFE,N,IFLAG)
          IF (IFLAG .GT. 0) RETURN
C
C CHECK THAT COMPUTED TANGENT  T  MAKES AN ANGLE NO LARGER THAN
C 60 DEGREES WITH CURRENT TANGENT  YP.  (I.E., COS OF ANGLE < .5)
C IF NOT, STEP SIZE WAS TOO LARGE, SO THROW AWAY Z0, AND TRY
C AGAIN WITH A SMALLER STEP.
C
          ALPHA = DOT_PRODUCT(T,YP)
          IF (ALPHA < 0.5_DP) EXIT NEWTON
C
C MAKE  DZ  ORTHOGONAL TO TANGENT DIRECTION  YP .
C
          SIGMA = -DOT_PRODUCT(DZ,YP)/DOT_PRODUCT(T,YP)
          DZ = DZ + SIGMA*T
C
C TAKE NEWTON STEP.
C
          Z0 = Z0 + DZ
C
C CHECK FOR CONVERGENCE.
C
          XSTEP=DNRM2(NP1,DZ,1)
          IF (XSTEP .LE. RELERR*DNRM2(NP1,Z0,1)+ABSERR) EXIT CONV
C
        END DO NEWTON   ! ***** END OF NEWTON LOOP. *****
C
C DIDN'T CONVERGE OR TANGENT AT NEW POINT DID NOT MAKE
C AN ANGLE SMALLER THAN 60 DEGREES WITH  YPOLD -- 
C TRY AGAIN WITH A SMALLER H.
C      
        FAILED = .TRUE.
        HFAIL = H
        IF (H .LE. FOURU*(ONE + S)) THEN
          IFLAG = 6
          RETURN
        ELSE
          H = H/2
        END IF
C
C END OF CONVERGENCE FAILURE SECTION.
C
      END DO CONV
C
C ***** CONVERGED -- MOP UP AND RETURN. *****
C
C COMPUTE TANGENT  T  AT  Z0 .
C
      CALL TANGNS(S,Z0,T,DZ,YP,A,MODE,LENQR,NFE,N,IFLAG)
      IF (IFLAG .GT. 0) RETURN
      ALPHA = DOT_PRODUCT(T,YP)
      ALPHA = ACOS(ALPHA)
C
C COMPUTE CORRECTOR DISTANCE.
C
      IF (START) THEN
        DZ = Y + H*YP
      ELSE
        DO I=1,NP1
          DZ(I)=QOFS(YOLD(I),YPOLD(I),Y(I),YP(I),HOLD,HOLD+H)
        END DO
      ENDIF
      DZ = DZ - Z0
      CORDIS = DNRM2(NP1,DZ,1)
C
C SET UP VARIABLES FOR NEXT CALL.
C
      YOLD = Y
      Y = Z0
      YPOLD = YP
      YP = T
C
C UPDATE ARCLENGTH   S = S + ||Y-YOLD||.
C
      HTEMP = HOLD
      Z0 = Z0 - YOLD
      HOLD = DNRM2(NP1,Z0,1)
      S = S+HOLD
C
C COMPUTE IDEAL ERROR FOR STEP SIZE ESTIMATION.
C
      IF (ITCNT .LE. 1) THEN
          THETA = 8.0_DP
      ELSE IF (ITCNT .EQ. 4) THEN
          THETA = ONE
      ELSE
          OMEGA=XSTEP/CORDIS
          IF (ITCNT .LT. 4) THEN
            LK = 4*ITCNT-7
            IF (OMEGA .GE. WRGE(LK)) THEN
              THETA = ONE
            ELSE IF (OMEGA .GE. WRGE(LK+1)) THEN
              THETA = ACOF(LK) + ACOF(LK+1)*LOG(OMEGA)
            ELSE IF (OMEGA .GE. WRGE(LK+2)) THEN
              THETA = ACOF(LK+2) + ACOF(LK+3)*LOG(OMEGA)
            ELSE 
              THETA = 8.0_DP
            END IF
          ELSE IF (ITCNT .GE. 7) THEN
            THETA = 0.125_DP
          ELSE
            LK = 4*ITCNT - 16
            IF (OMEGA .GT. WRGE(LK)) THEN
              LST = 2*ITCNT - 1
              THETA = ACOF(LST) + ACOF(LST+1)*LOG(OMEGA)
            ELSE
              THETA = 0.125_DP
            END IF
          END IF
      END IF
      IDLERR=THETA*IDLERR
C
C IDLERR SHOULD BE NO BIGGER THAN 1/2 PREVIOUS STEP.
C
      IDLERR = MIN(0.5_DP*HOLD,IDLERR)
C
C COMPUTE OPTIMAL STEP SIZE. 
C   WK = APPROXIMATE CURVATURE = 2*SIN(ALPHA/2)/HOLD  WHERE 
C        ALPHA = ARCCOS(YP*YPOLD).
C   GAMMA = EXPECTED CURVATURE FOR NEXT STEP, COMPUTED BY 
C        EXTRAPOLATING FROM CURRENT CURVATURE  WK, AND LAST 
C        CURVATURE  WKOLD.  GAMMA  IS FURTHER REQUIRED TO BE 
C        POSITIVE.
C
      IF (.NOT. START) WKOLD = WK
      WK = 2*ABS(SIN(0.5_DP*ALPHA))/HOLD
      IF (START) THEN
        GAMMA = WK
      ELSE 
        GAMMA = WK + HOLD/(HOLD+HTEMP)*(WK-WKOLD)
      END IF
      GAMMA = MAX(GAMMA, 0.01_DP)
      H = SQRT(2*IDLERR/GAMMA)
C
C     ENFORCE RESTRICTIONS ON STEP SIZE SO AS TO ENSURE STABILITY.
C        HMIN <= H <= HMAX, BMIN*HOLD <= H <= BMAX*HOLD.
C
      H = MIN(MAX(SSPAR(1),SSPAR(3)*HOLD,H),SSPAR(4)*HOLD,SSPAR(2))
      IF (FAILED) H = MIN(HFAIL,H)
      START = .FALSE.
C
C ***** END OF MOP UP SECTION. *****
C
      END SUBROUTINE STEPQS
C
      END MODULE HOMPACK_QS