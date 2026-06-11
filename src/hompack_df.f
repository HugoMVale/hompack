      MODULE HOMPACK_DF
C
      USE HOMPACK_KINDS, ONLY: DP
      IMPLICIT NONE
C
      CONTAINS
C
      SUBROUTINE FIXPDF(N,Y,IFLAG,ARCTOL,EPS,TRACE,A,NDIMA,
     &   NFE,ARCLEN)
C
C Subroutine  FIXPDF  finds a fixed point or zero of the
C N-dimensional vector function F(X), or tracks a zero curve
C of a general homotopy map RHO(A,LAMBDA,X).  For the fixed 
C point problem F(X) is assumed to be a C2 map of some ball 
C into itself.  The equation  X = F(X)  is solved by
C following the zero curve of the homotopy map
C
C  LAMBDA*(X - F(X)) + (1 - LAMBDA)*(X - A)  ,
C
C starting from LAMBDA = 0, X = A.  The curve is parameterized
C by arc length S, and is followed by solving the ordinary
C differential equation  D(HOMOTOPY MAP)/DS = 0  for
C Y(S) = (LAMBDA(S), X(S)).
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
C For the curve tracking problem RHO(A,LAMBDA,X) is assumed to
C be a C2 map from E**M X [0,1) X E**N into E**N, which for
C almost all parameter vectors A in some nonempty open subset
C of E**M satisfies
C
C  rank [D RHO(A,LAMBDA,X)/D LAMBDA , D RHO(A,LAMBDA,X)/DX] = N
C
C for all points (LAMBDA,X) such that RHO(A,LAMBDA,X)=0.  It is
C further assumed that
C
C           rank [ D RHO(A,0,X0)/DX ] = N  .
C
C With A fixed, the zero curve of RHO(A,LAMBDA,X) emanating
C from  LAMBDA = 0, X = X0  is tracked until  LAMBDA = 1  by
C solving the ordinary differential equation
C D RHO(A,LAMBDA(S),X(S))/DS = 0  for  Y(S) = (LAMBDA(S), X(S)),
C where S is arc length along the zero curve.  Also the homotopy
C map RHO(A,LAMBDA,X) is assumed to be constructed such that
C
C              D LAMBDA(0)/DS > 0  .
C
C This code is based on the algorithm in L. T. Watson, A
C globally convergent algorithm for computing fixed points of
C C2 maps, Appl. Math. Comput., 5 (1979) 297-311.
C
C
C For the fixed point and zero finding problems, the user
C must supply a subroutine  F(X,V)  which evaluates F(X) at X
C and returns the vector F(X) in V, and a subroutine  FJAC(X,V,K)
C which returns in V the Kth column of the Jacobian matrix of 
C F(X) evaluated at X.  For the curve tracking problem, the user must
C supply a subroutine  RHOA(V,LAMBDA,X)  which given 
C (LAMBDA,X) returns a parameter vector A in V such that 
C RHO(A,LAMBDA,X)=0, and a subroutine  RHOJAC(A,LAMBDA,X,V,K)
C which returns in V the Kth column of the N X (N+1) Jacobian 
C matrix [D RHO/D LAMBDA, D RHO/DX] evaluated at (A,LAMBDA,X).
C Whichever of the routines  F,  FJAC,  RHOA,  RHOJAC  are required
C should be supplied as external subroutines, conforming with the
C interfaces in the module  HOMOTOPY.
C FIXPDF  directly or indirectly uses the subroutines
C   F (or  RHOA ),  FJAC (or  RHOJAC ),  FODE,   ROOT,
C   SINTRP,  STEPS,  the LAPACK routine  DGEQPF,  auxiliary LAPACK 
C routines, and the BLAS functions  DCOPY,  DDOT,  DGEMV,  DGER,  
C   DNRM2,  DSCAL,  DSWAP,  IDAMAX.  
C The module  REAL_PRECISION  specifies 64-bit real arithmetic, which
C the user may want to change.
C
C ***Warning:  this subroutine is generally more robust than  FIXPNF
C and  FIXPQF, but may be slower than those subroutines by a
C factor of two.
C
C
C ON INPUT:
C
C N  is the dimension of X, F(X), and RHO(A,LAMBDA,X).
C
C Y  is an array of length  N + 1.  (Y(2),...,Y(N+1)) = A  is the
C    starting point for the zero curve for the fixed point and 
C    zero finding problems.  (Y(2),...,Y(N+1)) = X0  for the curve
C    tracking problem.
C
C IFLAG  can be -2, -1, 0, 2, or 3.  IFLAG  should be 0 on the 
C    first call to  FIXPDF  for the problem  X=F(X), -1 for the
C    problem  F(X)=0, and -2 for the problem  RHO(A,LAMBDA,X)=0.
C    In certain situations  IFLAG  is set to 2 or 3 by  FIXPDF,
C    and  FIXPDF  can be called again without changing  IFLAG.
C
C ARCTOL  is the local error allowed the ODE solver when
C    following the zero curve.  If  ARCTOL .LE. 0.0  on input
C    it is reset to  .5*SQRT(EPS).  Normally  ARCTOL  should
C    be considerably larger than  EPS.
C
C EPS  is the local error allowed the ODE solver when very
C    near the fixed point(zero).  EPS  is approximately the
C    mixed absolute and relative error in the computed fixed 
C    point(zero).
C
C TRACE  is an integer specifying the logical I/O unit for
C    intermediate output.  If  TRACE .GT. 0  the points computed on
C    the zero curve are written to I/O unit  TRACE .
C
C A(1:NDIMA) contains the parameter vector  A.  For the fixed point
C    and zero finding problems, A  need not be initialized by the
C    user, and is assumed to have length  N.  For the curve
C    tracking problem, A  has length  NDIMA  and must be initialized
C    by the user.
C
C NDIMA  is the dimension of  A, used for the curve tracking problem,
C    and must be N for the fixed point and zero finding problems.
C
C Y, ARCTOL, EPS, ARCLEN, NFE, and IFLAG should all be
C variables in the calling program.
C
C
C ON OUTPUT:
C
C N  and  TRACE  are unchanged.
C
C Y(1) = LAMBDA, (Y(2),...,Y(N+1)) = X, and Y is an approximate
C    zero of the homotopy map.  Normally LAMBDA = 1 and X is a
C    fixed point(zero) of F(X).  In abnormal situations LAMBDA
C    may only be near 1 and X is near a fixed point(zero).
C
C IFLAG =
C  -2   causes  FIXPDF  to initialize everything for the problem
C       RHO(A,LAMBDA,X) = 0 (use on first call).
C
C  -1   causes  FIXPDF  to initialize everything for the problem
C       F(X) = 0 (use on first call).
C
C   0   causes  FIXPDF  to initialize everything for the problem
C       X = F(X) (use on first call).
C
C   1   Normal return.
C
C   2   Specified error tolerance cannot be met.  EPS has been
C       increased to a suitable value.  To continue, just call
C       FIXPDF  again without changing any parameters.
C
C   3   STEPS  has been called 1000 times.  To continue, call
C       FIXPDF  again without changing any parameters.
C
C   4   Jacobian matrix does not have full rank.  The algorithm
C       has failed (the zero curve of the homotopy map cannot be
C       followed any further).
C
C   5   EPS  (or  ARCTOL) is too large.  The problem should be
C       restarted by calling  FIXPDF  with a smaller  EPS  (or
C       ARCTOL) and  IFLAG = 0 (-1, -2).
C
C   6   I - DF(X)  is nearly singular at the fixed point (DF(X) is
C       nearly singular at the zero, or  D RHO(A,LAMBDA,X)/DX  is
C       nearly singular at  LAMBDA = 1 ).  Answer may not be
C       accurate.
C
C   7   Illegal input parameters, a fatal error.
C
C   8   Memory allocation error, fatal.
C
C ARCTOL = EPS after a normal return (IFLAG = 1).
C
C EPS  is unchanged after a normal return (IFLAG = 1).  It is
C    increased to an appropriate value on the return IFLAG = 2.
C
C A  will (normally) have been modified.
C
C NFE  is the number of function evaluations (= number of
C    Jacobian matrix evaluations).
C
C ARCLEN  is the length of the path followed.
C
C
C Automatic work arrays:
C
C YP(1:N+1) is a work array containing the current tangent
C    vector to the zero curve.
C
C YPOLD(1:N+1) is a work array containing the previous tangent
C    vector to the zero curve.
C
C QR(1:N,1:N+1), ALPHA(1:3*N+3), TZ(1:N+1), and PIVOT(1:N+1) are 
C    all work arrays used by  FODE  to calculate the tangent
C    vector YP.
C
C WT(1:N+1), PHI(1:N+1,1:16), and P(1:N+1) are all work arrays
C    used by the ODE subroutine  STEPS  .
C
      USE HOMPACK_CORE_LEGACY, ONLY: F, RHOA
      USE HOMPACK_CORE_LEGACY, ONLY: SINTRP, ROOT
      IMPLICIT NONE
C
      INTEGER, INTENT(IN)::N,NDIMA,TRACE
      REAL (dp), DIMENSION(:), INTENT(IN OUT)::A,Y
      INTEGER, INTENT(IN OUT)::IFLAG
      REAL (dp), INTENT(IN OUT)::ARCTOL,EPS
      INTEGER, INTENT(OUT)::NFE
      REAL (dp), INTENT(OUT)::ARCLEN
C
C LOCAL VARIABLES.
      REAL (dp), SAVE:: CURSW,CURTOL,EPSSTP,EPST,H,HOLD,
     &  S,S99,SA,SB,SOUT,SQNP1,XOLD,Y1SOUT
      INTEGER, SAVE:: IFLAGC,ITER,IVC,JW,K,KGI,KOLD,
     &  KSTEPS,LCODE,LIMIT,NFEC,NP1
      LOGICAL, SAVE:: CRASH,START,ST99
C
C *****  ARRAY DECLARATIONS.  *****
C
C ARRAYS NEEDED BY THE ODE SUBROUTINE  STEPS .
      REAL (dp), ALLOCATABLE, SAVE:: P(:),PHI(:,:),WT(:),YP(:)
      REAL (dp), SAVE:: ALPHAS(12),G(13),GI(11),W(12)
      INTEGER, SAVE:: IV(10)
C
C ARRAYS NEEDED BY  FIXPDF , FODE , AND LAPACK ROUTINES.
      REAL (dp), DIMENSION(:), ALLOCATABLE, SAVE:: YPOLD
      REAL (dp):: ALPHA(3*N+3),AOLD(NDIMA),QR(N,N+1),TZ(N+1)
      INTEGER:: PIVOT(N+1)
C
C *****  END OF DIMENSIONAL INFORMATION.  *****
C
C LIMITD  IS AN UPPER BOUND ON THE NUMBER OF STEPS.  IT MAY BE
C CHANGED BY CHANGING THE FOLLOWING PARAMETER STATEMENT:
      INTEGER, PARAMETER:: LIMITD=1000
C
C :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :
      IF (N .LE. 0  .OR.  EPS .LE. 0.0  .OR.  N+1 .NE. SIZE(Y)
     &  .OR.  NDIMA .NE. SIZE(A)  .OR.
     &  ((IFLAG .EQ. -1  .OR.  IFLAG .EQ. 0) .AND.  N .NE. SIZE(A)))
     &  IFLAG=7
      IF (IFLAG .GE. -2  .AND.  IFLAG .LE. 0) GO TO 10
      IF (IFLAG .EQ. 2) GO TO 35
      IF (IFLAG .EQ. 3) GO TO 30
C ONLY VALID INPUT FOR  IFLAG  IS -2, -1, 0, 2, 3.
      IFLAG=7
      RETURN
C
C *****  INITIALIZATION BLOCK.  *****
C
10    ARCLEN=0.0
      S=0.0
      IF (ARCTOL .LE. 0.0) ARCTOL=.5*SQRT(EPS)
      NFEC=0
      IFLAGC=IFLAG
      NP1=N+1
      SQNP1=SQRT(DBLE(NP1))
      IF (ALLOCATED(P)) DEALLOCATE(P)
      IF (ALLOCATED(PHI)) DEALLOCATE(PHI)
      IF (ALLOCATED(WT)) DEALLOCATE(WT)
      IF (ALLOCATED(YP)) DEALLOCATE(YP)
      IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
      ALLOCATE(P(NP1),PHI(NP1,16),WT(NP1),YP(NP1),YPOLD(NP1),
     &  STAT=JW)
      IF (JW /= 0) THEN
        IFLAG=8
        RETURN
      END IF
C
C SWITCH FROM THE TOLERANCE  ARCTOL  TO THE (FINER) TOLERANCE  EPS  IF
C THE CURVATURE OF ANY COMPONENT OF  Y  EXCEEDS  CURSW.
C
      CURSW=10.0
C
      ST99=.FALSE.
      START=.TRUE.
      CRASH=.FALSE.
      HOLD=1.0
      H=.1
      EPSSTP=ARCTOL
      KSTEPS=0
C SET INITIAL CONDITIONS FOR ORDINARY DIFFERENTIAL EQUATION.
      YPOLD(1)=1.0
      YP(1)=1.0
      Y(1)=0.0
      YPOLD(2:NP1)=0.0
      YP(2:NP1)=0.0
C LOAD  A  FOR THE FIXED POINT AND ZERO FINDING PROBLEMS.
      IF (IFLAGC .GE. -1) THEN
        A=Y(2:NP1)
      ENDIF
30    LIMIT=LIMITD
C
C *****  END OF INITIALIZATION BLOCK.  *****
C
35    MAIN_LOOP: DO ITER=1,LIMIT  ! *****  MAIN LOOP.  *****
      IF (Y(1) .LT. 0.0) THEN
        ARCLEN=ARCLEN+S
        IFLAG=5
        CALL CLEANUP ; RETURN
      ENDIF
      IF (S .LE. 7.0*SQNP1) GO TO 80
C ARC LENGTH IS GETTING TOO LONG, THE PROBLEM WILL BE
C RESTARTED WITH A DIFFERENT  A  VECTOR.
      ARCLEN=ARCLEN+S
      S=0.0
60    START=.TRUE.
      CRASH=.FALSE.
C COMPUTE A NEW  A  VECTOR.
      IF (IFLAGC .EQ. -2) THEN
        AOLD=A
        CALL RHOA(A,Y(1),Y(2:NP1))
C IF NEW AND OLD  A  DIFFER BY TOO MUCH, TRACKING SHOULD NOT CONTINUE.
        IF (ANY(ABS(A-AOLD) .GT. 1.0+ABS(AOLD))) THEN
          ARCLEN=ARCLEN+S
          IFLAG=5
          CALL CLEANUP ; RETURN
        ENDIF
      ELSE
        CALL F(Y(2:NP1),YP(1:N))
        AOLD=A
        IF (IFLAGC .EQ. -1) THEN
          A=Y(1)*YP(1:N)/(1.0 - Y(1)) + Y(2:NP1)
        ELSE
          A=(Y(2:NP1) - Y(1)*YP(1:N))/(1.0 - Y(1))
        ENDIF
C IF NEW AND OLD  A  DIFFER BY TOO MUCH, TRACKING SHOULD NOT CONTINUE.
        IF (ANY(ABS(A-AOLD) .GT. 1.0+ABS(AOLD))) THEN
          ARCLEN=ARCLEN+S
          IFLAG=5
          CALL CLEANUP ; RETURN
        ENDIF
      ENDIF
      GO TO 100
80    IF (Y(1) .LE. .99  .OR. ST99) GO TO 100
C WHEN LAMBDA REACHES .99, THE PROBLEM WILL BE RESTARTED WITH
C A NEW  A  VECTOR.
90    ST99=.TRUE.
      EPSSTP=EPS
      ARCTOL=EPS
      GO TO 60
C
C SET DIFFERENT ERROR TOLERANCE FOR HIGH CURVATURE COMPONENTS OF THE
C TRAJECTORY Y(S).
100   CURTOL=CURSW*HOLD
      EPST=EPS/EPSSTP
      WHERE (ABS(YP-YPOLD) .LE. CURTOL)
        WT=(ABS(Y)+1.0)
      ELSEWHERE
        WT=(ABS(Y)+1.0)*EPST
      END WHERE
C
C TAKE A STEP ALONG THE CURVE.
      CALL STEPS(FODE,NP1,Y,S,H,EPSSTP,WT,START,HOLD,K,KOLD,CRASH,
     &     PHI,P,YP,ALPHAS,W,G,KSTEPS,XOLD,IVC,IV,KGI,GI,
     &     YPOLD,A,QR,ALPHA,TZ,PIVOT,NFEC,IFLAGC)
C PRINT LATEST POINT ON CURVE IF REQUESTED.
      IF (TRACE .GT. 0) THEN
        WRITE (TRACE,117) ITER,NFEC,S,Y(1),(Y(JW),JW=2,NP1)
117     FORMAT(/' STEP',I5,3X,'NFE =',I5,3X,'ARC LENGTH =',F9.4,3X,
     &  'LAMBDA =',F7.4,5X,'X VECTOR:'/(1X,6ES12.4))
      ENDIF
      NFE=NFEC
C CHECK IF THE STEP WAS SUCCESSFUL.
      IF (IFLAGC .EQ. 4) THEN
        ARCLEN=ARCLEN+S
        IFLAG=4
        CALL CLEANUP ; RETURN
      ENDIF
      IF (CRASH) THEN
C RETURN CODE FOR ERROR TOLERANCE TOO SMALL.
        IFLAG=2
C CHANGE ERROR TOLERANCES.
        EPS=EPSSTP
        IF (ARCTOL .LT. EPS) ARCTOL=EPS
C CHANGE LIMIT ON NUMBER OF ITERATIONS.
        LIMIT=LIMIT-ITER
        RETURN
      ENDIF
C
      IF (Y(1) .GE. 1.0) THEN
        IF (ST99) GO TO 160
C
C IF LAMBDA .GE. 1.0 BUT THE PROBLEM HAS NOT BEEN RESTARTED
C WITH A NEW  A  VECTOR, BACK UP AND RESTART.
C
        S99=S-.5*HOLD
C GET AN APPROXIMATE ZERO Y(S) WITH  Y(1)=LAMBDA .LT. 1.0  .
135     CALL SINTRP(S,Y,S99,WT,YP,NP1,KOLD,PHI,IVC,IV,KGI,GI,
     &     ALPHAS,G,W,XOLD,P)
        IF (WT(1) .LT. 1.0) GO TO 140
        S99=.5*(S-HOLD+S99)
        GO TO 135
C
140     Y=WT
        YPOLD=YP
        S=S99
        GO TO 90
      ENDIF
C
      END DO MAIN_LOOP  ! *****  END OF MAIN LOOP.  *****
C
C LAMBDA HAS NOT REACHED 1 IN 1000 STEPS.
      IFLAG=3
      RETURN
C
C
C USE INVERSE INTERPOLATION TO GET THE ANSWER AT LAMBDA = 1.0 .
C
160   SA=S-HOLD
      SB=S
      LCODE=1
170   CALL ROOT(SOUT,Y1SOUT,SA,SB,EPS,EPS,LCODE)
C ROOT  FINDS S SUCH THAT Y(1)(S) = LAMBDA = 1 .
      IF (LCODE .GT. 0) GO TO 190
      CALL SINTRP(S,Y,SOUT,WT,YP,NP1,KOLD,PHI,IVC,IV,KGI,GI,
     &     ALPHAS,G,W,XOLD,P)
      Y1SOUT=WT(1)-1.0
      GO TO 170
190   IFLAG=1
C SET IFLAG = 6 IF  ROOT  COULD NOT GET  LAMBDA = 1.0  .
      IF (LCODE .GT. 2) IFLAG=6
      ARCLEN=ARCLEN+SA
C LAMBDA(SA) = 1.0 .
      CALL SINTRP(S,Y,SA,WT,YP,NP1,KOLD,PHI,IVC,IV,KGI,GI,
     &     ALPHAS,G,W,XOLD,P)
C
      Y=WT
      CALL CLEANUP ; RETURN
C
      CONTAINS
        SUBROUTINE CLEANUP
        IF (ALLOCATED(P)) DEALLOCATE(P)
        IF (ALLOCATED(PHI)) DEALLOCATE(PHI)
        IF (ALLOCATED(WT)) DEALLOCATE(WT)
        IF (ALLOCATED(YP)) DEALLOCATE(YP)
        IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
        END SUBROUTINE CLEANUP
      END SUBROUTINE FIXPDF
C
C
      SUBROUTINE STEPS(FODE,NEQN,Y,X,H,EPS,WT,START,HOLD,K,KOLD,CRASH,
     &  PHI,P,YP,ALPHA,W,G,KSTEPS,XOLD,IVC,IV,KGI,GI, 
     &  FPWA1,FPWA2,FPWA3,FPWA4,FPWA5,IFPWA1,IFPC1,IFPC2)
C 
C   FORTRAN 90 MODIFICATION OF THE SUBROUTINE  STEP  WRITTEN BY
C   L. F. SHAMPINE AND M. K. GORDON
C 
C   ABSTRACT
C 
C   SUBROUTINE  STEPS  IS NORMALLY USED INDIRECTLY THROUGH SUBROUTINE 
C   DEABM .  BECAUSE  DEABM  SUFFICES FOR MOST PROBLEMS AND IS MUCH 
C   EASIER TO USE, USING IT SHOULD BE CONSIDERED BEFORE USING  STEPS
C   ALONE.
C 
C   SUBROUTINE STEPS INTEGRATES A SYSTEM OF  NEQN  FIRST ORDER ORDINARY 
C   DIFFERENTIAL EQUATIONS ONE STEP, NORMALLY FROM X TO X+H, USING A
C   MODIFIED DIVIDED DIFFERENCE FORM OF THE ADAMS PECE FORMULAS.  LOCAL 
C   EXTRAPOLATION IS USED TO IMPROVE ABSOLUTE STABILITY AND ACCURACY. 
C   THE CODE ADJUSTS ITS ORDER AND STEP SIZE TO CONTROL THE LOCAL ERROR 
C   PER UNIT STEP IN A GENERALIZED SENSE.  SPECIAL DEVICES ARE INCLUDED 
C   TO CONTROL ROUNDOFF ERROR AND TO DETECT WHEN THE USER IS REQUESTING 
C   TOO MUCH ACCURACY.
C 
C   THIS CODE IS COMPLETELY EXPLAINED AND DOCUMENTED IN THE TEXT, 
C   COMPUTER SOLUTION OF ORDINARY DIFFERENTIAL EQUATIONS, THE INITIAL 
C   VALUE PROBLEM  BY L. F. SHAMPINE AND M. K. GORDON.
C   FURTHER DETAILS ON USE OF THIS CODE ARE AVAILABLE IN *SOLVING 
C   ORDINARY DIFFERENTIAL EQUATIONS WITH ODE, STEP, AND INTRP*, 
C   BY L. F. SHAMPINE AND M. K. GORDON, SLA-73-1060.
C 
C 
C   THE PARAMETERS REPRESENT -- 
C      FODE -- SUBROUTINE TO EVALUATE DERIVATIVES
C      NEQN -- NUMBER OF EQUATIONS TO BE INTEGRATED 
C      Y(*) -- SOLUTION VECTOR AT X 
C      X -- INDEPENDENT VARIABLE
C      H -- APPROPRIATE STEP SIZE FOR NEXT STEP.  NORMALLY DETERMINED BY
C           CODE
C      EPS -- LOCAL ERROR TOLERANCE 
C      WT(*) -- VECTOR OF WEIGHTS FOR ERROR CRITERION 
C      START -- LOGICAL VARIABLE SET .TRUE. FOR FIRST STEP,  .FALSE.
C           OTHERWISE 
C      HOLD -- STEP SIZE USED FOR LAST SUCCESSFUL STEP
C      K -- APPROPRIATE ORDER FOR NEXT STEP (DETERMINED BY CODE)
C      KOLD -- ORDER USED FOR LAST SUCCESSFUL STEP
C      CRASH -- LOGICAL VARIABLE SET .TRUE. WHEN NO STEP CAN BE TAKEN,
C           .FALSE. OTHERWISE.
C      YP(*) -- DERIVATIVE OF SOLUTION VECTOR AT  X  AFTER SUCCESSFUL 
C           STEP
C      KSTEPS -- COUNTER ON ATTEMPTED STEPS 
C
C   THE VARIABLES X,XOLD,KOLD,KGI AND IVC AND THE ARRAYS Y,PHI,ALPHA,G, 
C   W,P,IV AND GI ARE REQUIRED FOR THE INTERPOLATION SUBROUTINE SINTRP. 
C   THE ARRAYS FPWA* AND IFPWA1 AND INTEGER CONSTANTS IFPC* ARE
C   WORKING STORAGE PASSED DIRECTLY THROUGH TO  FODE.
C 
C   INPUT TO STEPS
C 
C      FIRST CALL --
C 
C   THE USER MUST PROVIDE STORAGE IN HIS CALLING PROGRAM FOR ALL ARRAYS 
C   IN THE CALL LIST, NAMELY
C 
C     DIMENSION Y(NEQN),WT(NEQN),PHI(NEQN,16),P(NEQN),YP(NEQN), 
C    1  ALPHA(12),W(12),G(13),GI(11),IV(10),   FPWA1(NEQN),
C    2  FPWA2(NEQN-1),FPWA3(NEQN-1,NEQN),FPWA4(3*NEQN),
C    3  FPWA5(NEQN),IFPWA1(NEQN)
C                              --                --    **NOTE** 
C 
C   THE USER MUST ALSO DECLARE  START  AND  CRASH 
C   LOGICAL VARIABLES AND  FODE  AN EXTERNAL SUBROUTINE, SUPPLY THE
C   SUBROUTINE  FODE(X,Y,YP,FPWA1,FPWA2,FPWA3,FPWA4,FPWA5,IFPWA1,IFPC1,
C                 NEQN-1,IFPC2) TO EVALUATE
C      DY(I)/DX = YP(I) = F(X,Y(1),Y(2),...,Y(NEQN))
C   AND INITIALIZE ONLY THE FOLLOWING PARAMETERS. 
C      NEQN -- NUMBER OF EQUATIONS TO BE INTEGRATED 
C      Y(*) -- VECTOR OF INITIAL VALUES OF DEPENDENT VARIABLES
C      X -- INITIAL VALUE OF THE INDEPENDENT VARIABLE 
C      H -- NOMINAL STEP SIZE INDICATING DIRECTION OF INTEGRATION 
C           AND MAXIMUM SIZE OF STEP.  MUST BE VARIABLE 
C      EPS -- LOCAL ERROR TOLERANCE PER STEP.  MUST BE VARIABLE 
C      WT(*) -- VECTOR OF NON-ZERO WEIGHTS FOR ERROR CRITERION
C      START -- .TRUE.
C      KSTEPS -- SET KSTEPS TO ZERO 
C   DEFINE U TO BE THE MACHINE UNIT ROUNDOFF QUANTITY BY CALLING
C   THE INTRINSIC FUNCTION  EPSILON(1.0_dp), OR BY 
C   COMPUTING U SO THAT U IS THE SMALLEST POSITIVE NUMBER SUCH
C   THAT 1.0+U .GT. 1.0.
C 
C   STEPS  REQUIRES THAT THE L2 NORM OF THE VECTOR WITH COMPONENTS
C   LOCAL ERROR(L)/WT(L)  BE LESS THAN  EPS  FOR A SUCCESSFUL STEP.  THE
C   ARRAY  WT  ALLOWS THE USER TO SPECIFY AN ERROR TEST APPROPRIATE 
C   FOR HIS PROBLEM.  FOR EXAMPLE,
C      WT(L) = 1.0  SPECIFIES ABSOLUTE ERROR, 
C            = ABS(Y(L))  ERROR RELATIVE TO THE MOST RECENT VALUE OF THE
C                 L-TH COMPONENT OF THE SOLUTION, 
C            = ABS(YP(L))  ERROR RELATIVE TO THE MOST RECENT VALUE OF 
C                 THE L-TH COMPONENT OF THE DERIVATIVE, 
C            = MAX(WT(L),ABS(Y(L)))  ERROR RELATIVE TO THE LARGEST
C                 MAGNITUDE OF L-TH COMPONENT OBTAINED SO FAR,
C            = ABS(Y(L))*RELERR/EPS + ABSERR/EPS  SPECIFIES A MIXED 
C                 RELATIVE-ABSOLUTE TEST WHERE  RELERR  IS RELATIVE 
C                 ERROR,  ABSERR  IS ABSOLUTE ERROR AND  EPS =
C                 MAX(RELERR,ABSERR) .
C 
C      SUBSEQUENT CALLS --
C 
C   SUBROUTINE  STEPS  IS DESIGNED SO THAT ALL INFORMATION NEEDED TO
C   CONTINUE THE INTEGRATION, INCLUDING THE STEP SIZE  H  AND THE ORDER 
C   K , IS RETURNED WITH EACH STEP.  WITH THE EXCEPTION OF THE STEP 
C   SIZE, THE ERROR TOLERANCE, AND THE WEIGHTS, NONE OF THE PARAMETERS
C   SHOULD BE ALTERED.  THE ARRAY  WT  MUST BE UPDATED AFTER EACH STEP
C   TO MAINTAIN RELATIVE ERROR TESTS LIKE THOSE ABOVE.  NORMALLY THE
C   INTEGRATION IS CONTINUED JUST BEYOND THE DESIRED ENDPOINT AND THE 
C   SOLUTION INTERPOLATED THERE WITH SUBROUTINE  SINTRP .  IF IT IS 
C   IMPOSSIBLE TO INTEGRATE BEYOND THE ENDPOINT, THE STEP SIZE MAY BE 
C   REDUCED TO HIT THE ENDPOINT SINCE THE CODE WILL NOT TAKE A STEP 
C   LARGER THAN THE  H  INPUT.  CHANGING THE DIRECTION OF INTEGRATION,
C   I.E., THE SIGN OF  H , REQUIRES THE USER SET  START = .TRUE. BEFORE 
C   CALLING  STEPS  AGAIN.  THIS IS THE ONLY SITUATION IN WHICH  START
C   SHOULD BE ALTERED.
C 
C   OUTPUT FROM STEPS 
C 
C      SUCCESSFUL STEP -- 
C 
C   THE SUBROUTINE RETURNS AFTER EACH SUCCESSFUL STEP WITH  START  AND
C   CRASH  SET .FALSE. .  X  REPRESENTS THE INDEPENDENT VARIABLE
C   ADVANCED ONE STEP OF LENGTH  HOLD  FROM ITS VALUE ON INPUT AND  Y 
C   THE SOLUTION VECTOR AT THE NEW VALUE OF  X .  ALL OTHER PARAMETERS
C   REPRESENT INFORMATION CORRESPONDING TO THE NEW  X  NEEDED TO
C   CONTINUE THE INTEGRATION. 
C 
C      UNSUCCESSFUL STEP -- 
C 
C   WHEN THE ERROR TOLERANCE IS TOO SMALL FOR THE MACHINE PRECISION,
C   THE SUBROUTINE RETURNS WITHOUT TAKING A STEP AND  CRASH = .TRUE. .
C   AN APPROPRIATE STEP SIZE AND ERROR TOLERANCE FOR CONTINUING ARE 
C   ESTIMATED AND ALL OTHER INFORMATION IS RESTORED AS UPON INPUT 
C   BEFORE RETURNING.  TO CONTINUE WITH THE LARGER TOLERANCE, THE USER
C   JUST CALLS THE CODE AGAIN.  A RESTART IS NEITHER REQUIRED NOR 
C   DESIRABLE.
C***REFERENCES  SHAMPINE L.F., GORDON M.K., *SOLVING ORDINARY 
C                 DIFFERENTIAL EQUATIONS WITH ODE, STEP, AND INTRP*,
C                 SLA-73-1060, SANDIA LABORATORIES, 1973.
C
      USE HOMPACK_KINDS, ONLY: ZERO, ONE
      IMPLICIT NONE
C 
      REAL(DP):: ABSH,ALPHA,BETA,EPS,ERK,ERKM1,ERKM2,
     &  ERKP1,ERR,FOURU,FPWA1,FPWA2,FPWA3,FPWA4,FPWA5,G,GI,GSTR,H,
     &  HNEW,HOLD,P,P5EPS,PHI,PSI,R,REALI,REALNS,RHO,ROUND,SIG,
     &  SUM,TAU,TEMP1,TEMP2,TEMP3,TEMP4,TEMP5,TEMP6,TWO,TWOU,V,
     &  W,WT,X,XOLD,Y,YP
      INTEGER I,IFAIL,IFPC1,IFPC2,IFPWA1,IM1,IP1,IQ,IV,IVC,
     &  J,JV,K,KGI,KM1,KM2,KNEW,KOLD,KP1,KP2,KPREV,KSTEPS,
     &  L,LIMIT1,LIMIT2,NEQN,NS,NSM2,NSP1,NSP2
      LOGICAL START,CRASH,PHASE1,NORND
C
      DIMENSION Y(:),WT(NEQN),PHI(NEQN,16),P(NEQN),YP(NEQN),PSI(12), 
     &  ALPHA(12),BETA(12),SIG(13),V(12),W(12),G(13),GI(11),IV(10), 
     &  FPWA1(NEQN),FPWA2(:),FPWA3(NEQN-1,NEQN),FPWA4(3*NEQN),
     &  FPWA5(NEQN),IFPWA1(NEQN)
      DIMENSION TWO(13),GSTR(13)
C
      INTERFACE
        SUBROUTINE FODE(S,Y,YP,YPOLD,A,QR,ALPHA,TZ,PIVOT,NFE,N,IFLAG)
          IMPORT :: DP
          REAL(DP):: S
          INTEGER:: IFLAG,N,NFE
          REAL(DP):: A(:),Y(:),YP(N+1),YPOLD(N+1)
          REAL(DP):: ALPHA(3*N+3),QR(N,N+1),TZ(N+1)
          INTEGER, DIMENSION(N+1):: PIVOT
        END SUBROUTINE FODE
      END INTERFACE
C
C   ALL LOCAL VARIABLES ARE SAVED, RATHER THAN PASSED, IN THIS
C   SPECIALIZED VERSION OF STEPS.
C
      SAVE
C
      DATA TWO /2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp, 32.0_dp, 64.0_dp,
     &  128.0_dp, 256.0_dp, 512.0_dp, 1024.0_dp, 2048.0_dp,
     &  4096.0_dp, 8192.0_dp/
      DATA GSTR /0.500_dp, 0.0833_dp, 0.0417_dp, 0.0264_dp,
     &  0.0188_dp, 0.0143_dp, 0.0114_dp, 0.00936_dp, 0.00789_dp,
     &  0.00679_dp, 0.00592_dp, 0.00524_dp, 0.00468_dp/
C 
C 
C       ***     BEGIN BLOCK 0     *** 
C   CHECK IF STEP SIZE OR ERROR TOLERANCE IS TOO SMALL FOR MACHINE
C   PRECISION.  IF FIRST STEP, INITIALIZE PHI ARRAY AND ESTIMATE A
C   STARTING STEP SIZE. 
C                   *** 
C 
C   IF STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE
C 
C***FIRST EXECUTABLE STATEMENT
      TWOU = 2*EPSILON(ONE)
      FOURU = TWOU + TWOU
      CRASH = .TRUE.
      IF(ABS(H) .GE. FOURU*ABS(X)) GO TO 5
      H = SIGN(FOURU*ABS(X),H)
      RETURN
 5    P5EPS = 0.5_DP*EPS 
C 
C   IF ERROR TOLERANCE IS TOO SMALL, INCREASE IT TO AN ACCEPTABLE VALUE 
C 
      ROUND = ZERO 
      DO L = 1,NEQN
        ROUND = ROUND + (Y(L)/WT(L))**2 
      END DO
      ROUND = TWOU*SQRT(ROUND)
      IF(P5EPS .GE. ROUND) GO TO 15 
      EPS = 2*ROUND*(ONE + FOURU) 
      RETURN
 15   CRASH = .FALSE. 
      G(1) = ONE
      G(2) = 0.5_DP
      SIG(1) = ONE
      IF(.NOT.START) GO TO 99 
C 
C   INITIALIZE.  COMPUTE APPROPRIATE STEP SIZE FOR FIRST STEP 
C 
      CALL FODE(X,Y,YP,FPWA1,FPWA2,FPWA3,FPWA4,FPWA5,IFPWA1,
     &       IFPC1,NEQN-1,IFPC2)
      IF (IFPC2 .GT. 0) RETURN
      SUM = ZERO 
      DO L = 1,NEQN
        PHI(L,1) = YP(L)
        PHI(L,2) = ZERO
        SUM = SUM + (YP(L)/WT(L))**2
      END DO
      SUM = SQRT(SUM) 
      ABSH = ABS(H) 
      IF (EPS .LT. 16*SUM*H*H) ABSH = SQRT(EPS/SUM)/4 
      H = SIGN(MAX(ABSH,FOURU*ABS(X)),H)
C 
C*      U = D1MACH(3) 
C*      BIG = SQRT(D1MACH(2)) 
C*      CALL HSTART (F,NEQN,X,X+H,Y,YP,WT,1,U,BIG,
C*     1             PHI(1,3),PHI(1,4),PHI(1,5),PHI(1,6),RPAR,IPAR,H) 
C 
      HOLD = ZERO
      K = 1 
      KOLD = 0
      KPREV = 0 
      START = .FALSE. 
      PHASE1 = .TRUE. 
      NORND = .TRUE.
      IF (P5EPS .GT. 100*ROUND) GO TO 99 
      NORND = .FALSE. 
      DO L = 1,NEQN
        PHI(L,15) = ZERO 
      END DO
 99   IFAIL = 0 
C       ***     END BLOCK 0     *** 
C 
C       ***     BEGIN BLOCK 1     *** 
C   COMPUTE COEFFICIENTS OF FORMULAS FOR THIS STEP.  AVOID COMPUTING
C   THOSE QUANTITIES NOT CHANGED WHEN STEP SIZE IS NOT CHANGED. 
C                   *** 
C 
 100  KP1 = K+1 
      KP2 = K+2 
      KM1 = K-1 
      KM2 = K-2 
C 
C   NS IS THE NUMBER OF STEPS TAKEN WITH SIZE H, INCLUDING THE CURRENT
C   ONE.  WHEN K.LT.NS, NO COEFFICIENTS CHANGE
C 
      IF(H .NE. HOLD) NS = 0
      IF (NS.LE.KOLD) NS = NS+1 
      NSP1 = NS+1 
      IF (K .LT. NS) GO TO 199
C 
C   COMPUTE THOSE COMPONENTS OF ALPHA(*),BETA(*),PSI(*),SIG(*) WHICH
C   ARE CHANGED 
C 
      BETA(NS) = ONE
      REALNS = REAL(NS, DP) 
      ALPHA(NS) = ONE/REALNS
      TEMP1 = H*REALNS
      SIG(NSP1) = ONE 
      IF(K .LT. NSP1) GO TO 110 
      DO I = NSP1,K 
        IM1 = I-1 
        TEMP2 = PSI(IM1)
        PSI(IM1) = TEMP1
        BETA(I) = BETA(IM1)*PSI(IM1)/TEMP2
        TEMP1 = TEMP2 + H 
        ALPHA(I) = H/TEMP1
        REALI = REAL(I, DP) 
        SIG(I+1) = REALI*ALPHA(I)*SIG(I)
      END DO
 110  PSI(K) = TEMP1
C 
C   COMPUTE COEFFICIENTS G(*) 
C 
C   INITIALIZE V(*) AND SET W(*). 
C 
      IF(NS .GT. 1) GO TO 120 
      DO IQ = 1,K 
        TEMP3 = REAL(IQ*(IQ+1), DP) 
        V(IQ) = ONE/TEMP3 
        W(IQ) = V(IQ) 
      END DO
      IVC = 0 
      KGI = 0 
      IF (K .EQ. 1) GO TO 140 
      KGI = 1 
      GI(1) = W(2)
      GO TO 140 
C 
C   IF ORDER WAS RAISED, UPDATE DIAGONAL PART OF V(*) 
C 
 120  IF(K .LE. KPREV) GO TO 130
      IF (IVC .EQ. 0) GO TO 122 
      JV = KP1 - IV(IVC)
      IVC = IVC - 1 
      GO TO 123 
 122  JV = 1
      TEMP4 = REAL(K*KP1, DP) 
      V(K) = ONE/TEMP4
      W(K) = V(K) 
      IF (K .NE. 2) GO TO 123 
      KGI = 1 
      GI(1) = W(2)
 123  NSM2 = NS-2 
      IF(NSM2 .LT. JV) GO TO 130
      DO J = JV,NSM2
        I = K-J 
        V(I) = V(I) - ALPHA(J+1)*V(I+1) 
        W(I) = V(I) 
      END DO
      IF (I .NE. 2) GO TO 130 
      KGI = NS - 1
      GI(KGI) = W(2)
C 
C   UPDATE V(*) AND SET W(*)
C 
 130  LIMIT1 = KP1 - NS 
      TEMP5 = ALPHA(NS) 
      DO IQ = 1,LIMIT1
        V(IQ) = V(IQ) - TEMP5*V(IQ+1) 
        W(IQ) = V(IQ) 
      END DO
      G(NSP1) = W(1)
      IF (LIMIT1 .EQ. 1) GO TO 137
      KGI = NS
      GI(KGI) = W(2)
 137  W(LIMIT1+1) = V(LIMIT1+1) 
      IF (K .GE. KOLD) GO TO 140
      IVC = IVC + 1 
      IV(IVC) = LIMIT1 + 2
C 
C   COMPUTE THE G(*) IN THE WORK VECTOR W(*)
C 
 140  NSP2 = NS + 2 
      KPREV = K 
      IF(KP1 .LT. NSP2) GO TO 199 
      DO I = NSP2,KP1 
        LIMIT2 = KP2 - I
        TEMP6 = ALPHA(I-1)
        DO IQ = 1,LIMIT2
          W(IQ) = W(IQ) - TEMP6*W(IQ+1) 
        END DO
        G(I) = W(1) 
      END DO
 199  CONTINUE
C       ***     END BLOCK 1     *** 
C 
C       ***     BEGIN BLOCK 2     *** 
C   PREDICT A SOLUTION P(*), EVALUATE DERIVATIVES USING PREDICTED 
C   SOLUTION, ESTIMATE LOCAL ERROR AT ORDER K AND ERRORS AT ORDERS K, 
C   K-1, K-2 AS IF CONSTANT STEP SIZE WERE USED.
C                   *** 
C 
C   INCREMENT COUNTER ON ATTEMPTED STEPS
C 
      KSTEPS = KSTEPS + 1 
C 
C   CHANGE PHI TO PHI STAR
C 
      IF(K .LT. NSP1) GO TO 215 
      DO I = NSP1,K 
        TEMP1 = BETA(I) 
        DO L = 1,NEQN 
          PHI(L,I) = TEMP1*PHI(L,I) 
        END DO
      END DO
C 
C   PREDICT SOLUTION AND DIFFERENCES
C 
 215  DO L = 1,NEQN 
        PHI(L,KP2) = PHI(L,KP1) 
        PHI(L,KP1) = ZERO
        P(L) = ZERO
      END DO
      DO J = 1,K
        I = KP1 - J 
        IP1 = I+1 
        TEMP2 = G(I)
        DO L = 1,NEQN 
          P(L) = P(L) + TEMP2*PHI(L,I)
          PHI(L,I) = PHI(L,I) + PHI(L,IP1)
        END DO
      END DO
      IF(NORND) GO TO 240 
      DO L = 1,NEQN 
        TAU = H*P(L) - PHI(L,15)
        P(L) = Y(L) + TAU 
        PHI(L,16) = (P(L) - Y(L)) - TAU 
      END DO
      GO TO 250 
 240  DO L = 1,NEQN 
        P(L) = Y(L) + H*P(L)
      END DO
 250  XOLD = X
      X = X + H 
      ABSH = ABS(H) 
      CALL FODE(X,P,YP,FPWA1,FPWA2,FPWA3,FPWA4,FPWA5,IFPWA1,
     &       IFPC1,NEQN-1,IFPC2)
      IF (IFPC2 .GT. 0) RETURN
C 
C   ESTIMATE ERRORS AT ORDERS K,K-1,K-2 
C 
      ERKM2 = ZERO 
      ERKM1 = ZERO
      ERK = ZERO
      DO L = 1,NEQN 
        TEMP3 = ONE/WT(L) 
        TEMP4 = YP(L) - PHI(L,1)
        IF (KM2 .GT. 0) ERKM2 = ERKM2 + ((PHI(L,KM1)+TEMP4)*TEMP3)**2
        IF (KM2 .GE. 0) ERKM1 = ERKM1 + ((PHI(L,K)+TEMP4)*TEMP3)**2 
        ERK = ERK + (TEMP4*TEMP3)**2
      END DO
      IF (KM2 .GT. 0) ERKM2 = ABSH*SIG(KM1)*GSTR(KM2)*SQRT(ERKM2) 
      IF (KM2 .GE. 0) ERKM1 = ABSH*SIG(K)*GSTR(KM1)*SQRT(ERKM1) 
      TEMP5 = ABSH*SQRT(ERK)
      ERR = TEMP5*(G(K)-G(KP1)) 
      ERK = TEMP5*SIG(KP1)*GSTR(K)
      KNEW = K
C 
C   TEST IF ORDER SHOULD BE LOWERED 
C 
      IF (KM2 > 0) THEN
        IF(MAX(ERKM1,ERKM2) .LE. ERK) KNEW = KM1
      ELSE IF (KM2 == 0) THEN
        IF(ERKM1 .LE. ERK/2) KNEW = KM1 
      ENDIF
C 
C   TEST IF STEP SUCCESSFUL 
C 
      IF(ERR .LE. EPS) GO TO 400
C       ***     END BLOCK 2     *** 
C 
C       ***     BEGIN BLOCK 3     *** 
C   THE STEP IS UNSUCCESSFUL.  RESTORE  X, PHI(*,*), PSI(*) . 
C   IF THIRD CONSECUTIVE FAILURE, SET ORDER TO ONE.  IF STEP FAILS MORE 
C   THAN THREE TIMES, CONSIDER AN OPTIMAL STEP SIZE.  DOUBLE ERROR
C   TOLERANCE AND RETURN IF ESTIMATED STEP SIZE IS TOO SMALL FOR MACHINE
C   PRECISION.
C                   *** 
C 
C   RESTORE X, PHI(*,*) AND PSI(*)
C 
      PHASE1 = .FALSE.
      X = XOLD
      DO I = 1,K
        TEMP1 = ONE/BETA(I) 
        IP1 = I+1 
        DO L = 1,NEQN 
          PHI(L,I) = TEMP1*(PHI(L,I) - PHI(L,IP1))
        END DO
      END DO
      IF(K .LT. 2) GO TO 320
      DO I = 2,K
        PSI(I-1) = PSI(I) - H 
      END DO
C 
C   ON THIRD FAILURE, SET ORDER TO ONE.  THEREAFTER, USE OPTIMAL STEP 
C   SIZE
C 
 320  IFAIL = IFAIL + 1 
      TEMP2 = 0.5_DP 
      IF (IFAIL-3 .GT. 0) THEN
        IF(P5EPS .LT. ERK/4) TEMP2 = SQRT(P5EPS/ERK) 
      ENDIF
      IF (IFAIL-3 .GE. 0) KNEW = 1
      H = TEMP2*H 
      K = KNEW
      NS = 0
      IF(ABS(H) .GE. FOURU*ABS(X)) GO TO 340
      CRASH = .TRUE.
      H = SIGN(FOURU*ABS(X),H)
      EPS = EPS + EPS 
      RETURN
 340  GO TO 100 
C       ***     END BLOCK 3     *** 
C 
C       ***     BEGIN BLOCK 4     *** 
C   THE STEP IS SUCCESSFUL.  CORRECT THE PREDICTED SOLUTION, EVALUATE 
C   THE DERIVATIVES USING THE CORRECTED SOLUTION AND UPDATE THE 
C   DIFFERENCES.  DETERMINE BEST ORDER AND STEP SIZE FOR NEXT STEP. 
C                   *** 
 400  KOLD = K
      HOLD = H
C 
C   CORRECT AND EVALUATE
C 
      TEMP1 = H*G(KP1)
      IF (NORND) THEN
        DO L = 1,NEQN 
          TEMP3 = Y(L)
          Y(L) = P(L) + TEMP1*(YP(L) - PHI(L,1))
          P(L) = TEMP3
        END DO
      ELSE
        DO L = 1,NEQN 
          TEMP3 = Y(L)
          RHO = TEMP1*(YP(L) - PHI(L,1)) - PHI(L,16)
          Y(L) = P(L) + RHO 
          PHI(L,15) = (Y(L) - P(L)) - RHO 
          P(L) = TEMP3
        END DO
      ENDIF
      CALL FODE(X,Y,YP,FPWA1,FPWA2,FPWA3,FPWA4,FPWA5,IFPWA1,
     &       IFPC1,NEQN-1,IFPC2)
      IF (IFPC2 .GT. 0) RETURN
C 
C   UPDATE DIFFERENCES FOR NEXT STEP
C 
      DO L = 1,NEQN 
        PHI(L,KP1) = YP(L) - PHI(L,1) 
        PHI(L,KP2) = PHI(L,KP1) - PHI(L,KP2)
      END DO
      DO I = 1,K
        DO L = 1,NEQN 
          PHI(L,I) = PHI(L,I) + PHI(L,KP1)
        END DO
      END DO
C 
C   ESTIMATE ERROR AT ORDER K+1 UNLESS: 
C     IN FIRST PHASE WHEN ALWAYS RAISE ORDER, 
C     ALREADY DECIDED TO LOWER ORDER, 
C     STEP SIZE NOT CONSTANT SO ESTIMATE UNRELIABLE 
C 
      ERKP1 = ZERO 
      IF(KNEW .EQ. KM1  .OR.  K .EQ. 12) PHASE1 = .FALSE. 
      IF(PHASE1) GO TO 450
      IF(KNEW .EQ. KM1) GO TO 455 
      IF(KP1 .GT. NS) GO TO 460 
      DO L = 1,NEQN 
        ERKP1 = ERKP1 + (PHI(L,KP2)/WT(L))**2 
      END DO
      ERKP1 = ABSH*GSTR(KP1)*SQRT(ERKP1)
C 
C   USING ESTIMATED ERROR AT ORDER K+1, DETERMINE APPROPRIATE ORDER 
C   FOR NEXT STEP 
C 
      IF(K .GT. 1) GO TO 445
      IF(ERKP1 .GE. ERK/2) GO TO 460
      GO TO 450 
 445  IF(ERKM1 .LE. MIN(ERK,ERKP1)) GO TO 455 
      IF(ERKP1 .GE. ERK  .OR.  K .EQ. 12) GO TO 460 
C 
C   HERE ERKP1 .LT. ERK .LT. MAX(ERKM1,ERKM2) ELSE ORDER WOULD HAVE 
C   BEEN LOWERED IN BLOCK 2.  THUS ORDER IS TO BE RAISED
C 
C   RAISE ORDER 
C 
 450  K = KP1 
      ERK = ERKP1 
      GO TO 460 
C 
C   LOWER ORDER 
C 
 455  K = KM1 
      ERK = ERKM1 
C 
C   WITH NEW ORDER DETERMINE APPROPRIATE STEP SIZE FOR NEXT STEP
C 
 460  HNEW = H + H
      IF(PHASE1) GO TO 465
      IF(P5EPS .GE. ERK*TWO(K+1)) GO TO 465 
      HNEW = H
      IF(P5EPS .GE. ERK) GO TO 465
      TEMP2 = REAL(K+1, DP) 
      R = (P5EPS/ERK)**(ONE/TEMP2)
      HNEW = ABSH*MAX(0.5_dp, MIN(0.9_dp,R)) 
      HNEW = SIGN(MAX(HNEW,FOURU*ABS(X)),H) 
 465  H = HNEW
      RETURN
C       ***     END BLOCK 4     *** 
      END SUBROUTINE STEPS
C
      SUBROUTINE FODE(S,Y,YP,YPOLD,A,QR,ALPHA,TZ,PIVOT,NFE,N,IFLAG)
C
C SUBROUTINE  FODE  IS USED BY SUBROUTINE  STEPS  TO SPECIFY THE
C ORDINARY DIFFERENTIAL EQUATION  DY/DS = G(S,Y) , WHOSE SOLUTION
C IS THE ZERO CURVE OF THE HOMOTOPY MAP.  S = ARC LENGTH,
C YP = DY/DS, AND  Y(S) = (LAMBDA(S), X(S)) .
C
C CALLS  DGEQPF , DNRM2 .
C
      USE HOMPACK_KINDS, ONLY: ZERO, ONE
      USE HOMPACK_CORE_LEGACY, ONLY: F, FJAC, RHOJAC
      USE BLAS_INTERFACES, ONLY: DNRM2
      USE LAPACK_INTERFACES, ONLY: DGEQPF
      IMPLICIT NONE
C
      REAL(DP):: S,YPNORM
      INTEGER:: I,IFLAG,IK,K,KP1,LW,N,NFE,NP1
C
C *****  ARRAY DECLARATIONS.  *****
C
      REAL(DP):: A(:),Y(:),YP(N+1),YPOLD(N+1)
C
C ARRAYS FOR COMPUTING THE JACOBIAN MATRIX AND ITS KERNEL.
      REAL(DP):: ALPHA(3*N+3),QR(N,N+1),TZ(N+1)
      INTEGER, DIMENSION(N+1):: PIVOT
C
C *****  END OF DIMENSIONAL INFORMATION.  *****
C
C
      NP1=N+1
      NFE=NFE+1
C NFE CONTAINS THE NUMBER OF JACOBIAN EVALUATIONS.
C   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
C
C COMPUTE THE JACOBIAN MATRIX, STORE IT IN QR.
C
      IF (IFLAG .EQ. -2) THEN
C
C  QR = ( D RHO(A,LAMBDA,X)/D LAMBDA , D RHO(A,LAMBDA,X)/DX )  .
C
        DO K=1,NP1
          CALL RHOJAC(A,Y(1),Y(2:NP1),QR(:,K),K)
        END DO
      ELSE
        CALL F(Y(2:NP1),TZ(1:N))
        QR(:,1)=A
        IF (IFLAG .EQ. 0) THEN
C
C      QR = ( A - F(X), I - LAMBDA*DF(X) )  .
C
          QR(:,1)=QR(:,1)-TZ(1:N)
          DO K=1,N
            CALL FJAC(Y(2:NP1),TZ(1:N),K)
            KP1=K+1
            QR(:,KP1)=-Y(1)*TZ(1:N)
            QR(K,KP1)= ONE + QR(K,KP1)
          END DO
        ELSE
C
C   QR = ( F(X) - X + A, LAMBDA*DF(X) + (1 - LAMBDA)*I ) .
C
          QR(:,1)=QR(:,1)+TZ(1:N)-Y(2:NP1)
          DO K=1,N
            CALL FJAC(Y(2:NP1),TZ(1:N),K)
            KP1=K+1
            QR(:,KP1)=Y(1)*TZ(1:N)
            QR(K,KP1)=ONE - Y(1)+QR(K,KP1)
          END DO
        ENDIF
      ENDIF
C
C   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
C REDUCE THE JACOBIAN MATRIX TO UPPER TRIANGULAR FORM.
C
        PIVOT = 0
C
      CALL DGEQPF(N,NP1,QR,N,PIVOT,YP,ALPHA,K)
C
      IF (ABS(QR(N,N)) .LE. ABS(QR(1,1))*EPSILON(ONE)) THEN 
        IFLAG=4
        RETURN
      ENDIF 
C COMPUTE KERNEL OF JACOBIAN, WHICH SPECIFIES YP=DY/DS.
      TZ(NP1)=ONE
      DO LW=1,N
        I=NP1-LW
        IK=I+1
        TZ(I)=-DOT_PRODUCT(QR(I,IK:NP1),TZ(IK:NP1))/QR(I,I)
      END DO
      YPNORM=DNRM2(NP1,TZ,1)
      YP(PIVOT)=TZ/YPNORM
      IF (DOT_PRODUCT(YP,YPOLD) .LT. ZERO) YP=-YP
C
C SAVE CURRENT DERIVATIVE (= TANGENT VECTOR) IN  YPOLD .
      YPOLD=YP
      RETURN
      END SUBROUTINE FODE
C
      END MODULE HOMPACK_DF