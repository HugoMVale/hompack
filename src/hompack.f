      MODULE HOMPACK
!
!  This MODULE is an encapsulation of the HOMPACK90 drivers, and uses
!  the modules REAL_PRECISION (defines real precision for all
!  routines), HOMPACK_GLOBAL (allocatable global data structures for
!  sparse matrices), and HOMOTOPY (interfaces to user written routines
!  defining the problem).
!
!  The intended usage is that the calling program would include a
!  statement like, for example,
!     USE HOMPACK90, ONLY: FIXPNF
!
      USE HOMPACK_KINDS, ONLY: DP
      IMPLICIT NONE
C
      CONTAINS
!
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
      USE HOMPACK_INTERFACES, ONLY: F, RHOA
      USE HOMPACK_CORE, ONLY: FODE, STEPS, SINTRP, ROOT
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
!
      SUBROUTINE FIXPDS(N,Y,IFLAG,ARCTOL,EPS,TRACE,A,NDIMA,NFE,
     &     ARCLEN,MODE,LENQR)
C
C Subroutine  FIXPDS  finds a fixed point or zero of the
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
C Y(S) = (X(S), LAMBDA(S)).
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
C  rank [D RHO(A,X,LAMBDA)/D LAMBDA , D RHO(A,X,LAMBDA)/DX] = N
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
C This code is based on the algorithm in L. T. Watson, A
C globally convergent algorithm for computing fixed points of
C C2 maps, Appl. Math. Comput., 5 (1979) 297-311.
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
C  RHOA(V,LAMBDA,X)  which given (X,LAMBDA) returns a
C parameter vector A in V such that RHO(A,X,LAMBDA)=0, and a 
C subroutine  RHOJS(A,LAMBDA,X)  which, if
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
C Whichever of the routines  F,  FJACS,  RHOA,  RHOJS  are required
C should be supplied as external subroutines, conforming with the
C interfaces in the module  HOMOTOPY.
C
C
C Functions and subroutines directly or indirectly called by FIXPDS:
C DLAIC1  and  DLAMCH (LAPACK), F (or  RHOA ), FJACS (or  RHOJS ),
C FODEDS , GMFADS , GMRES , GMRILUDS , ILUFDS , ILUSOLVDS , MULTDS ,
C MULT2DS , PCGDS , ROOT , SINTRP , SOLVDS , STEPDS , and the BLAS
C functions  DDOT , DNRM2.  The module  REAL_PRECISION  specifies 64-bit
C real arithmetic, which the user may want to change.
C
C ***Warning:  this subroutine is generally more robust than  FIXPNS
C and  FIXPQS, but may be slower than those subroutines by a
C factor of two.
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
C    first call to  FIXPDS  for the problem  X=F(X), -1 for the
C    problem  F(X)=0, and -2 for the problem  RHO(A,X,LAMBDA)=0.
C    In certain situations  IFLAG  is set to 2 or 3 by  FIXPDS,
C    and  FIXPDS  can be called again without changing  IFLAG.
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
C A(1:NDIMA) contains the parameter vector  A .  For the fixed point
C    and zero finding problems, A  need not be initialized by the
C    user, and is assumed to have length  N.  For the curve
C    tracking problem, A  has length  NDIMA  and must be initialized
C    by the user.
C
C NDIMA  is the dimension of  A, used for the curve tracking problem,
C    and must be N for the fixed point and zero finding problems.
C
C MODE = 1 if the Jacobian matrix is symmetric and stored in a packed
C          skyline format;
C      = 2 if the Jacobian matrix is stored in a sparse row format.
C
C LENQR  is the number of nonzero entries in the sparse Jacobian
C    matrices, used to determine the sparse matrix data structures.
C
C A, Y, ARCTOL, EPS, ARCLEN, NFE, and IFLAG should all be
C variables in the calling program.
C
C
C ON OUTPUT:
C
C N  and  TRACE  are unchanged.
C
C (Y(1),...,Y(N)) = X, Y(N+1) = LAMBDA, and Y is an approximate
C    zero of the homotopy map.  Normally LAMBDA = 1 and X is a
C    fixed point(zero) of F(X).  In abnormal situations LAMBDA
C    may only be near 1 and X is near a fixed point(zero).
C
C IFLAG =
C  -2   causes  FIXPDS  to initialize everything for the problem
C       RHO(A,X,LAMBDA) = 0 (use on first call).
C
C  -1   causes  FIXPDS  to initialize everything for the problem
C       F(X) = 0 (use on first call).
C
C   0   causes  FIXPDS  to initialize everything for the problem
C       X = F(X) (use on first call).
C
C   1   Normal return.
C
C   2   Specified error tolerance cannot be met.  EPS has been
C       increased to a suitable value.  To continue, just call
C       FIXPDS  again without changing any parameters.
C
C   3   STEPDS  has been called 1000 times.  To continue, call
C       FIXPDS  again without changing any parameters.
C
C   4   Jacobian matrix does not have full rank or has a zero on the
C       diagonal, and/or the conjugate gradient iteration for the
C       kernel of the Jacobian matrix failed to converge.  The
C       algorithm has failed (the zero curve of the homotopy map
C       cannot be followed any further).
C
C   5   EPS  (or  ARCTOL ) is too large.  The problem should be
C       restarted by calling  FIXPDS  with a smaller  EPS  (or
C       ARCTOL ) and  IFLAG = 0 (-1, -2).
C
C   6   I - DF(X)  is nearly singular at the fixed point (DF(X) is
C       nearly singular at the zero, or  D RHO(A,X,LAMBDA)/DX  is
C       nearly singular at  LAMBDA = 1 ).  Answer may not be
C       accurate.
C
C   7   Illegal input parameters, a fatal error.
C
C ARCTOL = EPS after a normal return (IFLAG = 1).
C
C EPS  is unchanged after a normal return (IFLAG = 1).  It is
C    increased to an appropriate value on the return IFLAG = 2.
C
C A  will (normally) have been modified.
C
C NFE  is the number of function evaluations (= number of
C    Jacobian evaluations).
C
C ARCLEN  is the length of the path followed.
C
C
C Allocatable and automatic work arrays:
C
C YP(1:N+1) is a work array containing the current tangent
C    vector to the zero curve.
C
C YPOLD(1:N+1) is a work array containing the previous tangent
C    vector to the zero curve.
C
C QR(1:LENQR), PP(1:N), ROWPOS(1:N+2), COLPOS(1:LENQR) are all work
C    arrays used to define the sparse Jacobian matrices, allocated
C    here, and distributed via the module  HOMOTOPY .
C
C WT(1:N+1), PHI(1:N+1,1:16), and P(1:N+1) are all work arrays
C    used by the ODE subroutine  STEPDS  .
C
      USE HOMPACK_INTERFACES, ONLY: F, RHOA
      USE HOMPACK_CORE, ONLY: FODEDS, STEPDS, SINTRP, ROOT
      USE HOMPACK_GLOBAL, ONLY: QR => QRSPARSE, ROWPOS, COLPOS, PP,
     & PAR, IPAR
      IMPLICIT NONE
C
      INTEGER, INTENT(IN)::LENQR,MODE,N,NDIMA,TRACE
      REAL (dp), DIMENSION(:), INTENT(IN OUT)::A,Y
      INTEGER, INTENT(IN OUT)::IFLAG
      REAL (dp), INTENT(IN OUT)::ARCTOL,EPS
      INTEGER, INTENT(OUT)::NFE
      REAL (dp), INTENT(OUT)::ARCLEN
C
C *****  LOCAL VARIABLES.  *****
C
      REAL (dp), SAVE:: CURSW,CURTOL,EPSSTP,EPST,
     &  H,HOLD,S,S99,SA,SB,SOUT,SQNP1,XOLD,Y1SOUT
      INTEGER, SAVE:: IFLAGC,ITER,IVC,JW,K,KGI,KOLD,
     &  KSTEPS,LCODE,LIMIT,NFEC,NP1
      LOGICAL, SAVE:: CRASH,START,ST99
C
C ARRAYS NEEDED BY THE ODE SUBROUTINE  STEPDS .
      REAL (dp), SAVE:: ALPHAS(12),G(13),GI(11),W(12)
      REAL (dp), ALLOCATABLE, SAVE:: P(:),PHI(:,:),WT(:),YP(:)
      INTEGER, SAVE:: IV(10)
C
C ARRAYS NEEDED BY  FIXPDS , FODEDS , AND  PCGDS .
      REAL (dp), ALLOCATABLE, DIMENSION(:), SAVE:: AOLD,YPOLD
C
C *****  END OF DIMENSIONAL INFORMATION.  *****
C
C LIMITD  IS AN UPPER BOUND ON THE NUMBER OF STEPS.  IT MAY BE
C CHANGED BY CHANGING THE FOLLOWING PARAMETER STATEMENT:
      INTEGER, PARAMETER:: LIMITD=1000
C
C
C :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :
      IF (N .LE. 0  .OR.  EPS .LE. 0.0  .OR.  N+1 .NE. SIZE(Y)
     &  .OR.  NDIMA .NE. SIZE(A)  .OR.
     &  ((IFLAG .EQ. -1  .OR.  IFLAG .EQ. 0) .AND.  N .NE. SIZE(A))
     &  .OR.  MODE .LE. 0  .OR.  MODE .GE. 3)
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
      SQNP1=SQRT(REAL(NP1,kind=dp))
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
C ALLOCATE SAVED WORK ARRAYS.
      IF (ALLOCATED(AOLD)) DEALLOCATE(AOLD)
      IF (ALLOCATED(P)) DEALLOCATE(P)
      IF (ALLOCATED(PHI)) DEALLOCATE(PHI)
      IF (ALLOCATED(WT)) DEALLOCATE(WT)
      IF (ALLOCATED(YP)) DEALLOCATE(YP)
      IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
      ALLOCATE(AOLD(NDIMA),P(NP1),PHI(NP1,16),WT(NP1),YP(NP1),
     &  YPOLD(NP1))
C SET INITIAL CONDITIONS FOR ORDINARY DIFFERENTIAL EQUATION.
      YPOLD(NP1)=1.0
      YP(NP1)=1.0
      Y(NP1)=0.0
      YPOLD(1:N)=0.0
      YP(1:N)=0.0
C LOAD  A  FOR THE FIXED POINT AND ZERO FINDING PROBLEMS.
      IF (IFLAGC .GE. -1) THEN
        A=Y(1:N)
      ENDIF
30    LIMIT=LIMITD
C ALLOCATE ARRAYS FOR SPARSE JACOBIAN MATRIX DATA STRUCTURE.
35    SELECT CASE (MODE)
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
        ARCLEN=ARCLEN+S
        IFLAG=5
        CALL CLEANUPALL
        RETURN
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
        CALL RHOA(A,Y(NP1),Y(1:N))
C IF NEW AND OLD  A  DIFFER BY TOO MUCH, TRACKING SHOULD NOT CONTINUE.
        IF (ANY(ABS(A-AOLD) .GT. 1.0+ABS(AOLD))) THEN
            ARCLEN=ARCLEN+S
            IFLAG=5
            CALL CLEANUPALL
            RETURN
        ENDIF
      ELSE
        CALL F(Y(1:N),YP(1:N))
        AOLD=A
        IF (IFLAGC .EQ. -1) THEN
          A=Y(NP1)*YP(1:N)/(1.0 - Y(NP1)) + Y(1:N)
        ELSE
          A=(Y(1:N) - Y(NP1)*YP(1:N))/(1.0 - Y(NP1))
        ENDIF
C IF NEW AND OLD  A  DIFFER BY TOO MUCH, TRACKING SHOULD NOT CONTINUE.
        IF (ANY(ABS(A-AOLD) .GT. 1.0+ABS(AOLD))) THEN
            ARCLEN=ARCLEN+S
            IFLAG=5
            CALL CLEANUPALL
            RETURN
        ENDIF
      ENDIF
      GO TO 100
80    IF (Y(NP1) .LE. .99  .OR. ST99) GO TO 100
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
      CALL STEPDS(FODEDS,NP1,Y,S,H,EPSSTP,WT,START,HOLD,K,KOLD,CRASH,
     &     PHI,P,YP,ALPHAS,W,G,KSTEPS,XOLD,IVC,IV,KGI,GI,
     &     IFLAGC,YPOLD,A,NDIMA,LENQR,MODE,NFEC)
C PRINT LATEST POINT ON CURVE IF REQUESTED.
      IF (TRACE .GT. 0) THEN
        WRITE (TRACE,117) ITER,NFEC,S,Y(NP1),(Y(JW),JW=1,N)
117     FORMAT(/' STEP',I5,3X,'NFE =',I5,3X,'ARC LENGTH =',F9.4,3X,
     &  'LAMBDA =',F7.4,5X,'X vector:'/(1X,6ES12.4))
      ENDIF
      NFE=NFEC
C CHECK IF THE STEP WAS SUCCESSFUL.
      IF (IFLAGC .EQ. 4) THEN
        ARCLEN=ARCLEN+S
        IFLAG=4
        CALL CLEANUPALL
        RETURN
      ENDIF
      IF (CRASH) THEN
C RETURN CODE FOR ERROR TOLERANCE TOO SMALL.
        IFLAG=2
C CHANGE ERROR TOLERANCES.
        EPS=EPSSTP
        IF (ARCTOL .LT. EPS) ARCTOL=EPS
C CHANGE LIMIT ON NUMBER OF ITERATIONS.
        LIMIT=LIMIT-ITER
        CALL CLEANUP
        RETURN
      ENDIF
C
      IF (Y(NP1) .GE. 1.0) THEN
        IF (ST99) GO TO 160
C
C IF LAMBDA .GE. 1.0 BUT THE PROBLEM HAS NOT BEEN RESTARTED
C WITH A NEW  A  VECTOR, BACK UP AND RESTART.
C
        S99=S-.5*HOLD
C GET AN APPROXIMATE ZERO Y(S) WITH  Y(NP1)=LAMBDA .LT. 1.0  .
135     CALL SINTRP(S,Y,S99,WT,YP,NP1,KOLD,PHI,IVC,IV,KGI,GI,
     &     ALPHAS,G,W,XOLD,P)
        IF (WT(NP1) .LT. 1.0) GO TO 140
        S99=.5*(S-HOLD+S99)
        GO TO 135
C
140     Y=WT
        YPOLD=YP
        S=S99
        GO TO 90
      ENDIF
C
      END DO MAIN_LOOP ! *****  END OF MAIN LOOP.  *****
C
C LAMBDA HAS NOT REACHED 1 IN 1000 STEPS.
      IFLAG=3
      CALL CLEANUP
      RETURN
C
C USE INVERSE INTERPOLATION TO GET THE ANSWER AT LAMBDA = 1.0 .
C
160   SA=S-HOLD
      SB=S
      LCODE=1
170   CALL ROOT(SOUT,Y1SOUT,SA,SB,EPS,EPS,LCODE)
C ROOT  FINDS S SUCH THAT Y(NP1)(S) = LAMBDA = 1 .
      IF (LCODE .GT. 0) GO TO 190
      CALL SINTRP(S,Y,SOUT,WT,YP,NP1,KOLD,PHI,IVC,IV,KGI,GI,
     &     ALPHAS,G,W,XOLD,P)
      Y1SOUT=WT(NP1)-1.0
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
      CALL CLEANUPALL
      RETURN
C
      CONTAINS
C
      SUBROUTINE CLEANUPALL
      IF (ALLOCATED(AOLD)) DEALLOCATE(AOLD)
      IF (ALLOCATED(P)) DEALLOCATE(P)
      IF (ALLOCATED(PHI)) DEALLOCATE(PHI)
      IF (ALLOCATED(WT)) DEALLOCATE(WT)
      IF (ALLOCATED(YP)) DEALLOCATE(YP)
      IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
      CALL CLEANUP
      RETURN
      END SUBROUTINE CLEANUPALL
C      
      SUBROUTINE CLEANUP
      IF (ALLOCATED(QR)) DEALLOCATE(QR)
      IF (ALLOCATED(ROWPOS)) DEALLOCATE(ROWPOS)
      IF (ALLOCATED(COLPOS)) DEALLOCATE(COLPOS)
      IF (ALLOCATED(PP)) DEALLOCATE(PP)
      IF (ALLOCATED(PAR)) DEALLOCATE(PAR)
      IF (ALLOCATED(IPAR)) DEALLOCATE(IPAR)
      RETURN
      END SUBROUTINE CLEANUP
C
      END SUBROUTINE FIXPDS
C
      SUBROUTINE FIXPNF(N,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,A,
     &   SSPAR,NFE,ARCLEN,  POLY_SWITCH)
C
C Subroutine  FIXPNF  finds a fixed point or zero of the
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
C Y(S) = (LAMBDA(S), X(S)) using a Hermite cubic predictor and a
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
C
C For the fixed point and zero finding problems, the user must supply 
C a subroutine  F(X,V)  which evaluates F(X) at X and returns the 
C vector F(X) in V, and a subroutine  FJAC(X,V,K)  which returns in V 
C the Kth column of the Jacobian matrix of F(X) evaluated at X.  For 
C the curve tracking problem, the user must supply a subroutine  
C  RHO(A,LAMBDA,X,V)  which evaluates the homotopy map RHO at 
C (A,LAMBDA,X) and returns the vector RHO(A,LAMBDA,X) in V, and a
C subroutine  RHOJAC(A,LAMBDA,X,V,K)  which returns in V the Kth
C column of the N X (N+1) Jacobian matrix [D RHO/D LAMBDA, D RHO/DX]
C evaluated at (A,LAMBDA,X).  FIXPNF  directly or indirectly uses
C the subroutines  F (or  RHO ),  FJAC (or  RHOJAC ), 
C   ROOT,  ROOTNF,  STEPNF,  the LAPACK routines  DGEQPF,  DORMQR,  
C their auxiliary routines, and the BLAS routines  DCOPY,
C   DDOT,  DGEMM,  DGEMV,  DGER,  DNRM2,  DSCAL,  DSWAP,  DTRMM,  DTRMV, 
C   IDAMAX.  The module  REAL_PRECISION  specifies 64-bit
C real arithmetic, which the user may want to change.
C
C
C ON INPUT:
C
C N  is the dimension of X, F(X), and RHO(A,LAMBDA,X).
C
C Y(:)  is an array of length  N + 1.  (Y(2),...,Y(N+1)) = A  is the
C    starting point for the zero curve for the fixed point and 
C    zero finding problems.  (Y(2),...,Y(N+1)) = X0  for the curve
C    tracking problem.
C
C IFLAG  can be -2, -1, 0, 2, or 3.  IFLAG  should be 0 on the 
C    first call to  FIXPNF  for the problem  X=F(X), -1 for the
C    problem  F(X)=0, and -2 for the problem  RHO(A,LAMBDA,X)=0.
C    In certain situations  IFLAG  is set to 2 or 3 by  FIXPNF,
C    and  FIXPNF  can be called again without changing  IFLAG.
C
C ARCRE , ARCAE  are the relative and absolute errors, respectively,
C    allowed the normal flow iteration along the zero curve.  If
C    ARC?E .LE. 0.0  on input it is reset to  .5*SQRT(ANS?E) .
C    Normally  ARC?E should be considerably larger than  ANS?E .
C
C ANSRE , ANSAE  are the relative and absolute error values used for
C    the answer at LAMBDA = 1.  The accepted answer  Y = (LAMBDA, X)
C    satisfies
C
C       |Y(1) - 1|  .LE.  ANSRE + ANSAE           .AND.
C
C       ||Z||  .LE.  ANSRE*||X|| + ANSAE          where
C
C    (.,Z) is the Newton step to Y.
C
C TRACE  is an integer specifying the logical I/O unit for
C    intermediate output.  If  TRACE .GT. 0  the points computed on
C    the zero curve are written to I/O unit  TRACE .
C
C A(:)  contains the parameter vector  A .  For the fixed point
C    and zero finding problems, A  need not be initialized by the
C    user, and is assumed to have length  N.  For the curve
C    tracking problem, A  must be initialized by the user.
C
C SSPAR(1:8) = (LIDEAL, RIDEAL, DIDEAL, HMIN, HMAX, BMIN, BMAX, P)  is
C    a vector of parameters used for the optimal step size estimation.
C    If  SSPAR(J) .LE. 0.0  on input, it is reset to a default value
C    by  FIXPNF .  Otherwise the input value of  SSPAR(J)  is used.
C    See the comments below and in  STEPNF  for more information about
C    these constants.
C
C POLY_SWITCH  is an optional logical variable used only by the driver
C    POLSYS1H  for polynomial systems.
C
C
C ON OUTPUT:
C
C N , TRACE , A  are unchanged.
C
C Y(1) = LAMBDA, (Y(2),...,Y(N+1)) = X, and Y is an approximate
C    zero of the homotopy map.  Normally LAMBDA = 1 and X is a
C    fixed point(zero) of F(X).  In abnormal situations LAMBDA
C    may only be near 1 and X is near a fixed point(zero).
C
C IFLAG =
C  -2   causes  FIXPNF  to initialize everything for the problem
C       RHO(A,LAMBDA,X) = 0 (use on first call).
C
C  -1   causes  FIXPNF  to initialize everything for the problem
C       F(X) = 0 (use on first call).
C
C   0   causes  FIXPNF  to initialize everything for the problem
C       X = F(X) (use on first call).
C
C   1   Normal return.
C
C   2   Specified error tolerance cannot be met.  Some or all of
C       ARCRE , ARCAE , ANSRE , ANSAE  have been increased to 
C       suitable values.  To continue, just call  FIXPNF  again 
C       without changing any parameters.
C
C   3   STEPNF  has been called 1000 times.  To continue, call
C       FIXPNF  again without changing any parameters.
C
C   4   Jacobian matrix does not have full rank.  The algorithm
C       has failed (the zero curve of the homotopy map cannot be
C       followed any further).
C
C   5   The tracking algorithm has lost the zero curve of the
C       homotopy map and is not making progress.  The error tolerances
C       ARC?E  and  ANS?E  were too lenient.  The problem should be
C       restarted by calling  FIXPNF  with smaller error tolerances
C       and  IFLAG = 0 (-1, -2).
C
C   6   The normal flow Newton iteration in  STEPNF  or  ROOTNF
C       failed to converge.  The error tolerances  ANS?E  may be too
C       stringent.
C
C   7   Illegal input parameters, a fatal error.
C
C   8   Memory allocation error, fatal.
C
C ARCRE , ARCAE , ANSRE , ANSAE  are unchanged after a normal return 
C    (IFLAG = 1).  They are increased to appropriate values on the 
C    return  IFLAG = 2 .
C
C NFE  is the number of function evaluations (= number of
C    Jacobian matrix evaluations).
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
C QR(1:N,1:N+2), ALPHA(1:3*N+3), TZ(1:N+1), PIVOT(1:N+1), W(1:N+1),
C    WP(1:N+1), Z0(1:N+1), Z1(1:N+1)  are all work arrays used by
C    STEPNF  to calculate the tangent vectors and Newton steps.
C
C
      USE HOMPACK_CORE, ONLY: STEPNF, ROOTNF
      USE BLAS_INTERFACES, ONLY: DNRM2
      IMPLICIT NONE
C
      INTEGER, INTENT(IN)::N,TRACE
      REAL (dp), DIMENSION(:), INTENT(IN OUT)::A,Y
      INTEGER, INTENT(IN OUT)::IFLAG
      REAL (dp), INTENT(IN OUT)::ANSAE,ANSRE,ARCAE,ARCRE,
     &    SSPAR(8)
      INTEGER, INTENT(OUT)::NFE
      REAL (dp), INTENT(OUT)::ARCLEN
      LOGICAL, INTENT(IN), OPTIONAL::POLY_SWITCH
C
C LOCAL VARIABLES.
      REAL (dp), SAVE:: ABSERR,CURTOL,H,HOLD,RELERR,S
      INTEGER, SAVE:: IFLAGC,ITER,JW,LIMIT,NC,NFEC,NP1
      LOGICAL, SAVE:: CRASH,POLSYS,START
C
C ALLOCATABLE AND AUTOMATIC ARRAYS.
      REAL (dp), DIMENSION(:), ALLOCATABLE, SAVE:: YOLD,YP,YPOLD
      REAL (dp):: ALPHA(3*N+3),QR(N,N+2),TZ(N+1),
     &  W(N+1),WP(N+1),Z0(N+1),Z1(N+1)
      INTEGER:: PIVOT(N+1)
C
C ***** END OF DIMENSIONAL INFORMATION. *****
C
C LIMITD  IS AN UPPER BOUND ON THE NUMBER OF STEPS.  IT MAY BE
C CHANGED BY CHANGING THE FOLLOWING PARAMETER STATEMENT:
      INTEGER, PARAMETER:: LIMITD=1000
C
C SWITCH FROM THE TOLERANCE  ARC?E  TO THE (FINER) TOLERANCE  ANS?E  IF
C THE CURVATURE OF ANY COMPONENT OF  Y  EXCEEDS  CURSW.
      REAL (dp), PARAMETER:: CURSW=10.0
C
C :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :  :
C TEST LOGICAL SWITCH TO REFLECT INTENDED USAGE OF FIXPNF.
      IF (PRESENT(POLY_SWITCH)) THEN
        POLSYS=.TRUE.
      ELSE
        POLSYS=.FALSE.
      ENDIF
C
      IF (N .LE. 0  .OR.  ANSRE .LE. 0.0  .OR.  ANSAE .LT. 0.0
     &  .OR.  N+1 .NE. SIZE(Y)  .OR.
     &  ((IFLAG .EQ. -1  .OR.  IFLAG .EQ. 0) .AND.  N .NE. SIZE(A)))
     &  IFLAG=7
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
      IF (ALLOCATED(YP)) DEALLOCATE(YP)
      IF (ALLOCATED(YOLD)) DEALLOCATE(YOLD)
      IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
      ALLOCATE(YP(NP1),YOLD(NP1),YPOLD(NP1),STAT=JW)
      IF (JW /= 0) THEN
        IFLAG=8
        RETURN
      END IF
C SET INITIAL CONDITIONS FOR FIRST CALL TO  STEPNF .
      START=.TRUE.
      CRASH=.FALSE.
      HOLD=1.0
      H=.1
      S=0.0
      YPOLD(1)=1.0
      YP(1)=1.0
      Y(1)=0.0
      YPOLD(2:NP1)=0.0
      YP(2:NP1)=0.0
C SET OPTIMAL STEP SIZE ESTIMATION PARAMETERS.
C LET Z[K] DENOTE THE NEWTON ITERATES ALONG THE FLOW NORMAL TO THE
C DAVIDENKO FLOW AND Y THEIR LIMIT.
C IDEAL CONTRACTION FACTOR:  ||Z[2] - Z[1]|| / ||Z[1] - Z[0]||
      IF (SSPAR(1) .LE. 0.0) SSPAR(1)= .5
C IDEAL RESIDUAL FACTOR:  ||RHO(A, Z[1])|| / ||RHO(A, Z[0])||
      IF (SSPAR(2) .LE. 0.0) SSPAR(2)= .01
C IDEAL DISTANCE FACTOR:  ||Z[1] - Y|| / ||Z[0] - Y||
      IF (SSPAR(3) .LE. 0.0) SSPAR(3)= .5
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
      IF (IFLAGC .GE. -1) THEN
        A=Y(2:NP1)
      ENDIF
90    LIMIT=LIMITD
C
C *****  END OF INITIALIZATION BLOCK.  *****
C
120   MAIN_LOOP: DO ITER=1,LIMIT  ! *****  MAIN LOOP.  *****
      IF (Y(1) .LT. 0.0) THEN
        ARCLEN=S
        IFLAG=5
        CALL CLEANUP ; RETURN
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
      CALL STEPNF(NC,NFEC,IFLAGC,START,CRASH,HOLD,H,RELERR,ABSERR,
     &     S,Y,YP,YOLD,YPOLD,A,QR,ALPHA,TZ,PIVOT,W,WP,Z0,Z1,SSPAR)
C PRINT LATEST POINT ON CURVE IF REQUESTED.
      IF (TRACE .GT. 0) THEN
        WRITE (TRACE,217) ITER,NFEC,S,Y(1),(Y(JW),JW=2,NP1)
217     FORMAT(/' STEP',I5,3X,'NFE =',I5,3X,'ARC LENGTH =',F9.4,3X,
     &  'LAMBDA =',F7.4,5X,'X VECTOR:'/(1X,6ES12.4))
      ENDIF
      NFE=NFEC
C CHECK IF THE STEP WAS SUCCESSFUL.
      IF (IFLAGC .GT. 0) THEN
        ARCLEN=S
        IFLAG=IFLAGC
        CALL CLEANUP ; RETURN
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
        RETURN
      ENDIF
C
      IF (Y(1) .GE. 1.0) THEN
C
C USE HERMITE CUBIC INTERPOLATION AND NEWTON ITERATION TO GET THE 
C ANSWER AT LAMBDA = 1.0 .
C
C SAVE  YOLD  FOR ARC LENGTH CALCULATION LATER.
        Z0=YOLD
        CALL ROOTNF(NC,NFEC,IFLAGC,ANSRE,ANSAE,Y,YP,YOLD,YPOLD,
     &              A,QR,ALPHA,TZ,PIVOT,W,WP)
C
        NFE=NFEC
        IFLAG=1
C SET ERROR FLAG IF  ROOTNF  COULD NOT GET THE POINT ON THE ZERO
C CURVE AT  LAMBDA = 1.0  .
        IF (IFLAGC .GT. 0) IFLAG=IFLAGC
C CALCULATE FINAL ARC LENGTH.
        W=Y-Z0
        ARCLEN=S - HOLD + DNRM2(NP1,W,1)
        CALL CLEANUP ; RETURN
      ENDIF
C
C FOR POLYNOMIAL SYSTEMS AND THE  POLSYS1H  HOMOTOPY MAP,
C D LAMBDA/DS .GE. 0 NECESSARILY.  THIS CONDITION IS FORCED HERE IF
C THE  POLY_SWITCH  VARIABLE IS PRESENT.
C
      IF (POLSYS) THEN
        IF (YP(1) .LT. 0.0) THEN
C REVERSE TANGENT DIRECTION SO D LAMBDA/DS = YP(1) > 0 .
          YP=-YP
          YPOLD=YP
C FORCE  STEPNF  TO USE THE LINEAR PREDICTOR FOR THE NEXT STEP ONLY.
          START=.TRUE.
        ENDIF
      ENDIF
C
      END DO MAIN_LOOP   ! *****  END OF MAIN LOOP.  *****
C
C LAMBDA HAS NOT REACHED 1 IN 1000 STEPS.
      IFLAG=3
      ARCLEN=S
      RETURN
C
      CONTAINS
        SUBROUTINE CLEANUP
        IF (ALLOCATED(YP)) DEALLOCATE(YP)
        IF (ALLOCATED(YOLD)) DEALLOCATE(YOLD)
        IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
        END SUBROUTINE CLEANUP
      END SUBROUTINE FIXPNF
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
      USE HOMPACK_GLOBAL, ONLY: QR => QRSPARSE, ROWPOS, COLPOS, PP,
     & PAR, IPAR
      USE HOMPACK_CORE, ONLY: ROOTNS, STEPNS  
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
      SUBROUTINE FIXPQF(N,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,A,
     &     SSPAR,NFE,ARCLEN)
C
C Subroutine  FIXPQF  finds a fixed point or zero of the 
C N-dimensional vector function  F(X), or tracks a zero curve of a 
C general homotopy map  RHO(A,LAMBDA,X).  For the fixed point problem
C F(X) is assumed to be a C2 map of some ball into itself.  The 
C equation  X=F(X)  is solved by following the zero curve of the 
C homotopy map
C
C  LAMBDA*(X - F(X)) + (1 - LAMBDA)*(X - A) ,
C
C starting from  LAMBDA = 0, X = A.   The curve is parameterized
C by arc length  S, and is followed by solving the ordinary 
C differential equation  D(HOMOTOPY MAP)/DS = 0  for  
C Y(S) = (LAMBDA(S), X(S)).  This is done by using a Hermite cubic 
C predictor and a corrector which returns to the zero curve in a 
C hyperplane perpendicular to the tangent to the zero curve at the 
C most recent point.
C
C For the zero finding problem  F(X)  is assumed to be a C2 map
C such that for some  R > 0,  X*F(X) >= 0  whenever  NORM(X) = R.
C The equation  F(X) = 0  is solved by following the zero curve of
C the homotopy map
C
C  LAMBDA*F(X) + (1 - LAMBDA)*(X - A)
C
C emanating from  LAMBDA = 0, X = A.
C
C A  must be an interior point of the above mentioned balls.
C
C For the curve tracking problem RHO(A,LAMBDA,X) is assumed to 
C be a C2 map from  E**M X [0,1) X E**N  into  E**N, which for 
C almost all parameter vectors  A  in some nonempty open subset
C of E**M satisfies
C
C  rank [D RHO(A,LAMBDA,X)/D LAMBDA, D RHO(A,LAMBDA,X)/DX] = N
C
C for all points  (LAMBDA,X)  such that  RHO(A,LAMBDA,X) = 0.  It is
C further assumed that
C
C         rank [ D RHO(A,0,X0)/DX ] = N.
C
C With  A  fixed, the zero curve of  RHO(A,LAMBDA,X)  emanating from
C LAMBDA = 0, X = X0  is tracked until  LAMBDA = 1  by solving the 
C ordinary differential equation    D RHO(A,LAMBDA(S),X(S))/DS = 0
C for  Y(S) = (LAMBDA(S), X(S)), where  S  is arc length along the
C zero curve.  Also the homotopy map  RHO(A,LAMBDA,X)  is assumed to
C be constructed such that
C
C         D LAMBDA(0)/DS > 0.
C
C For the fixed point and zero finding problems, the user must supply
C a subroutine  F(X,V)  which evaluates  F(X)  at  X  and returns the
C vector F(X) in  V, and a subroutine  FJAC(X,V,K)  which returns in  V
C the Kth column of the Jacobian matrix of F(X) evaluated at X.  For
C the curve tracking problem, the user must supply a subroutine
C RHO(A,LAMBDA,X,V)  which evaluates the homotopy map  RHO at
C (A,LAMBDA,X)  and returns the vector  RHO(A,LAMBDA,X)  in  V, and
C a subroutine  RHOJAC(A,LAMBDA,X,V,K)  which returns in  V
C the Kth column of the  N X (N+1)  Jacobian matrix  
C [D RHO/D LAMBDA, D RHO/DX]  evaluated at  (A,LAMBDA,X).  FIXPQF
C directly or indirectly uses the subroutines  F (or RHO), 
C   FJAC (or RHOJAC),  ROOT,  ROOTQF,  STEPQF,  TANGQF,  UPQRQF,
C the LAPACK routines  DGEQRF,  DORGQR, their auxiliary routines,
C and the BLAS routines  DCOPY,  DDOT,  DGEMM,  DGEMV,
C   DGER,  DNRM2,  DSCAL,  DTPMV,  DTPSV,  DTRMM,  and   DTRMV.
C The module  REAL_PRECISION  specifies 64-bit real arithmetic,
C which the user may want to change.
C
C
C ON INPUT:
C
C N  is the dimension of X, F(X), and RHO(A,LAMBDA,X).
C
C Y(1:N+1)  contains the starting point for tracking the homotopy map.
C    (Y(2),...,Y(N+1)) = A  for the fixed point and zero finding 
C    problems.  (Y(2),...,Y(N+1)) = X0  for the curve tracking problem.
C    Y(1)  need not be defined by the user.
C
C IFLAG can be -2, -1, 0, 2, or 3.  IFLAG should be 0 on the first
C    call to  FIXPQF  for the problem  X=F(X), -1 for the problem
C    F(X)=0, and -2 for the problem  RHO(A,LAMBDA,X)=0.   In certain
C    situations  IFLAG  is set to 2 or 3 by  FIXPQF, and  FIXPQF  can
C    be called again without changing  IFLAG.
C
C ARCRE, ARCAE  are the relative and absolute errors, respectively,
C    allowed the quasi-Newton iteration along the zero curve.  If
C    ARC?E .LE. 0.0  on input, it is reset to  .5*SQRT(ANS?E).
C    Normally  ARC?E  should be considerably larger than  ANS?E.
C
C ANSRE, ANSAE  are the relative and absolute error values used for 
C    the answer at  LAMBDA = 1.  The accepted answer  Y = (LAMBDA, X)
C    satisfies
C
C      |Y(1) - 1| .LE. ANSRE + ANSAE      .AND.
C  
C      ||DZ|| .LE. ANSRE*||Y|| + ANSAE      where
C
C      DZ is the quasi-Newton step to Y.
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
C SSPAR(1:4) =  (HMIN, HMAX, BMIN, BMAX)  is a vector of parameters 
C    used for the optimal step size estimation.  A default value
C    can be specified for any of these four parameters by setting it
C    .LE. 0.0  on input.  See the comments in  STEPQF  for more
C    information about these parameters.
C
C
C ON OUTPUT:
C
C N , TRACE , A  are unchanged.
C
C Y(1) = LAMBDA, (Y(2),...,Y(N+1)) = X, and  Y  is an approximate
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
C       suitable values.  To continue, just call  FIXPQF  again
C       without changing any parameters.
C
C   3   STEPQF  has been called 1000 times.  To continue, call
C       FIXPQF  again without changing any parameters.
C
C   4   Jacobian matrix does not have full rank.  The algorithm
C       has failed (the zero curve of the homotopy map cannot be
C       followed any further).
C
C   5   The tracking algorithm has lost the zero curve of the 
C       homotopy map and is not making progress.  The error 
C       tolerances  ARC?E  and  ANS?E  were too lenient.  The problem 
C       should be restrarted by calling  FIXPQF  with smaller error 
C       tolerances and  IFLAG = 0 (-1, -2).
C
C   6   The quasi-Newton iteration in  STEPQF  or  ROOTQF  failed to
C       converge.  The error tolerances  ANS?E  may be too stringent.
C
C   7   Illegal input parameters, a fatal error.
C
C   8   Memory allocation error, fatal.
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
C YP(1:N+1)  is a work array containing the tangent vector to the
C    zero curve at the current point  Y.
C
C YOLD(1:N+1) is a work array containing the previous point found
C    on the zero curve.
C
C YPOLD(1:N+1) is a work array containing the tangent vector to
C    the zero curve at  YOLD.
C
C Q(1:N+1,1:N+1), R((N+1)*(N+2)/2), F0(1:N+1), F1(1:N+1), Z0(1:N+1),
C    DZ(1:N+1), W(1:N+1), T(1:N+1), YSAV(1:N+1)  are all work arrays 
C    used by  STEPQF, TANGQF and ROOTQF to calculate the tangent 
C    vectors and quasi-Newton steps.
C
C
C ***** DECLARATIONS *****
      USE HOMPACK_CORE, ONLY: ROOTQF, STEPQF
      USE BLAS_INTERFACES, ONLY: DNRM2
      IMPLICIT NONE
C
C     LOCAL VARIABLES 
C
      REAL (dp), SAVE:: ABSERR, H, HOLD, RELERR, S, WK 
      INTEGER, SAVE:: IFLAGC, ITER, JW, LIMIT, NP1
      LOGICAL, SAVE:: CRASH, START       
C
C     SCALAR ARGUMENTS 
C
      REAL (dp):: ARCRE, ARCAE, ANSRE, ANSAE, ARCLEN
      INTEGER:: N,IFLAG,TRACE,NFE
C
C     ARRAY DECLARATIONS 
C
      REAL (dp), DIMENSION(:), ALLOCATABLE, SAVE::
     &  R,YOLD,YP,YPOLD
      REAL (dp), DIMENSION(:,:), ALLOCATABLE, SAVE:: Q
      REAL (dp):: A(:), DZ(N+1), F0(N+1), F1(N+1),
     &    SSPAR(4), T(N+1), W(N+1), Y(:), YSAV(N+1), Z0(N+1)
C 
C ***** END OF DECLARATIONS *****
C
C LIMITD IS AN UPPER BOUND ON THE NUMBER OF STEPS.  IT MAY BE 
C CHANGED BY CHANGING THE FOLLOWING PARAMETER STATEMENT:
      INTEGER, PARAMETER:: LIMITD = 1000
C
C ***** FIRST EXECUTABLE STATEMENT *****
C
C CHECK IFLAG
C
      IF (N .LE. 0  .OR.  ANSRE .LE. 0.0  .OR.  ANSAE .LT. 0.0
     &  .OR.  N+1 .NE. SIZE(Y)  .OR.
     &  ((IFLAG .EQ. -1  .OR.  IFLAG .EQ. 0) .AND.  N .NE. SIZE(A)))
     &  IFLAG=7
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
      IF (ALLOCATED(Q)) DEALLOCATE(Q)
      IF (ALLOCATED(R)) DEALLOCATE(R)
      IF (ALLOCATED(YOLD)) DEALLOCATE(YOLD)
      IF (ALLOCATED(YP)) DEALLOCATE(YP)
      IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
      ALLOCATE(Q(NP1,NP1),R(NP1*(N+2)/2),YOLD(NP1),YP(NP1),YPOLD(NP1),
     &  STAT=JW)
      IF (JW /= 0) THEN
        IFLAG=8
        RETURN
      END IF
C 
C SET INITIAL CONDITIONS FOR FIRST CALL TO STEPQF.
C
      START=.TRUE.
      CRASH=.FALSE.
      RELERR = ARCRE
      ABSERR = ARCAE
      HOLD=1.0
      H=0.1
      S=0.0
      YPOLD(1) = 1.0
      Y(1) = 0.0
      YPOLD(2:NP1)=0.0
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
      IF (IFLAGC .GE. -1) THEN
        A=Y(2:NP1)
      ENDIF
C
40    LIMIT=LIMITD
C
C ***** END OF INITIALIZATION BLOCK. *****
C
50    DO ITER=1,LIMIT   ! ***** MAIN LOOP. *****
      IF (Y(1) .LT. 0.0) THEN
        ARCLEN = S
        IFLAG = 5
        CALL CLEANUP ; RETURN
      END IF
C
C TAKE A STEP ALONG THE CURVE.
C
      CALL STEPQF(N,NFE,IFLAGC,START,CRASH,HOLD,H,WK,
     &    RELERR,ABSERR,S,Y,YP,YOLD,YPOLD,A,Q,R,F0,F1,Z0,DZ,
     &    W,T,SSPAR) 
C
C PRINT LATEST POINT ON CURVE IF REQUESTED.
C
      IF (TRACE .GT. 0) THEN
         WRITE (TRACE,217) ITER,NFE,S,Y(1),(Y(JW),JW=2,NP1)
 217     FORMAT(/' STEP',I5,3X,'NFE =',I5,3X,'ARC LENGTH =',F9.4,3X,
     &   'LAMBDA =',F7.4,5X,'X VECTOR:'/(1X,6ES12.4))
      ENDIF
C
C CHECK IF THE STEP WAS SUCCESSFUL.
C
      IF (IFLAGC .GT. 0) THEN
        ARCLEN=S
        IFLAG=IFLAGC
        CALL CLEANUP ; RETURN
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
        END IF
        IF (ARCAE .LT. ABSERR) ARCAE=ABSERR
C
C         CHANGE LIMIT ON NUMBER OF ITERATIONS.
C
        LIMIT = LIMIT - ITER
        RETURN
      END IF
C
C IF LAMBDA >= 1.0, USE ROOTQF TO FIND SOLUTION.
C
      IF (Y(1) .GE. 1.0) GOTO 500
C
      END DO   ! ***** END OF MAIN LOOP *****
C
C DID NOT CONVERGE IN  LIMIT  ITERATIONS, SET  IFLAG  AND RETURN.
C
      ARCLEN = S
      IFLAG = 3
      RETURN
C
C ***** FINAL STEP -- FIND SOLUTION AT LAMBDA=1 *****
C
C SAVE  YOLD  FOR ARC LENGTH CALCULATION LATER.
C
 500  YSAV=YOLD
C
C FIND SOLUTION.
C
      CALL ROOTQF(N,NFE,IFLAGC,ANSRE,ANSAE,Y,YP,YOLD,
     &    YPOLD,A,Q,R,DZ,Z0,W,T,F0,F1)
C
C CHECK IF SOLUTION WAS FOUND AND SET  IFLAG  ACCORDINGLY.
C
      IFLAG=1
C
C     SET ERROR FLAG IF ROOTQF COULD NOT GET THE POINT ON THE ZERO
C     CURVE AT  LAMBDA = 1.0.
C
      IF (IFLAGC .GT. 0) IFLAG=IFLAGC
C
C CALCULATE FINAL ARC LENGTH.
C
      DZ = Y - YSAV
      ARCLEN = S - HOLD + DNRM2(NP1,DZ,1)
C
C ***** END OF FINAL STEP *****
C
      CALL CLEANUP ; RETURN
C
      CONTAINS
        SUBROUTINE CLEANUP
        IF (ALLOCATED(Q)) DEALLOCATE(Q)
        IF (ALLOCATED(R)) DEALLOCATE(R)
        IF (ALLOCATED(YOLD)) DEALLOCATE(YOLD)
        IF (ALLOCATED(YP)) DEALLOCATE(YP)
        IF (ALLOCATED(YPOLD)) DEALLOCATE(YPOLD)
        END SUBROUTINE CLEANUP
      END SUBROUTINE FIXPQF
!
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
      USE HOMPACK_GLOBAL, ONLY: QR => QRSPARSE, ROWPOS, COLPOS, PP,
     & PAR, IPAR
      USE HOMPACK_CORE, ONLY: ROOTNS, STEPQS
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
!
      SUBROUTINE POLSYS1H(N,NUMT,COEF,KDEG,IFLG1,IFLG2,EPSBIG,EPSSML,
     &     SSPAR,NUMRR,LAMBDA,ROOTS,ARCLEN,NFE)
C
C POLSYS1H finds all (complex) solutions to a system
C F(X)=0 of N polynomial equations in N unknowns
C with real coefficients. If IFLG=10 or IFLG=11, POLSYS1H
C returns the solutions at infinity also.
C
C The system F(X)=0 is described via the coefficents,
C "COEF", and the parameters "N, NUMT, KDEG", as follows.
C
C
C       NUMT(J)
C
C F(J) = SUM  COEF(J,K) * X(1)**KDEG(J,1,K)...X(N)**KDEG(J,N,K)
C
C        K=1
C
C FOR J=1, ..., N.
C
C
C POLSYS1H has two main run options:  automatic scaling and
C the projective transformation.  These are evoked via the
C flag "IFLG1", as described below.  The other input
C parameters are the same whether one or both of these options
C are specified, and the output is always returned unscaled
C and untransformed.
C
C If automatic scaling is specified, then the input
C coefficients are modified by subroutine  SCLGNP . The problem
C is solved with the scaled coefficients and scaled variables.
C The coefficients are returned scaled.
C
C If the projective transformation is specified, then
C essentially the system is reformulated in homogeneous
C coordinates, Z(1), ..., Z(N+1), and solved in complex
C projective space.  The resulting solutions are
C untransformed via
C
C X(J) = Z(J)/Z(N+1)   J=1, ..., N.
C
C On return,
C
C ROOTS(1,J,M) = real part of X(J) for the Mth path,
C
C ROOTS(2,J,M) = imaginary part of X(J) for the Mth path,
C
C for J=1, ..., N, and
C
C ROOTS(1,N+1,M) = real part of Z(N+1) for the Mth path,
C
C ROOTS(2,N+1,M) = imaginary part of Z(N+1) for the Mth path.
C
C If ROOTS(*,N+1,M) is small, then the associated solution
C should be regarded as being "near infinity".  Note that,
C when the projective transformation has been specified, the
C ROOTS values have been untransformed -- that is, divided
C through by Z(N+1) -- unless such division would have caused
C overflow.  In this latter case, the affected components of
C ROOTS are set to the largest floating point number (machine
C infinity).
C
C The code can be modified easily to solve systems with complex
C coefficients,  COEF .  Only the subroutines  INITP  and  FFUNP
C need be changed.
C
C The FORTRAN COMPLEX declaration is not used in POLSYS1H.
C Complex variables are represented by real arrays with first
C index dimensioned 2 and complex operations are evoked by
C subroutine calls.
C
C The total number of paths that will be tracked (if
C IFLG2(M)=-2 for all M) is equal to the "total degree" of the
C system, TOTDG.   TOTDG is equal to the products of the
C degrees of all the equations in the system.  The degree of
C an equation is the maximum of the degrees of its terms.  The
C degree of a term is the sum of the degrees of the variables.
C Thus, TOTDG = IDEG(1) * ... * IDEG(N) where IDEG(J) =
C MAX {JDEG(J,K) | K=1,...,NUMT(J)} where JDEG(J,K) = KDEG(J,1,K) +
C ... + KDEG(J,N,K).
C
C IFLG1  determines whether the system is to be automatically
C scaled by  POLSYS1H  and whether the projective transformation
C of the system is to be automatically evoked by POLSYS1H.  See
c "ON INPUT" below.
C
C IFLG2, EPSBIG, EPSSML, and  SSPAR  tell the path tracker
C FIXPNF  which paths to track and set parameters for the path
C tracker.
C
C NUMRR  tells  POLSYS1H  how many multiples of 1000 steps to try
C before abandoning a path.
C
C The output consists of  IFLG1, and of  LAMBDA, ROOTS, ARCLEN, and
C NFE  for each path.  IFLG1  returns input data error information.
C ROOTS  gives the solutions themselves, while  LAMBDA, ARCLEN,
C and  NFE  give information about the associated paths.
C
C
C The following subroutines are used directly or indirectly by
C POLSYS1H: 
C         Special for POLSYS1H:
C           INITP , STRPTP , OTPUTP , RHO , RHOJAC ,
C           HFUNP , HFUN1P , GFUNP , FFUNP ,
C           MULP , POWP , DIVP , SCLGNP .
C         From the general HOMPACK routines:
C           FIXPNF , ROOT , ROOTNF , STEPNF , TANGNF .
C         From LAPACK routines:
C           DGEQPF , DGEQRF , DORMQR .
C         From BLAS routines:
C           DCOPY ,  DDOT ,  DGEMM ,  DGEMV ,  DGER ,  
C           DNRM2 ,  DSCAL ,  DSWAP ,  DTRMM ,  DTRMV , DTRSV ,
C           IDAMAX ,  LSAME , XERBLA . 
C
C ON INPUT:
C
C N  is the number of equations and variables.
C
C NUMT(1:N)  is an integer array.  NUMT(J)  is the number of terms
C   in the Jth equation for J=1 to N.
C
C COEF(1:N,1:)  is a real array.  COEF(J,K)  is 
C   the Kth coefficient of the Jth equation for J=1 to N,
C   K=1 to NUMT(J).  The second dimension must be greater than or equal
C   to the maximum number of terms in each equation.  In other words,
C   SIZE(COEF,DIM=2) .GE. MAXT = MAX {NUMT(J) | J=1, ..., N} .
C
C KDEG(1:N,1:N+1,1:)  is an integer array.  
C   KDEG(J,L,K)  is the degree of the Lth variable in the Kth
C   term of the Jth equation for  J=1 to N, L=1 to N, K=1 to NUMT(J).
C   SIZE(KDEG,DIM=3) .GE. MAXT = MAX {NUMT(J) | J=1, ..., N} .
C
C IFLG1 =
C   00  if the problem is to be solved without
C       calling POLSYS1H' scaling routine, SCLGNP, and
C       without using the projective transformtion.
C
C   01  if scaling but no projective transformation is to be used.
C
C   10  if no scaling but projective transformation is to be used.
C
C   11  if both scaling and projective transformation are to be used.
C
C IFLG2(1:TOTDG)  is an integer array.  If IFLG2(M) = -2, then the 
C   Mth path is tracked.  Otherwise the Mth path is skipped.
C   Thus, to find all solutions set IFLG2(M) = -2 for M=1,...,TOTDG.
C   Selected paths can be rerun by setting IFLG2(M)=-2 for
C   the paths to be rerun and IFLG2(M).NE.-2 for the others.
C
C EPSBIG  is the local error tolerance allowed the path tracker along
C   the path.  ARCRE and ARCAE (in  FIXPNF ) are set to  EPSBIG.
C
C EPSSML  is the accuracy desired for the final solution.  ANSRE and
C   ANSAE (in  FIXPNF ) are set to  EPSSML.
C
C SSPAR(1:8) = (LIDEAL, RIDEAL, DIDEAL, HMIN, HMAX, BMIN, BMAX, P)  is
C    a vector of parameters used for the optimal step size estimation.
C    If  SSPAR(J) .LE. 0.0  on input, it is reset to a default value
C    by  FIXPNF .  Otherwise the input value of  SSPAR(J)  is used.
C    See the comments in  FIXPNF  and in  STEPNF  for more information
C    about these constants.
C
C NUMRR  is the number of multiples of 1000 steps that will be tried
C   before abandoning a path.
C
C
C ON OUTPUT:
C
C N, NUMT, COEF, KDEG, EPSBIG, EPSSML, and NUMRR are unchanged.
C
C IFLG1=
C   -1  if  NUMT  is incorrectly dimensioned or invalid.
C   -2  if  COEF  is incorrectly dimensioned.
C   -3  if  KDEG  is incorrectly dimensioned or invalid.
C   -4  if any of  IFLG2, LAMBDA, ROOTS, ARCLEN, or  NFE  are
C       incorrectly dimensioned.
C   -5  if the global work arrays  IPAR  and  PAR  could not be
C       allocated.
C   -6  if  IFLG1  on input is not 00 or 01 or 10 or 11.
C   Unchanged otherwise.
C
C IFLG2(1:TOTDG)  gives information about how the Mth path terminated:
C IFLG2(M) =
C   1   Normal return.
C
C   2   Specified error tolerance cannot be met.  Increase  EPSBIG
C       and  EPSSML  and rerun.
C
C   3   Maximum number of steps exceeded.  To track the path further,
C       increase  NUMRR  and rerun the path.  However, the path may
C       be diverging, if the  LAMBDA  value is near 1 and the  ROOTS 
C       values are large.
C
C   4   Jacobian matrix does not have full rank.  The algorithm
C       has failed (the zero curve of the homotopy map cannot be
C       followed any further).
C
C   5   The tracking algorithm has lost the zero curve of the
C       homotopy map and is not making progress.  The error tolerances
C       EPSBIG  and  EPSSML  were too lenient.  The problem should be
C       restarted with smaller error tolerances.
C
C   6   The normal flow Newton iteration in  STEPNF  or  ROOTNF
C       failed to converge.  The error tolerances  EPSBIG  or  EPSSML
C       may be too stringent.
C
C   7   Illegal input parameters, a fatal error.
C
C LAMBDA(M)  is the final LAMBDA value for the Mth path, M = 1, ...,
C   TOTDG, where LAMBDA is the continuation parameter.
C
C ROOTS(1,J,M), ROOTS(2,J,M)  are the real and imaginary parts
C   of the Jth variable respectively, for J = 1,...,N, for
C   the Mth path, for M = 1,...,TOTDG.  If  IFLG1 = 10 or 11, then
C   ROOTS(1,N+1,M)  and  ROOTS(2,N+1,M)  are the real and
C   imaginary parts respectively of the projective
C   coordinate of the solution.
C
C ARCLEN(M)  is the arc length of the Mth path for M = 1, ..., TOTDG.
C
C NFE(M)  is the number of Jacobian matrix evaluations required to 
C   track the Mth path for M =1, ..., TOTDG.
C
C ----------------------------------------------------------------------
      USE HOMPACK_GLOBAL, ONLY: IPAR, PAR
      USE HOMPACK_CORE, ONLY: INITP, STRPTP, OTPUTP
      IMPLICIT NONE
C
C TYPE DECLARATIONS FOR INPUT AND OUTPUT
C
      INTEGER, INTENT(IN):: N,NUMT(:),NUMRR
      REAL (dp), INTENT(IN OUT):: COEF(:,:),SSPAR(8)
      INTEGER, INTENT(IN OUT):: KDEG(:,:,:),IFLG1,IFLG2(:)
      REAL (dp), INTENT(IN):: EPSBIG,EPSSML
      REAL (dp), INTENT(OUT):: LAMBDA(:),ROOTS(:,:,:),ARCLEN(:)
      INTEGER, INTENT(OUT):: NFE(:)
C
C TYPE DECLARATIONS FOR LOCAL VARIABLES
C
      INTEGER:: I,ICOUNT(N),IDEG(N),IDUMMY,IFLAG,IJ,
     &  IPROFF(15),J,LIPAR(15),LPAR(25),MAXT,N2,N2P1,
     &  NNFE,NP1,NUMPAT,PROFF(25),TOTDG,TRACE
      REAL (dp):: AARCLN,ANSAE,ANSRE,ARCAE,ARCRE,CL(2,N+1),
     &  FACV(N),PDG(2,N),QDG(2*N),R(2,N),XNP1(2),Y(2*N+1)
C
C ----------------------------------------------------------------------
      N2=2*N
      NP1=N+1
      N2P1=N2+1
C
C CHECK THAT DIMENSIONS ARE VALID.
C
      IF ((SIZE(NUMT) /= N) .OR. ANY(NUMT .LE. 0)) THEN
        IFLG1=-1
        RETURN
      END IF
      MAXT = MAXVAL(NUMT)
      IF ((SIZE(COEF,DIM=1) /= N) .OR. (SIZE(COEF,DIM=2) < MAXT)) THEN
        IFLG1=-2
        RETURN
      END IF
      KDEG = ABS(KDEG)
      IF ((SIZE(KDEG,DIM=1) /= N) .OR. (SIZE(KDEG,DIM=2) /= NP1) .OR.
     &  (SIZE(KDEG,DIM=3) < MAXT) ) THEN
        IFLG1=-3
        RETURN
      END IF
      DO J=1,N
        IDEG(J)=MAXVAL(SUM(KDEG(J,1:N,1:NUMT(J)),DIM=1))
      END DO
      TOTDG = PRODUCT(IDEG)
      IF ((SIZE(IFLG2) < TOTDG) .OR. (SIZE(LAMBDA) < TOTDG) .OR.
     &  (SIZE(ROOTS,DIM=3) < TOTDG) .OR. (SIZE(ARCLEN) < TOTDG) .OR.
     &  (SIZE(NFE) < TOTDG) .OR. 
     &  (IFLG1 <= 1 .AND. SIZE(ROOTS,DIM=2) /= N) .OR.
     &  (IFLG1 >= 10 .AND. SIZE(ROOTS,DIM=2) /= NP1)) THEN
        IFLG1=-4
        RETURN
      END IF
      IF (IFLG1 /= 0 .AND. IFLG1 /= 1 .AND.
     &  IFLG1 /= 10 .AND. IFLG1 /= 11) THEN
        IFLG1=-6
        RETURN
      END IF
C
C ALLOCATE THE GLOBAL WORK ARRAYS  IPAR  AND  PAR, USED TO COMMUNICATE
C DATA BETWEEN SUBROUTINES VIA THE MODULE HOMOTOPY.
C
      ALLOCATE(IPAR(42 + 2*N + N*(N+1)*MAXT),
     &  PAR(2 + 28*N + 6*N**2 + 7*N*MAXT + 4*N**2*MAXT),STAT=IJ)
      IF (IJ .NE. 0) THEN
        IFLG1=-5
        RETURN
      END IF
C      
C INITIALIZATION
C
      CALL INITP(IFLG1,N,NUMT,KDEG,COEF,
     &                              IDEG,FACV,CL,PDG,QDG,R)
C
C INTEGER VARIABLES AND ARRAYS TO BE PASSED IN IPAR:
C
C    IPAR INDEX     VARIABLE NAME       LENGTH
C    ----------     -------------    -----------------
C          1                N               1
C          2             MAXT               1
C          3            PROFF               25
C          4           IPROFF               15
C          5             IDEG               N
C          6             NUMT               N
C          7             KDEG               N*(N+1)*MAXT
C
C
C DOUBLE PRECISION VARIABLES AND ARRAYS TO BE PASSED IN PAR:
C
C     PAR INDEX     VARIABLE NAME       LENGTH
C    ----------     -------------    -----------------
C          1              PDG               2*N
C          2               CL               2*(N+1)
C          3             COEF               N*MAXT
C          4                H               N2
C          5              DHX               N2*N2
C          6              DHT               N2
C          7            XDGM1               2*N
C          8              XDG               2*N
C          9              G                 2*N
C         10             DG                 2*N
C         11           PXDGM1               2*N
C         12             PXDG               2*N
C         13               F                2*N
C         14              DF                2*N*(N+1)
C         15               XX               2*N*(N+1)*MAXT
C         16              TRM               2*N*MAXT
C         17             DTRM               2*N*(N+1)*MAXT
C         18              CLX               2*N
C         19            DXNP1               2*N
C
C SET LENGTHS OF VARIABLES
      LIPAR(1)=1
      LIPAR(2)=1
      LIPAR(3)=25
      LIPAR(4)=15
      LIPAR(5)=N
      LIPAR(6)=N
      LIPAR(7)=N*(N+1)*MAXT
      LPAR( 1)=2*N
      LPAR( 2)=2*NP1
      LPAR( 3)=N*MAXT
      LPAR( 4)=N2
      LPAR( 5)=N2*N2
      LPAR( 6)=N2
      LPAR( 7)=2*N
      LPAR( 8)=2*N
      LPAR( 9)=2*N
      LPAR(10)=2*N
      LPAR(11)=2*N
      LPAR(12)=2*N
      LPAR(13)=2*N
      LPAR(14)=2*N*NP1
      LPAR(15)=2*N*NP1*MAXT
      LPAR(16)=2*N*MAXT
      LPAR(17)=2*N*NP1*MAXT
      LPAR(18)=2*N
      LPAR(19)=2*N
C
C PROFF AND IPROFF ARE OFFSETS THAT DEFINE THE VARIABLES LISTED ABOVE
      PROFF(1)=1
      DO I=2,19
          PROFF(I)=PROFF(I-1)+LPAR(I-1)
      END DO
      IPROFF(1)=1
      DO I=2,7
          IPROFF(I)=IPROFF(I-1)+LIPAR(I-1)
      END DO
C
C DEFINE VARIABLES
      IPAR(1)=N
      IPAR(2)=MAXT
      IPAR(IPROFF(3):IPROFF(3)+18) = PROFF(1:19)
      IPAR(IPROFF(4):IPROFF(4)+ 6) = IPROFF(1:7)
      IPAR(IPROFF(5):IPROFF(5)+N-1) = IDEG(1:N)
      IPAR(IPROFF(6):IPROFF(6)+N-1) = NUMT(1:N)
      IPAR(IPROFF(7):IPROFF(7)+LIPAR(7)-1) =
     &  PACK(KDEG(:,:,1:MAXT),.TRUE.)
      PAR(PROFF(1):PROFF(1)+LPAR(1)-1) = PACK(PDG,.TRUE.)
      PAR(PROFF(2):PROFF(2)+LPAR(2)-1) = PACK(CL,.TRUE.)
      PAR(PROFF(3):PROFF(3)+LPAR(3)-1) = PACK(COEF(:,1:MAXT),.TRUE.)
C
C ICOUNT IS A COUNTER USED BY "STRPTP"
      ICOUNT(1)=0
      ICOUNT(2:N)=1
C
C PATHS LOOP -- ITERATE THROUGH PATHS
C
      PATHS: DO NUMPAT = 1,TOTDG
C         GET A START POINT, Y, FOR THE PATH.
          Y(1) = 0.0
          CALL STRPTP(N,ICOUNT,IDEG,R,Y(2:N2P1))
C         CHECK WHETHER PATH IS TO BE FOLLOWED.
          IFLAG = IFLG2(NUMPAT)
          IF (IFLAG .NE. -2) CYCLE PATHS
          ARCRE = EPSBIG
          ARCAE = ARCRE
          ANSRE = EPSSML
          ANSAE = ANSRE
          TRACE = 0
C         TRACK A HOMOTOPY PATH.
          DO IDUMMY=1,MAX(NUMRR,1)
            CALL FIXPNF(N2,Y,IFLAG,ARCRE,ARCAE,ANSRE,ANSAE,TRACE,
     &        QDG,SSPAR,NNFE,AARCLN, POLY_SWITCH=.TRUE.)
            IF (IFLAG .NE. 2 .AND. IFLAG .NE. 3) EXIT
          END DO
C         UNSCALE AND UNTRANSFORM COMPUTED SOLUTION.
          CALL OTPUTP(N,NUMPAT,CL,FACV,
     &      PAR(PROFF(18):PROFF(18)+LPAR(18)-1),Y(2:N2P1),XNP1)
          LAMBDA(NUMPAT) = Y(1)
          ROOTS(1,1:N,NUMPAT) = Y(2:N2P1:2)
          ROOTS(2,1:N,NUMPAT) = Y(3:N2P1:2)
          ROOTS(1:2,NP1,NUMPAT) = XNP1
C
          ARCLEN(NUMPAT)= AARCLN
          NFE(NUMPAT)   = NNFE
          IFLG2(NUMPAT) = IFLAG
      END DO PATHS
C CLEAN UP WORK SPACE.
      IF (ALLOCATED(IPAR)) DEALLOCATE(IPAR)
      IF (ALLOCATED(PAR))  DEALLOCATE(PAR)
      RETURN
      END SUBROUTINE POLSYS1H
C
      END MODULE HOMPACK