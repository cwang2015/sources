cnea    PROGRAM SOLAVOF
C      PROGRAM SOLAVOF (INPUT,TAPE5=INPUT,OUTPUT,TAPE6=OUTPUT,TAPE7,TAPE8
C     1 )
C
C     *** SOLA-VOF  VOLUME OF FLUID METHOD ***
C
C     *** LIST OF PRIMARY VARIABLES ***
C
C     *** INPUT PARAMETERS (NAMELIST /XPUT/)
C
C     ALPHA     CONTROLS AMOUNT OF DONOR CELL FLUXING (=1.0 FOR FULL
C               DONOR CELL DIFFERENCING, =0.0 FOR CENTRAL DIFFERENCING)
C     AUTOT     AUTOMATIC TIME STEP FLAG (=1.0 FOR AUTOMATIC DELT
C               ADJUSTMENT, =0.0 FOR CONSTANT DELT)
C     CANGLE    CONTACT ANGLE, IN DEGRES, BETWEEN FLUID AND WALL
C     CSQ       MATERIAL SOUND SPEED SQUARED (=-1.0 FOR
C               INCOMPRESSIBLE MATERIAL)
C     DELT      TIME STEP
C     EPSI      PRESSURE ITERATION CONVERGENCE CRITERION
C     FLHT      FLUID HEIGHT, IN Y-DIRECTION
C     GX        BODY ACCELERATION IN POSITIVE X-DIRECTION
C     GY        BODY ACCELERATION IN POSITIVE Y-DIRECTION
C     ICYL      MESH GEOMETRY INDICATOR (=1 FOR CYLINDRICAL COORDINATES,
C               =0 FOR PLANE COORDINATES)
C     IMOVY     MOVIE INDICATOR (=1 FOR MOVIE FILM OUTPUT, =0 FOR
C               OTHER FILM OUTPUT)
C     ISURF10   SURFACE TENSION INDICATOR (=1 FOR SURFACE TENSION,
C               =0 FOR NO SURFACE TENSION)
C     ISYMPLT   SYMMETRY PLOT INDICATOR (=1 FOR SYMMETRY PLOT,
C               =0 FOR NO SYMMETRY PLOT)
C     NMAT      NUMBER OF MATERIALS
C     NPX       NUMBER OF PARTICLES IN X-DIRECTION IN RECTANGULAR SETUP
C     NPY       NUMBER OF PARTICLES IN Y-DIRECTION IN RECTANGULAR SETUP
C     NU        COEFFICIENT OF KINEMATIC VISCOSITY
C     OMG       OVER-RELAXATION FACTOR USED IN PRESSURE ITERATION
C     PLTDT     TIME INCREMENT BETWEEN PLOTS AND/OR PRINTS TO BE OUTPUT
C               ON FILM
C     PRTDT     TIME INCREMENT BETWEEN PRINTS ON PAPER
C     RHOF      FLUID DENSITY (FOR F=1.0 REGION)
C     RHOFC     FLUID DENSITY IN COMPLEMENT OF F REGION
C     SIGMA     SURFACE TENSION COEFFICIENT
C     TWFIN     PROBLEM TIME TO END CALCULATION
C     UI        X-DIRECTION VELOCITY USED FOR INITIALIZING MESH
C     VI        Y-DIRECTION VELOCITY USED FOR INITIALIZING MESH
C     VELMX     MAXIMUM VELOCITY EXPECTED IN PROBLEM, USED TO SCALE
C               VELOCITY VECTORS
C     WB        INDICATOR FOR BOUNDARY CONDITION TO BE USED ALONG THE
C               BOTTOM OF THE MESH (=1 FOR RIGID FREE-SLIP WALL,
C               =2 FOR RIGID NO-SLIP WALL, =3 FOR CONTINUATIVE
C               BOUNDARY, =4 FOR PERIODIC BOUNDARY, =5 FOR CONSTANT
C               PRESSURE BOUNDARY)
C     WL        INDICATOR FOR BOUNDARY CONDITION ALONG LEFT SIDE OF
C               MESH (SEE WB)
C     WR        INDICATOR FOR BOUNDARY CONDITION ALONG RIGHT SIDE OF
C               MESH (SEE WB)
C     WT        INDICATOR FOR BOUNDARY CONDITION ALONG TOP OF MESH
C               (SEE WB)
C     XPL       LOCATION OF LEFT SIDE OF RECTANGULAR PARTICLE REGION
C     XPR       LOCATION OF RIGHT SIDE OF RECTANGULAR PARTICLE REGION
C     YPB       LOCATION OF BOTTOM OF RECTANGULAR PARTICLE REGION
C     YPT       LOCATION OF TOP OF RECTANGULAR PARTICLE REGION
C
C     *** MESH SETUP INPUT (NAMELIST /MSHSET/)
C
C     DXMN(N)     MINIMUM SPACE INCREMENT IN X-DIRECTION IN SUBMESH N
C     DYMN(N)     MINIMUM SPACE INCREMENT IN Y-DIRECTION IN SUBMESH N
C     NKX         NUMBER OF SUBMESH REGIONS IN X-DIRECTION
C     NXL(N)      NUMBER OF CELLS BETWEEN LOCATIONS XL(N) AND XC(N) IN
C                 SUBMESH N
C     NXR(N)      NUMBER OF CELLS BETWEEN LOCATIONS XC(N) AND XL(N+1) IN
C                 SUBMESH N
C     NYL(N)      NUMBER OF CELLS BETWEEN LOCATIONS YL(N) AND YC(N) IN
C                 SUBMESH N
C     NYR(N)      NUMBER OF CELLS BETWEEN LOCATIONS YC(N) AND YL(N+1) IN
C                 SUBMESH N
C     XC(N)       X-COORDINATE OF THE CONVERGENCE POINT IN SUBMESH N
C     XL(N)       LOCATION OF THE LEFT EDGE OF SUBMESH N (NKX+1 VALUES
C                 OF XL(N) ARE NECESSARY BECAUSE THE RIGHT EDGE (XR) OF
C                 SUBMESH N IS DETERMINED BY THE LEFT EDGE OF
C                 SUBMESH N+1)
C     YC(N)       Y-COORDINATE OF THE CONVERGENCE POINT IN SUBMESH N
C     YL(N)       LOCATION OF THE BOTTOM OF SUBMESH N (NKY+1 VALUES OF
C                 YL(N) ARE NECESSARY BECAUSE THE TOP EDGE (YR) OF
C                 SUBMESH N IS DETERMINED BY THE BOTTOM EDGE OF
C                 SUBMESH N+1)
C
C     *** VARIABLES LISTED IN COMMON (EXCLUDING INPUT PARAMETERS)
C
C     CYCLE     CALCULATIONAL TIME CYCLE
C     CYL       MESH GEOMETRY INDICATOR (= ICYL)
C     DTSFT     MAXIMUM DELT VALUE ALLOWED BY THE SURFACE TENSION FORCES
C               STABILITY CRITERION (DELT IS AUTOMATICALLY ADJUSTED)
C     DTVIS     MAXIMUM DELT VALUE ALLOWED BY THE VISCOUS FORCES
C               STABILITY CRITERION (DELT IS AUTOMATICALLY ADJUSTED)
C     EMF       SMALL VALUE, TYPICALLY 1.0E-06, USED TO NEGATE ROUND-OFF
C               ERROR EFFECTS WHEN TESTING F=1.0 OR F=0.0
C     EMF1      =1.0-EMF
C     EM6       =1.0E-06
C     EM10      =1.0E-10
C     EP10      =1.0E+10
C     FLG       PRESSURE ITERATION CONVERGENCE TEST INDICATOR (=0.0 WHEN
C               THE CONVERGENCE TEST IS SATISFIED, =1.0 WHEN THE
C               CONVERGENCE TEST IS NOT SATISFIED)
C     FLGC      VOLUME OF FLUID FUNCTION CONVECTION LIMIT INDICATOR
C               (DELT REDUCED AND CYCLE STARTED OVER IF LIMIT
C               IS EXCEEDED)
C     FNOC      PRESSURE CONVERGENCE FAILURE INDICATOR (=1.0,
C               CONVERGENCE FAILED AND DELT IS REDUCED, =0.0 OTHERWISE)
C     IBAR      NUMBER OF REAL CELLS IN X-DIRECTION (EXCLUDES FICTICIOUS
C               CELLS)
C     IBAR2     =IBAR+2, SPECIFIED IN PARAMETER STATEMENT
C               (=IBAR+3, IF PERIODIC IN X-DIRECTION)
C     IMAX      TOTAL NUMBER OF MESH CELLS IN X-DIRECTION (=IBAR+2)
C               (=IBAR+3, IF PERIODIC IN X-DIRECTION)
C     IM1       VALUE OF THE INDEX I AT THE LAST REAL CELL IN THE
C               X-DIRECTION (=IMAX-1)
C     IM2       VALUE OF THE INDEX I AT THE NEXT TO THE LAST REAL CELL
C               IN THE X-DIRECTION (=IMAIN THE X-2)
C     IPL       LEFTMOST PRESSURE ITERATION INDEX IN X-DIRECTION
C               (=3 FOR CONSTANT PRESSURE BOUNDARY CONDITION, =2 FOR
C               ALL OTHER CASES)
C     IPR       RIGHTMOST PRESSURE ITERATION INDEX IN X-DIRECTION
C               (=IM2 FOR CONSTANT PRESSURE BOUNDARY CONDITION, =IM1 FOR
C               ALL OTHER CASES)
C     ITER      PRESSURE ITERATION COUNTER
C     JBAR      NUMBER OF REAL CELLS IN Y-DIRECTION (EXCLUDES FICTICIOUS
C               CELLS)
C     JBAR2     =JBAR+2, SPECIFIED IN PARAMETER STATEMENT
C               (=JBAR+3, IF PERIODIC IN Y-DIRECTION)
C     JMAX      TOTAL NUMBER OF MESH CELLS IN Y-DIRECTION (=JBAR+2)
C               (=JBAR+3, IF PERIODIC IN Y-DIRECTION)
C     JM1       VALUE OF THE INDEX J AT THE LAST REAL CELL IN THE
C               Y-DIRECTION (=JMAX-1)
C     JM2       VALUE OF THE INDEX J AT THE NEXT TO THE LAST REAL CELL
C               IN THE Y-DIRECTION (=JMAX-2)
C     JPB       BOTTOM PRESSURE ITERATION INDEX IN Y-DIRECTION
C               (=3 FOR CONSTANT PRESSURE BOUNDARY CONDITION, =2 FOR
C               ALL OTHER CASES)
C     JPT       TOP PRESSURE ITERATION INDEX IN Y-DIRECTION
C               (=JM2 FOR CONSTANT PRESSURE BOUNDARY CONDITION, =JM1 FOR
C               ALL OTHER CASES)
C     NFLGC     NUMBER OF CYCLES THE VOLUME OF FLUID FUNCTION CONVECTION
C               LIMIT (FLGC) IS EXCEEDED
C     NOCON     NUMBER OF CYCLES PRESSURE CONVEGENCE HAS FAILED
C     NP        TOTAL NUMBER OF PARTICLES COMPUTED TO BE IN MESH
C     NPRTS     NUMBER OF PARTICLES IN MESH, SPECFIED IN PARAMETER STATE
C               (USED TO SET ARRAY SIZE - MUST BE > 0)
C     NREG      NUMBER OF VOID REGIONS GENERATED IN CALCULATION
C     NVOR      MAXIMUM NUMBER OF VOID REGIONS ALLOWED, SPECIFIED IN
C               PARAMETER STATEMENT
C     NVRM      NUMBER OF VOID REGIONS
C     MESHX     NUMBER OF SUBMESH REGIONS IN X-DIRECTION, SPECIFIED
C               IN PARAMETER STATEMENT
C     MESHY     NUMBER OF SUBMESH REGIONS IN Y-DIRECTION, SPECIFIED
C               IN PARAMETER STATEMENT
C     PI        =3.141592654
C     RCSQ      RECIPROCAL OF CSQ
C     RHOD      DIFFERENCE IN FLUID DENSITIES (=RHOF-RHOFC)
C     RPD       DEGREES TO RADIANS CONVERSION FACTOR
C     SF        PLOT SCALING FACTOR
C     T         PROBLEM TIME
C     TANGLE    TANGENT OF CONTACT ANGLE, CANGLE
C     VCHGT     ACCUMULATED FLUID VOLUME CHANGE
C     VELMX1    VELMX NORMALIZED TO MINIMUM MESH CELL DIMENSION
C     XMAX      LOCATION OF RIGHT-HAND SIDE OF MESH
C     XMIN      LOCATION OF LEFT-HAND SIDE OF MESH
C     XSHFT     COMPUTED SHIFT ALONG THE PLOTTING ABSCISSA TO CENTER
C               THE PLOT FRAME ON FILM
C     YMAX      LOCATION OF THE TOP OF THE MESH
C     YMIN      LOCATION OF THE BOTTOM OF THE MESH
C     YSHFT     COMPUTED SHIFT ALONG THE PLOTTING ORDINATE TO CENTER THE
C               PLOT FRAME ON FILM
C
C     *** ARRAYS IN COMMON (EXCLUDING MESH SETUP PARAMETERS)
C
C     ACOM(1)     FIRST WORD IN COMMON
C     BETA(I,J)   PRESSURE ITERATION RELAXATION FACTOR IN CELL (I,J)
C     DELX(I)     MESH SPACING OF THE I-TH CELL ALONG THE X-AXIS
C     DELY(J)     MESH SPACING OF THE J-TH CELL ALONG THE Y-AXIS
C     F(I,J)      VOLUME OF FLUID PER UNIT VOLUME OF CELL (I,J) AT TIME
C                 LEVEL N+1
C     FN(I,J)     VOLUME OF FLUID PER UNIT VOLUME OF CELL (I,J) AT TIME
C                 LEVEL N
C     IP(K)       CELL INDEX FOR PARTICLE K ALONG X-AXIS
C     JP(K)       CELL INDEX FOR PARTICLE K ALONG Y-AXIS
C     NAME(10)    PROBLEM IDENTIFICATION LINE
C     NF(I,J)     FLAG OF SURFACE CELL (I,J) INDICATING THE LOCATION
C                 OF ITS NEIGHBORING PRESSURE INTERPOLATION CELL
C     NR(K)       LABEL OF VOID REGION, K > 5
C     P(I,J)      PRESSURE IN CELL (I,J) AT TIME LEVEL N+1
C     PETA(I,J)   PRESSURE INTERPOLATION FACTOR FOR CELL (I,J)
C     PN(I,J)     PRESSURE IN CELL (I,J) AT TIME LEVEL N
C     PR(K)       PRESSURE IN VOID REGION NR(K)
C     PS(I,J)     SURFACE PRESSURE IN CELL (I,J) COMPUTED FROM SURFACE
C                 TENSION FORCES
C     RDX(I)      RECIPROCAL OF DELX(I)
C     RDY(J)      RECIPROCAL OF DELY(J)
C     RX(I)       RECIPROCAL OF X(I)
C     RXI(I)      RECIPROCAL OF XI(I)
C     RYJ(J)      RECIPROCAL OF YJ(J)
C     TANTH(I,J)  SLOPE OF FLUID SURFACE IN CELL (I,J)
C     U(I,J)      X-DIRECTION VELOCITY COMPONENT IN CELL (I,J) AT TIME
C                 LEVEL N+1
C     UN(I,J)     X-DIRECTION VELOCITY COMPONENT IN CELL (I,J) AT TIME
C                 LEVEL N
C     V(I,J)      Y-DIRECTION VELOCITY COMPONENT IN CELL (I,J) AT TIME
C                 LEVEL N+1
C     VN(I,J)     Y-DIRECTION VELOCITY COMPONENT IN CELL (I,J) AT TIME
C                 LEVEL N
C     VOL(K)      VOLUME OF VOID REGION NR(K)
C     X(I)        LOCATION OF THE RIGHT-HAND BOUNDARY OF THE I-TH CELL
C                 ALONG THE X-AXIS
C     XI(I)       LOCATION OF THE CENTER OF THE I-TH CELL ALONG THE
C                 X-AXIS
C     XP(K)       X-COORDINATE OF PARTICLE K
C     Y(J)        LOCATION OF THE TOP BOUNDARY OF THE J-TH CELL ALONG TH
C                 Y-AXIS
C     YJ(J)       LOCATION OF THE CENTER OF THE J-TH CELL ALONG THE
C                 Y-AXIS
C     YP(K)       Y-COORDINATE OF PARTICLE K
C     ZCOM(1)     LAST WORD IN COMMON

      subroutine sola
c      implicit none
      include "include.f"
      NAMELIST /XPUT/ DELT,NU,ICYL,EPSI,GX,GY,UI,VI,VELMX,TWFIN,PRTDT
     1 ,PLTDT,OMG,ALPHA,WL,WR,WT,WB,IMOVY,AUTOT,FLHT,ISYMPLT,SIGMA
     2 ,ISURF10,CANGLE,CSQ,NMAT,RHOF,RHOFC,XPL,XPR,YPB,YPT,NPX,NPY
      NAMELIST /MSHSET/ NKX,XL,XC,NXL,NXR,DXMN,NKY,YL,YC,NYL,NYR,DYMN
C
      DATA EMF /1.0E-06/, EM6 /1.0E-06/, EM10 /1.0E-10/
      DATA EP10 /1.0E+10/
      DATA PI /3.141592654/, RPD /0.0174532925/
C
C     *** DEFAULT INPUT DATA
C     *** NOTE    USER MUST SUPPLY THE FOLLOWING REGARDLESS
C                      OF DEFAULTS: DELT,TWFIN,PRTDT,PLTDT
C
      DATA NU /0.0/, ICYL /0/, EPSI /1.0E-03/, GX /0.0/, GY /0.0/, UI /0
     1 .0/, VI /0.0/, VELMX /1.0/, IMOVY /0/, OMG /1.7/, ALPHA /1.0/, WL
     2 /1/, WR /1/, WT /1/, WB /1/, CSQ /-1.0/, AUTOT /1.0/, ISYMPLT /0/
     3 , ISURF10 /0/, SIGMA /0.0/, CANGLE /90.0/, NMAT /1/, RHOF /1.0/,
     4 RHOFC /1.0/, FLHT /0.0/, XPL /0.0/, YPB /0.0/, XPR /0.0/, YPT /0.
     5 0/, NPX /0/, NPY /0/
C
C     *** SETUP FILM INITIALIZATION TO SYSTEM
C     *** NOTE    FILMSET IS SYSTEM DEPENDANT
C
      open(101,file='D:/LiuSPH/data/vof_dam_break.dat')
      open(102,file='D:/LiuSPH/data/result.dat')
C      CALL FILMSET
C      write(*,*) 'filmset done!'
C
C     *** READ PROBLEM TITLE (NAME)
C
C      READ (5,110) NAME
C      write(*,*) 'read name done!'
C
C     *** READ INITIAL INPUT DATA
C
cwang      READ (5,XPUT)
      read(101,xput)
      write(*,*) 'delta=', delt
      write(*,*) 'read xput done!'
C
C     *** READ INPUT PARAMETERS FOR VARIABLE MESH
C
cwang      READ (5,MSHSET)
      read(101,mshset)
C
C     *** CALCULATE VARIABLE MESH DATA
C
      CALL MESHSET

C
C     *** PRINT INITIAL INPUT DATA
C
      CALL PRTPLT (1)
C
C     *** SET INITIAL CONDITIONS
C
      CALL SETUP
C
C     *** SET INITIAL BOUNDARY CONDITIONS
C
      CALL BC
C
      my_count = 0
      return
      end subroutine


      subroutine vof
      include "include.f"
      call sola
CC      GO TO 20
C
C     *** START TIME CYCLE
C
   10 CONTINUE

      if(my_count==0)then
         my_count = 1
         GOTO 20
      endif

      ITER=0             ! pressure iteration counter
      FLG=1.0            ! pressure convergence indicator
      FNOC=0.0           ! pressure convergence failure indicator
C
C     *** EXPLICITLY APPROXIMATE NEW TIME-LEVEL VELOCITIES
C
      CALL TILDE
C
      IF (NMAT.EQ.2.AND.ISURF10.EQ.1) CALL TMS10    ! surface tension
C
C     *** SET BOUNDARY CONDITIONS
C
      CALL BC
C
C     *** ITERATIVELY ADJUST CELL PRESSURE AND VELOCITY
C
      CALL PRESSIT
      write(*,*) 'adfafd',ITER
C
      IF (T.GT.EP10) GO TO 30  ! EP10= 1.0e+10 
C
C     *** UPDATE FLUID CONFIGURATION
C
   20 CALL VFCONV
C
      IF (FLGC.GT.0.5) GO TO 90    
C
C     *** SET BOUNDARY CONDITIONS
C
      CALL BC
C
C     *** MOVE MARKER PARTICLES
C
      CALL PARMOV
C
C     *** DETERMINE PRESSURE INTERPOLATION FACTOR AND NEIGHBOR
C     *** ALSO DETERMINE SURFACE TENSION PRESSURES AND
C     *** WALL ADHESION EFFECTS IN SURFACE CELLS
C
      CALL PETACAL
C
C     *** PRINT TIME AND CYCLE DATA ON PAPER AND/OR FILM
C
   30 CALL PRTPLT (2)
C
      IF (CYCLE.LE.0) GO TO 40              ! if(cycle.le.0)then
      IF (T+EM6.LT.TWPLT) GO TO 50          !    call prtplot(3)
      TWPLT=TWPLT+PLTDT                         
   40 CONTINUE                              ! elseif(T+EM6.lt.twplt)then
C                                           ! else   
C     *** PRINT FIELD VARIABLE DATA ON FILM !    twplt=twplt+pltdt
C                                           !    call prtplt(3)   
CC      CALL PRTPLT (3)                       ! endif
C                                            
C     *** PLOT VELOCITY VECTOR, FREE SURFACE, MESH,
C     *** AND MARKER PARTICLE ON FILM
C
cwang      CALL DRAW
C
   50 CONTINUE
      IF (CYCLE.LE.0) GO TO 60             ! if(cycle.le.0)then
      IF (T+EM6.LT.TWPRT) GO TO 70         !    call prtlt(4)
      TWPRT=TWPRT+PRTDT                    ! elseif(T+EM6.lt.twprt)then
   60 CONTINUE                             ! else
C                                          !    twprt = twprt+prtdt
C     *** PRINT FIELD VARIABLE DATA ON PAPER!   call prtplt(4)
C                                          ! endif
CC      CALL PRTPLT (4)
C
      
      write(102,*) 'VARIABLES="x","y","u","v","f","p","phi_f" '
      write(102,*) 'ZONE I=', IMAX*JMAX, ' F=BLOCK'
      write(102,*) ((xi(i),j=1,jmax),i=1,imax)
      write(102,*) ((yj(j),j=1,jmax),i=1,imax)
      write(102,*) ((u(i,j),j=1,jmax),i=1,imax)
      write(102,*) ((v(i,j),j=1,jmax),i=1,imax)
      write(102,*) ((f(i,j),j=1,jmax),i=1,imax)
      write(102,*) ((p(i,j),j=1,jmax),i=1,imax)
      write(102,*) ((phif(i,j),j=1,jmax),i=1,imax)
      
   70 CONTINUE

c out results (c.wang)
      write(*,*) imax,jmax,ibar,jbar
C      write(102,*) 'VARIABLES="x","y","u","v","f","p" '
C      write(102,*) 'ZONE I=', IMAX*JMAX, ' F=BLOCK'
C      write(102,*) ((xi(i),j=1,jmax),i=1,imax)
C      write(102,*) ((yj(j),j=1,jmax),i=1,imax)
C      write(102,*) ((u(i,j),j=1,jmax),i=1,imax)
C      write(102,*) ((v(i,j),j=1,jmax),i=1,imax)
C      write(102,*) ((f(i,j),j=1,jmax),i=1,imax)
C      write(102,*) ((p(i,j),j=1,jmax),i=1,imax)
      

C
C     *** SET THE ADVANCE TIME ARRAYS INTO THE TIME-N ARRAYS
C
      DO I=1,IMAX
      DO J=1,JMAX
      UN(I,J)=U(I,J)
      VN(I,J)=V(I,J)
      U(I,J)=0.0
      V(I,J)=0.0
      PN(I,J)=P(I,J)
      FN(I,J)=F(I,J)
      phifn(i,j)=phif(i,j)
      ENDDO
      ENDDO
      NREGN=NREG


C
C     *** ADJUST DELT
C
   90 CALL DELTADJ
C
C     *** ADVANCE TIME
C
      T=T+DELT
      IF (T.GT.TWFIN) GO TO 100  ! TWFIN: time to end calculation
C
C     *** ADVANCE CYCLE
C
      CYCLE=CYCLE+1
      IF (NFLGC.GE.25.OR.NOCON.GE.25) T=EP10  ! NOCON: f,p convergence failed
      GO TO 10
C
  100 CALL EXIT
C
  110 FORMAT (10A8)
      END

      subroutine vof_sph    !!! by CWANG
      
      include "include.f"
      save my_count
C
C      call sola
CC      GO TO 20
C
C     *** START TIME CYCLE
C

      
   10 CONTINUE

      if(my_count==0)then
         my_count = 1
         GOTO 20
      endif
      ITER=0             ! pressure iteration counter
      FLG=1.0            ! pressure convergence indicator
      FNOC=0.0           ! pressure convergence failure indicator
C
C     *** EXPLICITLY APPROXIMATE NEW TIME-LEVEL VELOCITIES
C
      CALL TILDE
C
      IF (NMAT.EQ.2.AND.ISURF10.EQ.1) CALL TMS10    ! surface tension
C
C     *** SET BOUNDARY CONDITIONS
C
      CALL BC
C
C     *** ITERATIVELY ADJUST CELL PRESSURE AND VELOCITY
C
      CALL PRESSIT

C
      IF (T.GT.EP10) GO TO 30  ! EP10= 1.0e+10 
C
C     *** UPDATE FLUID CONFIGURATION
C
   20 CALL VFCONV
C
      IF (FLGC.GT.0.5) GO TO 90    
C
C     *** SET BOUNDARY CONDITIONS
C
      CALL BC
C
C     *** MOVE MARKER PARTICLES
C
      CALL PARMOV
C
C     *** DETERMINE PRESSURE INTERPOLATION FACTOR AND NEIGHBOR
C     *** ALSO DETERMINE SURFACE TENSION PRESSURES AND
C     *** WALL ADHESION EFFECTS IN SURFACE CELLS
C
      CALL PETACAL
C
C     *** PRINT TIME AND CYCLE DATA ON PAPER AND/OR FILM
C
   30 CALL PRTPLT (2)
C
      IF (CYCLE.LE.0) GO TO 40              ! if(cycle.le.0)then
      IF (T+EM6.LT.TWPLT) GO TO 50          !    call prtplot(3)
      TWPLT=TWPLT+PLTDT                         
   40 CONTINUE                              ! elseif(T+EM6.lt.twplt)then
C                                           ! else   
C     *** PRINT FIELD VARIABLE DATA ON FILM !    twplt=twplt+pltdt
C                                           !    call prtplt(3)   
CC      CALL PRTPLT (3)                       ! endif
C                                            
C     *** PLOT VELOCITY VECTOR, FREE SURFACE, MESH,
C     *** AND MARKER PARTICLE ON FILM
C
cwang      CALL DRAW
C
   50 CONTINUE
      IF (CYCLE.LE.0) GO TO 60             ! if(cycle.le.0)then
      IF (T+EM6.LT.TWPRT) GO TO 70         !    call prtlt(4)
      TWPRT=TWPRT+PRTDT                    ! elseif(T+EM6.lt.twprt)then
   60 CONTINUE                             ! else
C                                          !    twprt = twprt+prtdt
C     *** PRINT FIELD VARIABLE DATA ON PAPER!   call prtplt(4)
C                                          ! endif
CC      CALL PRTPLT (4)
C
      
      write(102,*) 'VARIABLES="x","y","u","v","f","p","phi_f" '
      write(102,*) 'ZONE I=', IMAX*JMAX, ' F=BLOCK'
      write(102,*) ((xi(i),j=1,jmax),i=1,imax)
      write(102,*) ((yj(j),j=1,jmax),i=1,imax)
      write(102,*) ((u(i,j),j=1,jmax),i=1,imax)
      write(102,*) ((v(i,j),j=1,jmax),i=1,imax)
      write(102,*) ((f(i,j),j=1,jmax),i=1,imax)
      write(102,*) ((p(i,j),j=1,jmax),i=1,imax)
      write(102,*) ((phif(i,j),j=1,jmax),i=1,imax)
      
   70 CONTINUE

c out results (c.wang)
      write(*,*) imax,jmax,ibar,jbar
C      write(102,*) 'VARIABLES="x","y","u","v","f","p" '
C      write(102,*) 'ZONE I=', IMAX*JMAX, ' F=BLOCK'
C      write(102,*) ((xi(i),j=1,jmax),i=1,imax)
C      write(102,*) ((yj(j),j=1,jmax),i=1,imax)
C      write(102,*) ((u(i,j),j=1,jmax),i=1,imax)
C      write(102,*) ((v(i,j),j=1,jmax),i=1,imax)
C      write(102,*) ((f(i,j),j=1,jmax),i=1,imax)
C      write(102,*) ((p(i,j),j=1,jmax),i=1,imax)
      

C
C     *** SET THE ADVANCE TIME ARRAYS INTO THE TIME-N ARRAYS
C
      DO I=1,IMAX
      DO J=1,JMAX
      UN(I,J)=U(I,J)
      VN(I,J)=V(I,J)
      U(I,J)=0.0
      V(I,J)=0.0
      PN(I,J)=P(I,J)
      FN(I,J)=F(I,J)
      phifn(i,j)=phif(i,j)
      ENDDO
      ENDDO
      NREGN=NREG

      
C
C     *** ADJUST DELT
C
   90 CALL DELTADJ
C
C     *** ADVANCE TIME
C
      T=T+DELT
      IF (T.GT.TWFIN) GO TO 100  ! TWFIN: time to end calculation
C
C     *** ADVANCE CYCLE
C
      CYCLE=CYCLE+1
      IF (NFLGC.GE.25.OR.NOCON.GE.25) T=EP10  ! NOCON: f,p convergence failed
C      GO TO 10
C
      return
  100 CALL EXIT
C
  110 FORMAT (10A8)
      END
      
      
      
      SUBROUTINE TILDE
      include "include.f"
C
C     *** COMPUTE TEMPORARY U AND V EXPLICITLY
C
      DO 20 J=2,JM1
      DO 20 I=2,IM1
      U(I,J)=0.0
      RDELX=1.0/(DELX(I)+DELX(I+1))
      RDELY=1.0/(DELY(J)+DELY(J+1))
      IF (F(I,J)+F(I+1,J).LT.EMF.AND.NMAT.EQ.1) GO TO 10
      IF (BETA(I,J).LT.0.0.OR.BETA(I+1,J).LT.0.0) GO TO 10
      SGU=SIGN(1.0,UN(I,J))
      DUDR=(UN(I+1,J)*phifn(i+1,j)-UN(I,J)*phifn(i,j))*RDX(I+1)  ! by CWANG
      DUDL=(UN(I,J)*phifn(i,j)-UN(I-1,J)*phifn(i-1,j))*RDX(I)
      !DUDR=(UN(I+1,J)-UN(I,J))*RDX(I+1) 
      !DUDL=(UN(I,J)-UN(I-1,J))*RDX(I)      
      RDXA=DELX(I)+DELX(I+1)+ALPHA*SGU*(DELX(I+1)-DELX(I))
      RDXA=1.0/RDXA
      FUX=RDXA*UN(I,J)*(DELX(I)*DUDR+DELX(I+1)*DUDL+ALPHA*SGU*(DELX(I+1)
     1 *DUDL-DELX(I)*DUDR))
      VBT=(DELX(I)*VN(I+1,J)+DELX(I+1)*VN(I,J))*RDELX
      VBB=(DELX(I)*VN(I+1,J-1)+DELX(I+1)*VN(I,J-1))*RDELX
      VAV=0.5*(VBT+VBB)
      DYT=0.5*(DELY(J)+DELY(J+1))
      DYB=0.5*(DELY(J-1)+DELY(J))
      DUDT=(UN(I,J+1)*phifn(i,j+1)-UN(I,J)*phifn(i,j))/DYT   !!! by CWANG
      DUDB=(UN(I,J)*phifn(i,j)-UN(I,J-1)*phifn(i,j-1))/DYB
      !DUDT=(UN(I,J+1)-UN(I,J))/DYT   
      !DUDB=(UN(I,J)-UN(I,J-1))/DYB      
      SGV=SIGN(1.0,VAV)
      DYA=DYT+DYB+ALPHA*SGV*(DYT-DYB)
      FUY=(VAV/DYA)*(DYB*DUDT+DYT*DUDB+ALPHA*SGV*(DYT*DUDB-DYB*DUDT))
      UBDYT=(DELY(J)*UN(I,J+1)+DELY(J+1)*UN(I,J))/(DELY(J)+DELY(J+1))
      UBDYB=(DELY(J-1)*UN(I,J)+DELY(J)*UN(I,J-1))/(DELY(J)+DELY(J-1))
      DUDXSQ=2.0*(UN(I-1,J)*RDX(I)/(DELX(I)+DELX(I+1))+UN(I+1,J)*RDX(I+1
     1 )/(DELX(I)+DELX(I+1))-UN(I,J)*RDX(I)*RDX(I+1))
      DUDYT=(UN(I,J+1)*DELY(J)*RDY(J+1)-UN(I,J)*DELY(J+1)*RDY(J)-UBDYT*
     1 (DELY(J)*RDY(J+1)-DELY(J+1)*RDY(J)))/(0.5*(DELY(J)+DELY(J+1)))
      DUDYB=(UN(I,J)*DELY(J-1)*RDY(J)-UN(I,J-1)*DELY(J)*RDY(J-1)-UBDYB*
     1 (DELY(J-1)*RDY(J)-DELY(J)*RDY(J-1)))/(0.5*(DELY(J-1)+DELY(J)))
      DUDYSQ=(DUDYT-DUDYB)*RDY(J)
      DUDXL=(UN(I,J)-UN(I-1,J))*RDX(I)
      DUDXR=(UN(I+1,J)-UN(I,J))*RDX(I+1)
      RXDUDX=RX(I)*(DELX(I+1)*DUDXL+DELX(I)*DUDXR)/(DELX(I)+DELX(I+1))
      RXSQU=UN(I,J)*RX(I)**2
      VISX=NU*(DUDXSQ+DUDYSQ+CYL*RXDUDX-CYL*RXSQU)
      RHOX=(RHOFC+RHOD*F(I,J))*DELX(I+1)+(RHOFC+RHOD*F(I+1,J))*DELX(I)
      U(I,J)=UN(I,J)*phifn(i,j)+DELT*((P(I,J)-P(I+1,J))*2.0/RHOX       ! by CWANG
     1 *phifn(i,j)+GX*phifn(i,j)-FUX-FUY+VISX*phifn(i,j))
      u(i,j)=u(i,j)/phifn(i,j)      
      !U(I,J)=UN(I,J)+DELT*((P(I,J)-P(I+1,J))*2.0/RHOX+GX-FUX-FUY+VISX)      

   10 CONTINUE
      V(I,J)=0.0
      IF (F(I,J)+F(I,J+1).LT.EMF.AND.NMAT.EQ.1) GO TO 20
      IF (BETA(I,J).LT.0.0.OR.BETA(I,J+1).LT.0.0) GO TO 20
      UBR=(DELY(J+1)*UN(I,J)+DELY(J)*UN(I,J+1))*RDELY
      UBL=(DELY(J+1)*UN(I-1,J)+DELY(J)*UN(I-1,J+1))*RDELY
      UAV=0.5*(UBR+UBL)
      DXR=0.5*(DELX(I)+DELX(I+1))
      DXL=0.5*(DELX(I)+DELX(I-1))
      SGU=SIGN(1.0,UAV)
      DXA=DXR+DXL+ALPHA*SGU*(DXR-DXL)
      DVDR=(VN(I+1,J)*phifn(i+1,j)-VN(I,J)*phifn(i,j))/DXR    ! by CWang
      DVDL=(VN(I,J)*phifn(i,j)-VN(I-1,J)*phifn(i-1,j))/DXL    ! by CWang
!      DVDR=(VN(I+1,J)-VN(I,J))/DXR    
!      DVDL=(VN(I,J)-VN(I-1,J))/DXL    
      FVX=(UAV/DXA)*(DXL*DVDR+DXR*DVDL+ALPHA*SGU*(DXR*DVDL-DXL*DVDR))
      SGV=SIGN(1.0,VN(I,J))
      DYA=DELY(J+1)+DELY(J)+ALPHA*SGV*(DELY(J+1)-DELY(J))
      DVDT=(VN(I,J+1)*phifn(i,j+1)-VN(I,J)*phifn(i,j))*RDY(J+1) ! by CWang
      DVDB=(VN(I,J)*phifn(i,j)-VN(I,J-1)*phifn(i,j-1))*RDY(J)    ! by CWang
!      DVDT=(VN(I,J+1)-VN(I,J))*RDY(J+1) 
!      DVDB=(VN(I,J)-VN(I,J-1))*RDY(J)         
      FVY=(VN(I,J)/DYA)*(DELY(J)*DVDT+DELY(J+1)*DVDB+ALPHA*SGV*(DELY(J+1
     1 )*DVDB-DELY(J)*DVDT))
      VBDYR=(DELX(I+1)*VN(I,J)+DELX(I)*VN(I+1,J))/(DELX(I)+DELX(I+1))
      VBDYL=(DELX(I)*VN(I-1,J)+DELX(I-1)*VN(I,J))/(DELX(I)+DELX(I-1))
      DVDXR=(VN(I+1,J)*DELX(I)*RDX(I+1)-VN(I,J)*DELX(I+1)*RDX(I)-VBDYR*
     1 (DELX(I)*RDX(I+1)-DELX(I+1)*RDX(I)))/(0.5*(DELX(I+1)+DELX(I)))
      DVDXL=(VN(I,J)*DELX(I-1)*RDX(I)-VN(I-1,J)*DELX(I)*RDX(I-1)-VBDYL*
     1 (DELX(I-1)*RDX(I)-DELX(I)*RDX(I-1)))/(0.5*(DELX(I)+DELX(I-1)))
      DVDXSQ=(DVDXR-DVDXL)*RDX(I)
      DVDYSQ=2.0*(VN(I,J-1)*RDY(J)/(DELY(J+1)+DELY(J))-VN(I,J)*RDY(J+1)
     1 *RDY(J)+VN(I,J+1)*RDY(J+1)/(DELY(J+1)+DELY(J)))
      DVDXRX=(VBDYR-VBDYL)*RDX(I)*RXI(I)
      VISY=NU*(DVDXSQ+DVDYSQ+CYL*DVDXRX)
      RHOY=(RHOFC+RHOD*F(I,J))*DELY(J+1)+(RHOFC+RHOD*F(I,J+1))*DELY(J)
!      V(I,J)=VN(I,J)+DELT*((P(I,J)-P(I,J+1))*2.0/RHOY+GY-FVX-FVY+VISY)

      V(I,J)=VN(I,J)*phifn(i,j)+DELT*((P(I,J)-P(I,J+1))*2.0/RHOY      ! by CWANG
     1 *phifn(i,j)
     1 +phifn(i,j)*GY-FVX-FVY+VISY*phifn(i,j))
      v(i,j)=v(i,j)/phifn(i,j)      
   20 CONTINUE
      RETURN
      END
      SUBROUTINE TMS10
      include "include.f"
C
C     *** TWO MATERIAL SURFACE TENSION
C
C     *** NOTE: THIS ROUTINE INTRODUCES SOME NUMERICAL NOISE
C     ***       AND MAY BE REPLACED IN THE FUTURE
C
      DO 30 I=2,IM1
      DO 30 J=2,JM1
      IF (NF(I,J).EQ.0.OR.NF(I,J).GE.5.OR.BETA(I,J).LT.0.0) GO TO 30
      WHTL=0.5
      WHTR=0.5
      WHTT=0.5
      WHTB=0.5
      IF (NF(I,J).GT.2) GO TO 10
      WHTL=1.0-F(I,J)
      IF (NF(I,J).EQ.2) WHTL=1.0-WHTL
      WHTR=1.0-WHTL
      STFX=PS(I,J)*DELY(J)
      IF (NF(I,J).EQ.1) STFX=-STFX
      STFY=STFX*TANTH(I,J)
      GO TO 20
   10 WHTB=1.0-F(I,J)
      IF (NF(I,J).EQ.4) WHTB=1.0-WHTB
      WHTT=1.0-WHTB
      STFY=PS(I,J)*DELX(I)
      IF (NF(I,J).EQ.3) STFY=-STFY
      STFX=-STFY*TANTH(I,J)
   20 CONTINUE
      RHOXR=(RHOFC+RHOD*F(I,J))*DELX(I+1)+(RHOFC+RHOD*F(I+1,J))*DELX(I)
      U(I,J)=U(I,J)+2.0*DELT*WHTR*STFX/(RHOXR*DELY(J))
      RHOXL=(RHOFC+RHOD*F(I-1,J))*DELX(I)+(RHOFC+RHOD*F(I,J))*DELX(I-1)
      U(I-1,J)=U(I-1,J)+2.0*DELT*WHTL*STFX/(RHOXL*DELY(J))
      RHOYT=(RHOFC+RHOD*F(I,J))*DELY(J+1)+(RHOFC+RHOD*F(I,J+1))*DELY(J)
      V(I,J)=V(I,J)+2.0*DELT*WHTT*STFY/(RHOYT*DELX(I))
      RHOYB=(RHOFC+RHOD*F(I,J-1))*DELY(J)+(RHOFC+RHOD*F(I,J))*DELY(J-1)
      V(I,J-1)=V(I,J-1)+2.0*DELT*WHTB*STFY/(RHOYB*DELX(I))
   30 CONTINUE
      RETURN
      END
      SUBROUTINE VFCONV
      include "include.f"
C
C     *** CONVECT THE VOLUME OF FLUID FUNCTION F
C
      IF (CYCLE.LT.1) GO TO 40
      FLGC=0.0
      DO 30 J=1,JM1
      DO 30 I=1,IM1
      VX=U(I,J)*DELT
      VY=V(I,J)*DELT
      ABVX=ABS(VX)
      ABVY=ABS(VY)
      IF (ABVX.GT.0.5*DELX(I).OR.ABVY.GT.0.5*DELY(J)) FLGC=1.0
      IA=I+1
      ID=I
      IDM=MAX0(I-1,1)
      RB=X(I)
      RA=XI(I+1)
      RD=XI(I)
      IF (VX.GE.0.0) GO TO 10
      IA=I
      ID=I+1
      IDM=MIN0(I+2,IMAX)
      RA=XI(I)
      RD=XI(I+1)
   10 CONTINUE
      IAD=IA
      IF (NF(ID,J).EQ.3.OR.NF(ID,J).EQ.4) IAD=ID
      IF (FN(IA,J).LT.EMF.OR.FN(IDM,J).LT.EMF) IAD=IA
      FDM=AMAX1(FN(IDM,J),FN(ID,J))
      FX1=FN(IAD,J)*ABS(VX)+AMAX1((FDM-FN(IAD,J))*ABS(VX)-(FDM-FN(ID,J))
     1 *DELX(ID),0.0)
      FX=AMIN1(FX1,FN(ID,J)*DELX(ID))
      F(ID,J)=F(ID,J)-FX*RDX(ID)*((ABS(RB/RD))*CYL+(1.0-CYL))
      F(IA,J)=F(IA,J)+FX*RDX(IA)*((ABS(RB/RA))*CYL+(1.0-CYL))
      JA=J+1
      JD=J
      JDM=MAX0(J-1,1)
      IF (VY.GE.0.0) GO TO 20
      JA=J
      JD=J+1
      JDM=MIN0(J+2,JMAX)
   20 CONTINUE
      JAD=JA
      IF (NF(I,JD).EQ.1.OR.NF(I,JD).EQ.2) JAD=JD
      IF (FN(I,JA).LT.EMF.OR.FN(I,JDM).LT.EMF) JAD=JA
      FDM=AMAX1(FN(I,JDM),FN(I,JD))
      FY1=FN(I,JAD)*ABS(VY)+AMAX1((FDM-FN(I,JAD))*ABS(VY)-(FDM-FN(I,JD))
     1 *DELY(JD),0.0)
      FY=AMIN1(FY1,FN(I,JD)*DELY(JD))
      F(I,JD)=F(I,JD)-FY*RDY(JD)
      F(I,JA)=F(I,JA)+FY*RDY(JA)
   30 CONTINUE
   40 CONTINUE
      DO 80 J=2,JM1
      DO 80 I=2,IM1
      IF (BETA(I,J).LT.0.0) GO TO 80
      VCHG=0.0
      IF (F(I,J).GT.EMF.AND.F(I,J).LT.EMF1) GO TO 60
      IF (F(I,J).GE.EMF1) GO TO 50
      VCHG=F(I,J)
      F(I,J)=0.0
      GO TO 60
   50 CONTINUE
      VCHG=-(1.0-F(I,J))
      F(I,J)=1.0
   60 CONTINUE
      VCHGT=VCHGT+VCHG*DELX(I)*DELY(J)*(XI(I)*2.0*PI*CYL+(1.0-CYL))
      IF (F(I,J).LT.EMF1) GO TO 80
      IF (F(I+1,J).LT.EMF) GO TO 70
      IF (F(I-1,J).LT.EMF) GO TO 70
      IF (F(I,J+1).LT.EMF) GO TO 70
      IF (F(I,J-1).LT.EMF) GO TO 70
      GO TO 80
   70 F(I,J)=F(I,J)-1.1*EMF
      VCHG=1.1*EMF
      VCHGT=VCHGT+VCHG*DELX(I)*DELY(J)*(XI(I)*2.0*PI*CYL+(1.0-CYL))
   80 CONTINUE
C
C     *** SPECIAL BOUNDARY CONDITIONS FOR F
C
      RETURN
      END


      SUBROUTINE BC
      include "include.f"
C
C     *** SET BOUNDARY CONDITIONS
C
      DO 100 J=1,JMAX
      F(1,J)=F(2,J)
      F(IMAX,J)=F(IM1,J)
      P(1,J)=P(2,J)
      P(IMAX,J)=P(IM1,J)
      GO TO (10,20,30,40,30), WL
   10 U(1,J)=0.0
      V(1,J)=V(2,J)
      GO TO 50
   20 U(1,J)=0.0
      V(1,J)=-V(2,J)*DELX(1)/DELX(2)
      GO TO 50
   30 IF (ITER.GT.0) GO TO 50
      U(1,J)=U(2,J)*(X(2)*RX(1)*CYL+1.0-CYL)
      V(1,J)=V(2,J)
      GO TO 50
   40 U(1,J)=U(IM2,J)
      V(1,J)=V(IM2,J)
      F(1,J)=F(IM2,J)
   50 GO TO (60,70,80,90,80), WR
   60 U(IM1,J)=0.0
      V(IMAX,J)=V(IM1,J)
      GO TO 100
   70 U(IM1,J)=0.0
      V(IMAX,J)=-V(IM1,J)*DELX(IMAX)/DELX(IM1)
      GO TO 100
   80 IF (ITER.GT.0) GO TO 100
      U(IM1,J)=U(IM2,J)*(X(IM2)*RX(IM1)*CYL+1.0-CYL)
      V(IMAX,J)=V(IM1,J)
      GO TO 100
   90 U(IM1,J)=U(2,J)
      V(IM1,J)=V(2,J)
      P(IM1,J)=P(2,J)
      PS(IM1,J)=PS(2,J)
      F(IM1,J)=F(2,J)
      V(IMAX,J)=V(3,J)
      F(IMAX,J)=F(3,J)
  100 CONTINUE
      DO 200 I=1,IMAX
      F(I,1)=F(I,2)
      F(I,JMAX)=F(I,JM1)
      P(I,1)=P(I,2)
      P(I,JMAX)=P(I,JM1)
      GO TO (110,120,130,140,130), WT
  110 V(I,JM1)=0.0
      U(I,JMAX)=U(I,JM1)
      GO TO 150
  120 V(I,JM1)=0.0
      U(I,JMAX)=-U(I,JM1)*DELY(JMAX)/DELY(JM1)
      GO TO 150
  130 IF (ITER.GT.0) GO TO 150
      V(I,JM1)=V(I,JM2)
      U(I,JMAX)=U(I,JM1)
      GO TO 150
  140 V(I,JM1)=V(I,2)
      U(I,JM1)=U(I,2)
      P(I,JM1)=P(I,2)
      PS(I,JM1)=PS(I,2)
      F(I,JM1)=F(I,2)
      U(I,JMAX)=U(I,3)
      F(I,JMAX)=F(I,3)
  150 GO TO (160,170,180,190,180), WB
  160 V(I,1)=0.0
      U(I,1)=U(I,2)
      GO TO 200
  170 V(I,1)=0.0
      U(I,1)=-U(I,2)*DELY(1)/DELY(2)
      GO TO 200
  180 IF (ITER.GT.0) GO TO 200
      V(I,1)=V(I,2)
      U(I,1)=U(I,2)
      GO TO 200
  190 V(I,1)=V(I,JM2)
      U(I,1)=U(I,JM2)
      F(I,1)=F(I,JM2)
  200 CONTINUE
C
C     *** FREE SURFACE AND SLOPED BOUNDARY CONDITIONS
C
      DO 450 I=2,IM1
      XRP=RDX(I)+0.5*RXI(I)
      RXRP=1./XRP
      XRM=RDX(I)-0.5*RXI(I)
      IF (XRM.GT.0.0) GO TO 210
      RXRM=0.0
      GO TO 220
  210 CONTINUE
      RXRM=1./XRM
  220 CONTINUE
      DO 450 J=2,JM1
      IF (BETA(I,J).GT.0.0) GO TO 230
      BMR=0.0
      BMT=0.0
      BML=0.0
      BMB=0.0
      F(I,J)=0.0
      P(I,J)=0.0
      IF (BETA(I+1,J).GT.0.0) BMR=1.0
      IF (BETA(I,J+1).GT.0.0) BMT=1.0
      IF (BETA(I-1,J).GT.0.0) BML=1.0
      IF (BETA(I,J-1).GT.0.0) BMB=1.0
      BMTOT=BMR+BMT+BML+BMB
      IF (BMTOT.LE.0.0) GO TO 450
      F(I,J)=(BMR*F(I+1,J)+BMT*F(I,J+1)+BML*F(I-1,J)+BMB*F(I,J-1))/BMTOT
      P(I,J)=(BMR*P(I+1,J)+BMT*P(I,J+1)+BML*P(I-1,J)+BMB*P(I,J-1))/BMTOT
      GO TO 450
  230 CONTINUE
      IF (NMAT.EQ.2) GO TO 450
      IF (F(I,J).LT.EMF.OR.F(I,J).GT.EMF1) GO TO 450
      NFSB=0
      IF (F(I+1,J).LT.EMF) NFSB=NFSB+1
      IF (F(I,J+1).LT.EMF) NFSB=NFSB+2
      IF (F(I-1,J).LT.EMF) NFSB=NFSB+4
      IF (F(I,J-1).LT.EMF) NFSB=NFSB+8
      IF (NFSB.EQ.0) GO TO 450
      IF (NFSB.GT.8) GO TO 240
      GO TO (250,260,270,280,290,300,310,320), NFSB
  240 NFSB1=NFSB-8
      GO TO (330,340,350,360,370,380,390), NFSB1
  250 U(I,J)=(U(I-1,J)-DELX(I)*RDY(J)*(V(I,J)-V(I,J-1)))*(1.0-CYL)+CYL*
     1 (U(I-1,J)*XRM*RXRP-RDY(J)*RXRP*(V(I,J)-V(I,J-1)))
      GO TO 410
  260 V(I,J)=(V(I,J-1)-DELY(J)*RDX(I)*(U(I,J)-U(I-1,J)))*(1.0-CYL)+CYL*
     1 (V(I,J-1)-DELY(J)*(XRP*U(I,J)-XRM*U(I-1,J)))
      GO TO 410
  270 U(I,J)=U(I-1,J)*(1.0-CYL)+CYL*U(I-1,J)
      GO TO 260
  280 U(I-1,J)=(U(I,J)+DELX(I)*RDY(J)*(V(I,J)-V(I,J-1)))*(1.0-CYL)+CYL*
     1 (U(I,J)*XRP*RXRM+RDY(J)*RXRM*(V(I,J)-V(I,J-1)))
      GO TO 410
  290 U(I-1,J)=U(I-1,J-1)
      GO TO 250
  300 U(I-1,J)=U(I,J)*(1.0-CYL)+CYL*U(I,J)
      GO TO 260
  310 U(I-1,J)=U(I-1,J-1)
      U(I,J)=U(I,J-1)
      GO TO 260
  320 V(I,J-1)=(V(I,J)+DELY(J)*RDX(I)*(U(I,J)-U(I-1,J)))*(1.0-CYL)+CYL*
     1 (V(I,J)+DELY(J)*(XRP*U(I,J)-XRM*U(I-1,J)))
      GO TO 410
  330 U(I,J)=U(I-1,J)*(1.0-CYL)+CYL*U(I-1,J)
      GO TO 320
  340 V(I,J)=V(I-1,J)
      GO TO 320
  350 V(I,J)=V(I-1,J)
      V(I,J-1)=V(I-1,J-1)
      GO TO 250
  360 U(I-1,J)=U(I,J)*(1.0-CYL)+CYL*U(I,J)
      GO TO 320
  370 U(I,J)=U(I,J+1)
      U(I-1,J)=U(I-1,J+1)
      GO TO 320
  380 V(I,J)=V(I+1,J)
      V(I,J-1)=V(I+1,J-1)
      GO TO 280
  390 U(I,J)=U(I-1,J)*(1.0-CYL)+CYL*U(I-1,J)*XRM*RXRP
      V(I,J-1)=V(I,J)
      V(I,J+1)=V(I,J)
      GO TO 410
C
C     *** SET VELOCITIES IN EMPTY CELLS ADJACENT TO PARTIAL FLUID CELLS
C
  410 CONTINUE
      IF (FLG.GT.0.5.AND.ITER.GT.0) GO TO 450
      IF (F(I+1,J).GT.EMF) GO TO 420
      IF (F(I+1,J+1).LT.EMF) V(I+1,J)=V(I,J)
      IF (F(I+1,J-1).LT.EMF) V(I+1,J-1)=V(I,J-1)
  420 IF (F(I,J+1).GT.EMF) GO TO 430
      IF (F(I+1,J+1).LT.EMF) U(I,J+1)=U(I,J)
      IF (F(I-1,J+1).LT.EMF) U(I-1,J+1)=U(I-1,J)
  430 IF (F(I-1,J).GT.EMF) GO TO 440
      IF (F(I-1,J+1).LT.EMF) V(I-1,J)=V(I,J)
      IF (F(I-1,J-1).LT.EMF) V(I-1,J-1)=V(I,J-1)
  440 IF (F(I,J-1).GT.EMF) GO TO 450
      IF (F(I+1,J-1).LT.EMF) U(I,J-1)=U(I,J)
      IF (F(I-1,J-1).LT.EMF) U(I-1,J-1)=U(I-1,J)
  450 CONTINUE
C
C     *** SPECIAL VELOCITY BOUNDARY CONDITIONS
C
      RETURN
      END
      SUBROUTINE CAVOVO
      include "include.f"
C
C     *** CALCULATE VOID VOLUMES
C
C     *** INITIALIZE VOID VOLUMES
C
      DO 10 K=1,NVRM
   10 VOL(K)=0.0
C
C     *** COMPUTE VOID REGION VOLUMES
C
      DO 30 J=2,JM1
      DO 30 I=2,IM1
      INF=NF(I,J)
      IF (INF.EQ.0.OR.BETA(I,J).LT.0.0) GO TO 30
      VOLA=(1.0-F(I,J))*DELX(I)*DELY(J)*(1.0-CYL+CYL*2.0*3.14159*XI(I))
      IF (INF.GT.5) GO TO 20
      INFR=NF(I+1,J)
      INFT=NF(I,J+1)
      INFL=NF(I-1,J)
      INFB=NF(I,J-1)
      INF=MAX0(INFR,INFT,INFL,INFB)
   20 VOL(INF)=VOL(INF)+VOLA
   30 CONTINUE
      RETURN
      END
      SUBROUTINE DELTADJ
      include "include.f"
C
C     *** DELT (TIME STEP) ADJUSTMENT
C
      DELTN=DELT
      IF (FLGC.LT.0.5) GO TO 20
      T=T-DELT
      CYCLE=CYCLE-1
      DELT=0.5*DELT
      DO 10 I=1,IMAX
      DO 10 J=1,JMAX
      P(I,J)=PN(I,J)
      F(I,J)=FN(I,J)
      U(I,J)=0.0
      V(I,J)=0.0
   10 CONTINUE
      NFLGC=NFLGC+1
   20 CONTINUE
      IF (AUTOT.LT.0.5.AND.FNOC.LT.0.5) GO TO 35
      DUMX=EM10
      DVMX=EM10
      IF (FNOC.GT.0.5) DELT=0.5*DELT
      DO 30 I=2,IM1
      DO 30 J=2,JM1
      UDM=ABS(UN(I,J))/(XI(I+1)-XI(I))
      VDM=ABS(VN(I,J))/(YJ(J+1)-YJ(J))
      DUMX=AMAX1(DUMX,UDM)
      DVMX=AMAX1(DVMX,VDM)
   30 CONTINUE
      DTMP=1.01
      IF (ITER.GT.25) DTMP=0.99
      DELTO=DELT*DTMP
      CON=0.25
      DELT=AMIN1(DELTO,CON/DUMX,CON/DVMX,DTVIS,DTSFT)
      IF (IMOVY.GT.0) DELT=AMIN1(DELT,PLTDT)
   35 IF(DELT.EQ.DELTN .AND. NMAT.EQ.1) GO TO 50
      CTOS=DELT*RDTEXP
      COMG=AMIN1(CTOS**2,1.0)
      OMG1=(OMG-1.0)*COMG+1.0
      DO 40 I=1,IMAX
      DO 40 J=1,JMAX
      IF (BETA(I,J).LT.0.0) GO TO 40
      RHXR=(RHOFC+RHOD*F(I,J))*DELX(I+1)+(RHOFC+RHOD*F(I+1,J))*DELX(I)
      RHXL=(RHOFC+RHOD*F(I-1,J))*DELX(I)+(RHOFC+RHOD*F(I,J))*DELX(I-1)
      RHYT=(RHOFC+RHOD*F(I,J))*DELY(J+1)+(RHOFC+RHOD*F(I,J+1))*DELY(J)
      RHYB=(RHOFC+RHOD*F(I,J-1))*DELY(J)+(RHOFC+RHOD*F(I,J))*DELY(J-1)
      XX=DELT*RDX(I)*(2.0/RHXL+2.0/RHXR)+DELT*RDY(J)*(2.0/RHYT+2.0/RHYB)
      RHOR=RHOF/(RHOFC+RHOD*F(I,J))
      BETA(I,J)=OMG1/(XX*COMG+RCSQ*RHOR/DELT)
   40 CONTINUE
   50 CONTINUE
      RETURN
      END
      SUBROUTINE DRAW
      include "include.f"
      integer*2 CNN
C
C     *** PLOT VELOCITY VECTOR, FREE SURFACE, MESH,
C     *** AND MARKER PARTICLE ON FILM
C
      IF (IMOVY.EQ.1.AND.NP.GT.0) GO TO 80
C
C     *** VELOCITY VECTOR PLOT
C
cwang      CALL ADV (1)
      CALL FRAME (XMIN,XMAX,YMAX,YMIN)
      IF (IMOVY.EQ.1) GO TO 10
cwang      CALL LINCNT (62)
      WRITE (12,150) NAME
cwang      CALL LINCNT (64)
      WRITE (12,140) T,CYCLE
   10 CONTINUE
      NCOLR =10
      CALL DRWOBS
      AROFAC = 0.25*VELMX1
      DO 200 J=2,JM1
      DO 200 I=2,IM1
      IF (F(I,J).LT.0.5.AND.NMAT.EQ.1) GO TO 200
      IF (BETA(I,J).LT.0.0) GO TO 200
      UVEC=(U(I-1,J)+U(I,J))*0.5
      VVEC=(V(I,J-1)+V(I,J))*0.5
      VEC = SQRT(UVEC*UVEC+VVEC*VVEC)
      vmax   = MAX (VMAX,VEC)
      VMIN   = MIN (VMIN,VEC)
      RANGE  = VMAX - VMIN
      FCOLR  =  NCOLR
      DINT   = RANGE / FCOLR
200   CONTINUE
C
      DO 20 J=2,JM1
      DO 20 I=2,IM1
      IF (F(I,J).LT.0.5.AND.NMAT.EQ.1) GO TO 20
      IF (BETA(I,J).LT.0.0) GO TO 20
      XCC=XI(I)
      YCC=0.5*(Y(J)+Y(J-1))
      UVEC=(U(I-1,J)+U(I,J))*0.5
      VVEC=(V(I,J-1)+V(I,J))*0.5
      VEC = SQRT(UVEC*UVEC+VVEC*VVEC)
      IF (VEC .GT.VMAX ) THEN
           NCN = NCOLR
      ELSE IF (VEC .LT. VMIN) THEN
           NCN = 2
      ELSE
           CNN = (VEC -VMIN)/DINT
           CNN = MOD (CNN,50) + 1.5
           NCN = CNN
      ENDIF
cwang      CALL GRGRAL (0,0,NCN)
      ARROW = VEC*AROFAC
      COSAL = UVEC /VEC
      SINAL = VVEC / VEC
      EA = ARROW *COSAL
      BA = ARROW *SINAL
      CA = .17*EA
      DA = .17*BA
      VECU = UVEC*VELMX1 + XCC
      VECV = VVEC*VELMX1 + YCC
      XL3 = VECU -EA - DA
      YL3 = VECV -BA + CA
      XL4 = VECU -EA + DA
      YL4 = VECV -BA - CA
      CALL DRWVEC (XCC,YCC,VECU,VECV,1)
      CALL DRWVEC (VECU,VECV,XL3,YL3,1)
      CALL DRWVEC (XL4,YL4,VECU,VECV,1)
C      CALL PLTPT (XCC,YCC,53B,1)
C      CALL PLTPT (XCC,YCC,53 ,1)
   20 CONTINUE
C
C     *** DRAW FREE SURFACE
C
      FPL=0.5
      DO 40 I=2,IM1
      DO 40 J=2,JM1
      IF (BETA(I,J).LT.0.0) GO TO 40
      FATR=0.25*(F(I,J)+F(I+1,J)+F(I,J+1)+F(I+1,J+1))
      FXTR=0.5*(F(I+1,J+1)+F(I+1,J)-F(I,J+1)-F(I,J))/(XI(I+1)-XI(I))
      FYTR=0.5*(F(I,J+1)+F(I+1,J+1)-F(I,J)-F(I+1,J))/(YJ(J+1)-YJ(J))
      FTRS=FXTR**2+FYTR**2
      IF (FTRS.EQ.0.0) FTRS=EP10
      XTR=0.5*(XI(I+1)+XI(I))+(FPL-FATR)*FXTR/FTRS
      XTR=AMAX1(XTR,XI(I))
      XTR=AMIN1(XTR,XI(I+1))
      XTRM=-XTR
      YTR=0.5*(YJ(J)+YJ(J+1))+(FPL-FATR)*FYTR/FTRS
      YTR=AMAX1(YTR,YJ(J))
      YTR=AMIN1(YTR,YJ(J+1))
      IF (F(I,J).GT.0.5.AND.F(I+1,J).GT.0.5) GO TO 30
      IF (F(I,J).LT.0.5.AND.F(I+1,J).LT.0.5) GO TO 30
      FABR=0.25*(F(I,J)+F(I+1,J)+F(I,J-1)+F(I+1,J-1))
      FXBR=0.5*(F(I+1,J)+F(I+1,J-1)-F(I,J)-F(I,J-1))/(XI(I+1)-XI(I))
      FYBR=0.5*(F(I,J)+F(I+1,J)-F(I,J-1)-F(I+1,J-1))/(YJ(J)-YJ(J-1))
      FBRS=FXBR**2+FYBR**2
      IF (FBRS.EQ.0.0) FBRS=EP10
      XBR=0.5*(XI(I+1)+XI(I))+(FPL-FABR)*FXBR/FBRS
      XBR=AMAX1(XBR,XI(I))
      XBR=AMIN1(XBR,XI(I+1))
      YBR=0.5*(YJ(J)+YJ(J-1))+(FPL-FABR)*FYBR/FBRS
      YBR=AMAX1(YBR,YJ(J-1))
      YBR=AMIN1(YBR,YJ(J))
      CALL DRWVEC (XBR,YBR,XTR,YTR,1)
   30 CONTINUE
      IF (F(I,J).GT.0.5.AND.F(I,J+1).GT.0.5) GO TO 40
      IF (F(I,J).LT.0.5.AND.F(I,J+1).LT.0.5) GO TO 40
      FATL=0.25*(F(I,J)+F(I,J+1)+F(I-1,J)+F(I-1,J+1))
      FXTL=0.5*(F(I,J+1)+F(I,J)-F(I-1,J+1)-F(I-1,J))/(XI(I)-XI(I-1))
      FYTL=0.5*(F(I-1,J+1)+F(I,J+1)-F(I-1,J)-F(I,J))/(YJ(J+1)-YJ(J))
      FTLS=FXTL**2+FYTL**2
      IF (FTLS.EQ.0.0) FTLS=EP10
      XTL=0.5*(XI(I-1)+XI(I))+(FPL-FATL)*FXTL/FTLS
      XTL=AMAX1(XTL,XI(I-1))
      XTL=AMIN1(XTL,XI(I))
      YTL=0.5*(YJ(J)+YJ(J+1))+(FPL-FATL)*FYTL/FTLS
      YTL=AMAX1(YTL,YJ(J))
      YTL=AMIN1(YTL,YJ(J+1))
      CALL DRWVEC (XTL,YTL,XTR,YTR,1)
   40 CONTINUE
C
C     *** MESH PLOT
C
      IF (IMOVY.EQ.1) GO TO 130
      IF (T.GT.0.0) GO TO 70
cwang      CALL ADV (1)
      CALL DRWOBS
      DO 50 J=1,JM1
      YCC=Y(J)
      CALL DRWVEC (XMIN,YCC,XMAX,YCC,0)
   50 CONTINUE
      DO 60 I=1,IM1
      XCC=X(I)
      CALL DRWVEC (XCC,YMIN,XCC,YMAX,1)
   60 CONTINUE
   70 CONTINUE
C
C     *** PLOT PARTICLES
C
      IF (NP.EQ.0) GO TO 110
   80 CONTINUE
cwang      CALL ADV (1)
      CALL DRWOBS
      CALL FRAME (XMIN,XMAX,YMAX,YMIN)
      IF (IMOVY.EQ.1) GO TO 90
cwang      CALL LINCNT (62)
      WRITE (12,150) NAME
cwang      CALL LINCNT (64)
      WRITE (12,140) T,CYCLE
   90 CONTINUE
      DO 100 N=1,NP
C     CALL PLTPT (XP(N),YP(N),300B,1)
      CALL PLTPT (XP(N),YP(N),300 ,1)
  100 CONTINUE
  110 CONTINUE
      IF (IMOVY.EQ.1) GO TO 130
  130 CONTINUE
cwang      call grvide
      RETURN
C
  140 FORMAT (1H+,80X,2HT=,1PE10.3,4X,6HCYCLE=,I4)
  150 FORMAT (1H ,18X,10A8,1X,A10,2(1X,A8))
      END
       SUBROUTINE DRV(IX1,IY1,IX2,IY2)
C
       X1 = IX1
       X2 = IX2
       Y1 = IY1
       Y2 = IY2
cwang       CALL GRDEPA(X1,Y1)
cwang       CALL GRDESA(X2,Y2)
       RETURN
       END
      SUBROUTINE DRWOBS
      include "include.f"
C
C     *** DRAW AROUND ALL OBSTACLES
C
      DO 20 I=2,IM1
      DO 20 J=2,JM1
      IF ((BETA(I,J).LT.0.0.AND.BETA(I+1,J).LT.0.0).OR.(BETA(I,J).GT.0.0
     1 .AND.BETA(I+1,J).GT.0.0)) GO TO 10
      XONE=X(I)
      XTWO=XONE
      YONE=Y(J-1)
      YTWO=Y(J)
      CALL DRWVEC (XONE,YONE,XTWO,YTWO,1)
   10 IF ((BETA(I,J).LT.0.0.AND.BETA(I,J+1).LT.0.0).OR.(BETA(I,J).GT.0.0
     1 .AND.BETA(I,J+1).GT.0.0)) GO TO 20
      XONE=X(I-1)
      XTWO=X(I)
      YONE=Y(J)
      YTWO=YONE
      CALL DRWVEC (XONE,YONE,XTWO,YTWO,1)
   20 CONTINUE
      RETURN
      END
      SUBROUTINE DRWVEC (XONE,YONE,XTWO,YTWO,ISYM)
      include "include.f"
C
C     *** DRAW A VECTOR
C     *** PROVIDES A SYSTEM DEPENDANT CALL
C
      IC=0
      X1=XONE
      Y1=YONE
      X2=XTWO
      Y2=YTWO
   10 X01=(X1-XMIN)*SF+XSHFT
      Y01=(Y1-YMIN)*SF+YSHFT
      X02=(X2-XMIN)*SF+XSHFT
      Y02=(Y2-YMIN)*SF+YSHFT
      IX1=50.+920.0*X01
      IX2=50.+920.0*X02
      IY1=50.+920.0*Y01
      IY2=50.+920.0*Y02
      CALL DRV (IX1,IY1,IX2,IY2)
      IF (ISYMPLT.EQ.0.OR.ISYM.EQ.0) GO TO 20
      IC=IC+1
      IF (IC.GT.1) GO TO 20
      X1=-X1
      X2=-X2
      GO TO 10
   20 RETURN
      END

      SUBROUTINE FRAME (XXL,XXR,YYT,YYB)
C
C     *** DRAW A FRAME AROUND THE PLOT
C
      CALL DRWVEC (XXL,YYT,XXR,YYT,0)
      CALL DRWVEC (XXL,YYT,XXL,YYB,0)
      CALL DRWVEC (XXL,YYB,XXR,YYB,0)
      CALL DRWVEC (XXR,YYB,XXR,YYT,0)
      RETURN
      END
      SUBROUTINE LAVORE
      include "include.f"
C
C     *** LABEL VOID REGIONS - - VOID REGIONS ARE NF.EQ.6 AND ABOVE
C
      NNR=6
      NVR=6
      DO 30 J=2,JM1
      DO 30 I=2,IM1
      IF (NF(I,J).LT.6) GO TO 30
      INFB=NF(I,J-1)
      INFL=NF(I-1,J)
      IF (INFB.LT.6.AND.INFL.LT.6) GO TO 20
      IF (INFB.LT.6.OR.INFL.LT.6) GO TO 10
      NF(I,J)=MIN0(INFB,INFL)
      INRB=NR(INFB)
      INRL=NR(INFL)
      INRMN=MIN0(INRB,INRL)
      NR(INFB)=INRMN
      NR(INFL)=INRMN
      GO TO 30
   10 NF(I,J)=INFB
      IF (INFB.LT.6) NF(I,J)=INFL
      GO TO 30
   20 NF(I,J)=NVR
      NR(NVR)=NNR
      NVR=NVR+1
      NNR=NNR+1
   30 CONTINUE
C
C     *** REDEFINE REGION NUMBERS TO BE CONSECUTIVE
C
      NVR1=NVR-1
      NNR1=NNR-1
      KKN=7
      DO 50 KK=7,NNR1
      KFLG=0
      DO 40 K=7,NVR1
      IF (NR(K).NE.KK) GO TO 40
      NR(K)=KKN
      KFLG=1
   40 CONTINUE
      IF (KFLG.EQ.1) KKN=KKN+1
   50 CONTINUE
      NREG=KKN-6
C
C     *** REDEFINE VOID NUMBERS TO BE CONSECUTIVE IF NREG.GT.1
C
      DO 60 J=2,JM1
      DO 60 I=2,IM1
      INF=NF(I,J)
      IF (INF.LT.6) GO TO 60
      NF(I,J)=NR(INF)
   60 CONTINUE
      RETURN
      END
       subroutine lib4020
c
cwang       call grdebu(0)
c
       ig =1
cwang       call grstde(ig,4125,480,itest)
c
cwang       call grgrac(0,0,1,0,4,0.,0.)
c
cwang       call grcade(1,itest)
cwang       call grclot(-33.,100.,0.,100.)
cwang       call grfene(0.,1024.,0.,1024.)
       return
       end
cwang      SUBROUTINE LINCNT
c      ENTRY DRV
cwang      ENTRY GFR80
cwang      ENTRY GRPHFTN
cwang      ENTRY GRPHLUN
c      ENTRY LIB4020
c      ENTRY LINCNT
c      ENTRY PLT
cwang      ENTRY SETFLSH
cwang      RETURN
cwang      END
      SUBROUTINE MESHSET
      include "include.f"
C
C     *** MESH SETUP  (GENERATION)
C
      I=1
      J=1
      X(1)=XL(1)
      Y(1)=YL(1)
      DO 30 K=1,NKX
      DXML=(XC(K)-XL(K))/NXL(K)
      DXMR=(XL(K+1)-XC(K))/NXR(K)
      DXMN1=DXMN(K)
      NT=NXL(K)
      TN=NT
      TN=AMAX1(TN,1.0+EM6)
      DXMN(K)=AMIN1(DXMN1,DXML)
      CMC=(XC(K)-XL(K)-TN*DXMN(K))*TN/(TN-1.0)
      IF (NT.EQ.1) CMC=0.0
      BMC=XC(K)-XL(K)-CMC
      DO 10 L=1,NT
      I=I+1
      RLN=(FLOAT(L)-TN)/TN
   10 X(I)=XC(K)+BMC*RLN-CMC*RLN*RLN
      NT=NXR(K)
      TN=NT
      TN=AMAX1(TN,1.0+EM6)
      DXMN(K)=AMIN1(DXMN1,DXMR)
      CMC=(XL(K+1)-XC(K)-TN*DXMN(K))*TN/(TN-1.0)
      IF (NT.EQ.1) CMC=0.0
      BMC=XL(K+1)-XC(K)-CMC
      DO 20 L=1,NT
      I=I+1
      RLN=FLOAT(L)/TN
   20 X(I)=XC(K)+BMC*RLN+CMC*RLN*RLN
   30 CONTINUE
      IF (WR.NE.4) GO TO 40
      I=I+1
      X(I)=X(I-1)+X(2)-X(1)
   40 CONTINUE
      DO 70 K=1,NKY
      DYML=(YC(K)-YL(K))/NYL(K)
      DYMR=(YL(K+1)-YC(K))/NYR(K)
      DYMN1=DYMN(K)
      NT=NYL(K)
      TN=NT
      TN=AMAX1(TN,1.0+EM6)
      DYMN(K)=AMIN1(DYMN1,DYML)
      CMC=(YC(K)-YL(K)-TN*DYMN(K))*TN/(TN-1.0)
      IF (NT.EQ.1) CMC=0.0
      BMC=YC(K)-YL(K)-CMC
      DO 50 L=1,NT
      J=J+1
      RLN=(FLOAT(L)-TN)/TN
   50 Y(J)=YC(K)+BMC*RLN-CMC*RLN*RLN
      NT=NYR(K)
      TN=NT
      TN=AMAX1(TN,1.0+EM6)
      DYMN(K)=AMIN1(DYMN1,DYMR)
      CMC=(YL(K+1)-YC(K)-TN*DYMN(K))*TN/(TN-1.0)
      IF (NT.EQ.1) CMC=0.0
      BMC=YL(K+1)-YC(K)-CMC
      DO 60 L=1,NT
      J=J+1
      RLN=FLOAT(L)/TN
   60 Y(J)=YC(K)+BMC*RLN+CMC*RLN*RLN
   70 CONTINUE
      IF (WT.NE.4) GO TO 80
      J=J+1
      Y(J)=Y(J-1)+Y(2)-Y(1)
   80 CONTINUE
      NUMX=I
      NUMY=J
      NUMXM1=NUMX-1
      NUMYM1=NUMY-1
      NUMXP1=NUMX+1
      NUMYP1=NUMY+1
      IBAR=NUMX-1
      JBAR=NUMY-1
      IMAX=IBAR+2
      JMAX=JBAR+2
      IM1=IMAX-1
      JM1=JMAX-1
      IM2=IMAX-2
      JM2=JMAX-2
C
C     *** CALCULATE VALUES NEEDED FOR VARIABLE MESH
C
      DO 100 I=1,NUMX
      IF (X(I).EQ.0.0) GO TO 90
      RX(I)=1.0/X(I)
      GO TO 100
   90 RX(I)=0.0
  100 CONTINUE
      DO 110 I=2,NUMX
      XI(I)=0.5*(X(I-1)+X(I))
      DELX(I)=X(I)-X(I-1)
      RXI(I)=1.0/XI(I)
  110 RDX(I)=1.0/DELX(I)
      DELX(1)=DELX(2)
      XI(1)=XI(2)-DELX(2)
      RXI(1)=1.0/XI(1)
      RDX(1)=1.0/DELX(1)
      DELXA=DELX(NUMX)
      IF(WR.EQ.4) DELXA=DELX(3)
      DELX(NUMXP1)=DELXA
      XI(NUMXP1)=XI(NUMX)+DELXA
      X(NUMXP1)=XI(NUMXP1)+0.5*DELX(NUMXP1)
      RXI(NUMXP1)=1.0/XI(NUMXP1)
      RDX(NUMXP1)=1.0/DELX(NUMXP1)
      DO 120 I=2,NUMY
      YJ(I)=0.5*(Y(I-1)+Y(I))
      RYJ(I)=1.0/YJ(I)
      DELY(I)=Y(I)-Y(I-1)
      RDY(I)=1.0/DELY(I)
  120 CONTINUE
      DELY(1)=DELY(2)
      RDY(1)=1.0/DELY(1)
      YJ(1)=YJ(2)-DELY(2)
      RYJ(1)=1.0/YJ(1)
      DELYA=DELY(NUMY)
      IF(WT.EQ.4) DELYA=DELY(3)
      DELY(NUMYP1)=DELYA
      YJ(NUMYP1)=YJ(NUMY)+DELYA
      Y(NUMYP1)=YJ(NUMYP1)+0.5*DELY(NUMYP1)
      RYJ(NUMYP1)=1.0/YJ(NUMYP1)
      RDY(NUMYP1)=1.0/DELY(NUMYP1)
      WRITE (6,190)
      DO 130 I=1,NUMXP1
      WRITE (6,200) I,X(I),I,RX(I),I,DELX(I),I,RDX(I),I,XI(I),I,RXI(I)
  130 CONTINUE
      WRITE (6,190)
      DO 140 I=1,NUMYP1
      WRITE (6,210) I,Y(I),I,DELY(I),I,RDY(I),I,YJ(I),I,RYJ(I)
  140 CONTINUE
      IF (IMOVY.EQ.1) GO TO 170
      WRITE (12,190)
      DO 150 I=1,NUMXP1
      WRITE (12,200) I,X(I),I,RX(I),I,DELX(I),I,RDX(I),I,XI(I),I,RXI(I)
  150 CONTINUE
      WRITE (12,190)
      DO 160 I=1,NUMYP1
      WRITE (12,210) I,Y(I),I,DELY(I),I,RDY(I),I,YJ(I),I,RYJ(I)
  160 CONTINUE
  170 CONTINUE
C
C     *** TEST ARRAY SIZE
C
      IF (IMAX.LE.IBAR2.AND.JMAX.LE.JBAR2) GO TO 180
      WRITE (6,220)
C
      CALL EXIT
C
  180 CONTINUE
      RETURN
C
  190 FORMAT (1H1)
  200 FORMAT (1X,2HX(,I2,2H)=,1PE12.5,2X,3HRX(,I2,2H)=,1PE12.5,2X,5HDELX
     1(,I2,2H)=,1PE12.5,1X,4HRDX(,I2,2H)=,1PE12.5,2X,3HXI(,I2,2H)=,1PE12
     2 .5,2X,4HRXI(,I2,2H)=,1PE12.5)
  210 FORMAT (1X,2HY(,I2,2H)=,1PE12.5,3X,5HDELY(,I2,2H)=,1PE12.5,3X,4HRD
     1Y(,I2,2H)=,1PE12.5,3X,3HYJ(,I2,2H)=,1PE12.5,3X,4HRYJ(,I2,2H)=,1PE1
     2 2.5)
  220 FORMAT (41H  MESH SIZE GREATER THAN ARRAY DIMENSIONS)
      END
      SUBROUTINE PARMOV
      include "include.f"
C
C     *** MARKER PARTICLE MOVEMENT SECTION
C
      NPT=0
      NPN=0
      K=1
      KN=1
      PFLG=1.0
      IPER=IM1
      IF(WR.EQ.4) IPER=IM2
      JPER=JM1
      IF(WT.EQ.4) JPER=JM2
   10 IF (NP.EQ.NPT) GO TO 150
C
C     *** CLACULATE U WEIGHTED VELOCITY OF PARTICLE
C
      I=IP(K)
      J=JP(K)
      IF (YP(K).GT.YJ(J)) GO TO 20
      HPX=X(I)-XP(K)
      HMX=DELX(I)-HPX
      HPY=YJ(J)-YP(K)
      NORMY=(DELY(J)+DELY(J-1))*0.5
      HMY=NORMY-HPY
      UTOP=(U(I-1,J)*HPX+U(I,J)*HMX)*RDX(I)
      UBOT=(U(I-1,J-1)*HPX+U(I,J-1)*HMX)*RDX(I)
      UPART=(UTOP*HMY+UBOT*HPY)/NORMY
      GO TO 30
   20 HPX=X(I)-XP(K)
      HMX=DELX(I)-HPX
      HPY=YJ(J+1)-YP(K)
      NORMY=(DELY(J+1)+DELY(J))*0.5
      HMY=NORMY-HPY
      UTOP=(U(I-1,J+1)*HPX+U(I,J+1)*HMX)*RDX(I)
      UBOT=(U(I-1,J)*HPX+U(I,J)*HMX)*RDX(I)
      UPART=(UTOP*HMY+UBOT*HPY)/NORMY
C
C     *** CALCULATE V WEIGHTED VELOCITY OF PARTICLE
C
   30 IF (XP(K).GT.XI(I)) GO TO 40
      NORMX=(DELX(I)+DELX(I-1))*0.5
      RNORMX=1.0/NORMX
      HPX=XI(I)-XP(K)
      HMX=NORMX-HPX
      HPY=Y(J)-YP(K)
      HMY=DELY(J)-HPY
      VTOP=(V(I-1,J)*HPX+V(I,J)*HMX)*RNORMX
      VBOT=(V(I-1,J-1)*HPX+V(I,J-1)*HMX)*RNORMX
      VPART=(VTOP*HMY+VBOT*HPY)*RDY(J)
      GO TO 50
   40 NORMX=(DELX(I)+DELX(I+1))*0.5
      RNORMX=1.0/NORMX
      HPX=XI(I+1)-XP(K)
      HMX=NORMX-HPX
      HPY=Y(J)-YP(K)
      HMY=DELY(J)-HPY
      VTOP=(V(I,J)*HPX+V(I+1,J)*HMX)*RNORMX
      VBOT=(V(I,J-1)*HPX+V(I+1,J-1)*HMX)*RNORMX
      VPART=(VTOP*HMY+VBOT*HPY)*RDY(J)
   50 XPART=XP(K)+UPART*DELT
      YPART=YP(K)+VPART*DELT
      IF (XPART.GT.X(I)) IP(KN)=IP(K)+1
      IF (XPART.LT.X(I-1)) IP(KN)=IP(K)-1
      IF (YPART.GT.Y(J)) JP(KN)=JP(K)+1
      IF (YPART.LT.Y(J-1)) JP(KN)=JP(K)-1
      XP(KN)=XPART
      YP(KN)=YPART
      IF(XP(KN).LT.X(1)) GO TO 90
      IF(YP(KN).LT.Y(1)) GO TO 100
      IF(XP(KN).GT.X(IPER)) GO TO 110
      IF(YP(KN).GT.Y(JPER)) GO TO 120
      GO TO 130
   90 IF(WL.LE.2) GO TO 130
      IF(WL.NE.4) GO TO 140
      XP(KN)=XP(KN)+X(IM2)-X(1)
      IP(KN)=IP(KN)+IM2-1
      GO TO 130
  100 IF(WB.LE.2) GO TO 130
      IF(WB.NE.4) GO TO 140
      YP(KN)=YP(KN)+Y(JM2)-Y(1)
      JP(KN)=JP(KN)+JM2-1
      GO TO 130
  110 IF(WR.LE.2) GO TO 130
      IF(WR.NE.4) GO TO 140
      XP(KN)=XP(KN)-X(IM2)+X(1)
      IP(KN)=IP(KN)-IM2+1
      GO TO 130
  120 IF(WT.LE.2) GO TO 130
      IF(WT.NE.4) GO TO 140
      YP(KN)=YP(KN)-Y(JM2)+Y(1)
      JP(KN)=JP(KN)-JM2+1
  130 KN=KN+1
      NPN=NPN+1
  140 K=K+1
      NPT=NPT+1
      PFLG=1.0
      GO TO 10
  150 NP=NPN
      RETURN
      END
      SUBROUTINE PETACAL
      include "include.f"
C
C     *** DETERMINE THE PRESSURE INTERPOLATION FACTOR PETA
C     *** DETERMINE THE SURFACE TENSION PRESSURE AND
C     *** WALL ADHESION EFFECTS IN SURFACE CELLS
C
      DO 10 I=1,IMAX
      DO 10 J=1,JMAX
      NF(I,J)=0
      PS(I,J)=0.0
   10 PETA(I,J)=1.0
      IPASS=0
      DO 150 I=2,IM1
      DO 150 J=2,JM1
      TANTH(I,J)=EP10
      IF (BETA(I,J).LT.0.0) GO TO 150
      IF (F(I,J).LT.EMF) NF(I,J)=6
      IF (F(I,J).LT.EMF.OR.F(I,J).GT.EMF1) GO TO 150
      IF (F(I+1,J).LT.EMF) GO TO 20
      IF (F(I,J+1).LT.EMF) GO TO 20
      IF (F(I-1,J).LT.EMF) GO TO 20
      IF (F(I,J-1).LT.EMF) GO TO 20
      GO TO 150
   20 CONTINUE
C
C     *** CALCULATE THE PARTIAL DERIVATIVES OF F
C
      DXR=0.5*(DELX(I)+DELX(I+1))
      DXL=0.5*(DELX(I)+DELX(I-1))
      DYT=0.5*(DELY(J)+DELY(J+1))
      DYB=0.5*(DELY(J)+DELY(J-1))
      RXDEN=1.0/(DXR*DXL*(DXR+DXL))
      RYDEN=1.0/(DYT*DYB*(DYT+DYB))
      FL=F(I-1,J+1)
      IF (BETA(I-1,J+1).LT.0.0.OR.(I.EQ.2.AND.WL.LT.3)) FL=1.0
      FC=F(I,J+1)
      IF (BETA(I,J+1).LT.0.0) FC=1.0
      FR=F(I+1,J+1)
      IF (BETA(I+1,J+1).LT.0.0.OR.(I.EQ.IM1.AND.WR.LT.3)) FR=1.0
      AVFT=FL*DELX(I-1)+FC*DELX(I)+FR*DELX(I+1)
      FL=F(I-1,J-1)
      IF (BETA(I-1,J-1).LT.0.0.OR.(I.EQ.2.AND.WL.LT.3)) FL=1.0
      FC=F(I,J-1)
      IF (BETA(I,J-1).LT.0.0) FC=1.0
      FR=F(I+1,J-1)
      IF (BETA(I+1,J-1).LT.0.0.OR.(I.EQ.IM1.AND.WR.LT.3)) FR=1.0
      AVFB=FL*DELX(I-1)+FC*DELX(I)+FR*DELX(I+1)
      FL=F(I-1,J)
      IF (BETA(I-1,J).LT.0.0.OR.(I.EQ.2.AND.WL.LT.3)) FL=1.0
      FR=F(I+1,J)
      IF (BETA(I+1,J).LT.0.0.OR.(I.EQ.IM1.AND.WR.LT.3)) FR=1.0
      AVFCY=FL*DELX(I-1)+F(I,J)*DELX(I)+FR*DELX(I+1)
      FB=F(I,J-1)
      IF (BETA(I,J-1).LT.0.0.OR.(J.EQ.2.AND.WB.LT.3)) FB=1.0
      FT=F(I,J+1)
      IF (BETA(I,J+1).LT.0.0.OR.(J.EQ.JM1.AND.WT.LT.3)) FT=1.0
      AVFCX=FB*DELY(J-1)+F(I,J)*DELY(J)+FT*DELY(J+1)
      FB=F(I-1,J-1)
      IF (BETA(I-1,J-1).LT.0.0.OR.(J.EQ.2.AND.WB.LT.3)) FB=1.0
      FC=F(I-1,J)
      IF (BETA(I-1,J).LT.0.0) FC=1.0
      FT=F(I-1,J+1)
      IF (BETA(I-1,J+1).LT.0.0.OR.(J.EQ.JM1.AND.WT.LT.3)) FT=1.0
      AVFL=FB*DELY(J-1)+FC*DELY(J)+FT*DELY(J+1)
      FB=F(I+1,J-1)
      IF (BETA(I+1,J-1).LT.0.0.OR.(J.EQ.2.AND.WB.LT.3)) FB=1.0
      FC=F(I+1,J)
      IF (BETA(I+1,J).LT.0.0) FC=1.0
      FT=F(I+1,J+1)
      IF (BETA(I+1,J+1).LT.0.0.OR.(J.EQ.JM1.AND.WT.LT.3)) FT=1.0
      AVFR=FB*DELY(J-1)+FC*DELY(J)+FT*DELY(J+1)
C
C     *** BOUNDARY CONDITIONS FOR WALL ADHESION
C
      IF (ISURF10.EQ.0.OR.CANGLE.EQ.0.0) GO TO 60
      IF (BETA(I+1,J).GE.0.0.AND.I.NE.IM1) GO TO 30
      AVFR=AVFCX+0.5*(DELX(I)+DELX(I+1))/TANCA
      IF (F(I,J+1).LT.EMF.AND.F(I,J-1).GE.EMF) AVFT=AVFCY-0.5*(DELY(J)
     1 +DELY(J+1))*TANCA
      IF (F(I,J-1).LT.EMF.AND.F(I,J+1).GE.EMF) AVFB=AVFCY-0.5*(DELY(J)
     1 +DELY(J-1))*TANCA
   30 IF (BETA(I,J+1).GE.0.0.AND.J.NE.JM1) GO TO 40
      AVFT=AVFCY+0.5*(DELY(J)+DELY(J+1))/TANCA
      IF (F(I+1,J).LT.EMF.AND.F(I-1,J).GE.EMF) AVFR=AVFCX-0.5*(DELX(I)
     1 +DELX(I+1))*TANCA
      IF (F(I-1,J).LT.EMF.AND.F(I+1,J).GE.EMF) AVFL=AVFCX-0.5*(DELX(I)
     1 +DELX(I-1))*TANCA
   40 IF (BETA(I,J-1).GE.0.0.AND.J.NE.2) GO TO 50
      AVFB=AVFCY+0.5*(DELY(J)+DELY(J-1))/TANCA
      IF (F(I+1,J).LT.EMF.AND.F(I-1,J).GE.EMF) AVFR=AVFCX-0.5*(DELX(I)
     1 +DELX(I+1))*TANCA
      IF (F(I-1,J).LT.EMF.AND.F(I+1,J).GE.EMF) AVFL=AVFCX-0.5*(DELX(I)
     1 +DELX(I-1))*TANCA
   50 IF (BETA(I-1,J).GE.0.0.AND.I.NE.2) GO TO 60
      IF (CYL.GT.0.5.AND.X(1).EQ.0.0) GO TO 60
      AVFL=AVFCX+0.5*(DELX(I)+DELX(I-1))/TANCA
      IF (F(I,J+1).LT.EMF.AND.F(I,J-1).GE.EMF) AVFT=AVFCY-0.5*(DELY(J)
     1 +DELY(J+1))*TANCA
      IF (F(I,J-1).LT.EMF.AND.F(I,J+1).GE.EMF) AVFB=AVFCY-0.5*(DELY(J)
     1 +DELY(J-1))*TANCA
   60 CONTINUE
      XTHM=3.0*AMAX1(AVFT,AVFCY,AVFB)/(DELX(I-1)+DELX(I)+DELX(I+1))
      YTHM=3.0*AMAX1(AVFL,AVFCX,AVFR)/(DELY(J-1)+DELY(J)+DELY(J+1))
      PFX=RXDEN*((AVFR-AVFCX)*DXL**2+(AVFCX-AVFL)*DXR**2)
      PFY=RYDEN*((AVFT-AVFCY)*DYB**2+(AVFCY-AVFB)*DYT**2)
      PF=PFX**2+PFY**2
      IF (PF.GT.EM10) GO TO 70
      NF(I,J)=5
      P(I,J)=0.25*(P(I+1,J)+P(I,J+1)+P(I-1,J)+P(I,J-1))
      GO TO 150
   70 CONTINUE
C
C     *** DETERMINE THE PRESSURE INTERPOLATION CELL NF
C
      ABPFX=ABS(PFX)
      ABPFY=ABS(PFY)
      L=I
      M=J
      IF (ABPFY.GE.ABPFX) GO TO 80
      DXDYR=DELY(J)*RDX(I)
      PFMN=PFY
      NF(I,J)=2
      L=I+1
      DMX=DELX(I)
      DMIN=0.5*(DMX+DELX(I+1))
      IF (PFX.GT.0.0) GO TO 90
      NF(I,J)=1
      PFMN=-PFY
      L=I-1
      DMX=DELX(I)
      DMIN=0.5*(DMX+DELX(I-1))
      GO TO 90
   80 CONTINUE
      DXDYR=DELX(I)*RDY(J)
      PFMN=-PFX
      NF(I,J)=4
      M=J+1
      DMX=DELY(J)
      DMIN=0.5*(DMX+DELY(J+1))
      IF (PFY.GT.0.0) GO TO 90
      NF(I,J)=3
      PFMN=PFX
      M=J-1
      DMX=DELY(J)
      DMIN=0.5*(DMX+DELY(J-1))
   90 CONTINUE
      TANTH(I,J)=PFMN
      ABTAN=ABS(TANTH(I,J))
C
C     *** DETERMINE THE CURVATURE AND SURFACE PRESSURE
C
      DFS=(0.5-F(I,J))*DMX
      IF (F(I,J).LT.0.5*ABTAN*DXDYR) DFS=0.5*DMX*(1.0+DXDYR*ABTAN-SQRT(8
     1 .0*F(I,J)*DXDYR*ABTAN))
      IF (ISURF10.LT.1) GO TO 140
      NFC=NF(I,J)
      PXR=(AVFR-AVFCX)/DXR
      PXL=(AVFCX-AVFL)/DXL
      PYT=(AVFT-AVFCY)/DYT
      PYB=(AVFCY-AVFB)/DYB
      YDFS=-DFS
      IF(NFC.EQ.2 .OR. NFC.EQ.4) YDFS=DFS
      IF(NFC.GT.2) GO TO 100
      DXDN=DELY(J)
      XINB=YDFS+0.5*TANTH(I,J)*DXDN
      XINT=2.0*YDFS-XINB
      GP1=PYT
      PX1=PXL
      IF(XINT.GT.0.0) PX1=PXR
      IF(ABS(PX1).LT.ABS(GP1)) GP1=SIGN(1.0,GP1)/(ABS(PX1)+EM10)
      GP2=PYB
      PX2=PXR
      IF(XINB.LT.0.0) PX2=PXL
      IF(ABS(PX2).LT.ABS(GP2)) GP2=SIGN(1.0,GP2)/(ABS(PX2)+EM10)
      GO TO 110
  100 DXDN=DELX(I)
      YINR=YDFS+0.5*TANTH(I,J)*DXDN
      YINL=2.0*YDFS-YINR
      GP1=PXR
      PY1=PYT
      IF(YINR.LT.0.0) PY1=PYB
      IF(ABS(PY1).LT.ABS(GP1)) GP1=SIGN(1.0,GP1)/(ABS(PY1)+EM10)
      GP2=PXL
      PY2=PYB
      IF(YINL.GT.0.0) PY2=PYT
      IF(ABS(PY2).LT.ABS(GP2)) GP2=SIGN(1.0,GP2)/(ABS(PY2)+EM10)
  110 GP1D=1.0+GP1*GP1
      GP2D=1.0+GP2*GP2
      CURVXY=(GP2/SQRT(GP2D)-GP1/SQRT(GP1D))/DXDN
      CURVCYL=0.0
      IF (CYL.LT.1.0) GO TO 120
      XLITLR=XI(I)
      IF (NFC.EQ.1) XLITLR=X(I-1)+F(I,J)*DELX(I)
      IF (NFC.EQ.2) XLITLR=X(I)-F(I,J)*DELX(I)
      RLITLR=AMIN1(1.0/XLITLR,RXI(2))
      TRIG=ABS(SIN(ATAN(ABTAN)))
      IF (NFC.LE.2) TRIG=ABS(COS(ATAN(ABTAN)))
      CURVCYL=-CYL*TRIG*SIGN(1.0,PFX)*RLITLR
  120 CURV=CURVXY+CURVCYL
      PS(I,J)=SIGMA*CURV
      IF (XTHM.LT.1.0.OR.YTHM.LT.1.0) PS(I,J)=0.0
  140 CONTINUE
C
C     *** CALCULATE PETA
C
      NFSB=0
      IF (F(I+1,J).LT.EMF.OR.I.EQ.IM1.OR.BETA(I+1,J).LT.0.0) NFSB=NFSB+1
      IF (F(I,J+1).LT.EMF.OR.BETA(I,J+1).LT.0.0) NFSB=NFSB+2
      IF (F(I-1,J).LT.EMF.OR.BETA(I-1,J).LT.0.0) NFSB=NFSB+4
      IF (F(I,J-1).LT.EMF.OR.BETA(I,J-1).LT.0.0) NFSB=NFSB+8
      IF (NFSB.EQ.15) PS(I,J)=0.0
      IF (NMAT.EQ.2) GO TO 150
      PETA(I,J)=1.0/(1.0-DFS/DMIN)
      IF (L.EQ.1.OR.L.EQ.IMAX) PETA(I,J)=1.0
      IF (M.EQ.1.OR.M.EQ.JMAX) PETA(I,J)=1.0
      IF (BETA(L,M).LT.0.0) PETA(I,J)=1.0
  150 CONTINUE
C
      CALL LAVORE
C
      CALL CAVOVO
C
C     IF NECESSARY, DETERMINE PRESSURES PR FOR VOID REGIONS NF
C
      IF (NMAT.EQ.2) GO TO 300
C
C     *** SET PETA IN ADJACENT FULL CELL
C
      DO 290 J=1,JMAX
      DO 290 I=1,IMAX
      NFF=NF(I,J)
      IF (NFF.EQ.0.OR.BETA(I,J).LT.0.0) GO TO 290
      IF (NFF.GT.5) GO TO 280
      L=I
      M=J
      GO TO (230,240,250,260,290), NFF
  230 L=I-1
      DMX=DELX(L)
      DMIN=0.5*(DMX+DELX(I))
      GO TO 270
  240 L=I+1
      DMX=DELX(L)
      DMIN=0.5*(DMX+DELX(I))
      GO TO 270
  250 M=J-1
      DMX=DELY(M)
      DMIN=0.5*(DMX+DELY(J))
      GO TO 270
  260 M=J+1
      DMX=DELY(M)
      DMIN=0.5*(DMX+DELY(J))
  270 CONTINUE
      IF (NF(L,M).GT.0) GO TO 290
      CTOS=DELT*RDTEXP
      COMG=AMIN1(CTOS**2,1.0)
      BPD=1.0/PETA(L,M)-BETA(L,M)*(1.0-PETA(I,J))
     1*DELT/(DMIN*DMX)*(COMG/RHOF)
      PETA(L,M)= 1.0/BPD
      GO TO 290
  280 CONTINUE
      P(I,J)=PR(NFF)
  290 CONTINUE
  300 CONTINUE
      RETURN
      END
cwang      SUBROUTINE PLT (IX1,IY1,ICHAR)
cwang      INTEGER ICHAR(*)
cwang      x1 = ix1
cwang      y1 = iy1
cwang      call grsymb(3)
cwang      call grsyma(x1,y1)
cwang      RETURN
cwang      END
      SUBROUTINE PLTPT (XONE,YONE,ICHAR,ISYM)
      include "include.f"
C
C     *** PLOT (DRAW) A POINT
C     *** PROVIDES A SYSTEM DEPENDANT CALL
C
      IC=0
      X1=XONE
      Y1=YONE
   10 X01=(X1-XMIN)*SF+XSHFT
      Y01=(Y1-YMIN)*SF+YSHFT
      IX1=50.+920.0*X01
      IY1=50.+920.0*Y01
C      CALL PLT (IX1,IY1,42)
      IF (ABS(X1).LE.EM6) GO TO 20
      IF (ISYMPLT.EQ.0.OR.ISYM.EQ.0) GO TO 20
      IC=IC+1
      IF (IC.GT.1) GO TO 20
      X1=-X1
      GO TO 10
   20 RETURN
      END
      SUBROUTINE PRESSIT
      include "include.f"
C
C     *** PRESSURE ITERATION
C
C     *** TEST FOR CONVERGENCE
C
   10 IF (FLG.EQ.0.) GO TO 140
      ITER=ITER+1
      ITMAX=1000
      IF (ITER.LT.ITMAX) GO TO 20
      FNOC=1.0
      NOCON=NOCON+1
      GO TO 140
   20 FLG=0.0
C
C     *** COMPUTE UPDATED CELL PRESSURE AND VELOCITIES
C
      DO 130 J=JPB,JPT
      DO 130 I=IPL,IPR
      IF (BETA(I,J).LT.0.0) GO TO 130
      IF (NMAT.EQ.2) GO TO 80
      IF (F(I,J).LT.EMF) GO TO 130
      IF (NF(I,J).EQ.0) GO TO 80
C
C     *** CALCULATE PRESSURE FOR SURFACE CELLS
C
      NFF=NF(I,J)   ! pressure interploation cell
      L=I
      M=J
      GO TO (30,40,50,60,130), NFF
   30 L=I-1
      GO TO 70
   40 L=I+1
      GO TO 70
   50 M=J-1
      GO TO 70
   60 M=J+1
   70 CONTINUE
      NFEL=NF(I-1,J)
      NFER=NF(I+1,J)
      NFEB=NF(I,J-1)
      NFET=NF(I,J+1)
      NFE=MAX0(NFEL,NFER,NFEB,NFET)
      PSURF=PS(I,J)+PR(NFE)
      PLM=P(L,M)
      IF (NF(L,M).NE.0.AND.BETA(I,J).GT.0.0) PLM=PSURF
      DELP=(1.0-PETA(I,J))*PLM+PETA(I,J)*PSURF-P(I,J)
      GO TO 90
   80 CONTINUE
      DIJ=RDX(I)*(U(I,J)-U(I-1,J))+RDY(J)*(V(I,J)-V(I,J-1))+CYL*0.5*RXI
     1 (I)*(U(I,J)+U(I-1,J))
      RHOR=RHOF/(RHOFC+RHOD*F(I,J))
      DFUN=DIJ+RHOR*RCSQ*(P(I,J)-PN(I,J))/DELT
C
C     *** SET FLAG INDICATING CONVERGENCE
C
      IF (ABS(DFUN).GE.EPSI) FLG=1.0
      DELP=-BETA(I,J)*DFUN*PETA(I,J)
   90 CONTINUE
      P(I,J)=P(I,J)+DELP
      CTOS=DELT*RDTEXP
      COMG=AMIN1(CTOS**2,1.0)
      DPTC=2.0*DELT*DELP*COMG
      IF (BETA(I+1,J).LT.0.0) GO TO 100
      RHOXR=(RHOFC+RHOD*F(I,J))*DELX(I+1)+(RHOFC+RHOD*F(I+1,J))*DELX(I)
      U(I,J)=U(I,J)+DPTC/RHOXR
  100 IF (BETA(I-1,J).LT.0.0) GO TO 110
      RHOXL=(RHOFC+RHOD*F(I-1,J))*DELX(I)+(RHOFC+RHOD*F(I,J))*DELX(I-1)
      U(I-1,J)=U(I-1,J)-DPTC/RHOXL
  110 IF (BETA(I,J+1).LT.0.0) GO TO 120
      RHOYT=(RHOFC+RHOD*F(I,J))*DELY(J+1)+(RHOFC+RHOD*F(I,J+1))*DELY(J)
      V(I,J)=V(I,J)+DPTC/RHOYT
  120 IF (BETA(I,J-1).LT.0.0) GO TO 130
      RHOYB=(RHOFC+RHOD*F(I,J-1))*DELY(J)+(RHOFC+RHOD*F(I,J))*DELY(J-1)
      V(I,J-1)=V(I,J-1)-DPTC/RHOYB
  130 CONTINUE
      CALL BC
      GO TO 10
  140 CONTINUE
      RETURN
      END
      SUBROUTINE PRTPLT (N)
      include "include.f"
C
C     *** PRINT AND PLOT
C     *** PROVIDES FORMATTED WRITES TO PAPER AND FILM
C
      GO TO (10,70,90,130), N
C
C     *** PRTPLT (1) WRITE OUT INITIAL DATA AND MESH DATA
C
   10 WRITE (6,170)
      WRITE (6,180) NAME
      WRITE (6,220) IBAR,JBAR,DELT,NU,ICYL,EPSI,GX,GY,UI,VI,VELMX,TWFIN
     1 ,PRTDT,PLTDT,OMG,ALPHA,WL,WR,WT,WB,IMOVY,AUTOT,FLHT,ISYMPLT,SIGMA
     2 ,ISURF10,CANGLE,CSQ,NMAT,RHOF,RHOFC
      IF (IMOVY.GT.0) GO TO 40
      WRITE (12,170)
      WRITE (12,180) NAME
      WRITE (12,220) IBAR,JBAR,DELT,NU,ICYL,EPSI,GX,GY,UI,VI,VELMX,TWFIN
     1 ,PRTDT,PLTDT,OMG,ALPHA,WL,WR,WT,WB,IMOVY,AUTOT,FLHT,ISYMPLT,SIGMA
     2 ,ISURF10,CANGLE,CSQ,NMAT,RHOF,RHOFC
C
C     *** WRITE ON FILM VARIABLE MESH INPUT DATA
C
      WRITE (12,260) NKX
      DO 20 I=1,NKX
      WRITE (12,270) I,XL(I),XC(I),XL(I+1),NXL(I),NXR(I),DXMN(I)
   20 CONTINUE
      WRITE (12,280) NKY
      DO 30 I=1,NKY
      WRITE (12,275) I,YL(I),YC(I),YL(I+1),NYL(I),NYR(I),DYMN(I)
   30 CONTINUE
   40 CONTINUE
C
C     *** PRINT VARIABLE MESH INPUT DATA
C
      WRITE (6,260) NKX
      DO 50 I=1,NKX
      WRITE (6,270) I,XL(I),XC(I),XL(I+1),NXL(I),NXR(I),DXMN(I)
   50 CONTINUE
      WRITE (6,280) NKY
      DO 60 I=1,NKY
      WRITE (6,275) I,YL(I),YC(I),YL(I+1),NYL(I),NYR(I),DYMN(I)
   60 CONTINUE
      GO TO 160
C
C     *** PRTPLT (2)  WRITE TIME STEP, CYCLE INFORMATION
C
   70 CONTINUE
      WRITE (6,210) ITER,T,DELT,CYCLE,VCHGT
      IF (IMOVY.EQ.1) GO TO 80
      IF (T.GT.0.) GO TO 80
      WRITE (12,210) ITER,T,DELT,CYCLE,VCHGT
   80 CONTINUE
      GO TO 160
C
C     *** PRTPLT (3)  WRITE FIELD VARIABLES TO FILM
C
   90 IF (IMOVY.EQ.1) GO TO 120
cwang      CALL ADV (1)
      WRITE (12,250) NAME
      WRITE (12,210) ITER,T,DELT,CYCLE,VCHGT
      WRITE (12,240)
      WRITE (12,290) NREG
      WRITE (12,300)
      KNR=NREG+5
      DO 100 K=6,KNR
      WRITE (12,310) K,VOL(K),PR(K)
  100 CONTINUE
      WRITE (12,190)
      DO 110 I=1,IMAX
      DO 110 J=1,JMAX
      DIJ=RDX(I)*(U(I,J)-U(I-1,J))+RDY(J)*(V(I,J)-V(I,J-1))+CYL*0.5*RXI
     1 (I)*(U(I,J)+U(I-1,J))
      WRITE (12,200) I,J,U(I,J),V(I,J),P(I,J),DIJ,PS(I,J),F(I,J),NF(I,J)
     1 ,PETA(I,J)
  110 CONTINUE
  120 CONTINUE
      GO TO 160
C
C     *** PRTPLT (4)  WRITE FIELD VARIABLES TO PAPER
C
  130 WRITE (6,170)
      WRITE (6,250) NAME
      WRITE (6,210) ITER,T,DELT,CYCLE,VCHGT
      WRITE (6,240)
      WRITE (6,290) NREG
      WRITE (6,300)
      KNR=NREG+5
      DO 140 K=6,KNR
      WRITE (6,310) K,VOL(K),PR(K)
  140 CONTINUE
      WRITE (6,240)
      WRITE (6,190)
      DO 150 I=1,IMAX
      DO 150 J=1,JMAX
      DIJ=RDX(I)*(U(I,J)-U(I-1,J))+RDY(J)*(V(I,J)-V(I,J-1))+CYL*0.5*RXI
     1 (I)*(U(I,J)+U(I-1,J))
      WRITE (6,200) I,J,U(I,J),V(I,J),P(I,J),DIJ,PS(I,J),F(I,J),NF(I,J)
     1 ,PETA(I,J)
  150 CONTINUE
  160 RETURN
C
  170 FORMAT (1H1)
  180 FORMAT (10A8)
  190 FORMAT (4X,1HI,5X,1HJ,9X,1HU,14X,1HV,15X,1HP,15X,1HD,12X,2HPS,13X,
     1 1HF,11X,2HNF,9X,4HPETA)
  200 FORMAT (2X,I3,3X,I3,6(3X,1PE12.5),3X,I3,3X,E12.5)
  210 FORMAT (6X,6HITER= ,I5,5X,6HTIME= ,1PE12.5,5X,6HDELT= ,1PE12.5,5X,
     1 7HCYCLE= ,I4,5X,7HVCHGT= ,1PE12.5)
  220 FORMAT (1H ,5X,6HIBAR= ,I3/6X,6HJBAR= ,I3/6X,6HDELT= ,1PE12.5/8X,4
     1 HNU= ,E12.5/6X,6HICYL= ,I2/6X,6HEPSI= ,E12.5/8X,4HGX= ,E12.5/8X,4
     2 HGY= ,E12.5/8X,4HUI= ,E12.5/8X,4HVI= ,E12.5/5X,7HVELMX= ,E12.5/5X
     3 ,7HTWFIN= ,E12.5/5X,7HPRTDT= ,E12.5/5X,7HPLTDT= ,E12.5/7X,5HOMG=
     4 ,E12.5/5X,7HALPHA= ,E12.5/8X,4HWL= ,I2/8X,4HWR= ,I2/8X,4HWT= ,I2/
     5 8X,4HWB= ,I2/5X,7HIMOVY= ,E12.5/5X,7HAUTOT= ,E12.5/6X,6HFLHT=
     6 ,E12.5/3X,9HISYMPLT= ,I2/5X,7HSIGMA= ,E12.5/3X,9HISURF10= ,I2/4X,
     7 8HCANGLE= ,E12.5/7X,5HCSQ= ,E12.5/6X,6HNMAT= ,I2/6X,6HRHOF= ,E12.
     8 5/,5X,7HRHOFC= ,E12.5/)
  230 FORMAT (6X,6HCWTD= ,E12.5/6X,6HTRST= ,E12.5/)
  240 FORMAT (1H0)
  250 FORMAT (1H ,18X,10A8,1X,A10,2(1X,A8))
  260 FORMAT (2X,5HNKX= ,I4)
  270 FORMAT(2X,8HMESH-X= ,I4,3X,4HXL= ,1PE12.5,3X,4HXC= ,E12.5,3X,
     1 4HXR= ,E12.5,3X,5HNXL= ,I4,3X,5HNXR= ,I4,3X,6HDXMN= ,E12.5)
  275 FORMAT(2X,8HMESH-Y= ,I4,3X,4HYL= ,1PE12.5,3X,4HYC= ,E12.5,3X,
     1 4HYR= ,E12.5,3X,5HNYL= ,I4,3X,5HNYR= ,I4,3X,6HDYMN= ,E12.5)
  280 FORMAT (2X,5HNKY= ,I4)
  290 FORMAT (2X,6HNREG= ,I4)
  300 FORMAT (15X,1HK,6X,6HVOL(K),9X,5HPR(K))
  310 FORMAT (13X,I3,2X,1PE12.5,3X,E12.5)
      END
      SUBROUTINE SETUP
      include "include.f"
C
C     *** COMPUTE CONSTANT TERMS AND INITIALIZE NECESSARY VARIABLES
C
C     *** SET PARAMETER STATEMENT VALUE INTO CONSTANT
C
      NVRM=NVOR         ! max number of void regions
C
      CYL=FLOAT(ICYL)
      EMF1=1.0-EMF
      T=0.0
      ITER=0
      CYCLE=0
      TWPRT=0.0
      TWPLT=0.0
      VCHGT=0.0  !Accumulated fluid volume change
      NOCON=0
      NFLGC=0
      FNOC=0.0              ! pressure convergence failure indicator
      RCSQ=1.0/(RHOF*CSQ)
      IF (CSQ.LT.0.0) RCSQ=0.0   ! For imcompressible fluid
      IF (NMAT.EQ.1) RHOFC=RHOF
      RHOD=RHOF-RHOFC
      IF (CANGLE.EQ.90.0) CANGLE=CANGLE-EM6
      CANGLE=CANGLE*RPD
      TANCA=TAN(CANGLE)
      IPL=2
      IF (WL.EQ.5) IPL=3   ! 5: constant pressure
      IPR=IM1
      IF (WR.EQ.5) IPR=IM2
      JPB=2
      IF (WB.EQ.5) JPB=3
      JPT=JM1
      IF (WT.EQ.5) JPT=JM2
C
C     *** SET CONSTANT TERMS FOR PLOTTING
C
      XMIN=X(1)
      XMAX=X(IM1)
      IF (ISYMPLT.GT.0) XMIN=-XMAX
      YMIN=Y(1)
      YMAX=Y(JM1)
      D1=XMAX-XMIN
      D2=YMAX-YMIN
      D3=AMAX1(D1,D2)
      SF=1.0/D3                ! scaling factor
      XSHFT=0.5*(1.0-D1*SF)
      YSHFT=0.5*(1.0-D2*SF)
      DXMIN=EP10
      DO 10 I=2,IM1
   10 DXMIN=AMIN1(DELX(I),DXMIN)
      DYMIN=EP10
      DO 20 I=2,JM1
   20 DYMIN=AMIN1(DELY(I),DYMIN)
      VELMX1=AMIN1(DXMIN,DYMIN)/VELMX
C
C     *** DETERMINE SLOPED BOUNDARY LOCATION
C
C     *** COMPUTE INITIAL VOLUME FRACTION FUNCTION F IN CELLS
C
      DO 40 I=1,IMAX
      DO 30 J=2,JMAX
      F(I,J)=1.0
      IF (FLHT.GT.Y(J-1).AND.FLHT.LT.Y(J)) F(I,J)=RDY(J)*(FLHT-Y(J-1))
      IF (Y(J-1).GE.FLHT) F(I,J)=0.0
   30 CONTINUE
      F(I,1)=F(I,2)
   40 CONTINUE
      
      do i = 1, imax                   !!! Add by cwang!
          do j = 2,jmax
              if(x(i).gt.0.05)f(i,j)=0.0
          enddo
      enddo
      
 
      do i = 1, imax
          do j = 2, jmax
              phif(i,j) = 1.0; phifn(i,j)=1.0
              if(y(j).lt.0.05.and.x(i).lt.0.06)then
                  phif(i,j)=0.45; phifn(i,j)=0.45
              endif
          enddo
      enddo
      
      
C
C     *** GENERATE SPECIAL F-FUNCTION (FLUID) CONFIGURATION
C
C     *** CALCULATE DTVIS AND DTSFT
C
      DS=1.0E+10
      DTVIS=1.0E+10
      DTSFT=1.0E+10
      DO 50 I=2,IM1
      DO 50 J=2,JM1
      DXSQ=DELX(I)**2
      DYSQ=DELY(J)**2
      RDSQ=DXSQ*DYSQ/(DXSQ+DYSQ)
      RDSQ=RDSQ/(3.0*NU+1.0E-10)
      DTVIS=AMIN1(DTVIS,RDSQ)
      DS=AMIN1(DELX(I),DELY(J),DS)
   50 CONTINUE
      SIGX=SIGMA
      RHOMN=AMIN1(RHOF,RHOFC)
      IF(SIGX.EQ.0.0) SIGX=EM10
      DTM=SQRT(RHOMN*DS**3/(SIGX*4.0*(1.0+CYL)))
      DTSFT=AMIN1(DTSFT,DTM)
C
C     *** CALCULATE BETA(I,J) FOR MESH    ! pressure relaxation factor
C
      RDTEXP= 2.0*SQRT(ABS(CSQ))/DS
      IF(CSQ.LT.0.0) RDTEXP= 1.0E+10
      CTOS=DELT*RDTEXP
      COMG= AMIN1(CTOS**2,1.0)
      OMG1=(OMG-1.0)*COMG+1.0
      DO 55 I=2,IM1
      DO 55 J=2,JM1
      RHXR=(RHOFC+RHOD*F(I,J))*DELX(I+1)+(RHOFC+RHOD*F(I+1,J))*DELX(I)
      RHXL=(RHOFC+RHOD*F(I-1,J))*DELX(I)+(RHOFC+RHOD*F(I,J))*DELX(I-1)
      RHYT=(RHOFC+RHOD*F(I,J))*DELY(J+1)+(RHOFC+RHOD*F(I,J+1))*DELY(J)
      RHYB=(RHOFC+RHOD*F(I,J-1))*DELY(J)+(RHOFC+RHOD*F(I,J))*DELY(J-1)
      XX=DELT*RDX(I)*(2.0/RHXL+2.0/RHXR)+DELT*RDY(J)*(2.0/RHYT+2.0/RHYB)
      RHOR=RHOF/(RHOFC+RHOD*F(I,J))
      BETA(I,J)=OMG1/(XX*COMG+RCSQ*RHOR/DELT)
   55 CONTINUE
C
C     *** SET BETA(I,J)= -1.0 IN OBSTACLE CELLS
C         MUST BE DONE BY HAND IN GENERAL
C
C     *** PRINT BETA(I,J) ON FILM AND PAPER
C
      IF (IMOVY.EQ.1) GO TO 70
      WRITE (12,210)
      DO 60 J=1,JM1
      DO 60 I=1,IM1
      WRITE (12,220) I,J,BETA(I,J)
   60 CONTINUE
   70 CONTINUE
      WRITE (6,210)
      DO 80 J=1,JM1
      DO 80 I=1,IM1
      WRITE (6,220) I,J,BETA(I,J)
   80 CONTINUE
C
C     *** CALCULATE HYDROSTATIC PRESSURE
C
      DO 90 I=2,IM1
      P(I,JMAX)=0.0
      DO 90 J=2,JM1
      JJ=JM1-J+2
      RHOYA=(RHOFC+RHOD*F(I,JJ))*DELY(JJ)*0.5+(RHOFC+RHOD*F(I,JJ+1))
     1 *DELY(JJ+1)*0.5
      IF (NMAT.EQ.1) RHOYA=(AMIN1(F(I,JJ+1),0.5)*DELY(JJ+1)+AMAX1(0.0,F
     1 (I,JJ)-0.5)*DELY(JJ))*RHOF
      P(I,JJ)=P(I,JJ+1)-GY*RHOYA
   90 CONTINUE
C
C     *** PARTICLE SET UP
C
      NP=NPY*(1+NPX)
      IF (NP.EQ.0) GO TO 160
      DXP=(XPR-XPL)/FLOAT(NPX)
      DYP=(YPT-YPB)/FLOAT(NPY)
      K=0
      DO 100 JN=1,NPY,2
      DO 100 IN=1,NPX
      K=K+1
      XP(K)=XPL+(FLOAT(IN)-0.5)*DXP
      YP(K)=YPB+(FLOAT(JN)-1.0)*DYP
      IF (YP(K).GT.YPT) YP(K)=YPT
  100 CONTINUE
      DO 110 JN=2,NPY,2
      K=K+1
      XP(K)=XPL
      YP(K)=YPB+(FLOAT(JN)-1.0)*DYP
      IF (YP(K).GT.YPT) YP(K)=YPT
      DO 110 IN=1,NPX
      K=K+1
      XP(K)=XPL+FLOAT(IN)*DXP
      YP(K)=YPB+(FLOAT(JN)-1.0)*DYP
      IF (YP(K).GT.YPT) YP(K)=YPT
  110 CONTINUE
      NP=K
      DO 150 K=1,NP
      DO 120 I=2,IM1
      IF (XP(K).GE.X(I-1).AND.XP(K).LE.X(I)) IP(K)=I
      IF (X(I-1).GT.XPR) GO TO 130
  120 CONTINUE
  130 DO 140 J=2,JM1
      IF (YP(K).GE.Y(J-1).AND.YP(K).LE.Y(J)) JP(K)=J
      IF (Y(J-1).GT.YPT) GO TO 150
  140 CONTINUE
  150 CONTINUE
  160 CONTINUE
C
C     *** SET INITIAL SURFACE PRESSURE
C
      DO 170 J=2,JM1
      DO 170 I=2,IM1
      PS(I,J)=0.0
  170 CONTINUE
C
C     *** SET INITIAL VELOCITY FIELD INTO U AND V ARRAYS
C
      DO 180 I=2,IM1
      DO 180 J=2,JM1
      V(I,J)=VI
      U(I,J)=UI
      IF (F(I,J).GT.EMF.OR.NMAT.EQ.2) GO TO 180
      U(I,J)=0.0
      V(I,J)=0.0
  180 CONTINUE
C
C     *** SET INITIAL VOID REGION QUANTITIES
C
      DO 190 K=1,NVRM
      NR(K)=0
      PR(K)=0.0
  190 VOL(K)=0.0
  200 RETURN
C
  210 FORMAT (1H1)
  220 FORMAT (2X,5HBETA(,I2,1H,,I2,2H)=,1PE14.7)
      END
