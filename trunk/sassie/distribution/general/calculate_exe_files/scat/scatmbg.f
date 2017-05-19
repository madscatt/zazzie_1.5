      PROGRAM SCATMBG
C*******************************************************************
C     This program reads a set of coordinates (x,y,z,b) from one or
C     more specified files, where (x,y,z) are the cartesian
C     coordinates for the scatterer and b is the corresponding 
C     scattering length density. When the program has calculated the
C     distance distribution function, it is written to a file
C     (of extension *.PR) and Fourier transformed to give the
C     scattering profile, which is written to the specified file.
C     Original Program SCAT by Steen Hansen (1991)
C     Modified by Susan Krueger to accommodate larger arrays and to
C     allow multiple input (x,y,z,b) files (1997)
C     Modified by Susan Krueger to run under gfortran on a PC (2008)      
C********************************************************************
      INTEGER PNTOT,PNDAT
	PARAMETER (PNTOT=1000000)
      PARAMETER (PNDAT=1000001)
      DIMENSION B(PNTOT), X(PNTOT), Y(PNTOT), Z(PNTOT)
      DIMENSION COR(0:PNTOT),FINT(0:PNTOT)
      CHARACTER*40 ANAME(10),BNAME,TNAME,DNAME
      DATA COR/PNDAT*0./,FINT/PNDAT*0./
      REAL INTS
      PI=3.14159
C********************************************************************
C      Input
C********************************************************************
      I=-1
      WRITE(6,1)
    1 FORMAT(1X,'No. of input files         => ',$)
      READ(5,2)NFILES
    2 FORMAT(I2)
      DO I=1,NFILES
        WRITE(6,4)I
    4   FORMAT(1X,'Inputfile ',I2,'               => ',$)
        READ(5,5)ANAME(I)
    5   FORMAT(A)
      END DO
      WRITE(6,7)
    7 FORMAT(1X,'Outputfile                 => ',$)
      READ(5,6)BNAME
    6 FORMAT(A)
      TNAME = BNAME(1:INDEX(BNAME,'.')-1)
      CALL GETLENGTH(TNAME,L,40)
      DNAME = TNAME(1:L)//'.pr'
      WRITE(6,8)
    8 FORMAT(1X,'Approx no points in output => ',$)
      READ(5,*)NPOINTS
      WRITE(6,9)
    9 FORMAT(1X,'Qmax                       => ',$)
      READ(5,*)QMAX
      XMAX=-1.E30
      XMIN=-XMAX
      YMIN=XMIN
      YMAX=XMAX
      ZMIN=XMIN
      ZMAX=XMAX
      NTOT=0
      DO K=1,NFILES
        OPEN(UNIT=2,FILE=ANAME(K),STATUS='OLD')
   10   DO 15 I=NTOT+1,PNTOT
          READ(2,*,ERR=990,END=20)X(I),Y(I),Z(I),B(I)
          XMIN=MIN(XMIN,X(I))
          XMAX=MAX(XMAX,X(I))
          YMIN=MIN(YMIN,Y(I))
          YMAX=MAX(YMAX,Y(I))
          ZMIN=MIN(ZMIN,Z(I))
          ZMAX=MAX(ZMAX,Z(I))
          NTOT=NTOT+1
   15   CONTINUE
   20   WRITE(6,25)NTOT
      END DO
      CLOSE(UNIT=2)
      write(6,24)xmin,xmax,ymin,ymax,zmin,zmax
   24 FORMAT(1X,'Dimensions : ',6F8.1)
   25 FORMAT(1X,'No of scatterers              =   ',I5)
      RMAX=SQRT((XMAX-XMIN)**2+(YMAX-YMIN)**2+(ZMAX-ZMIN)**2)
      IF(NTOT .GT. 1000)THEN
        RSTEP=RMAX/1000.
      ELSE
        RSTEP=RMAX/NTOT
      ENDIF  
      WRITE(6,30)RMAX
   30 FORMAT(1X,'Radial Dimension = ',F8.2)
      WRITE(6,35)RSTEP
   35 FORMAT(1X,'Steplength for dd-function    = ',F10.4)
C*********************************************************************
C     Guinier-radius
C*********************************************************************
      CMX=0
      CMY=0
      CMZ=0
      BTOT=0
      DO 210 I=1,NTOT
        CMX=CMX+B(I)*X(I)
        CMY=CMY+B(I)*Y(I)
        CMZ=CMZ+B(I)*Z(I)
        BTOT=BTOT+B(I)
  210 CONTINUE
      CMX=CMX/BTOT
      CMY=CMY/BTOT
      CMZ=CMZ/BTOT
      WRITE(6,211)CMX,CMY,CMZ
  211 FORMAT(1X,'Center of (scattering-)mass   = ',3(F8.2))
      RG=0
      DO 220 I=1,NTOT
      RG=RG+B(I)*((X(I)-CMX)**2+(Y(I)-CMY)**2+(Z(I)-CMZ)**2)
  220 CONTINUE
      RG=SQRT(RG/BTOT)
      WRITE(6,221)RG
  221 FORMAT(1X,'Guinier-radius                = ',F8.2)
C*********************************************************************
C     Distance distribution function
C*********************************************************************
C      NDIV=JMOD(NTOT,2)
C      NTOT=NTOT-NDIV
      DO 101 J=1,NTOT-1
        DO 100 K=2,NTOT
          DXX=X(J)-X(K)
          DYY=Y(J)-Y(K)
          DZZ=Z(J)-Z(K)
          D=SQRT(DXX*DXX+DYY*DYY+DZZ*DZZ)
          L=NINT(D/RSTEP)
          COR(L)=COR(L)+B(J)*B(K)
  100   CONTINUE
  101 CONTINUE
      OPEN(UNIT=3,FILE=DNAME,STATUS='UNKNOWN')
      CORMAX=0.0
      NMAX=0
      DO 105 I=1,NTOT
        IF(COR(I).NE.0) NMAX=I
        IF(COR(I).GE.CORMAX)CORMAX=COR(I)
  105 CONTINUE
C      TYPE *,'CORMAX = ',CORMAX
      WRITE(6,130)NMAX*RSTEP
  130 FORMAT(1X,'Maximum length of molecule    = ',F10.3)
      COR(0)=0.0
      NTOT=NMAX
      ZERO=0.0
      DO 200 J=0,NMAX  
        WRITE(3,*)J*RSTEP,COR(J)/CORMAX,ZERO
  200 CONTINUE
      CLOSE(UNIT=3)
C*********************************************************************
C     Fouriertransformation => scattering profile
C*********************************************************************
      WRITE(6,222)QMAX
  222 FORMAT(1X,'Qmax at Fouriertransformation =   ',F7.4)
      QSTEP=QMAX/NPOINTS
      CALL FT(COR,RSTEP,FINT,QSTEP,NMAX,NPOINTS)
C*********************************************************************
C     Output
C*********************************************************************
      OPEN(UNIT=3,FILE=BNAME,STATUS='UNKNOWN')
      DO 500 I=0,NPOINTS
        WRITE(3,*)I*QSTEP,FINT(I)/FINT(0),ZERO
  500 CONTINUE
      CLOSE(UNIT=3)
C**********************************************************************
C     Calculation of Guinier radius from scattering profile 
C**********************************************************************
      NGUI=0
      DO 510 I=1,NTOT
        IF(I*QSTEP*RMAX.GE.1.) GOTO 511
        NGUI=NGUI+1
  510 CONTINUE
  511 IF(NGUI .LT. 3) GOTO 999
      SUM=0
      DO 515 I=0,NGUI-1
        SUM1=LOG(FINT(I+1)/FINT(I))/((I+1)**2-I**2)
        SUM=SUM+SUM1
  515 CONTINUE
      RGC=0
      ALPHA=-ABS(SUM)/(QSTEP**2*NGUI)
      RGC=SQRT(-ALPHA*3)
      RQ=RGC*NGUI*QSTEP
      WRITE(6,521)RGC
  521 FORMAT(1X,'Guinier-radius from scattering profile = ',F8.2)
      WRITE(6,522)NGUI+1,RQ
  522 FORMAT(1X,'Calculated from',I5,'  points with R*Q <',F7.2)
      GOTO 999
 
  990 WRITE(6,*)'ERROR AT I=',I
 
  999 STOP
      END
 
**********************************************************************
*     Fourier transformation
**********************************************************************
      SUBROUTINE FT(COR,RSTEP,FINT,QSTEP,NTOT,NPOINTS)
      REAL COR(0:NTOT),FINT(0:NTOT)
      DO 400 K=1,NPOINTS                              ! K>1
        Q=K*QSTEP
        FADD=0.
        DO 300 I=1,NTOT                               ! I>1
          R=I*RSTEP
          DEBEYE=SIN(Q*R)/(Q*R)
          FADD=COR(I)*DEBEYE
          FINT(K)=FINT(K)+FADD
  300   CONTINUE
  400 CONTINUE
      DO 500 I=1,NTOT
        FINT(0)=FINT(0)+COR(I)
  500 CONTINUE
      RETURN
      END
c      options/extend_source
      subroutine getlength(string,length,maxlength)

      implicit none

      character*(*) string

      integer i,j,length,maxlength

      do i=1,maxlength
       j=ichar(string(i:i))
c        type *,' i,char: ',i,j
       if(j .ne. 32 .and. j .ne. 0)then
         length=i
       end if
      end do

c      type *,' length: ',length

      RETURN

      END




