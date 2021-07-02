!================================================================================
!Name        : module_CAMS
!Version     : Beta 1.0
!Date        : agosto 03, 2020
!Authors     : Alfredo Blanco Serrano <alfredoblancoserrano@gmail.com>
!              Alfonso Allen-Perkins  <alfonso.allen.perkins@gmail.com>
!              Roberto F. S. Andrade  <randrade@ufba.br>
!================================================================================
!About this module: This compendium of subroutines are used for the purpose
!                   of estimating the mean square displacement or MSD as described
!                   in the main_CAMS.f90 file. The subroutines contained in this
!                   module are listed  below:
!
! Initialization  : This code allows the assignment of values to variables that
!                   do not depend on main loop.
! PV              : Calculation of the probability vector considering that the
!                   HT is a circulating matrix and therefore the number of loops
!                   is considerably reduced.
! Line_TH         : Calculate the distance between the first node and the
!                   remaining N-1 nodes.This calculation is performed with
!                   the Eq(15) of the manuscript.
! Heaviside       : Heaviside step function, or the unit step function.
! Num_Deri        : Calculation of the numerical derivative.
! Data_Reading    : Subroutine that reads the initial data from the
!                   file "in_CAME_data.dat".
! Write_Data      : Writing the results obtained.
!================================================================================
MODULE module_CAMS
!==========================universal variables====================================
INTEGER, PARAMETER :: sp=4,dp=8

REAL(dp),DIMENSION(:,:),ALLOCATABLE :: MSDT
REAL(dp),DIMENSION(:),ALLOCATABLE   :: MSD
REAL(dp),DIMENSION(:),ALLOCATABLE   :: Deri
REAL(dp),DIMENSION(:),ALLOCATABLE   :: alp
REAL(dp),DIMENSION(:),ALLOCATABLE   :: NV
REAL(dp),DIMENSION(:),ALLOCATABLE   :: SNV
REAL(dp),DIMENSION(:),ALLOCATABLE   :: P
REAL(dp),DIMENSION(:),ALLOCATABLE   :: x
REAL(dp),DIMENSION(:),ALLOCATABLE   :: TMPD

REAL(dp)      :: flag_alpha,N
INTEGER(sp)   :: ALLOCATESTATUS,tf,NAS,inN,dmax
CHARACTER(50) :: input,result,data_alpha
!===================================================================================
CONTAINS
!===================================================================================
SUBROUTINE Initialization
  IMPLICIT NONE
  INTEGER(sp) :: i
  REAL(dp)    :: cont
  inN=INT(N)
!==============================Variable initialization==============================
   ALLOCATE(MSDT(tf,2), STAT=ALLOCATESTATUS)
   IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
   ALLOCATE(MSD(tf), STAT=ALLOCATESTATUS)
   IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
   ALLOCATE(Deri(tf), STAT=ALLOCATESTATUS)
   IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
   ALLOCATE(alp(tf), STAT=ALLOCATESTATUS)
   IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
   ALLOCATE(NV(inN), STAT=ALLOCATESTATUS)
   IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
   ALLOCATE(SNV(inN), STAT=ALLOCATESTATUS)
   IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
   ALLOCATE(P(inN), STAT=ALLOCATESTATUS)
   IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
   ALLOCATE(x(inN), STAT=ALLOCATESTATUS)
   IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
   ALLOCATE(TMPD(inN), STAT=ALLOCATESTATUS)
   IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"

   MSD=0.d0;MSDT=1.d0
   dmax=INT((n+1.0)/2.0)
!===================================================================================
!The vector NV has the distances between node i=1 and the remaining N-1 nodes
!in the TH.
!===================================================================================
   CALL LINE_TH(NV)
   SNV=NV*NV
!========================Reading alpha data=========================================
  IF(flag_alpha .LT. 0)THEN
    OPEN(UNIT=10,FILE=data_alpha,STATUS='unknown')
    READ(10,*)
    READ(10,*)
    DO i=1,NAS
      READ(10,*)alp(i+1)
    END DO
    alp(1)=1000
    CLOSE(10)
  ELSE
    alp(1:tf)=flag_alpha
  END IF
!===================================================================================

  cont=1.d0
  DO i=1,tf
    MSDT(i,1)=cont
    cont=cont+1.d0
  END DO

  TMPD(2:inN)=(NV(2:inN))**(-1.d0*alp(1))
  TMPD(2:inN)=(TMPD(2:inN))/(SUM(TMPD(2:inN)))
  P=TMPD
   RETURN
END SUBROUTINE Initialization
!======================Probability vector============================================
SUBROUTINE PV(t,P)
IMPLICIT NONE
INTEGER(sp),          INTENT(IN)     :: t
REAL(dp),DIMENSION(:),INTENT(IN OUT) :: P

INTEGER(sp):: i,j

x=0.d0

TMPD(2:inN)=(NV(2:inN))**(-1.d0*alp(t))
TMPD(2:inN)=(TMPD(2:inN))/(SUM(TMPD(2:inN)))

DO i=1,dmax
  DO j=1,i
    x(i)=x(i)+P(j)*TMPD(i-j+1)
  END DO
END DO

DO i=1,dmax
  DO j=1,inN-i
    x(i)=x(i)+P(i+j)*TMPD(inN-j+1)
  END DO
END DO

p(1:dmax)=x(1:dmax)
DO i=1,dmax-1
  p(dmax+i)=x(dmax-i+1)
END DO
RETURN
END SUBROUTINE PV
!===================================================================================
SUBROUTINE Line_TH(NV)
!===Estimation of the distance between the first node and the other N-1 nodes.======
!=====================See Eq(15) in the manuscript==================================
  IMPLICIT NONE
  REAL(dp),DIMENSION(:),INTENT(OUT) :: NV

  REAL(sp)    :: x,jj1
  INTEGER(SP) :: JO2,JO,j,h,h2,sj,xj
  INTEGER(SP) :: JJ0,X1,INSQRT
  INSQRT=INT(SQRT(N))

  DO j=2,inN

    JJ0=j-((inN+3)/2)
    jj1=((Real(inN)+2.0)/2.0)-real(j)
    jj1=sign(1.0,jj1)
    CALL Heaviside(JJ0,X1)
    xj=(inN+2)*X1+(j)*INT(jj1)-1

    sj=inN-j+1
    JO=MOD(xj+1,INSQRT)
    JO2=JO-NINT(((REAL(INSQRT)+3.0)/2.0))
    CALL Heaviside(-JO,H);CALL Heaviside(JO2,h2)

    x=(REAL(xj)/REAL(INSQRT))+((REAL(JO)-1.0)*((REAL(INSQRT)-1.0)/REAL(INSQRT)))

    NV(j)=(NINT(x)+2*H-2*h2*JO2)
  END DO
  RETURN
END SUBROUTINE Line_TH
!==================================================================================
SUBROUTINE Heaviside(x,y)
!===============================heaviside function=================================
  IMPLICIT NONE
  INTEGER(sp),INTENT(IN)  :: x
  INTEGER(sp),INTENT(OUT) :: y
  y=NINT(0.5*(sign(1,x)+1))
  return
END SUBROUTINE Heaviside
!=================================================================================
SUBROUTINE Num_Deri(A,B)
!===========================Log10 Numerical Derivative============================
  IMPLICIT NONE
  REAL(dp),DIMENSION(:,:),INTENT(IN)  :: A
  REAL(dp),DIMENSION(:),  INTENT(OUT) :: B

  REAL(dp)    :: yo,y1,xo,x1,m,d,c
  INTEGER(sp) :: i

  DO i=1,tf-1
    yo=dlog10(A(i,2))
    y1=dlog10(A(i+1,2))

    xo=dlog10(A(i,1))
    x1=dlog10(A(i+1,1))

    d=(y1-yo);c=(x1-xo)
    m=d/c
    B(i)=m
  END DO
  B(tf)=B(tf-1)
  RETURN
END SUBROUTINE Num_Deri
!====================================================================================
SUBROUTINE Data_Reading
!====================================================================================
  IMPLICIT NONE
  input="in_CAMS_data.dat"
  OPEN(UNIT=10,FILE=input,STATUS='unknown')
  READ(10,*)N,NAS,tf
  READ(10,*)flag_alpha
  READ(10,*)data_alpha
  READ(10,*)result
  CLOSE(10)
  RETURN
END SUBROUTINE Data_Reading
!===================================================================================
SUBROUTINE Write_Data(secs)
!================================================================================
  IMPLICIT NONE
  REAL(sp), INTENT(IN) :: secs
  INTEGER(sp) :: i

  OPEN(UNIT=10,FILE=result,STATUS='unknown')
  WRITE(10,*)'#==================================================#'
  WRITE(10,502) "#",'time= ',secs, '[s]','input= ',input,'#'
  WRITE(10,*)'#==================================================#'
  WRITE(10,*)
  WRITE(10,*)"#",'  Log10(t)','       Log10(MSD)', '      Derivative'
  DO i=1,tf
   WRITE(10,501) dlog10(MSDT(i,1)),dlog10(MSDT(i,2)),deri(i)
  END DO
  CLOSE(10)
  501 FORMAT(2x,f12.10,4x,f12.10,4x,f12.10)
  502 FORMAT(1x,a,a6,f11.7,a3,2x,a7,a20,1x,a)
  RETURN
END SUBROUTINE Write_Data
END MODULE module_CAMS
