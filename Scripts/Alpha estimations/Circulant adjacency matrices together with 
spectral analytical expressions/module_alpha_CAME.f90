!================================================================================
!Name        : module_alpha_CAME
!Version     : Beta 1.0
!Date        : agosto 03, 2020
!Authors     : Alfredo Blanco Serrano <alfredoblancoserrano@gmail.com>
!              Alfonso Allen-Perkins <alfonso.allen.perkins@gmail.com>
!              Roberto F. S. Andrade <randrade@ufba.br>
!================================================================================
!About this module: This compedium of subroutines are used for the purpose of
!                   estimating the set of alpha values that generate a
!                   superdiffusive spectrum by using Lévi random walks (LRW)
!                   in the helical tori (TH). The subroutines contained in this
!                   module are listed  below:
!
! Initialization  : This code allows the assignment of values to variables that
!                   do not depend on main loop.
! Eq35            : The MSD value is estimated following the Eq(35) derived
!                   in the supplementary material of the manuscript.
! Bisection       : root-finding method.
! Lambda          : Calculation of the eigenvalues of a helical tori following
!                   the Eq(15) derived in the supplementary material of the
!                   manuscript.
! Line_TH         : Calculate the distance between the first node and the
!                   remaining N-1 nodes.This calculation is performed with
!                   the Eq(15) of the manuscript.
! Heaviside       : Heaviside step function, or the unit step function.
! Data_Reading    : Subroutine that reads the initial data from the
!                   file "in_CAME_data.dat".
! Write_Data      : Writing the results obtained.
!================================================================================
!===================================================================================
MODULE module_alpha_CAME
!==========================universal variables====================================
INTEGER,PARAMETER :: sp=4,dp=8
REAL(dp),PARAMETER::pi=4.d0*DATAN(1.d0)

REAL(dp),DIMENSION(:),ALLOCATABLE   :: EMSD
REAL(dp),DIMENSION(:),ALLOCATABLE   :: alp
REAL(dp),DIMENSION(:),ALLOCATABLE   :: x1
REAL(dp),DIMENSION(:),ALLOCATABLE   :: x2
REAL(dp),DIMENSION(:),ALLOCATABLE   :: x3
REAL(dp),DIMENSION(:),ALLOCATABLE   :: d
REAL(dp),DIMENSION(:),ALLOCATABLE   :: NV
REAL(dp),DIMENSION(:),ALLOCATABLE   :: la

INTEGER(sp)   :: inN,ALLOCATESTATUS,tf,nm
REAL(dp)      :: satura,N,xt,alp2,gam,tol
CHARACTER(50) :: input,result

!===================================================================================
CONTAINS
!===================================================================================
SUBROUTINE Initialization
  IMPLICIT NONE

  REAL(dp),DIMENSION(:),ALLOCATABLE :: lambda_ini
  REAL(dp),DIMENSION(:),ALLOCATABLE :: Theta

  REAL(dp)    :: s,x
  INTEGER(sp) :: i,k

!==============================Variable initialization==============================

  inN=INT(N)
  nm=INT((N-1.d0)/2.d0)

  ALLOCATE(EMSD(500), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
  ALLOCATE(alp(500), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(x1(inN), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(x2(inN), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(x3(inN), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(d(inN), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(NV(inN), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(la(inN), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(lambda_ini(inN), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(Theta(inN), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"

  x=2.d0;alP2=1000;tf=500
  DO i=2,500
    EMSD(i)=x**(gam)
    x=x+1.d0
  END DO
!===================================================================================
!The vector NV has the distances between node i=1 and the remaining N-1 nodes
!for a TH.
!===================================================================================
  CALL Line_TH(NV)
  d(2:inN)=NV(2:inN)**(-alP2)
!=======Initialization the Eq(35) of the supplementary material=====================
  satura=((N-1)*(7*N-3))/(24*N)
  x=2.d0

  DO i=2,inN
    theta(i)=(2*pi*(x-1.d0))/N
    x=x+1.d0
  END DO

  s=SUM(d(2:inN))

  DO i=2,inN
    CALL Lambda(i-1,x)
    lambda_ini(i)=x/s
  END DO
  x1=lambda_ini

  DO k=2,inN
    X=0.d0
    DO i=2,inN
      x=x+NV(i)*NV(i)*DCOS(theta(K)*DBLE(i-1))
    END DO
    X2(k)=x
  END DO

  DEALLOCATE(Theta, STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  DEALLOCATE(lambda_ini, STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
RETURN
END SUBROUTINE Initialization
!============Calculation of the Eq(35) of the supplementary material===============
SUBROUTINE Eq35(x)
  IMPLICIT NONE
  REAL(dp), INTENT(OUT) :: x

  REAL(dp)    :: s
  INTEGER(sp) :: k,i

  x=0.d0
  d(2:inN)=NV(2:inN)**(-alp2)
  s=SUM(d(2:inN))
  x3=x1
  la=0.d0

  DO k=1,inN-1
    DO i=1,nm
      la(k)=la(k)+d(i+1)*DCOS((2.d0*pi*DBLE(k*i))/N)
    END DO
  END DO
  la=2.d0*la

  DO k=2,inN
    x3(k)=x3(k)*(la(k-1)/s)
  END DO

  DO k=2,inN
    x=x+x3(k)*x2(k)
  END DO
  x=satura+(x/N)
RETURN
END SUBROUTINE Eq35
!===================================================================================
!=======================================================================================
SUBROUTINE Bisection(control,min,max,me)
  IMPLICIT NONE
  REAL(dp),INTENT(IN OUT) :: min,max,me
  REAL(dp), INTENT(IN)    :: control

  IF (control .GT. 0)THEN
    max=me
    me=max-ABS(max-min)/2.d0
  ELSE
    min=me
    me=min+ABS(max-min)/2.d0
  END IF
  RETURN
END SUBROUTINE Bisection
!===================================================================================
SUBROUTINE Lambda(l,x)
!==========================Eigenvalues of the TH====================================
!=====================See Eq(27) in the supplementary material======================
  IMPLICIT NONE
  REAL(dp),INTENT(OUT)   :: x
  INTEGER(sp),INTENT(IN) :: l

  INTEGER(sp) i
  x=0.d0

  DO i=1,nm
    x=x+d(i+1)*DCOS((2.d0*pi*DBLE(l*i))/N)
  END DO
  x=2.D0*x
RETURN
END SUBROUTINE Lambda
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
    jj1=((REAL(inN)+2.0)/2.0)-REAL(j)
    jj1=SIGN(1.0,jj1)
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
!================================================================================
SUBROUTINE Heaviside(l,m)
!===============================heaviside function===============================
  IMPLICIT NONE
  INTEGER(sp),INTENT(IN)   :: l
  INTEGER(sp), INTENT(OUT) :: m

  m=NINT(0.5*(SIGN(1,l)+1))
  RETURN
END SUBROUTINE Heaviside
!================================================================================
SUBROUTINE Data_Reading
!================================================================================
  IMPLICIT NONE
  input="in_CAME_data.dat"
  OPEN(UNIT=10,FILE=input,STATUS='unknown')
  READ(10,*)N,gam,tol
  READ(10,*)result
  CLOSE(10)
  RETURN
END SUBROUTINE Data_Reading
!================================================================================
SUBROUTINE Write_Data(secs)
!================================================================================
  IMPLICIT NONE
  REAL(sp), INTENT(IN) :: secs

  INTEGER(sp) :: i

  OPEN(UNIT=10,FILE=result,STATUS='unknown')
  WRITE(10,*)'#================================================================'
  WRITE(10,502) "#",'time= ',secs, '[s]','input= ',input
  WRITE(10,*)'#================================================================'
  WRITE(10,*)
  WRITE(10,*)'#   Alpha(t)'
  WRITE(10,*)
  DO i=1,tf
   WRITE(10,501) alp(i)
  END DO
  CLOSE(10)
  501 FORMAT(2x,f12.7)
  502 FORMAT(1x,a,a6,f12.7,a3,2x,a7,a20)
  RETURN
END SUBROUTINE Write_Data
END MODULE module_alpha_CAME
