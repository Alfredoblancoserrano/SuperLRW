!================================================================================
!Name        : module_GAM
!Version     : Beta 1.0
!Date        : agosto 03, 2020
!Authors     : Alfredo Blanco Serrano <alfredoblancoserrano@gmail.com>
!              Alfonso Allen-Perkins <alfonso.allen.perkins@gmail.com>
!              Roberto F. S. Andrade <randrade@ufba.br>
!================================================================================
!About this module: This compendium of subroutines are used for the purpose
!                   of estimating the mean square displacement or MSD as described
!                   in the main_CAME.f90 file. The subroutines contained in this
!                   module are listed  below:
!
! initialization  : This code allows the assignment of values to variables that
!                   do not depend on main loop.
! Eq35            : The MSD value is estimated following the Eq(35) derived
!                   in the supplementary material of the manuscript.
! Lambda          : Calculation of the eigenvalues of a helical tori following
!                   the Eq(15) derived in the supplementary material of the
!                   manuscript.
! Line_TH         : Calculate the distance between the first node and the
!                   remaining N-1 nodes.This calculation is performed with
!                   the Eq(15) of the manuscript.
! Heaviside       : Heaviside step function, or the unit step function.
! Num_Deri        : Calculation of the numerical derivative.
! Data_Reading    : Subroutine that reads the initial data from the
!                   file "in_CAME_data.dat".
! Write_Data      : Writing the results obtained.
!================================================================================
!===================================================================================
MODULE module_CAME
!==========================universal variables====================================
INTEGER,PARAMETER :: sp=4,dp=8
REAL(dp),PARAMETER::pi=4.d0*DATAN(1.d0)

REAL(dp),DIMENSION(:,:),ALLOCATABLE :: MSDT
REAL(dp),DIMENSION(:),ALLOCATABLE   :: MSD
REAL(dp),DIMENSION(:),ALLOCATABLE   :: Deri
REAL(dp),DIMENSION(:),ALLOCATABLE   :: alp
REAL(dp),DIMENSION(:),ALLOCATABLE   :: x1
REAL(dp),DIMENSION(:),ALLOCATABLE   :: x2
REAL(dp),DIMENSION(:),ALLOCATABLE   :: d
REAL(dp),DIMENSION(:),ALLOCATABLE   :: NV
REAL(dp),DIMENSION(:),ALLOCATABLE   :: la

REAL(dp)      :: satura,N,xt,flag_alpha
INTEGER(sp)   :: inN,ALLOCATESTATUS,tf,nm,NAS
CHARACTER(50) :: input,result,data_alpha

!===================================================================================
CONTAINS
!===================================================================================
SUBROUTINE Initialization
  IMPLICIT NONE
  INTEGER(sp)::i,k
  REAL(dp)::x
  REAL(dp) :: s
  REAL(dp),DIMENSION(:),ALLOCATABLE :: lambda_ini
  REAL(dp),DIMENSION(:),ALLOCATABLE :: Theta

!==============================Variable initialization==============================

  inN=INT(N)
  nm=INT((N-1.d0)/2.d0)

  ALLOCATE(MSDT(tf,2), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(MSD(tf), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(Deri(tf), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(alp(tf), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(x1(inN), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(x2(inN), STAT=ALLOCATESTATUS)
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

  x=2.d0;MSDT=1.d0
  DO i=2,tf
    MSDT(i,1)=x
    x=x+1.d0
  END DO
  !========================Reading alpha data=======================================
    IF(flag_alpha .LT. 0)THEN
        OPEN(UNIT=10,FILE=data_alpha,STATUS='unknown')
        DO i=1,6
          READ(10,*)
        END DO
        DO i=1,Nas
          READ(10,*)alp(i)
        END DO
        CLOSE(10)
    ELSE
        alp(1:tf)=flag_alpha
    END IF
!===================================================================================
!The vector NV has the distances between node i=1 and the remaining N-1 nodes
!in the TH.
!===================================================================================
  CALL Line_TH(NV)
  d(2:inN)=NV(2:inN)**(-alP(1))
!=======Initialization the Eq(30) of the supplementary material=====================
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
      x=x+NV(i)*NV(i)*dcos(theta(K)*dble(i-1))
    END DO
    X2(k)=x
  END DO
  DEALLOCATE(Theta, STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  DEALLOCATE(lambda_ini, STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
RETURN
END SUBROUTINE Initialization
!============Calculation of the Eq(30) of the supplementary material===============
SUBROUTINE Eq35(t,xt)
  IMPLICIT NONE
  INTEGER(sp),INTENT(IN) :: t
  REAL(dp), INTENT(OUT)  :: xt

  REAL(dp)    :: s
  INTEGER(sp) :: k,i

  xt=0.d0
  d(2:inN)=NV(2:inN)**(-alp(t))
  s=SUM(d(2:inN))
  la=0.d0

  DO k=1,inN-1
    DO i=1,nm
      la(k)=la(k)+d(i+1)*dcos((2.d0*pi*dble(k*i))/N)
    END DO
  END DO
  la=2.d0*la

  DO k=2,inN
    x1(k)=x1(k)*(la(k-1)/s)
  END DO

  DO k=2,inN
    xt=xt+x1(k)*x2(k)
  END DO
  xt=xt/N
RETURN
END SUBROUTINE Eq35
!===================================================================================
SUBROUTINE Lambda(l,x)
!==========================Eigenvalues of the TH====================================
!=====================See Eq(22) in the supplementary material======================
  IMPLICIT NONE
  INTEGER(sp),INTENT(IN)  :: l
  REAL(dp),   INTENT(OUT) :: x

  INTEGER(sp) :: i
  x=0.d0

  DO i=1,nm
    x=x+d(i+1)*dcos((2.d0*pi*dble(l*i))/N)
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
!================================================================================
SUBROUTINE Heaviside(x,y)
!===============================heaviside function===============================
  IMPLICIT NONE
  INTEGER(sp),INTENT(IN)  :: x
  INTEGER(sp),INTENT(OUT) :: y
  y=NINT(0.5*(sign(1,x)+1))
  RETURN
END SUBROUTINE Heaviside
!================================================================================
SUBROUTINE Num_Deri(A,B)
!===========================Log10 Numerical Derivative============================
  IMPLICIT NONE
  REAL(dp),DIMENSION(:,:),INTENT(IN)  :: A
  REAL(dp),DIMENSION(:),  INTENT(OUT) :: B

  INTEGER(sp) :: i
  REAL(dp)    :: yo,y1,xo,x1,m,d,c

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
!=================================================================================
SUBROUTINE Data_Reading
!================================================================================
  IMPLICIT NONE
  input="in_CAME_data.dat"
  OPEN(UNIT=10,FILE=input,STATUS='unknown')
  READ(10,*)N,NAS,tf
  READ(10,*)flag_alpha
  READ(10,*)data_alpha
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
END MODULE module_CAME
