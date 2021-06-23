!================================================================================
!Name        : module_GAM
!Version     : Beta 1.0
!Date        : agosto 03, 2020
!Authors     : Alfredo Blanco Serrano <alfredoblancoserrano@gmail.com>
!              Alfonso Allen-Perkins <alfonso.allen.perkins@gmail.com>
!              Roberto F. S. Andrade <randrade@ufba.br>
!================================================================================
!About this module  : This module contains all the necessary subroutines
!                     to solve the inverse problem as described in the
!                     Main_alpha GAM.f90. In the following, the subroutines
!                     contained in this module are described.
!
! Initialization    : This code allows the assignment of values to variables that
!                     do not depend on main loop.
! GAM_MSD           : Calculation of the mean square displacement.
! Transision_Matrix : The transission matrix(S) is formed from the Mellin
!                     transformed d-path adjacency matrix(TM).
! Bisection         : root-finding method.
! MELLI             : Mellin transformed d-path adjacency matrix(TM)
! Heaviside         : Heaviside step function, or the unit step function.
! TN                : Calculation of the neighborhood matrix and the diameter for
!                     a TN of size N odd.
! LINE_TN & dis_TN  : These codes are used by the TN algorithm to estimate the
!                     neighborhood matrix.
! Data_Reading      : Subroutine that reads the initial data from the
!                     file "in_GAM_data.dat".
! Write_Data        : Writing the results obtained.
!================================================================================
MODULE module_alpha_GAM
!==========================universal variables====================================
INTEGER,PARAMETER :: sp=4,dp=8

REAL(dp),   DIMENSION(:,:),ALLOCATABLE :: TM
REAL(dp),   DIMENSION(:,:),ALLOCATABLE :: S
REAL(dp),   DIMENSION(:,:),ALLOCATABLE :: PTR
REAL(dp),   DIMENSION  (:),ALLOCATABLE :: alp
REAL(dp),   DIMENSION  (:),ALLOCATABLE :: EMSD
INTEGER(sp),DIMENSION(:,:),ALLOCATABLE :: ADJA
INTEGER(sp),DIMENSION(:,:),ALLOCATABLE :: MC

REAL(dp)      :: alp2,gam,tol
INTEGER(sp)   :: N,tf,dmax,ALLOCATESTATUS,flag_T
CHARACTER(30) :: input,result

!================================================================================
CONTAINS
!================================================================================
SUBROUTINE Initialization
  IMPLICIT NONE
  REAL(dp)    :: x
  INTEGER(sp) :: i
!========================Initialization of required variables=====================
  ALLOCATE(ADJA(N,N), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
  ALLOCATE(S(N,N), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
  ALLOCATE(TM(N,N), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
  ALLOCATE(PTR(N,N), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
  ALLOCATE(MC(N,N), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
  ALLOCATE(EMSD(1000), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
  ALLOCATE(alp(1000), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
  x=2.d0
!========================Initialization of required variables=====================

  IF(flag_T .GT. 0)THEN
    CALL TH(N,MC,ADJA,dmax)
  ELSE
    CALL TN(N,MC,ADJA,dmax)
  END IF
  PTR=REAL(ADJA)

  DO i=2,1000
    EMSD(i)=x**(gam)
    x=x+1.d0
  END DO
  RETURN
END SUBROUTINE Initialization
!================================================================================
SUBROUTINE GAM_MSD(PTI,PTRR,c)
  IMPLICIT NONE
  REAL(dp),DIMENSION(:,:),INTENT(IN OUT) :: PTI
  REAL(dp),DIMENSION(:,:),INTENT(IN OUT) :: PTRR
  REAL(dp)               ,INTENT   (OUT) :: c

  REAL(dp)    :: h
  INTEGER(sp) :: i,j

  CALL MELLI(ADJA,MC,dmax,alp2,TM)
  CALL Transision_matrix(TM,PTI)
  PTRR=MATMUL(PTR,PTI)
  h=0.d0
  DO i=1,n
    DO j=1,n
      h=h+(MC(i,j)*MC(i,j)*PTRR(j,i))
    END DO
  END DO
  c=h/REAL(N)
  RETURN
END SUBROUTINE GAM_MSD
!================================================================================
SUBROUTINE Transision_Matrix(TM,S)
!======================= stochastic transition matrix============================
  IMPLICIT NONE
  REAL(dp), DIMENSION(:,:), INTENT(IN) :: TM
  REAL(dp), DIMENSION(:,:), INTENT(OUT) :: S

  INTEGER(sp) :: i
  REAL(dp) :: k
  REAL(dP), DIMENSION(:,:), ALLOCATABLE :: B

  ALLOCATE(B(n,n), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"

  B=0.d0;S=0.d0

  DO i=1,n
    k=SUM(TM(i,1:n))
    S(i,1:n)=TM(i,1:n)/k
  END DO

  DO i=1,n
    B(1:n,i)=S(i,1:n)
  END DO
  S(1:n,1:n)=B(1:n,1:n)

  DEALLOCATE(B, STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  RETURN
END SUBROUTINE Transision_Matrix
!================================================================================]
SUBROUTINE Bisection(control,min,max,me)
  IMPLICIT NONE
  REAL(dp),INTENT(IN)     :: control
  REAL(dp),INTENT(IN OUT) :: min,max,me

  IF (control .GT. 0)THEN
    max=me
    me=max-ABS(max-min)/2.d0
  ELSE
    min=me
    me=min+ABS(max-min)/2.d0
  END IF
  RETURN
END SUBROUTINE Bisection
!================================================================================
!================================================================================
SUBROUTINE MELLI(A,MC,dmax,Alpha,TM)
!================Mellin transformed d−path adjacency matrix======================
  IMPLICIT NONE
  INTEGER(sp),DIMENSION(:,:),INTENT(IN)  :: A
  INTEGER(sp),DIMENSION(:,:),INTENT(IN)  :: MC
  INTEGER(sp)               ,INTENT(IN)  :: dmax
  REAL(dp)                  ,INTENT(IN)  :: Alpha
  REAL(dp),DIMENSION(:,:)   ,INTENT(OUT) :: TM

  REAL(dp),ALLOCATABLE,DIMENSION(:,:) :: B,DUMMY
  REAL(dp)    :: d
  INTEGER(sp) :: Sample,i,j,k

  ALLOCATE(B(n,n), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(DUMMY(n,n), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"

  B(1:n,1:n)=REAL(A(1:n,1:n))
  DO i=2,dmax
    DUMMY(1:n,1:n)=0.d0
    DO j=1,n
      DO k=1,n
        Sample=MC(j,k)
        IF(Sample .EQ. i )DUMMY(j,k)=1.d0
      END DO
    END DO
    d=1.d0/REAL(i)**(alpha)

    B(1:n,1:n)=B(1:n,1:n)+d*DUMMY(1:n,1:n)
  END DO

  TM(1:n,1:n)=B(1:n,1:n)
  DEALLOCATE(B, STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  DEALLOCATE(DUMMY, STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  RETURN
END SUBROUTINE MELLI
!=================================================================================
SUBROUTINE TH(n,MC,A,dmax)
!===============Adjacency and neighborhood matrix of the TH=======================
  IMPLICIT NONE

  INTEGER(sp)               ,INTENT(IN)  :: n
  INTEGER(sp),DIMENSION(:,:),INTENT(OUT) :: MC
  INTEGER(sp),DIMENSION(:,:),INTENT(OUT) :: A
  INTEGER(sp)               ,INTENT(OUT) :: dmax

  INTEGER(sp),DIMENSION(:,:),ALLOCATABLE :: B
  INTEGER(sp),DIMENSION(:,:),ALLOCATABLE :: E
  INTEGER(sp),DIMENSION(:)  ,ALLOCATABLE :: dist

  REAL(sp)    :: x,jj1
  INTEGER(SP) :: JO2,JO,j,h,h2,sj,xj,i,d1,d2
  INTEGER(SP) :: JJ0,X1,nq,INSQRT,inN

  ALLOCATE(B(n,n), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(E(n,n), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(dist(n), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  inN=N;nq=inN;B=0;dist=0;E=0;d2=1
  INSQRT=INT(SQRT(REAL(nq)))

  DO j=1,nq

    JJ0=j-((NQ+3)/2)
    jj1=((Real(nq)+2.0)/2.0)-real(j)
    jj1=sign(1.0,jj1)
    CALL Heaviside(JJ0,X1)
    xj=(nq+2)*X1+(j)*INT(jj1)-1

    sj=inN-j+1
    JO=MOD(xj+1,INSQRT)
    JO2=JO-NINT(((REAL(INSQRT)+3.0)/2.0))
    CALL Heaviside(-JO,H);CALL Heaviside(JO2,h2)

    x=(REAL(xj)/REAL(INSQRT))+((REAL(JO)-1.0)*((REAL(INSQRT)-1.0)/REAL(INSQRT)))

    dist(j)=(NINT(x)+2*H-2*h2*JO2)
  END DO

  DO i=1,nq
    DO j=i,nq
      B(i,j)=dist(j-i+1)
    END DO
  END DO

  DO i=1,nq
    E(i:nq,i)=B(i,i:nq)
  END DO
  E=B+E
  B=0
  DO i=1,nq
    DO j=1,nq
      d1=E(i,j)
      IF(d1 .EQ. 1)B(i,j)=1
      IF(d1 .GT. d2)d2=d1
    END DO
  END DO
  A=B;MC=E;dmax=d2
  DEALLOCATE(B, STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  DEALLOCATE(E, STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  DEALLOCATE(dist, STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  RETURN
END SUBROUTINE TH
!================================================================================
SUBROUTINE Heaviside(x,y)
!===============================heaviside function===============================
  IMPLICIT NONE
  INTEGER(sp),INTENT(IN)  :: x
  INTEGER(sp),INTENT(OUT) :: y
  y=NINT(0.5*(sign(1,x)+1))
  return
END SUBROUTINE Heaviside
!================================================================================
SUBROUTINE TN(n,A,C,dmax)
!===============Adjacency and neighborhood matrix of the TN=======================
  IMPLICIT NONE
  INTEGER(sp)                ,INTENT(IN)  :: n
  INTEGER(sp), DIMENSION(:,:),INTENT(OUT) :: A
  INTEGER(sp), DIMENSION(:,:),INTENT(OUT) :: c
  INTEGER(sp)                ,INTENT(OUT) :: dmax

  INTEGER(sp),DIMENSION(:,:),ALLOCATABLE :: B
  INTEGER(sp),DIMENSION(:,:),ALLOCATABLE :: E

  INTEGER(sp) :: i,j,INn,d,d1

  INn=INT(SQRT(REAL(n)))

  ALLOCATE(B(n,n), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(E(n,n), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"

  B=0;E=0;d1=1
  DO i=1,n
    DO j=1,n
      CALL LINE_TN(i,j,n,d)
      B(i,j)=d
      IF(d .GT. d1)d1=d
    END DO
  END DO

  DO i=1,n
    DO j=1,n
      d=B(i,j)
      IF (d .EQ. 1)E(i,j)=1
    END DO
  END DO

  A=B;C=E;dmax=d1
  DEALLOCATE(B, STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  DEALLOCATE(E, STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  RETURN
END SUBROUTINE TN
!================================================================================
SUBROUTINE Line_TN(a,b,n,D)
!========Calculation of the first row of the neighborhood matrix=================
  IMPLICIT NONE
  INTEGER(sp),INTENT(IN)  :: a,b,n
  INTEGER(sp),INTENT(OUT) :: D

  INTEGER(sp) :: x1,x2,xmax1,xmin1,xmax2,xmin2
  INTEGER(sp) :: a1,i,j1,j2,a2
  INTEGER(sp) :: D1,D2,Rn


  Rn=INT(SQRT(REAL(n)))
  x1=CEILING(REAL(a)/REAL(Rn));x2=CEILING(REAL(b)/REAL(Rn))
  xmax1=x1*Rn;xmin1=Rn*(x1-1)+1
  xmax2=x2*Rn;xmin2=Rn*(x2-1)+1

  j1=xmin1
  j2=xmin2

  DO i=1,Rn
    IF(a .EQ. j1)THEN
      a1=i
      EXIT
    END IF
    j1=xmin1+i
  END DO

  DO i=1,Rn
    IF(b .EQ. j2)THEN
      a2=i
      EXIT
    END IF
    j2=xmin2+i
  END DO
  CALL dis_TN(a1,a2,n,D1)
  CALL dis_TN(x1,x2,n,D2)
  D=D1+D2
  RETURN
END SUBROUTINE Line_TN
!================================================================================
SUBROUTINE dis_TN(a,b,n,D)
!=====================Distance between two nodes in the TN=======================
IMPLICIT NONE
INTEGER(sp), INTENT(IN) :: a,b
INTEGER(sp), INTENT(IN) :: n
INTEGER(sp),INTENT(OUT) :: D

REAL(sp) :: D1
INTEGER(sp) :: fun,x1,Rn

Rn=INT(SQRT(REAL(n)))

x1=abs(a-b)
D1=(REAL(Rn)/2)
IF(x1 .LT. D1)THEN
  fun=x1
ELSE
  fun=Rn-x1
END IF

D=fun
RETURN
END SUBROUTINE dis_TN
!================================================================================
SUBROUTINE Data_Reading
!================================================================================
  IMPLICIT NONE
  input="in_GAM_data.dat"
  OPEN(UNIT=10,FILE=input,STATUS='unknown')
  READ(10,*)N,tol,gam
  READ(10,*)flag_T
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
  502 FORMAT(1x,a,a6,f14.7,a3,2x,a7,a20)
  RETURN
END SUBROUTINE Write_Data
END MODULE module_alpha_GAM
!================================================================================
