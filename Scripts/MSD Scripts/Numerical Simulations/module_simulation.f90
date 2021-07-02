!================================================================================
!Name        : Module_simulation.f90
!Version     : Beta 1.0
!Authors     : Alfredo Blanco Serrano <alfredoblancoserrano@gmail.com>
!              Alfonso Allen-Perkins  <alfonso.allen.perkins@gmail.com>
!              Roberto F. S. Andrade  <randrade@ufba.br>
!================================================================================
!About this module: This compendium of subroutines are used for the purpose
!                   of estimating the mean square displacement or MSD as described
!                   in the main_simulation.f90 file. The subroutines contained in
!                   this module are listed  below:
!
! Initialization    : This code allows the assignment of values to variables that
!                     do not depend on main loop.
! Hops              : Subroutine that decides to which node the random walker will
!                     jump to.
! Prob              : Calculation of the probability normalization factor.
!                     For more details see Eq(4) of the manuscript.
! TH                : Calculation of the neighborhood matrix and the diameter for
!                     a TH of size N odd.!
! Heaviside         : Heaviside step function, or the unit step function.
! TN                : Calculation of the neighborhood matrix and the diameter for
!                     a TN of size N odd.
! Line_TN & Dis_TN  : These codes are used by the TN algorithm to estimate the
!                     neighborhood matrix.
! Rseed             : Random seed is established. If the input value is zero the
!                     seed is random, otherwise the seed will be determined by
!                     the input value.
! Num_Deri          : Calculation of the numerical derivative.
! Data_Reading      : Subroutine that reads the initial data from the
!                     file "in_CAME_data.dat".
! Write_Data        : Writing the results obtained.
!================================================================================
!================================================================================
MODULE module_simulation
!==========================universal variables====================================
INTEGER, PARAMETER :: sp=4,dp=8

REAL(dp),   DIMENSION(:,:),ALLOCATABLE :: MSDT
REAL(dp),   DIMENSION(:),  ALLOCATABLE :: MSD
REAL(dp),   DIMENSION(:),  ALLOCATABLE :: Deri
REAL(dp),   DIMENSION(:),  ALLOCATABLE :: alp
REAL(dp),   DIMENSION(:),  ALLOCATABLE :: distt
REAL(dp),   DIMENSION(:),  ALLOCATABLE :: PRO
REAL(dp),   DIMENSION(:),  ALLOCATABLE :: multiplicity
INTEGER(sp),DIMENSION(:,:),ALLOCATABLE :: MC

REAL(dp)      :: N
INTEGER(sp)   :: Rn,inN,flag_T,flag_alpha,NAS,walkers
INTEGER(sp)   :: ALLOCATESTATUS,dmax,tf
CHARACTER(50) :: data_alpha,input,result

!===================================================================================
CONTAINS
!===================================================================================
SUBROUTINE Initialization
  IMPLICIT NONE
  INTEGER(sp)  :: i
  REAL(dp)     :: x
!==============================Variable initialization==============================

  Rn=INT(SQRT(N));inN=INT(N)

  ALLOCATE(MSDT(tf,2), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(MSD(tf), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(Deri(tf), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(alp(tf), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(MC(inN,inN), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
  ALLOCATE(multiplicity(rn), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(distt(inN), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(PRO(tf), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"

  distt=0;PRO=0.d0;MSD=0.d0
  MSDT=0.d0;Deri=0.d0;multiplicity=0
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
!========================Random seed configuration================================
  CALL RSEED(0)
  x=2.d0;MSDT=1.d0
  DO i=2,tf
    MSDT(i,1)=x
    x=x+1.d0
  END DO
!========================Choose between HT or TN==================================
  IF(flag_T .GT. 0)THEN
    CALL TH(inN,MC,dmax)
  ELSE
    CALL TN(inN,MC,dmax)
  END IF

  CALL PROB(pro)
  RETURN
END SUBROUTINE Initialization
!=================================================================================
SUBROUTINE Hops(ini,dis,new)
  IMPLICIT NONE
  INTEGER(sp),INTENT(IN)  :: ini
  INTEGER(sp),INTENT(IN)  :: dis
  INTEGER(sp),INTENT(OUT) :: new

  INTEGER(sp) :: j,R,h,o
  REAL(dp) :: x2

  o=FLOOR(multiplicity(dis))
  CALL RANDOM_NUMBER(x2)
  r=FLOOR(x2*o)+1
  h=0

  DO j=1,inN
    IF(MC(ini,j) .EQ. dis)THEN
      h=h+1
      IF(h .EQ. R)new=j
    END IF
  END DO

  RETURN
END SUBROUTINE Hops
!=======================================================================================#
SUBROUTINE Prob(p)
  IMPLICIT NONE
  REAL(dp),DIMENSION(:),INTENT(OUT) :: p

  INTEGER(sp) :: i,j
  REAL(dp) :: x
  x=0.d0
  DO i=2,tf
    DO j=2,inN
      x=x+distt(j)**(-alp(i))
    END DO
    P(i)=x
    x=0.d0
  END DO
  RETURN
END SUBROUTINE Prob
!=======================================================================================#
SUBROUTINE TH(nn,MC,dmax)
!===============Adjacency and neighborhood matrix of the TH=======================
  IMPLICIT NONE

  INTEGER(sp),                INTENT(IN)  :: nn
  INTEGER(sp), DIMENSION(:,:),INTENT(OUT) :: MC
  INTEGER(sp),                INTENT(OUT) :: dmax

  INTEGER(sp), DIMENSION(:,:),ALLOCATABLE :: B
  INTEGER(sp), DIMENSION(:,:),ALLOCATABLE :: E
  INTEGER(sp), DIMENSION(:),  ALLOCATABLE :: dist

  REAL(sp)    :: x,jj1
  INTEGER(SP) :: JO2,JO,j,h,h2,sj,xj,i,d1,d2
  INTEGER(SP) :: JJ0,X1,INSQRT,mul

  ALLOCATE(B(nn,nn), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(E(nn,nn), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  ALLOCATE(dist(nn), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"

  B=0;dist=0;E=0;d2=1
  INSQRT=INT(SQRT(N))

  DO j=1,inN

    JJ0=j-((inN+3)/2)
    jj1=((Real(N)+2.0)/2.0)-real(j)
    jj1=sign(1.0,jj1)
    CALL Heaviside(JJ0,X1)
    xj=(inN+2)*X1+(j)*INT(jj1)-1

    sj=inN-j+1
    JO=MOD(xj+1,INSQRT)
    JO2=JO-NINT(((REAL(INSQRT)+3.0)/2.0))
    CALL Heaviside(-JO,H);CALL Heaviside(JO2,h2)

    x=(REAL(xj)/REAL(INSQRT))+((REAL(JO)-1.0)*((REAL(INSQRT)-1.0)/REAL(INSQRT)))

    dist(j)=(NINT(x)+2*H-2*h2*JO2)
  END DO

  DO i=1,inN
    DO j=i,inN
      B(i,j)=dist(j-i+1)
    END DO
  END DO

  DO i=1,inN
    E(i:inN,i)=B(i,i:inN)
  END DO
  E=B+E
  B=0
  DO i=1,inN
    DO j=1,inN
      d1=E(i,j)
      IF(d1 .EQ. 1)B(i,j)=1
      IF(d1 .GT. d2)d2=d1
    END DO
  END DO

  distt(1:nn)=DBLE(E(1,1:nn))
  DO i=1,d2
    DO j=1,nn
      MUL=INT(distt(j))
      IF(mul .EQ. i)multiplicity(i)=multiplicity(i)+1.d0
    END DO
  END DO

  MC=E;dmax=d2
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
  INTEGER(sp),INTENT(IN)   :: x
  INTEGER(sp),INTENT(OUT)  :: y
  y=NINT(0.5*(sign(1,x)+1))
  return
END SUBROUTINE Heaviside
!================================================================================
SUBROUTINE TN(nn,A,dmax)
!===============Adjacency and neighborhood matrix of the TN=======================
  IMPLICIT NONE
  INTEGER(sp),                INTENT(IN)  :: nn
  INTEGER(sp), DIMENSION(:,:),INTENT(OUT) :: A
  INTEGER(sp),                INTENT(OUT) :: dmax

  INTEGER(sp), DIMENSION(:,:), ALLOCATABLE :: B
  INTEGER(sp) :: i,j,d,d1,mul

  ALLOCATE(B(nn,nn), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"

  B=0;d1=1
  DO i=1,nn
    DO j=1,nn
      CALL LINE_TN(i,j,d)
      B(i,j)=d
      IF(d .GT. d1)d1=d
    END DO
  END DO

  distt(1:nn)=DBLE(B(1,1:nn))

  DO i=1,d1
    DO j=1,nn
      MUL=FLOOR(distt(j))
      IF(mul .EQ. i)multiplicity(i)=multiplicity(i)+1.d0
    END DO
  END DO

  A=B;dmax=d1
  DEALLOCATE(B, STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0 )STOP "*** NOT ENOUGH MEMORY ***"
  RETURN
END SUBROUTINE TN
!================================================================================
SUBROUTINE Line_TN(a,b,D)
!========Calculation of the first row of the neighborhood matrix=================
  IMPLICIT NONE
  INTEGER(sp),INTENT(IN)  :: a,b
  INTEGER(sp),INTENT(OUT) :: D

  INTEGER(sp) :: x1,x2,xmax1,xmin1,xmax2,xmin2
  INTEGER(sp) :: a1,i,j1,j2,a2
  INTEGER(sp) :: D1,D2

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
  CALL dis_TN(a1,a2,D1)
  CALL dis_TN(x1,x2,D2)
  D=D1+D2
  RETURN
END SUBROUTINE Line_TN
!=================================================================================
SUBROUTINE Dis_TN(a,b,D)
!=====================Distance between two nodes in the TN=======================
IMPLICIT NONE
INTEGER(sp),INTENT(IN)  :: a,b
INTEGER(sp),INTENT(OUT) :: D

REAL(sp)    :: D1
INTEGER(sp) :: fun,x1

x1=abs(a-b)
D1=(REAL(Rn)/2)
IF(x1 .LT. D1)THEN
  fun=x1
ELSE
  fun=Rn-x1
END IF

D=fun
RETURN
END SUBROUTINE Dis_TN
!=================================================================================
SUBROUTINE Rseed(FLAG)
  IMPLICIT NONE
  INTEGER(sp), INTENT(IN) :: FLAG
  INTEGER(sp), ALLOCATABLE :: SEED(:)
  INTEGER(sp) :: i,j,N

  CALL RANDOM_SEED(size=N)
  ALLOCATE(SEED(N))
  IF(FLAG .EQ. 0)THEN
    DO i=1,N
      CALL SYSTEM_CLOCK(j)
      SEED(i)=j
    END DO
    CALL RANDOM_SEED(PUT = SEED)
  ELSE
    DO i=1,N
      SEED(i)=FLAG
    END DO
    CALL RANDOM_SEED(PUT = SEED)
  END IF
  RETURN
END SUBROUTINE Rseed
!=================================================================================
SUBROUTINE Num_Deri(A,B)
!===========================Log10 Numerical Derivative============================
  IMPLICIT NONE
  REAL(dp),DIMENSION(:,:),INTENT(IN) :: A
  REAL(dp),DIMENSION(:),INTENT(OUT) :: B

  REAL(dp)    :: yo,y1,xo,x1,m,d,c
  INTEGER(sp) :: i

  DO i=1,tf-1
    yo=LOG10(A(i,2))
    y1=LOG10(A(i+1,2))

    xo=LOG10(A(i,1))
    x1=LOG10(A(i+1,1))

    d=(y1-yo);c=(x1-xo)
    m=d/c
    B(i)=m
  END DO
  B(tf)=B(tf-1)
  RETURN
END SUBROUTINE Num_Deri
!=================================================================================
!================================================================================
SUBROUTINE Data_Reading
!================================================================================
  IMPLICIT NONE
  input="in_simulation.dat"
  OPEN(UNIT=10,FILE=input,STATUS='unknown')
  READ(10,*)N,NAS,tf,walkers
  READ(10,*)flag_T,flag_alpha
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
END MODULE module_simulation
