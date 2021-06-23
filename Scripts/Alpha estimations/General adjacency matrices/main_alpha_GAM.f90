!================================================================================
!Name        : Main_GAM
!Version     : Beta 1.0
!Date        : agosto 03, 2020
!Authors     : Alfredo Blanco Serrano <alfredoblancoserrano@gmail.com>
!              Alfonso Allen-Perkins <alfonso.allen.perkins@gmail.com>
!              Roberto F. S. Andrade <randrade@ufba.br>
!================================================================================
!About this code: The implementation of this code allows to solve the inverse
!                 problem for a normal tori or for a helical tori. This code is
!                 based on the Markov framework presented in detail in the
!                 sections 2.1 and 2.2 of the manuscript,
!
! "Efficient approach to time-dependent super-diffusive Lévy random walks
!            on finite 2D-tori using circulant analogues".
!================================================================================
PROGRAM alpha_GAM
  USE module_alpha_GAM
  IMPLICIT NONE
!============================Variable declaration ===============================
  REAL(dp),DIMENSION(:,:),ALLOCATABLE :: PTI
  REAL(dp),DIMENSION(:,:),ALLOCATABLE :: PTRR

  REAL(sp)    :: secs
  REAL(dp)    :: max,min,control,MSD
  INTEGER(sp) :: start, finish,computime,count_rate, count_max
  INTEGER(sp) :: t,i
!============================data_reading========================================
  CALL data_reading

  ALLOCATE(PTI(N,N), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
  ALLOCATE(PTRR(N,N), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"

!===================subroutines for computing time===============================
  CALL system_clock(count_max=count_max, count_rate=count_rate)
  CALL SYSTEM_CLOCK(start)
!==========================Estimation of MSD=====================================

  CALL Initialization
  CALL Transision_matrix(PTR,PTI)
  PTR=PTI;alp(1)=1000.d0
  alP2=100.d0

  DO t=2,1000
    max=alP2;min=0.d0
    DO i=1,1000
      CALL GAM_MSD(PTI,PTRR,MSD)
      control=EMSD(t)-MSD
      IF(ABS(control) .LT. tol)THEN
        alp(t)=alp2
        PTR=PTRR
        write(*,*)alp(t)
        EXIT
      END IF
      CALL Bisection(control,min,max,alp2)
    END DO
    IF(alp(t) .LE. 1.d-200)THEN
      tf=t
      EXIT
    END IF
  END DO
!===================subroutines for computing time===============================
  CALL SYSTEM_CLOCK(finish)
  computime=finish-start
  secs=REAL(computime)/real(count_rate)
!=============================Results============================================
  CALL write_data(secs)
END PROGRAM alpha_GAM
!===========================================================================
