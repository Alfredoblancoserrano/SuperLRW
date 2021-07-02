!================================================================================
!Name        : Main_GAM
!Version     : Beta 1.0
!Authors     : Alfredo Blanco Serrano <alfredoblancoserrano@gmail.com>
!              Alfonso Allen-Perkins  <alfonso.allen.perkins@gmail.com>
!              Roberto F. S. Andrade  <randrade@ufba.br>
!================================================================================
!About this code: We use a Markov chain formalism that considered random walkers
!                 in a 2D tori with two different jumping strategies. The first
!                 is the classical where the walker has the an equal probability
!                 of jumping between adjacent nodes and the second one consideres
!                 a time-dependent probability distribution of long-distance jumps.
!                 This code can be used either to the normal tori(NT) or for
!                 the helical tori(HT). A theoretical description for this
!                 formalism and the results can be found in the  manuscript:
!
! "Efficient approach to time-dependent super-diffusive Lévy random walks
!            on finite 2D-tori using circulant analogues".
!================================================================================
PROGRAM MSD_GAM
  USE module_GAM
  IMPLICIT NONE
!============================Variable declaration ===============================
  REAL(sp) :: secs
  REAL(dp),DIMENSION(:,:),ALLOCATABLE :: PTI
  REAL(dp),DIMENSION(:),ALLOCATABLE :: MSD
  INTEGER(sp) :: start, finish,computime,count_rate, count_max
  INTEGER(sp) :: t,i,j
!============================data_reading========================================
  CALL Data_Reading

  ALLOCATE(PTI(N,N), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
  ALLOCATE(MSD(tf), STAT=ALLOCATESTATUS)
  IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
!===================subroutines for computing time===============================
  CALL System_Clock(count_max=count_max, count_rate=count_rate)
  CALL System_Clock(start)
!==========================Estimation of MSD=====================================

  MSD=0.0;MSD(1)=1.d0
  CALL Initialization
  CALL ST(PTR,PTI)
  PTR=PTI


  DO t=2,tf
    CALL Melli(ADJA,MC,dmax,alp,t,TM)
    CALL ST(TM,PTI)
    PTR=MATMUL(PTR,PTI)
    DO i=1,n
      DO j=1,n
        MSD(t)=MSD(t)+(MC(i,j)*MC(i,j)*PTR(j,i))
      END DO
    END DO
  END DO
  MSDT(2:tf,2)=MSD(2:tf)/REAL(N)
!==================Estimation of the numerical derivite for the MSD==============
  CALL Num_Deri(MSDT,deri)
!===================subroutines for computing time===============================
  CALL System_Clock(finish)
  computime=finish-start
  secs=REAL(computime)/real(count_rate)
!=============================Results============================================
  CALL Write_Data(secs)
END PROGRAM MSD_GAM
!===========================================================================
