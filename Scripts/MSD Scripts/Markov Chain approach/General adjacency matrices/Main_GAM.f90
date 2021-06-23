!================================================================================
!Name        : Main_GAM
!Version     : Beta 1.0
!Date        : agosto 03, 2020
!Authors     : Alfredo Blanco Serrano <alfredoblancoserrano@gmail.com>
!              Alfonso Allen-Perkins <alfonso.allen.perkins@gmail.com>
!              Roberto F. S. Andrade <randrade@ufba.br>
!================================================================================
!About this code: This code, written in fortran, computes the mean square distan-
!                 ce or MSD using the Markov chain formalism. In this code
!                 random walkers are considered in a 2D tori with the classical
!                 strategy of jumping with equal probability between adjacent no-
!                 des, and also, the strategy of long distance jumps using a
!                 set of values alpha(t). This 2D tori can be either the normal
!                 tori (NT) or the helical tori (HT). A theoretical description
!                 of this formalism and the results can be found in the
!                 manuscript:

! "Efficient approach to time-dependent super-diffusive Lévy random walks
!            on finite 2D-tori using circulant analogues".
!================================================================================
PROGRAM caminantenormal

USE module_GAM
IMPLICIT NONE
!============================Variable declaration ===============================
REAL(sp) :: secs
REAL(dp),DIMENSION(:,:),ALLOCATABLE :: PTI
REAL(dp),DIMENSION(:),ALLOCATABLE :: MSD
INTEGER(sp) :: start, finish,computime,count_rate, count_max
INTEGER(sp) :: t,i,j
!============================data_reading=======================================
CALL data_reading

ALLOCATE(PTI(N,N), STAT=ALLOCATESTATUS)
IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
ALLOCATE(MSD(tf), STAT=ALLOCATESTATUS)
IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
!===================subroutines for computing time===============================
call system_clock(count_max=count_max, count_rate=count_rate)
CALL SYSTEM_CLOCK(start)
!==========================Estimation of MSD=====================================

MSD=0.0;MSD(1)=1.d0
CALL initialization
CALL S_GRADE(PTR,PTI)
PTR=PTI


DO t=2,tf
  CALL MELLI(ADJA,MC,dmax,alp,t,TM)
  CALL S_GRADE(TM,PTI)
  PTR=MATMUL(PTR,PTI)

  DO i=1,n
    DO j=1,n
      MSD(t)=MSD(t)+(MC(i,j)*MC(i,j)*PTR(j,i))
    END DO
  END DO
END DO

MSDT(2:tf,2)=MSD(2:tf)/REAL(N)
!==================Estimation of the numerical derivite for the MSD==============
CALL Num_DERI(MSDT,deri)
DEALLOCATE(PTI, STAT=ALLOCATESTATUS)
IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
DEALLOCATE(MSD, STAT=ALLOCATESTATUS)
IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
!=============================Results============================================
CALL SYSTEM_CLOCK(finish)
computime=finish-start
secs=REAL(computime)/real(count_rate)
CALL write_data(secs)
END PROGRAM caminantenormal
!===========================================================================
