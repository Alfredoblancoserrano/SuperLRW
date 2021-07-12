!==============================================================================#
!================================================================================
!Name        : Main_CAMS
!Version     : Beta 1.0
!Authors     : Alfredo Blanco Serrano <alfredoblancoserrano@gmail.com>
!              Alfonso Allen-Perkins  <alfonso.allen.perkins@gmail.com>
!              Roberto F. S. Andrade  <randrade@ufba.br>
!================================================================================
!About this code: In this code we use a Markov chain formalism that considered
!                 random walkers in a 2D tori with two different jumping strategies.
!                 The first is the classical where the walker has the an equal
!                 probability of jumping between adjacent nodes and the second
!                 one consideres a time-dependent probability distribution of
!                 long-distance jumps. The code considers a helical tori(HT)
!                 with N odd nodes. The estimation of the MSD is performed
!                 considering the fact that the adjacency matrix of the tori is
!                 circulant and therefore the number of loops is considerably
!                 reduced.his approach significantly reduces the normal
!                 computational complexity of the MSD calculation. The results
!                 of this coding can be found in the manuscript:
!
! "Efficient approach to time-dependent super-diffusive Lévy random walks
!            on finite 2D-tori using circulant analogues".
!================================================================================
PROGRAM MSD_CAMS
USE module_CAMS
IMPLICIT NONE
!============================Variable declaration ===============================
  REAL(sp)    :: secs
  INTEGER(sp) :: start, finish,computime,count_rate, count_max,j,t
!============================data_reading========================================
  CALL Data_Reading
!===================subroutines for computing time===============================
  CALL System_Clock(count_max=count_max, count_rate=count_rate)
  CALL System_Clock(start)
!==========================Estimation of MSD=====================================
  CALL Initialization

  DO t=2,tf
    CALL PV(t,P)
    DO j=2,dmax
     MSD(t)=MSD(t)+(SNV(j)*p(j))
    END DO
  END DO

  MSDT(2:tf,2)=2.0*MSD(2:tf)
!==================Estimation of the numerical derivite for the MSD==============
  CALL Num_Deri(MSDT,Deri)
!===================subroutines for computing time===============================
  CALL System_Clock(finish)
  computime=finish-start
  secs=REAL(computime)/real(count_rate)
!=============================Results============================================
  CALL Write_Data(secs)
!================================================================================
END PROGRAM MSD_CAMS
