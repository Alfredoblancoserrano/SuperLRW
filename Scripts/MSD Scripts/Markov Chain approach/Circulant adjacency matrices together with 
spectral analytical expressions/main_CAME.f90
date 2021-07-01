!================================================================================
!Name        : Main_GAM
!Version     : Beta 1.0
!Date        : agosto 03, 2020
!Authors     : Alfredo Blanco Serrano <alfredoblancoserrano@gmail.com>
!              Alfonso Allen-Perkins <alfonso.allen.perkins@gmail.com>
!              Roberto F. S. Andrade <randrade@ufba.br>
!================================================================================
!About this code: We use a Markov chain formalism that considered random walkers
!                 in a 2D tori with two different jumping strategies. The first
!                 is the classical where the walker has the an equal probability
!                 of jumping between adjacent nodes and the second one consideres
!                 a time-dependent probability distribution of long-distance jumps.
!                 The code considers a helical tori(HT) with N odd nodes. The
!                 estimation of the MSD is performed considering the fact that
!                 the adjacency matrix of the tori is circulant and then its
!                 spectral properties are used to derive an expression.
!                 This theoretical development can be found in the supplementary
!                 material of the manuscript.
!
! "Efficient approach to time-dependent super-diffusive Lévy random walks
!            on finite 2D-tori using circulant analogues".
!================================================================================
PROGRAM HELI_MSD 
  USE module_CAME
  IMPLICIT NONE
!============================Variable declaration ===============================
  INTEGER(sp) :: t
  REAL(sp) :: secs
  INTEGER(sp) :: start, finish,computime,count_rate, count_max
!============================data_reading========================================
  CALL data_reading
!===================subroutines for computing time===============================
  CALL system_clock(count_max=count_max, count_rate=count_rate)
  CALL SYSTEM_CLOCK(start)
!==========================Estimation of MSD=====================================
  CALL initialization

  DO t=2,tf
    CALL Eq30(t,xt)
    MSD(t)=satura+xt
  END DO

  MSDT(2:tf,2)=MSD(2:tf)
!==================Estimation of the numerical derivite for the MSD==============
  CALL Num_DERI(MSDT,Deri)
!===================subroutines for computing time===============================
  CALL SYSTEM_CLOCK(finish)
  computime=finish-start
  secs=REAL(computime)/REAL(count_rate)
!=============================Results============================================
  CALL write_data(secs)
END PROGRAM HELI_MSD
!==============================================================================#
