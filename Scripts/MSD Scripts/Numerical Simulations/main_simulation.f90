!==============================================================================#
!================================================================================
!Name        : main_simulation
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
!                 This code performs a simulation of random walkers on a helical
!                 tori(HT) and a normal tori (NT) with N odd nodes.
!                 Note that if you want the MSD to follow a super-diffusion
!                 state you must give to the program the values of alpha(t)
!                 in a file. For more details see the manuscript,
!
! "Efficient approach to time-dependent super-diffusive Lévy random walks
!            on finite 2D-tori using circulant analogues".
!================================================================================
PROGRAM MSD_simulation
USE module_simulation
IMPLICIT NONE
!============================Variable declaration ===============================
  REAL(sp)    :: secs
  REAL(dp)    :: x,test,p
  INTEGER(sp) :: t,new,k,jini,i,d
  INTEGER(sp) :: start, finish,computime,count_rate, count_max
!============================data_reading========================================
  CALL Data_Reading
!===================subroutines for computing time===============================
  CALL System_Clock(count_max=count_max, count_rate=count_rate)
  CALL System_Clock(start)
!==========================Estimation of MSD=====================================
  CALL Initialization
  DO i=1,walkers
    jini=1
    new=jini
    MSD(1)=MSD(1)+1.d0
    CALL Hops(jini,1,new)

    DO t=2,tf
      k=new
  1000 CONTINUE
      CALL RANDOM_NUMBER(x)
      d=int(x*Dmax)+1
      p=multiplicity(d)/(DBLE(d)**alp(t))
      p=p/PRO(t)

      CALL RANDOM_NUMBER(X)
      IF (X .LT. p)THEN
        CALL Hops(k,d,new)
        d=mc(jini,new)
        MSD(t)=MSD(t)+DBLE(d*d)
        test=MSD(t)
      ELSE
        GO TO 1000
      END IF

      IF(abs(test-MSD(t-1)) .LE. 1.d-4)THEN
        MSD(t+1:tf)=test
        EXIT
      END IF
    END DO
  END DO

  MSD(1:tf)=MSD(1:tf)/DBLE(walkers)
  MSDT(1:tf,2)=MSD(1:tf)
!==================Estimation of the numerical derivite for the MSD==============
  CALL Num_Deri(MSDT,Deri)
!===================subroutines for computing time===============================
  CALL System_Clock(finish)
  computime=finish-start
  secs=REAL(computime)/real(count_rate)
!=============================Results============================================
  CALL Write_Data(secs)
!==============================================================================#
END PROGRAM MSD_simulation
