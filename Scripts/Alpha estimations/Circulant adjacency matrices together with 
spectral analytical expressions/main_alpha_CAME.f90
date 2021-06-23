!================================================================================
!Name        : main_alpha_CAME
!Version     : Beta 1.0
!Date        : agosto 03, 2020
!Authors     : Alfredo Blanco Serrano <alfredoblancoserrano@gmail.com>
!              Alfonso Allen-Perkins <alfonso.allen.perkins@gmail.com>
!              Roberto F. S. Andrade <randrade@ufba.br>
!================================================================================
!About this code: This code implements a solution of the inverse problem
!                 by considering the spectral properties of a circulant matrix
!                 as the adjacency of a helical tori. Thus this program finds
!                 the set of alpha(t) values that allow to generate a
!                 superdiffusive regime in a helical tori using Lévy random
!                 walks (LRW). The results and the analysis of this procedure
!                 are presented in detail in the manuscript,
!
! "Efficient approach to time-dependent super-diffusive Lévy random walks
!            on finite 2D-tori using circulant analogues".
!================================================================================
PROGRAM alpha_CAME
  USE module_alpha_CAME
  IMPLICIT NONE
!============================Variable declaration ===============================
  REAL(dp)    :: max,min,control,MSD
  REAL(sp)    :: secs
  INTEGER(sp) :: start, finish,computime,count_rate, count_max
  INTEGER(sp) :: t,i
!============================data_reading========================================
  CALL Data_Reading
!===================subroutines for computing time===============================
  CALL System_Clock(count_max=count_max, count_rate=count_rate)
  CALL System_Clock(start)
!==========================Estimation of alpha(t)=================================
  CALL Initialization

  alp2=100.d0;alp(1)=1000.d0
  DO t=2,tf
    max=alp2;min=0.d0
    DO i=1,1000
      CALL Eq35(MSD)
      control=EMSD(t)-MSD
      IF(ABS(control) .LT. tol)THEN
        x1=x3
        alp(t)=alp2
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
  CALL System_Clock(finish)
  computime=finish-start
  secs=REAL(computime)/REAL(count_rate)
!=============================Results============================================
  CALL Write_Data(secs)
END PROGRAM alpha_CAME
!==============================================================================#
