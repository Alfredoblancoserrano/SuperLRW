!===========================================================================
!Name        :
!Description :
!Version     : Beta1.0
!Date        : agosto 03, 2020
!Authors     : Alfredo Blanco Serrano <alfredoblancoserrano@gmail.com>
!              Alfonso Allen-Perkins <alfonso.allen.perkins@gmail.com>
!              Roberto F. S. Andrade <randrade@ufba.br>
!===========================================================================
!About this code:
!===========================================================================
PROGRAM caminantenormal

USE module_caminante_normal

IMPLICIT NONE
!===================Declaración de variables de entrada=====================
REAL(sp) :: secs
REAL(dp),DIMENSION(:,:),ALLOCATABLE :: PTI
REAL(dp),DIMENSION(:),ALLOCATABLE :: MSD
INTEGER(sp) :: start, finish,computime,count_rate, count_max
INTEGER(sp) :: t,i,j

CALL data_reading

ALLOCATE(PTI(N,N), STAT=ALLOCATESTATUS)
IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
ALLOCATE(MSD(tf), STAT=ALLOCATESTATUS)
IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
!===================Llamado de subrutinas en el modulo======================
call system_clock(count_max=count_max, count_rate=count_rate)
CALL SYSTEM_CLOCK(start)
!===================Llamado de subrutinas en el modulo======================

MSD=0.0;MSD(1)=1.d0
CALL initialization
CALL S_GRAU(PTR,PTI)
PTR=PTI


DO t=2,tf
  CALL MELLI(ADJA,MC,dmax,alp,t,TM)
  CALL S_GRAU(TM,PTI)
  PTR=MATMUL(PTR,PTI)

  DO i=1,n
    DO j=1,n
      MSD(t)=MSD(t)+(MC(i,j)*MC(i,j)*PTR(j,i))
    END DO
  END DO
END DO

MSDT(2:tf,2)=MSD(2:tf)/REAL(N)
CALL DERIVADA(MSDT,deri)
DEALLOCATE(PTI, STAT=ALLOCATESTATUS)
IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
DEALLOCATE(MSD, STAT=ALLOCATESTATUS)
IF(ALLOCATESTATUS .NE. 0)STOP "***NOT ENOUGH MEMORY ***"
!===================Llamado de subrutinas en el modulo======================
CALL SYSTEM_CLOCK(finish)
computime=finish-start
secs=REAL(computime)/real(count_rate)
CALL write_data(secs)
END PROGRAM caminantenormal
!===========================================================================
