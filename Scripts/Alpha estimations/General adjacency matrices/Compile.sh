#!/usr/bin/env bash
#===========================================================================
gfortran -O3 -funroll-loops module_caminante_normal.f90 main_caminante_normal.f90 -o NMC_49.out
