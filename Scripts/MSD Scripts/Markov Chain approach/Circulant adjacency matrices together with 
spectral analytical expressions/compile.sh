#!/bin/bash
gfortran -O -fopenmp -funroll-loops module_CAME.f90 main_CAME.f90 -o TH_cir_49.out
  
