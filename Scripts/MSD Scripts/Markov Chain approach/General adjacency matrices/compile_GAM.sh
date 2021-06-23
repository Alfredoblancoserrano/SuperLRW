#!/usr/bin/env bash
#===========================================================================
gfortran -O3 -funroll-loops module_GAM.f90 Main_GAM.f90 -o TH1225.out
