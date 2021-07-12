#!/usr/bin/env bash
#===========================================================================
gfortran -O -funroll-loops module_GAM.f90 Main_GAM.f90 -o MSD_GAM.out
