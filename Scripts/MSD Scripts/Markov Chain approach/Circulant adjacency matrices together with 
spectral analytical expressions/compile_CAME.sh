#!/bin/bash
gfortran -O -funroll-loops module_CAME.f90 main_CAME.f90 -o MSD_CAME.out
