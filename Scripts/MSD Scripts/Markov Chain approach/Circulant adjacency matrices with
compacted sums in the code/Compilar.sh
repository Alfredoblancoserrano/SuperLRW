#!/bin/bash
gfortran -O -funroll-loops module_CAMS.f90 main_CAMS.f90 -o a.out
