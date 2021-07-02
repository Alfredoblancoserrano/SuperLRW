#!/bin/bash
gfortran -O  -funroll-loops Module.f90 main.f90 -o MSD_simulation.out
