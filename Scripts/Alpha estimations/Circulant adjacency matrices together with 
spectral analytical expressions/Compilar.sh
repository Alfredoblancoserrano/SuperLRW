#!/bin/bash
gfortran -O -fopenmp -funroll-loops helicoidal_module.f90 helicoidal_main.f90 -o Ver1.out
