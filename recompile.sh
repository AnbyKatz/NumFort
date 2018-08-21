#!/bin/bash

ifort -O2 -c kinds.f90
ifort -O2 -mkl -c lapack.f90
ifort -O2 -c PLplots.f90 -I/home/anthony/bin/PLplot/install_dir/include/plplot -I/home/anthony/bin/PLplot/install_dir/lib/fortran/modules/plplot -L/home/anthony/bin/PLplot/install_dir/lib -lplplotfortran -lplplot
ifort -O2 -mkl  -c numFort.f90 -lmkl_lapack95_lp64 
