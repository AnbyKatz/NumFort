#!/bin/bash

ifort -O2 -c kinds.f90
ifort -O2 -mkl -c lapack.f90
ifort -O2 -xHost -mkl -I/home/anthony/bin/PLplot/install_dir/include/plplot -I/home/anthony/bin/PLplot/install_dir/lib/fortran/modules/plplot -c numFort.f90 -lmkl_lapack95_lp64 -L/home/anthony/bin/PLplot/install_dir/lib -lplplotfortran -lplplot
