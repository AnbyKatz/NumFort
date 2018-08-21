#!/bin/bash

# Example recompile file
# type chmod +x recompileEx.sh and then ./recompileEX.sh to compile

ifort -O2 -c kinds.f90
ifort -O2 -mkl -c lapack.f90

# Comment out this line if you wish to not use PLplot/haven't installed it
ifort -O2 -c PLplots.f90 -I/path/to/PLplotInstall/include/plplot -I/path/to/PLplotInstall/lib/fortran/modules/plplot -L/path/to/PLplotInstall/lib -lplplotfortran -lplplot

ifort -O2 -mkl  -c numFort.f90 -lmkl_lapack95_lp64 
