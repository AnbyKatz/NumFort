#!/bin/bash

rm -rf build
mkdir build

cwd=$(pwd)/src
cd build
ifort -O2 -c $cwd/kinds.f90
ifort -O2 -mkl -c $cwd/lapack.f90
ifort -O2 -c $cwd/minf.f90
ifort -O2 -mkl -c $cwd/numFort.f90 -lmkl_lapack95_lp64
cd ..
