#!/bin/bash

rm -rf Build
mkdir Build

cwd=$(pwd)/src
cd Build
ifort -O2 -c $cwd/kinds.f90
ifort -O2 -mkl -c $cwd/lapack.f90
ifort -O2 -c $cwd/minf.f90
ifort -O2 -mkl -c $cwd/numFort.f90 -lmkl_lapack95_lp64
cd ..
